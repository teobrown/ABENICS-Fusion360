"""
ABENICS gear set generator for Autodesk Fusion

Builds, in a new design:
  - the CS gear (cross spherical gear, the ball)
  - MP gear A and MP gear B (monopole gears), each cut by its own SH cutter
    and placed in its meshing position around the ball

Math follows Abe, Tadakuma and Tadakuma, "ABENICS: Active Ball Joint Mechanism
With Three-DoF Based on Spherical Gear Meshings", IEEE T-RO 37(5), 2021:
  d_sph = m * z_sph                      (1)
  z_mpl = z_sph / 2                      (5)
  R_cr  = m * (z_sph + z_mpl) / 2        (6)  carrier length = center distance
  psi_hb = psi_cr * z_mpl / z_sph        (7)  cutter spin relative to the carrier
  Constraint 1: z_sph is even
  Constraint 2: the x and y revolve axes sit on the center of a tooth peak or valley

Frame used here (same as the paper):
  - the ball is centered on the origin
  - tooth family A is the spur profile in the XY plane revolved about X
  - tooth family B is the same profile revolved about Y
  - Z is the ball's free third axis
  - MP gear A sits on +X, MP gear B sits MODULE_SPACING_DEG further around Z.
    Both MP axes are parallel to Z. At this pose family A's axis points at
    MP gear A and family B's axis points at MP gear B (pitch = 0 for both).

How the MP gear is cut (paper, section on MP gear generation):
  the blank stays still, the SH cutter (the ball profile revolved about ONE
  axis) rides a carrier of length R_cr around it. At carrier angle psi the
  cutter has turned psi + psi_hb = 1.5 psi in total, and is subtracted from
  the blank. One full carrier lap cuts every tooth, including the pole.

All lengths in this file are millimetres. Fusion's API uses centimetres, so
values are divided by 10 only where they are handed to Fusion.
"""

import math
import time
import traceback

try:
    import adsk.core
    import adsk.fusion
except ImportError:          # lets the math below be imported and tested outside Fusion
    adsk = None


# ============================ SETTINGS ============================
# Module and tooth count are asked for when the script runs. These are the
# defaults shown in those prompts.
DEFAULT_MODULE_MM = 2.0
DEFAULT_TEETH = 32           # CS gear teeth, must be even. MP gears get half.
DEFAULT_STEPS = 120          # cutter positions per MP gear over one carrier lap

# 3D printing clearance, also asked for at run time. Every tooth flank on the
# ball and on the MP gears is pulled back by this much, measured normal to the
# flank, so mating flanks end up 2x this far apart. Common FFF guidance is
# 0.1 to 0.2 mm per gear. Enter 0 for exact CAD gears (no print clearance).
DEFAULT_PRINT_CLEARANCE_MM = 0.15

PRESSURE_ANGLE_DEG = 20.0
POLE_FEATURE = 'tooth'       # what sits on the ball's X axis: 'tooth' or 'valley' (both meet constraint 2)
ROOT_FILLET = 0.38           # CS root fillet radius, x module. Clipped to what fits below the base circle.

# The paper's footnote says the SH cutter's addendum and dedendum "need to be
# adjusted". Making the cutter teeth longer than the ball's teeth by this much
# (x module) gives the MP gear a standard 0.25 m root clearance.
CUTTER_EXTRA_ADDENDUM = 0.25

# Only used when the print clearance is 0. Negative values cut the MP gear
# deeper (the cutter rides closer), which adds backlash and trims the MP tips by
# the same amount. -0.1 mm is the rack shift listed for the paper's machined
# PEEK prototype and gives about 0.07 mm of backlash at the MP pitch circle.
MP_PROFILE_SHIFT_MM = -0.1

MP_FACE_WIDTH = 0.75         # MP gear width along its axis, x MP pitch diameter
MP_BLANK = 'throat'          # 'throat': tips follow the ball at constant clearance
                             # 'cylinder': plain cylinder blank, teeth get shallower toward the faces
MP_BORE_DIAMETER_MM = 0.0    # 0 = no bore
MODULE_SPACING_DEG = 90.0    # angle between the MP gears around the ball. 90 = paper's perpendicular type.
# ==================================================================


# ------------------------------------------------------------------
# Plain gear math (no Fusion calls)
# ------------------------------------------------------------------

def inv(a):
    """Involute function."""
    return math.tan(a) - a


def polar(r, a):
    return (r * math.cos(a), r * math.sin(a))


class Profile2D:
    """One spur gear outline in the XY plane: involute flanks, tip arcs,
    radial lines below the base circle, root fillets, root arcs."""

    def __init__(self, module, teeth, alpha_deg, tip_addendum, root_dedendum, fillet, flank_offset=0.0,
                 axis_side_extra=0.0, spline_tips=False):
        # flank_offset (mm): moves every involute flank along its own normal.
        # Negative thins the teeth, positive thickens them. For an involute,
        # a normal offset d is exactly a rotation of the flank by d / r_base.
        # axis_side_extra (mm): moves only the flanks that face the nearer end
        # of the revolve axis (angle 0 or 180 deg in half_outline) out by this
        # much more. Used for the CS gear, see make_cs_body.
        self.flank_offset = float(flank_offset)
        self.extra = float(axis_side_extra)
        # spline_tips: draw the tip lands as splines through points on the tip
        # circle instead of arcs. Used for the SH cutter, see SH_CUTTER note.
        self.spline_tips = bool(spline_tips)
        self.m = float(module)
        self.z = int(teeth)
        self.alpha = math.radians(alpha_deg)
        self.p = 2.0 * math.pi / self.z
        self.r = self.m * self.z / 2.0
        self.rb = self.r * math.cos(self.alpha)
        self.ra = self.r + tip_addendum * self.m
        self.rf = self.r - root_dedendum * self.m
        self.notes = []

        # half tooth angle at the base circle
        self.th_b = math.pi / (2.0 * self.z) + inv(self.alpha)

        # Pointed-tip guard: keep at least a small land on the tip.
        min_land = 0.02 * self.p
        if self.half_angle(self.ra) < min_land:
            lo, hi = max(self.rb, self.rf), self.ra
            for _ in range(80):
                mid = 0.5 * (lo + hi)
                if self.half_angle(mid) > min_land:
                    lo = mid
                else:
                    hi = mid
            self.notes.append('tip radius reduced from %.3f to %.3f mm (tooth came to a point)'
                              % (self.ra, lo))
            self.ra = lo
        self.ha = self.half_angle(self.ra)

        # Lowest point of the involute and the flank angle there
        self.r_inv0 = max(self.rb, self.rf)
        self.th_root = self.half_angle(self.r_inv0)
        self.gap = self.p / 2.0 - self.th_root - self.extra / self.rb   # flank bottom to space center
        if self.gap <= 0:
            raise ValueError('Teeth overlap at the root circle. Use more teeth or a smaller pressure angle.')

        # Root fillet between the radial line and the root circle
        rho = fillet * self.m
        if self.rf < self.rb and rho > 0:
            # biggest fillet that fits below the base circle: its top end then
            # sits on the base circle, tangent to the start of the involute
            rho_line = (self.rb ** 2 - self.rf ** 2) / (2.0 * self.rf)
            s = math.sin(0.9 * self.gap)
            rho_gap = self.rf * s / (1.0 - s)
            rho_fit = min(rho, rho_line, rho_gap)
            # A very short radial line between the fillet and the involute turns
            # into a ring face a few hundredths of a mm wide, and with those the
            # CS intersection failed in Fusion (z = 32, clearance 0.15). Grow the
            # fillet to the base circle instead.
            if (self.rb - math.sqrt(self.rf ** 2 + 2.0 * self.rf * rho_fit) < 0.05 * self.m
                    and rho_line <= rho_gap):
                rho_fit = rho_line
            if rho_fit < 0.05 * self.m:
                self.notes.append('no root fillet (only %.3f mm fits below the base circle)' % rho_fit)
                rho_fit = 0.0
            elif rho_fit < rho - 1e-9:
                self.notes.append('root fillet clipped from %.3f to %.3f mm to fit below the base circle'
                                  % (rho, rho_fit))
            elif rho_fit > rho + 1e-9:
                self.notes.append('root fillet grown from %.3f to %.3f mm to reach the base circle'
                                  % (rho, rho_fit))
            rho = rho_fit
        else:
            if rho > 0:
                self.notes.append('no root fillet (root circle is above the base circle)')
            rho = 0.0
        self.rho = rho
        if rho > 0:
            self.delta = math.asin(rho / (self.rf + rho))
            self.r_lt = math.sqrt(self.rf ** 2 + 2.0 * self.rf * rho)   # fillet meets the radial line here
            if rho == rho_line:
                self.r_lt = self.rb
        else:
            self.delta = 0.0
            self.r_lt = self.rf

    def half_angle(self, rho, extra=0.0):
        """Half the angular tooth thickness at radius rho (rho >= base radius),
        for a flank moved out by `extra` mm beyond flank_offset."""
        rho = max(rho, self.rb)
        phi = math.acos(min(1.0, self.rb / rho))
        return math.pi / (2.0 * self.z) + inv(self.alpha) - inv(phi) + (self.flank_offset + extra) / self.rb

    # Which flanks get axis_side_extra. Tooth centers are never closer than half
    # a pitch to 90 deg unless they sit exactly on it.
    def _extra_ccw(self, c):
        """Extra for the counterclockwise-facing flank of the tooth centered at c."""
        return self.extra if c >= 0.5 * math.pi - 1e-9 else 0.0

    def _extra_cw(self, c):
        """Extra for the clockwise-facing flank of the tooth centered at c."""
        return self.extra if c <= 0.5 * math.pi + 1e-9 else 0.0

    # ---- segments ----
    # ('spline', [pts]), ('line', p0, p1), ('arc', p0, pmid, p1)

    def tip_arc(self, a0, a1):
        if self.spline_tips:
            n = 6
            return ('spline', [polar(self.ra, a0 + (a1 - a0) * k / float(n)) for k in range(n + 1)])
        return ('arc', polar(self.ra, a0), polar(self.ra, 0.5 * (a0 + a1)), polar(self.ra, a1))

    def tooth_tip(self, c):
        return self.tip_arc(c - self.half_angle(self.ra, self._extra_cw(c)),
                            c + self.half_angle(self.ra, self._extra_ccw(c)))

    def space_first_half(self, c, extra=0.0, n_inv=12):
        """From the left tip corner of the tooth centered at angle c, down the
        flank, around the fillet, along the root to the space center."""
        segs = []
        pts = []
        for i in range(n_inv + 1):
            # denser near the base circle where curvature is highest
            t = i / float(n_inv)
            rr = self.ra - (self.ra - self.r_inv0) * (1.0 - (1.0 - t) ** 2)
            pts.append(polar(rr, c + self.half_angle(rr, extra)))
        segs.append(('spline', pts))
        a_flank = c + self.half_angle(self.r_inv0, extra)
        if self.rf < self.rb:
            if self.rb - self.r_lt > 1e-6:
                segs.append(('line', polar(self.rb, a_flank), polar(self.r_lt, a_flank)))
            if self.rho > 0:
                ac = a_flank + self.delta
                cx, cy = polar(self.rf + self.rho, ac)
                t1 = polar(self.r_lt, a_flank)
                t2 = polar(self.rf, ac)
                u1 = ((t1[0] - cx) / self.rho, (t1[1] - cy) / self.rho)
                u2 = ((t2[0] - cx) / self.rho, (t2[1] - cy) / self.rho)
                bx, by = u1[0] + u2[0], u1[1] + u2[1]
                bl = math.hypot(bx, by)
                pm = (cx + self.rho * bx / bl, cy + self.rho * by / bl)
                segs.append(('arc', t1, pm, t2))
            a_root0 = a_flank + self.delta
        else:
            a_root0 = a_flank
        a_mid = c + self.p / 2.0
        if a_mid - a_root0 > 1e-9:
            segs.append(('arc', polar(self.rf, a_root0), polar(self.rf, 0.5 * (a_root0 + a_mid)),
                         polar(self.rf, a_mid)))
        return segs

    def space_second_half(self, c, extra=0.0):
        """Mirror of space_first_half about the space center: from the space
        center up the right flank of the next tooth to its tip corner."""
        beta = c + self.p / 2.0
        c2, s2 = math.cos(2 * beta), math.sin(2 * beta)

        def refl(q):
            return (q[0] * c2 + q[1] * s2, q[0] * s2 - q[1] * c2)

        out = []
        for seg in reversed(self.space_first_half(c, extra)):
            if seg[0] == 'spline':
                out.append(('spline', [refl(q) for q in reversed(seg[1])]))
            elif seg[0] == 'line':
                out.append(('line', refl(seg[2]), refl(seg[1])))
            else:
                out.append(('arc', refl(seg[3]), refl(seg[2]), refl(seg[1])))
        return out

    def half_outline(self, feature_at_zero):
        """Closed outline of the half of the gear between angle 0 and pi
        (y >= 0), with a tooth or a valley centered on angle 0. Returns the
        segment list (CCW) and the closing line along the X axis."""
        p, zh = self.p, self.z // 2
        segs = []
        if feature_at_zero == 'tooth':
            segs.append(self.tip_arc(0.0, self.half_angle(self.ra, self._extra_ccw(0.0))))
            for k in range(zh):
                c = k * p
                segs += self.space_first_half(c, self._extra_ccw(c))
                segs += self.space_second_half(c, self._extra_cw(c + p))
                if k + 1 < zh:
                    segs.append(self.tooth_tip((k + 1) * p))
                else:
                    segs.append(self.tip_arc(math.pi - self.half_angle(self.ra, self._extra_cw(math.pi)), math.pi))
            r_axis = self.ra
        else:
            segs += self.space_second_half(-p / 2.0, self._extra_cw(p / 2.0))
            for k in range(zh):
                c = p / 2.0 + k * p
                segs.append(self.tooth_tip(c))
                segs += self.space_first_half(c, self._extra_ccw(c))
                if k + 1 < zh:
                    segs += self.space_second_half(c, self._extra_cw(c + p))
            r_axis = self.rf
        segs = snap_ends(segs, (r_axis, 0.0), (-r_axis, 0.0))
        closing = ('line', (-r_axis, 0.0), (r_axis, 0.0))
        return segs, closing


def _replace_point(seg, which, q):
    if seg[0] == 'spline':
        pts = list(seg[1])
        pts[0 if which == 0 else -1] = q
        return ('spline', pts)
    lst = list(seg)
    lst[1 if which == 0 else len(lst) - 1] = q
    return tuple(lst)


def snap_ends(segs, first, last):
    """Put the first and last points exactly on the revolve axis."""
    segs = list(segs)
    segs[0] = _replace_point(segs[0], 0, first)
    segs[-1] = _replace_point(segs[-1], 1, last)
    return segs


def rotate_seg(seg, ang):
    c, s = math.cos(ang), math.sin(ang)

    def rot(q):
        return (q[0] * c - q[1] * s, q[0] * s + q[1] * c)
    if seg[0] == 'spline':
        return ('spline', [rot(q) for q in seg[1]])
    return tuple([seg[0]] + [rot(q) for q in seg[1:]])


def other_feature(f):
    return 'valley' if f == 'tooth' else 'tooth'


def family_outline(prof, family, pole_feature):
    """Half outline to revolve for tooth family 'A' (about X, uses y >= 0) or
    'B' (about Y, uses x >= 0). Both families are the SAME XY-plane profile;
    for B the feature on +Y is a tooth when z is a multiple of 4 and the
    opposite feature when z = 2 mod 4."""
    if family == 'A':
        return prof.half_outline(pole_feature)
    feat = pole_feature if prof.z % 4 == 0 else other_feature(pole_feature)
    segs, closing = prof.half_outline(feat)
    # local angle 0 -> global -90 deg, local 180 -> global +90 deg
    r_axis = closing[2][0]
    segs = [rotate_seg(s, -math.pi / 2.0) for s in segs]
    segs = snap_ends(segs, (0.0, -r_axis), (0.0, r_axis))
    closing = ('line', (0.0, r_axis), (0.0, -r_axis))
    return segs, closing


class GearSet:
    """Every number the build needs, in mm."""

    def __init__(self, module, teeth, steps=DEFAULT_STEPS, print_clearance=DEFAULT_PRINT_CLEARANCE_MM,
                 pressure_angle=PRESSURE_ANGLE_DEG,
                 pole_feature=POLE_FEATURE, root_fillet=ROOT_FILLET,
                 cutter_extra=CUTTER_EXTRA_ADDENDUM, mp_shift=MP_PROFILE_SHIFT_MM,
                 face_width=MP_FACE_WIDTH, blank=MP_BLANK, bore=MP_BORE_DIAMETER_MM,
                 spacing_deg=MODULE_SPACING_DEG):
        m = float(module)
        z = int(teeth)
        if m <= 0:
            raise ValueError('Module must be positive.')
        if z % 2 != 0:
            raise ValueError('Tooth count must be even (paper constraint 1). Got %d.' % z)
        if z < 12:
            raise ValueError('Use at least 12 teeth on the ball (6 on each MP gear).')
        if pole_feature not in ('tooth', 'valley'):
            raise ValueError("POLE_FEATURE must be 'tooth' or 'valley'.")
        if blank not in ('throat', 'cylinder'):
            raise ValueError("MP_BLANK must be 'throat' or 'cylinder'.")
        if print_clearance < 0 or print_clearance > 0.25 * m:
            raise ValueError('Print clearance must be between 0 and %.2f mm (0.25 x module).' % (0.25 * m))
        self.m, self.z, self.steps = m, z, int(steps)
        self.clearance = float(print_clearance)
        self.z_mpl = z // 2
        self.pole = pole_feature
        self.alpha_deg = pressure_angle
        self.notes = []

        # CS gear (eq. 1-3) and its SH cutter. With a print clearance the ball's
        # flanks are pulled back by c and the cutter's are pushed out by c, which
        # pulls the generated MP flanks back by c as well.
        c = self.clearance
        self._cs_args = (m, z, pressure_angle, 1.0, 1.25, root_fillet)
        self.cs = Profile2D(m, z, pressure_angle, 1.0, 1.25, root_fillet, flank_offset=-c)
        # The cutter's tip lands are splines, not arcs (see the SH_CUTTER note
        # above cut_mp for why).
        self.cutter = Profile2D(m, z, pressure_angle, 1.0 + cutter_extra, 1.25, root_fillet, flank_offset=c,
                                spline_tips=True)
        self.notes += ['CS gear: ' + n for n in self.cs.notes]
        self.notes += ['SH cutter: ' + n for n in self.cutter.notes]

        self.d_sph = m * z
        self.r_sph = self.d_sph / 2.0
        self.d_mpl = m * self.z_mpl
        self.r_mpl = self.d_mpl / 2.0
        self.R = m * (z + self.z_mpl) / 2.0          # eq. 6, = 0.75 d_sph
        self.shift = 0.0 if self.clearance > 0 else float(mp_shift)
        self.R_cut = self.R + self.shift             # carrier length while cutting

        # MP blank
        self.blank = blank
        self.W = face_width * self.d_mpl
        self.r_tip_mid = self.r_mpl + m + self.shift
        self.rho_t = self.R - self.r_tip_mid         # tip clearance sphere radius about the ball center
        if blank == 'throat':
            if self.W / 2.0 >= 0.95 * self.rho_t:
                raise ValueError('MP face width is too large for a throat blank.')
            gap = math.sqrt(self.rho_t ** 2 - (self.W / 2.0) ** 2)
            self.r_tip_face = self.R - gap
            self.r_blank = self.r_tip_face + min(m, 0.5 * gap)   # oversize; the torus trims it
        else:
            self.r_tip_face = self.r_tip_mid
            self.r_blank = self.r_tip_mid
        self.bore = float(bore)
        if self.bore > 0 and self.bore / 2.0 >= self.r_mpl - 1.25 * m - 2 * m:
            raise ValueError('MP bore is too large for the MP root.')

        # Placement
        self.spacing = math.radians(spacing_deg)
        self.u = {'A': (1.0, 0.0, 0.0),
                  'B': (math.cos(self.spacing), math.sin(self.spacing), 0.0)}
        self.C = {k: (self.R * v[0], self.R * v[1], 0.0) for k, v in self.u.items()}
        center_gap = 2.0 * self.R * math.sin(self.spacing / 2.0)
        if center_gap < 2.0 * self.r_blank + 0.5 * m:
            raise ValueError('MP gears would hit each other at %.1f deg spacing. '
                             'Increase MODULE_SPACING_DEG.' % spacing_deg)

        # features on each revolve axis of the ball
        self.feature = {'A': self.pole,
                        'B': self.pole if z % 4 == 0 else other_feature(self.pole)}
        self.identical_mp = (z % 4 == 0)
        self.reuse_A_for_B = self.identical_mp and abs(spacing_deg - 90.0) < 1e-9

    def cs_profile(self, nudge=0.0):
        """The ball profile with the flanks that face the revolve axis moved
        out by `nudge` mm (see make_cs_body)."""
        if nudge == 0.0:
            return self.cs
        return Profile2D(*self._cs_args, flank_offset=-self.clearance, axis_side_extra=nudge)

    def tool_pose(self, family, psi):
        """Rotation angle about Z and translation (mm) of the SH cutter at
        carrier angle psi, for cutting MP gear `family` in place.
        Cutter body is modeled centered on the origin in the ball's home pose."""
        ux, uy, _ = self.u[family]
        cx, cy, _ = self.C[family]
        # cutter center at psi = 0 sits R_cut from the MP axis, toward the ball
        ox, oy = cx - self.R_cut * ux, cy - self.R_cut * uy
        vx, vy = ox - cx, oy - cy
        cp, sp = math.cos(psi), math.sin(psi)
        qx = cx + cp * vx - sp * vy
        qy = cy + sp * vx + cp * vy
        spin = psi * (1.0 + float(self.z_mpl) / self.z)      # psi_cr + psi_hb (eq. 7)
        return spin, (qx, qy, 0.0)

    def blank_volume(self):
        """Exact MP blank volume (mm^3), used to check the Fusion blank."""
        h = self.W / 2.0
        rb2 = (self.bore / 2.0) ** 2 if self.bore > 0 else 0.0
        if self.blank == 'cylinder':
            return math.pi * (self.r_tip_mid ** 2 - rb2) * self.W
        n = 400
        dz = 2.0 * h / n
        tot = 0.0
        for i in range(n):
            zz = -h + (i + 0.5) * dz
            r = self.R - math.sqrt(self.rho_t ** 2 - zz ** 2)
            tot += math.pi * (r * r - rb2) * dz
        return tot

    def summary(self):
        lines = [
            'Module m = %.3f mm, pressure angle %.1f deg' % (self.m, self.alpha_deg),
            'CS gear: %d teeth, pitch dia %.2f, tip dia %.2f, root dia %.2f mm'
            % (self.z, self.d_sph, 2 * self.cs.ra, 2 * self.cs.rf),
            'MP gears: %d teeth, pitch dia %.2f mm, face width %.2f mm (%s blank)'
            % (self.z_mpl, self.d_mpl, self.W, self.blank),
            'Center distance R_cr = %.3f mm' % self.R,
            'Ball X axis: %s, Y axis: %s' % (self.feature['A'], self.feature['B']),
        ]
        if not self.identical_mp:
            lines.append('z is not a multiple of 4, so MP gear A and MP gear B are different parts.')
        a = math.radians(self.alpha_deg)
        if self.clearance > 0:
            # both flanks of both gears pulled back by c: total backlash 4c along the flank normal
            j = 4.0 * self.clearance / math.cos(a)
            lines.append('Print clearance %.3f mm per gear per flank: %.3f mm gap on each side of a '
                         'tooth when centered' % (self.clearance, 2.0 * self.clearance))
        else:
            j = 2.0 * abs(self.shift) * math.tan(a)
            lines.append('Exact CAD gears, MP profile shift %.3f mm' % self.shift)
        lines.append('Total backlash at the pitch circles %.3f mm: %.2f deg of MP rotation, %.2f deg at the ball'
                     % (j, math.degrees(j / self.r_mpl), math.degrees(j / self.r_sph)))
        lines += self.notes
        return '\n'.join(lines)


# ------------------------------------------------------------------
# Fusion side
# ------------------------------------------------------------------

def P(x, y, z=0.0):
    """mm -> Fusion Point3D in cm"""
    return adsk.core.Point3D.create(x / 10.0, y / 10.0, z / 10.0)


def _dist2(a, b):
    return (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2


def draw_outline(sketch, segs, closing):
    """Draws the segments end to end, sharing sketch points, then the closing line."""
    curves = sketch.sketchCurves
    sketch.isComputeDeferred = True
    first_sp = None
    last_sp = None
    tol2 = 1e-10   # mm^2
    for seg in segs:
        kind = seg[0]
        if kind == 'spline':
            pts = seg[1]
            coll = adsk.core.ObjectCollection.create()
            coll.add(last_sp if last_sp is not None else P(*pts[0]))
            for q in pts[1:]:
                coll.add(P(*q))
            crv = curves.sketchFittedSplines.add(coll)
            if first_sp is None:
                first_sp = crv.startSketchPoint
            last_sp = crv.endSketchPoint
        elif kind == 'line':
            p0, p1 = seg[1], seg[2]
            if _dist2(p0, p1) < tol2:
                continue
            crv = curves.sketchLines.addByTwoPoints(last_sp if last_sp is not None else P(*p0), P(*p1))
            if first_sp is None:
                first_sp = crv.startSketchPoint
            last_sp = crv.endSketchPoint
        else:
            p0, pm, p1 = seg[1], seg[2], seg[3]
            if _dist2(p0, p1) < tol2:
                continue
            crv = curves.sketchArcs.addByThreePoints(last_sp if last_sp is not None else P(*p0), P(*pm), P(*p1))
            s, e = crv.startSketchPoint, crv.endSketchPoint
            target = P(*p1)
            if s.geometry.distanceTo(target) < e.geometry.distanceTo(target):
                near, far = s, e
            else:
                near, far = e, s
            if first_sp is None:
                first_sp = far
            last_sp = near
    curves.sketchLines.addByTwoPoints(last_sp, first_sp)
    sketch.isComputeDeferred = False


def revolve_sketch(comp, sketch, axis, name):
    profs = adsk.core.ObjectCollection.create()
    for i in range(sketch.profiles.count):
        profs.add(sketch.profiles.item(i))
    if profs.count == 0:
        raise RuntimeError('Sketch "%s" has no closed profile.' % sketch.name)
    revs = comp.features.revolveFeatures
    ri = revs.createInput(profs, axis, adsk.fusion.FeatureOperations.NewBodyFeatureOperation)
    ri.setAngleExtent(False, adsk.core.ValueInput.createByReal(2.0 * math.pi))
    rev = revs.add(ri)
    body = rev.bodies.item(0)
    if rev.bodies.count > 1:
        tools = adsk.core.ObjectCollection.create()
        for i in range(1, rev.bodies.count):
            tools.add(rev.bodies.item(i))
        ci = comp.features.combineFeatures.createInput(body, tools)
        ci.operation = adsk.fusion.FeatureOperations.JoinFeatureOperation
        body = comp.features.combineFeatures.add(ci).bodies.item(0)
    body.name = name
    return body


def build_family_body(comp, prof, family, pole, name):
    segs, closing = family_outline(prof, family, pole)
    sk = comp.sketches.add(comp.xYConstructionPlane)
    sk.name = name + ' profile'
    draw_outline(sk, segs, closing)
    axis = comp.xConstructionAxis if family == 'A' else comp.yConstructionAxis
    return revolve_sketch(comp, sk, axis, name)


def new_component(root, name):
    occ = root.occurrences.addNewComponent(adsk.core.Matrix3D.create())
    occ.component.name = name
    return occ


def z_pose_matrix(angle, t_mm):
    """Rotation about Z by angle, then translation t (mm)."""
    c, s = math.cos(angle), math.sin(angle)
    mx = adsk.core.Matrix3D.create()
    mx.setWithCoordinateSystem(P(*t_mm),
                               adsk.core.Vector3D.create(c, s, 0),
                               adsk.core.Vector3D.create(-s, c, 0),
                               adsk.core.Vector3D.create(0, 0, 1))
    return mx


def make_mp_blank(tbm, gs, family):
    cx, cy, _ = gs.C[family]
    h = gs.W / 2.0
    blank = tbm.createCylinderOrCone(P(cx, cy, -h), gs.r_blank / 10.0, P(cx, cy, h), gs.r_blank / 10.0)
    if gs.blank == 'throat':
        # Fusion's temporary torus can come back inside-out depending on the
        # axis direction (seen with +Z in Fusion 2705). Keep the one whose tube
        # contains the ball center and not the MP axis.
        inside = adsk.fusion.PointContainment.PointInsidePointContainment
        outside = adsk.fusion.PointContainment.PointOutsidePointContainment
        torus = None
        for sgn in (-1.0, 1.0):
            t = tbm.createTorus(P(cx, cy, 0), adsk.core.Vector3D.create(0, 0, sgn), gs.R / 10.0, gs.rho_t / 10.0)
            if t.pointContainment(P(0, 0, 0)) == inside and t.pointContainment(P(cx, cy, 0)) == outside:
                torus = t
                break
        if torus is None:
            raise RuntimeError('Could not build a valid torus for the MP throat. Set MP_BLANK = "cylinder".')
        if not tbm.booleanOperation(blank, torus, adsk.fusion.BooleanTypes.DifferenceBooleanType):
            raise RuntimeError('Could not trim the MP blank to the throat shape.')
    if gs.bore > 0:
        hole = tbm.createCylinderOrCone(P(cx, cy, -h - 1.0), gs.bore / 20.0, P(cx, cy, h + 1.0), gs.bore / 20.0)
        tbm.booleanOperation(blank, hole, adsk.fusion.BooleanTypes.DifferenceBooleanType)
    got = blank.volume * 1000.0 if blank.faces.count else 0.0
    want = gs.blank_volume()
    if abs(got - want) > 0.01 * want:
        raise RuntimeError('MP blank came out wrong: %.0f mm^3, expected %.0f mm^3.' % (got, want))
    return blank


# SH_CUTTER note: as the cutter steps around the blank, the tip land of each
# cutter position survives on the MP only as a thin strip between its
# neighbours (a few hundred strips per MP gear, many under 0.5 mm wide, some
# under 0.01 mm). With arc tips every strip is a piece of a sphere, and
# Fusion's mesher fills some thin sphere strips in as the whole sphere (R =
# cutter tip radius, 34.5 mm at m2 z32). Which strips break depends on the
# mesh quality, so the bubbles moved around between the screen, STL and 3MF
# exports. With spline tips the strips are NURBS faces, which Fusion meshes
# correctly at every quality.
def cut_mp(tbm, ui, gs, family, cutter_tmp, label):
    work = make_mp_blank(tbm, gs, family)
    n = gs.steps
    prog = ui.createProgressDialog() if ui else None
    if prog:
        prog.isCancelButtonShown = True
        prog.show('ABENICS gears', label + ': cut %v of %m', 0, n, 0)
    failed = []
    cancelled = False
    try:
        for i in range(n):
            psi = 2.0 * math.pi * i / n
            ok = False
            for nudge in (0.0, 1e-4, -1e-4):
                spin, t = gs.tool_pose(family, psi + nudge)
                tool = tbm.copy(cutter_tmp)
                tbm.transform(tool, z_pose_matrix(spin, t))
                if tbm.booleanOperation(work, tool, adsk.fusion.BooleanTypes.DifferenceBooleanType):
                    ok = True
                    break
            if not ok:
                failed.append(i)
            adsk.doEvents()
            if prog:
                prog.progressValue = i + 1
                if prog.wasCancelled:
                    cancelled = True
                    break
    finally:
        if prog:
            prog.hide()
    vol = work.volume * 1000.0 if work.faces.count else 0.0
    if not cancelled and not (0.3 * gs.blank_volume() < vol < gs.blank_volume()):
        raise RuntimeError('%s came out wrong (%.0f mm^3 from a %.0f mm^3 blank).'
                           % (label, vol, gs.blank_volume()))
    return work, failed, cancelled


def add_body(comp, tmp_body, name):
    bf = comp.features.baseFeatures.add()
    bf.startEdit()
    comp.bRepBodies.add(tmp_body, bf)
    bf.finishEdit()
    body = bf.bodies.item(0)
    try:
        bf.name = name
        body.name = name
    except Exception:
        pass
    return body


def _check_points(gs, n=500):
    """Fixed set of test points (mm) spread through the toothed shell of the
    ball, half of them packed close to the XY plane where the two tooth
    families touch."""
    pts = []
    golden = math.pi * (3.0 - math.sqrt(5.0))
    r0 = gs.cs.rf - 0.5 * gs.m
    r1 = gs.cs.ra + 0.2 * gs.m
    for i in range(n):
        frac = (i * 0.6180339887) % 1.0
        rad = r0 + (r1 - r0) * frac
        th = golden * i
        if i % 2:
            zz = gs.m * (2.0 * ((i * 0.7548776662) % 1.0) - 1.0)    # |z| < m
        else:
            zz = rad * (1.0 - 2.0 * (i + 0.5) / n)
        rr = math.sqrt(max(0.0, rad * rad - zz * zz))
        pts.append((rr * math.cos(th), rr * math.sin(th), zz))
    return pts


def _bad_intersection(res, a, b, pts):
    """Compares a boolean result with point tests on its two inputs.
    Returns (mismatches, points used)."""
    inside = adsk.fusion.PointContainment.PointInsidePointContainment
    outside = adsk.fusion.PointContainment.PointOutsidePointContainment
    bad = used = 0
    for q in pts:
        p = P(*q)
        ca, cb, cr = a.pointContainment(p), b.pointContainment(p), res.pointContainment(p)
        if any(c not in (inside, outside) for c in (ca, cb, cr)):
            continue
        used += 1
        if (cr == inside) != (ca == inside and cb == inside):
            bad += 1
    return bad, used


# Why the CS gear needs a nudge:
# family A and family B both contain the tooth outline in the XY plane, and both
# surfaces are tangent to Z there, so the two bodies touch tangentially along
# every flank of that outline. Fusion's solid modeler can grind on that for a
# minute or more and then fail with ASM_INCONS_REL (z = 32 with print
# clearance did). On each flank near the XY plane only one family forms the
# ball's surface; the other just touches it. Which one depends on which way the
# flank faces: the flanks that face a family's own revolve axis are the ones
# where that family only touches. Moving just those flanks out by CS_NUDGE_MM
# opens a real gap, so nothing touches, and the ball near the XY plane comes
# out unchanged. Away from the XY plane some of those flanks do form the
# surface, and there the ball's teeth end up at most CS_NUDGE_MM thicker.
# (The nudge alone was not enough for z = 32. The root fillet change in
# Profile2D, which removes a 0.03 mm wide ring face, was needed as well.)
CS_NUDGE_MM = (0.005, 0.02)      # tried in order


def make_cs_body(comp, gs, ui=None):
    """CS gear = family A intersected with family B. Returns (body, note)."""
    design = comp.parentDesign
    tbm = adsk.fusion.TemporaryBRepManager.get()
    inter = adsk.fusion.BooleanTypes.IntersectionBooleanType
    pts = _check_points(gs, 200)
    tried = []
    prog = ui.createProgressDialog() if ui else None
    if prog:
        prog.isCancelButtonShown = False
        prog.show('ABENICS gears', 'Intersecting the CS gear tooth families. This can take a minute '
                  'and Fusion will not respond meanwhile.', 0, 1, 0)
        adsk.doEvents()
    try:
        for nudge in CS_NUDGE_MM:
            mark = design.timeline.count
            prof = gs.cs_profile(nudge)
            try:
                body_a = build_family_body(comp, prof, 'A', gs.pole, 'Family A (about X)')
                body_b = build_family_body(comp, prof, 'B', gs.pole, 'Family B (about Y)')
            except Exception as e:
                tried.append('%g mm: tooth family failed (%s)' % (nudge, str(e).splitlines()[0]))
                design.timeline.markerPosition = mark
                design.timeline.deleteAllAfterMarker()
                continue
            adsk.doEvents()
            tmp_a, tmp_b = tbm.copy(body_a), tbm.copy(body_b)
            v_max = min(tmp_a.volume, tmp_b.volume) * 1000.0
            res = tbm.copy(tmp_a)
            t0 = time.time()
            try:
                ok = tbm.booleanOperation(res, tbm.copy(tmp_b), inter)
                why = 'failed'
            except Exception as e:
                ok, why = False, str(e).split(' - ')[0]
            t_bool = time.time() - t0
            vol = res.volume * 1000.0 if (ok and res.faces.count) else 0.0
            if ok and not (0.5 * v_max < vol < 1.0001 * v_max):
                ok, why = False, 'bad volume %.0f mm^3' % vol
            if ok:
                bad, used = _bad_intersection(res, tmp_a, tmp_b, pts)
                if bad or used < 0.8 * len(pts):
                    ok, why = False, '%d of %d test points wrong' % (bad, used)
            if ok:
                for body in (body_a, body_b):
                    comp.features.removeFeatures.add(body)
                cs_body = add_body(comp, res, 'CS gear')
                note = ('CS gear: intersection took %.0f s, built with a %g mm flank nudge '
                        '(teeth at most that much thicker in places)' % (t_bool, nudge))
                if tried:
                    note += '; earlier tries: %s' % '; '.join(tried)
                return cs_body, note
            tried.append('%g mm nudge: %s after %.0f s' % (nudge, why, t_bool))
            # didn't work: take these sketches and revolves back out of the timeline
            design.timeline.markerPosition = mark
            design.timeline.deleteAllAfterMarker()
    finally:
        if prog:
            prog.hide()
    raise RuntimeError('Could not intersect the two tooth families of the CS gear.\n' + '\n'.join(tried))


# ---- mesh check ------------------------------------------------------------
# STL and 3MF exports, the 3D Print command and the screen all mesh the B-rep
# with Fusion's mesher, each at its own quality. A face the mesher can't handle
# comes out as junk far bigger than the face (the bubbles). The check below
# meshes each finished body at all four qualities and compares the mesh with
# the body: junk sticks out of the body's box or inflates the surface area.

def _mesh_qualities():
    Q = adsk.fusion.TriangleMeshQualityOptions
    return (('Low', Q.LowQualityTriangleMesh), ('Normal', Q.NormalQualityTriangleMesh),
            ('High', Q.HighQualityTriangleMesh), ('Very high', Q.VeryHighQualityTriangleMesh))


def _mesh_vs_entity(entity, quality):
    """Meshes a body or face. Returns (how far the mesh sticks out of the
    entity's bounding box in mm, mesh area / true area)."""
    mc = entity.meshManager.createMeshCalculator()
    mc.setQuality(quality)
    mesh = mc.calculate()
    c = mesh.nodeCoordinatesAsDouble
    idx = mesh.nodeIndices
    if not idx:
        return 0.0, 0.0
    xs, ys, zs = c[0::3], c[1::3], c[2::3]
    bb = entity.boundingBox
    out = max(bb.minPoint.x - min(xs), bb.minPoint.y - min(ys), bb.minPoint.z - min(zs),
              max(xs) - bb.maxPoint.x, max(ys) - bb.maxPoint.y, max(zs) - bb.maxPoint.z, 0.0) * 10.0
    area = 0.0
    for i in range(0, len(idx), 3):
        a, b, d = idx[i] * 3, idx[i + 1] * 3, idx[i + 2] * 3
        ux, uy, uz = c[b] - c[a], c[b + 1] - c[a + 1], c[b + 2] - c[a + 2]
        vx, vy, vz = c[d] - c[a], c[d + 1] - c[a + 1], c[d + 2] - c[a + 2]
        area += 0.5 * math.sqrt((uy * vz - uz * vy) ** 2 + (uz * vx - ux * vz) ** 2 + (ux * vy - uy * vx) ** 2)
    return out, area / entity.area


def _junk(out_mm, ratio):
    # junk spheres stick out by tens of mm; a clean coarse mesh can run a few
    # percent over the true area on thin strips
    return out_mm > 0.05 or ratio > 1.1


def mesh_report(body):
    """{quality name: (stick-out mm, area ratio)} for the qualities that make junk."""
    bad = {}
    for nm, q in _mesh_qualities():
        out, ratio = _mesh_vs_entity(body, q)
        if _junk(out, ratio):
            bad[nm] = (out, ratio)
    return bad


def check_mesh(comp, name):
    """Checks a finished body at every mesh quality. If some quality makes junk,
    finds the faces responsible and removes them with a Delete Face feature
    (Fusion extends the neighbouring faces to close the gap), then checks again.
    With the spline-tip cutter this should never have anything to do."""
    adsk.doEvents()
    body = comp.bRepBodies.itemByName(name)
    bad = mesh_report(body)
    if not bad:
        return '%s: meshes clean at Low, Normal, High and Very high' % name
    quals = [(nm, q) for nm, q in _mesh_qualities() if nm in bad]
    faces = adsk.core.ObjectCollection.create()
    for f in body.faces:
        for nm, q in quals:
            if _junk(*_mesh_vs_entity(f, q)):
                faces.add(f)
                break
    first = ', '.join('%s (%.1f mm out, area x%.2f)' % (k, v[0], v[1]) for k, v in bad.items())
    if faces.count == 0:
        return 'WARNING %s: mesh junk at %s, but no single face to blame. Check the export.' % (name, first)
    try:
        feat = comp.features.deleteFaceFeatures.add(faces)
        feat.name = 'Remove unmeshable faces (%s)' % name
    except Exception:
        return 'WARNING %s: mesh junk at %s from %d faces that could not be removed.' % (name, first, faces.count)
    adsk.doEvents()
    left = mesh_report(comp.bRepBodies.itemByName(name))
    return '%s: mesh junk at %s, removed %d faces%s' % (
        name, first, faces.count, ', now clean' if not left else ', STILL BAD at %s' % ', '.join(left))


def ask(ui, prompt, title, default):
    res = ui.inputBox(prompt, title, default)
    if isinstance(res, (list, tuple)):
        val, cancelled = res[0], res[1]
    else:
        val, cancelled = res, False
    if cancelled:
        return None
    return val.strip()


def get_inputs(ui):
    title = 'ABENICS gears'
    while True:
        m = ask(ui, 'Module (mm):', title, str(DEFAULT_MODULE_MM))
        if m is None:
            return None
        try:
            m = float(m)
            if m <= 0:
                raise ValueError
            break
        except ValueError:
            ui.messageBox('Module must be a positive number.', title)
    while True:
        z = ask(ui, 'Teeth on the CS gear (even number, 12 or more).\n'
                    'Each MP gear gets half.', title, str(DEFAULT_TEETH))
        if z is None:
            return None
        try:
            z = int(z)
        except ValueError:
            ui.messageBox('Tooth count must be a whole number.', title)
            continue
        if z % 2 != 0:
            ui.messageBox('%d is odd. The paper requires an even tooth count (constraint 1).' % z, title)
            continue
        if z < 12:
            ui.messageBox('Use at least 12 teeth.', title)
            continue
        break
    while True:
        c = ask(ui, 'Print clearance per gear in mm.\n'
                    'Each tooth flank on the ball and the MP gears is pulled back by this much.\n'
                    'Common FFF guidance is 0.1 to 0.2. Enter 0 for exact CAD gears.',
                title, str(DEFAULT_PRINT_CLEARANCE_MM))
        if c is None:
            return None
        try:
            c = float(c)
        except ValueError:
            ui.messageBox('Clearance must be a number.', title)
            continue
        if c < 0 or c > 0.25 * m:
            ui.messageBox('Use a clearance between 0 and %.2f mm (0.25 x module).' % (0.25 * m), title)
            continue
        break
    while True:
        n = ask(ui, 'Cutter positions per MP gear (more = smoother flanks, slower).\n'
                    'Try 24 for a quick check.', title, str(DEFAULT_STEPS))
        if n is None:
            return None
        try:
            n = int(n)
            if n < 4:
                raise ValueError
            break
        except ValueError:
            ui.messageBox('Use a whole number of 4 or more.', title)
    return m, z, n, c


def run(context):
    ui = None
    try:
        app = adsk.core.Application.get()
        ui = app.userInterface
        inputs = get_inputs(ui)
        if inputs is None:
            return
        text = build(app, ui, *inputs)
        ui.messageBox(text, 'ABENICS gears')
    except Exception:
        if ui:
            ui.messageBox('Failed:\n%s' % traceback.format_exc(), 'ABENICS gears')


def build(app, ui, m, z, steps, clearance=DEFAULT_PRINT_CLEARANCE_MM):
    """Creates a new design with the CS gear, both SH cutters (hidden) and
    both MP gears. Returns a summary string."""
    gs = GearSet(m, z, steps, print_clearance=clearance)
    t0 = time.time()

    app.documents.add(adsk.core.DocumentTypes.FusionDesignDocumentType)
    design = adsk.fusion.Design.cast(app.activeProduct)
    design.designType = adsk.fusion.DesignTypes.ParametricDesignType
    root = design.rootComponent
    tag = 'm%g z%d' % (gs.m, gs.z) + (' print c%g' % gs.clearance if gs.clearance > 0 else '')

    # ---- CS gear: family A (about X) intersected with family B (about Y)
    cs_occ = new_component(root, 'CS gear %s' % tag)
    cs = cs_occ.component
    cs_body, cs_note = make_cs_body(cs, gs, ui)
    report = [cs_note]

    # ---- SH cutters: ball profile with longer teeth, revolved about one axis only
    cut_occ = new_component(root, 'SH cutters (generation tools)')
    cut = cut_occ.component
    cutter_a = build_family_body(cut, gs.cutter, 'A', gs.pole, 'SH cutter A')
    cutter_b = None if gs.identical_mp else build_family_body(cut, gs.cutter, 'B', gs.pole, 'SH cutter B')
    cut_occ.isLightBulbOn = False
    adsk.doEvents()

    # ---- MP gears
    tbm = adsk.fusion.TemporaryBRepManager.get()
    mp_a_occ = new_component(root, 'MP gear A %s' % tag)
    t1 = time.time()
    tmp_a, failed_a, cancelled = cut_mp(tbm, ui, gs, 'A', tbm.copy(cutter_a), 'MP gear A')
    add_body(mp_a_occ.component, tmp_a, 'MP gear A')
    report.append('MP gear A: %d cuts in %.0f s, %d failed%s'
                  % (steps, time.time() - t1, len(failed_a), ' (cancelled)' if cancelled else ''))
    report.append(check_mesh(mp_a_occ.component, 'MP gear A'))

    mp_b_occ = new_component(root, 'MP gear B %s' % tag)
    if cancelled:
        pass
    elif gs.reuse_A_for_B:
        # copy the finished A
        tmp_b = tbm.copy(mp_a_occ.component.bRepBodies.itemByName('MP gear A'))
        tbm.transform(tmp_b, z_pose_matrix(math.pi / 2.0, (0.0, 0.0, 0.0)))
        add_body(mp_b_occ.component, tmp_b, 'MP gear B')
        report.append('MP gear B: same part as A, rotated 90 deg about the ball')
        report.append(check_mesh(mp_b_occ.component, 'MP gear B'))
    else:
        if gs.identical_mp:
            cut_tool_b = tbm.copy(cutter_a)
            tbm.transform(cut_tool_b, z_pose_matrix(math.pi / 2.0, (0.0, 0.0, 0.0)))
        else:
            cut_tool_b = tbm.copy(cutter_b)
        t1 = time.time()
        tmp_b, failed_b, cancelled = cut_mp(tbm, ui, gs, 'B', cut_tool_b, 'MP gear B')
        add_body(mp_b_occ.component, tmp_b, 'MP gear B')
        report.append('MP gear B: %d cuts in %.0f s, %d failed%s'
                      % (steps, time.time() - t1, len(failed_b), ' (cancelled)' if cancelled else ''))
        report.append(check_mesh(mp_b_occ.component, 'MP gear B'))
    report.append(check_mesh(cs, 'CS gear'))

    app.activeViewport.fit()
    return '%s\n\n%s\nTime: %.0f s' % (gs.summary(), '\n'.join(report), time.time() - t0)
