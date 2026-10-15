# ABENICS gear generator for Fusion 360

A Fusion 360 script that builds your own ABENICS spherical gear set: the cross spherical (CS) gear, which is the ball, and both monopole (MP) gears, placed in mesh around it. The math follows the ABENICS paper ([Abe, Tadakuma and Tadakuma, IEEE T-RO 2021](https://doi.org/10.1109/TRO.2021.3070124)). The MP gears are cut the same way the paper does it: a cutter with the ball's tooth profile steps around the MP blank and gets subtracted at each step.

## Install and run

1. Download `ABENICS_Gears.zip` from [Releases](https://github.com/teobrown/ABENICS-Fusion360/releases) and unzip it. You can also use Code > Download ZIP; the script is in the `ABENICS_Gears` folder.
2. In Fusion, go to Utilities > Add-Ins > Scripts and Add-Ins.
3. On the Scripts tab, click + and select the `ABENICS_Gears` folder.
4. Select ABENICS_Gears and click Run.

It opens a new design with the CS gear, MP gear A and MP gear B. The CS gear step can take a minute and Fusion won't respond while it runs. When it finishes, a message box lists the gear sizes and the backlash.

## Parameters

The script asks for these when it runs:

| Prompt | Default | Notes |
|---|---|---|
| Module (mm) | 2 | Sets the gear size. Ball pitch diameter = module x teeth. |
| CS gear teeth | 32 | Even, 12 or more. Each MP gear gets half. A multiple of 4 makes both MP gears the same part. |
| Print clearance (mm) | 0.15 | Taken off every tooth flank so printed gears don't bind. 0 gives exact CAD gears. |
| Cutter positions | 120 | Cuts per MP gear. More is smoother and slower. 24 is a quick test. |

The defaults give a 64 mm pitch diameter ball and two 16 tooth MP gears.

## Changing the script

Everything else is in the SETTINGS block at the top of `ABENICS_Gears.py`: pressure angle, whether a tooth or a valley sits on the ball's axis, root fillet, MP face width, MP blank shape, bore diameter, and the angle between the two MP gears. Change a value, save, and run it again. To open the file from Fusion, select the script and click Edit.

The gear math (the `Profile2D` and `GearSet` classes) doesn't call Fusion, so you can import the file in plain Python to check sizes before building anything.

## License

MIT
