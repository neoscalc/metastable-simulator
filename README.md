![neoscalc_metastable-simulator_2022](https://user-images.githubusercontent.com/54409312/169957104-0dc632c0-912a-429e-9f1b-b00b207a4b86.png)


# metastable-simulator
Collection of programs for energy calculation of metastable assemblages using theriak

# Input file MetastabilitySimulatorIN.txt

An example of input file is provided bellow. Each entry is described in the table

| Short code      | Description                                                   |
| ----------------| --------------------------------------------------------------|
| Version         | Version of metastable simulator for compatibility check. Changing this value requires updating the format of the input file... |
| TheriakPath     | Absolute path to Theriak       |
| Database        | Name of the database for Theriak        |
| Metacalc        | Name of the database with the excluded solution for mode 1 (nucleation))        |
| Bulk            | Input for the original bulk rock composition in the Theriak format        |
| --------------- |                                                               |
| Mode            | Mode to run the program for this job        |
| --- Options (1) | Options for mode 1 (nucleation)                                                              |
| OverstepMin     | Define the mineral which is prevented to be stable    |
| FracMin         | Mineral kept metastable (default: NONE)      |
| FracMolFrac     | Molar fraction of the mineral kept metastable in the interval [0,1] (default: 1)      |
| --- Options (2) | Options for mode 2 (persistence)                                                              |
| Frac2Meta       | ON/OFF (default: ON) If activated, the program keep adding the newly formed mineral in the fractionated list to the metastable system      |
| EquiMin         | Mineral names (default: NONE) separated by space/tabulation. These minerals if stable will re-equilibrate at each stage!       |
| EquiMolFrac     | Molar fractions in the interval [0,1] (default: 1) separated by space/tabulation. The mole fraction of the mineral that re-equilibrate        |
| --------------- |                                                               |
| GenerateSeeds   | ON/OFF (default: OFF) if activated the program only generates seeds and stop         |
| SaveOutput      | ON/OFF (default: ON)       |
| Print           | ON/OFF (default: OFF) if activated print all details in the command window              |
| Pause           | ON/OFF (default: OFF) if activated the program pauses after each loop                     |
| --------------- |                                                                |
| >               | PT path + system behaviour (0=metastable; 1=equilibrated). Each row bellow this keyword defines a P–T step with the entries separated by space or tabulation: temperature (°C); pressure (bar); key for system behaviour |




```
Version:			1.7
TheriakPath:			/Users/pierrelanari/Geologie/Programs/TheriakDominoCompiled/theriak
Database:			JUN92.bs
Metacalc:			JUN92_ExclPl.bs
Bulk:				AL(0.30866)CA(0.02801)FE(0.098918)MG(0.071276)MN(0.0027789)NA(0.086481)SI(1.0716)TI(0.011115)K(0.090249)H(0.03744)O(?)   *  LB_3  
-----------------------
Mode:				2	[1]NUCLEATION | [2]PERSISTENCE
------ Options (1)
OverstepMin:			FSP2
FracMin:			GARNET
FracMolFrac:			0.2		
------ Options (2)
Frac2Meta:			YES
EquiMin:			BIOTITE
EquiMolFrac:			1
-----------------------
GenerateSeeds:			OFF
SaveOutput:			ON
Print:				OFF		
Pause:				OFF


> PT path + system behaviour (0=metastable; 1=equilibrated)
600	6000	1
610	6200	0
620	6400	0
630	6600	0
640	6800	0
650	7000	0
660	7200	0
670	7400	0
680	7600	0
```


