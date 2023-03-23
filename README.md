![neoscalc_metastable-simulator_2022](https://user-images.githubusercontent.com/54409312/169957104-0dc632c0-912a-429e-9f1b-b00b207a4b86.png)


# metastable-simulator
Collection of programs for energy calculation of metastable assemblages using theriak

# Input file MetastabilitySimulatorIN.txt

An example of input file is provided bellow. Each entry is described in the table

| Short code      | Description                                     |
| ----------------| ------------------------------------------------|
| Version         |         |
| TheriakPath     |         |
| Database        |         |
| Metacalc        |         |
| Bulk            |         |
| --------------- |                                                 |
| Mode            |         |
| --- Options (1) |                                                 |
| FracMin         |         |
| FracMolFrac     |         |
| --- Options (2) |                                                 |
| EquiMin         |         |
| EquiMolFrac     |         |
| --------------- |                                                 |
| GenerateSeeds   |         |
| SaveOutput      |         |
| Print           |         |
| Pause           |         |
| --------------- |                                                 |
| >               | PT path + system behaviour (0=metastable; 1=equilibrated) |




```
Version:			1.6
TheriakPath:		/Users/pierrelanari/Geologie/Programs/TheriakDominoCompiled/theriak
Database:			JUN92.bs
Metacalc:			JUN92_ExclPl.bs
Bulk:				AL(0.30866)CA(0.02801)FE(0.098918)MG(0.071276)MN(0.0027789)NA(0.086481)SI(1.0716)TI(0.011115)K(0.090249)H(0.03744)O(?)   *  LB_3  
-----------------------
Mode:				2		[1]NUCLEATION | [2]PERSISTENCE
------ Options (1)
FracMin:			GARNET
FracMolFrac:		0.2		
------ Options (2)
EquiMin:			BIOTITE
EquiMolFrac:		1
-----------------------
GenerateSeeds:		OFF
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


