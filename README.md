
******Readme for FebRNA1 package******  by Tan-group at Wuhan University

FebRNA1 is a package for building RNA 3D structures with input their sequence and secondary structures 
based on coarse-grained fragment ensembles. The program of FebRNA1 is run in Python,
and numpy, biopython and scipy modules are required.


Please run FebRNA1 as follows:


## Run in the example dir 
python ./Run.py (or python3 ./Run.py)
(in the file directory depending on the installed Python version) .

## According to corresponding instructions from FebRNA1, please input :
- (a) sequence information, 
- (b) secondary structure in  dot-bracket form, 
- (c) number of structures required (n), and
- (d) whether all-atom construction is required accordingly.
 
## Wait for a while (usually within several minutes) to obtain the results.
- (a) The results are placed in the './RESULT'; 
- (b) './RESULT/CG_Result' contains all the predicted coarse-grained conformations;
- (c) './RESULT/Select_Result' contains a selection of TOP-n coarse-grained conformations;
- (d) './RESULT/AA_Result' contains the rebuilt all-atom structures of selected coarse-grained structures.

An example is:
```
python Run.py 
Sequence:GCGGCACCGUCCGCUCAAACAAACGG
Secondary Structure:((((..[[[.)))).........]]]
Seleted Num(0=all):5
All-atom rebuilding?(y/n):y
Finish in folder ./RESULT
Running time :37.020s
```

If you have any questions about FebRNA1, please contact us by the email: zjtan@whu.edu.cn .
