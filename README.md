# Biodiff
A bioinformatics program to compare two data files and output
### Program:
    This program can find differences between two files containing bioinformatics data. 

    If provided two files A (from-file) and B (to-file), program will generate all lines
    in A-B, A&B_A, A&B_B and B-A in terms of the criteria given. The file format of 
    file A and B can be different.

    There will be two styles for comparison: one is coordinate based (option –c ) and 
    the other is name based (option –n). The default style is option -c. 

    The two styles were described as follows. 

    1) Coordinate-based diff. Two or more columns from file A and B will be selected and
    compared to check if the two regions overlap. If two regions from the two files overlap,
    then these two regions will be put into to A&B_A and A&B_B; those regions in A but not 
    in A&B will be put into A-B; and those in B but not in A&B will be put into B-A. 

    2) Name-based diff. Two columns from file A and B will be selected and compared in terms
    of string comparison. Users need to specify the column numbers in two files to be compared. 
    If their names “overlap”, it should generate 4 result files corresponding to A&B_A, A&B_B, 
    A-B, and B-A, where A&B_A contains those lines from file A and overlapping with some entries
    in file B; A&B_B contains lines from file B and overlapping with entries in file A; A-B 
    contains those lines from file A and with no overlapping entries in B; and B-A stands for 
    those lines from file B but with no overlapping entries in A. 

### Usage:
    This program can be used by following code when execute:
```
        ./Biodiff [options] from-file to-file
```
    In which, 
        option             --- the option you choose. See also OPTION in Reference
        from-file, to-file --- the absolute or relative path to the two files
### Reference:
    OPTION:
        -c output according to coordinate-based diff
        -n output according to name-based diff
        -a specify the criteria of from-file(following by the criteria)
        -b specify the criteria of to-file(following by the criteria)
        -h get the help of this program
### Example:
```
   ./Biodiff -a 3,4 -b 3,4 geneA.gtf geneB.gtf
   ./Biodiff -c -a 3,4 -b 3,4 geneA.gtf geneB.gtf
   ./Biodiff -n -a 0 -b 8 geneA.gtf geneB.gtf
```
### Author:
    Zhang
