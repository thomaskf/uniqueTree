# Software: uniqueTree 

This program reports:
1. The set of unique trees from the input trees;
2. The Maximum Agreement SubTrees (MASTs) of the input trees; and
3. The maximum agreement subtrees of MASTs (KASTs).

## Installation

The software was written in C++, and it has been tested under linux and MacOS platform. You need
to have C++ compiler installed in the machine in order to compile the source codes. The compilation
steps are shown as follows:

```
tar -zxvf uniqueTree-1.0.tar.gz
cd uniqueTree-1.0
make
```

Then an executable file named *uniqueTree* will appear

## Usage

Syntax:
```
./uniqueTree [input trees] [options]
```
```
[input trees] : trees in newick format
```

Options: 
```
-m : output all the MAST trees
     (this option may significantly slow down the program
      and consume much memory when there are too many MASTs)
```

Output files:

1. [input trees w/o ext].summary.txt :
                the summary of the result

2. [input trees w/o ext].unique_binary.nwk :
                the unique binary trees in newick format

3. [input trees w/o ext].unique_binary_weights.nwk :
                the unique binary trees with frequencies

4. [input trees w/o ext].unique_binary_linenums.nwk :
                the unique binary trees with the line
                numbers of the same trees

5. [input trees w/o ext].MAST.nwk :
                the maximum agreement subtrees (MASTs)
                (exists if -m option is used)

6. [input trees w/o ext].KAST.nwk :
                the maximum agreement subtrees of MASTs

7. [input trees w/o ext].seqs.list.txt :
                for each sequence, showing the number
                of MASTs or KASTs containing the sequence

8. [input trees w/o ext].polytomic_trees.nwk :
                the polytomic trees which are discarded
