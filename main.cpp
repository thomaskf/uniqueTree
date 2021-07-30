//
//  main.cpp
//  unique_tree
//
//  Created by kfwong on 18/4/13.
//  Copyright (c) 2013 CSIRO. All rights reserved.
//

#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <set>
#include "Tree.h"
#include "MAST.h"

#define VERSION "1.10"

using namespace std;

string removeExt(string str) {
    // remove the extension part of the file name
    unsigned int found = (unsigned int) str.find_last_of(".");
    if (found < str.length()) {
        if (found > 0)
            return str.substr(0,found);
        else
            return "";
    }
    return str;
}

// return true if c is space or non-printable character
static bool isSpaceOrNonprintableChar(char c) {
    return (c <= 32 || c >= 127);
}

static void trim(string& str) {
    int i;
    for (i=(int)str.length()-1; i>=0; i--) {
        if (!isSpaceOrNonprintableChar(str[i]))
            break;
    }
    if (i==(int)str.length()-1) {
        // do nothing
    } else if (i>=0) {
        // resize the string
        str.resize(i+1);
    } else {
        // empty the string
        str = "";
    }
}

class TreeInfo {
public:
	string topStr;
	double weight;
	string lineNums;
	TreeInfo(string topStr, double weight, string lineNums) {
		this->topStr = topStr;
		this->weight = weight;
		this->lineNums = lineNums;
	}
};

bool myfunction (pair<double,TreeInfo*> i, pair<double,TreeInfo*> j) { return (i.first > j.first); }

int main(int argc, char** argv) {

    // tree manipulation program
    
    // =================================
    // return the unique trees
    // =================================

    if (argc!=2 && argc!=3) {
        cerr << "=========================================================" << endl;
        cerr << "                      Unique Tree " << VERSION << endl;
        cerr << "=========================================================" << endl;

        cerr << "Usage: " << argv[0] << " [input trees] [options]" << endl << endl;
        cerr << "[input trees] : trees in newick format" << endl << endl << endl;
        
		cerr << "Options: " << endl;
		cerr << "-m : output all the MAST trees" << endl;
		cerr << "     (this option may significantly slow down the program" << endl;
		cerr << "      and consume much memory when there are too many MASTs)" << endl << endl;
		
        cerr << "Output files:" << endl << endl;
        cerr << "1. [input trees w/o ext].summary.txt :" << endl;
        cerr << "                the summary of the result" << endl;
        cerr << endl;
        cerr << "2. [input trees w/o ext].unique_binary.nwk :" << endl;
        cerr << "                the unique binary trees in newick format" << endl;
        cerr << endl;
        cerr << "3. [input trees w/o ext].unique_binary_weights.nwk :" << endl;
        cerr << "                the unique binary trees with frequencies" << endl;
        cerr << endl;
        cerr << "4. [input trees w/o ext].unique_binary_linenums.nwk :" << endl;
        cerr << "                the unique binary trees with the line" << endl;
        cerr << "                numbers of the same trees" << endl;
        cerr << endl;
        cerr << "5. [input trees w/o ext].MAST.nwk :" << endl;
        cerr << "                the maximum agreement subtrees (MASTs)" << endl;
		cerr << "                (exists if -m option is used)" << endl;
        cerr << endl;
        cerr << "6. [input trees w/o ext].KAST.nwk :" << endl;
        cerr << "                the maximum agreement subtrees of MASTs" << endl;
        cerr << endl;
        cerr << "7. [input trees w/o ext].seqs.list.txt :" << endl;
        cerr << "                for each sequence, showing the number" << endl;
        cerr << "                of MASTs or KASTs containing the sequence" << endl;
        cerr << endl;
        cerr << "8. [input trees w/o ext].polytomic_trees.nwk :" << endl;
        cerr << "                the polytomic trees which are discarded" << endl;
        cerr << endl;
        cerr << "Contact: Thomas Wong <Thomas.Wong@anu.edu.au>" << endl;
        cerr << "=========================================================" << endl;
        exit(1);
    }
    
    ifstream fin;
    ofstream fout0, fout1, fout2, fout3, fout4, fout5;
    ofstream* fout6 = NULL;
    string prefix = removeExt(argv[1]);
	
	int outputMAST = 0;
	if (argc >= 3 && strcmp(argv[2],"-m")==0)
		outputMAST = 1;

    fin.open(argv[1]);
    if (!fin.is_open()) {
        cerr << "Error in opening the file :" << argv[1] << endl;
        exit(1);

    }

    string outFile0 = prefix + ".summary.txt";
    fout0.open(outFile0.c_str());
    if (!fout0.is_open()) {
        cerr << "Error in opening the file :" << outFile0 << endl;
        exit(1);
    }

    string outFile1 = prefix + ".unique_binary.nwk";
    fout1.open(outFile1.c_str());
    if (!fout1.is_open()) {
        cerr << "Error in opening the file :" << outFile1 << endl;
        exit(1);
    }
    
    string outFile2 = prefix + ".unique_binary_weights.nwk";
    fout2.open(outFile2.c_str());
    if (!fout2.is_open()) {
        cerr << "Error in opening the file :" << outFile2 << endl;
        exit(1);
    }

    string outFile3 = prefix + ".MAST.nwk";
	if (outputMAST) {
		fout3.open(outFile3.c_str());
		if (!fout3.is_open()) {
			cerr << "Error in opening the file :" << outFile3 << endl;
			exit(1);
		}
	}

    string outFile4 = prefix + ".KAST.nwk";
    fout4.open(outFile4.c_str());
    if (!fout4.is_open()) {
        cerr << "Error in opening the file :" << outFile4 << endl;
        exit(1);
    }
	
	string outFile5 = prefix + ".seqs.list.txt";
	fout5.open(outFile5.c_str());
	if (!fout5.is_open()) {
		cerr << "Error in opening the file :" << outFile5 << endl;
		exit(1);
	}

    string outFile6 = prefix + ".polytomic_trees.nwk";

    cout << "=================" << endl;
    cout << " Unique Tree " << VERSION << endl;
    cout << "=================" << endl;

    //===============================================================================//
    // For computation of unique set of trees
    //===============================================================================//

	cout << "computing unique set of trees..." << endl;
    map<string,TreeInfo* > allTrees;
    map<string,TreeInfo* >::iterator itr;
    string aline;
    string aTreeline = "";
    Tree tree;
    string orderedTopology;
    int i;
    while (getline(fin,aline)) {
        trim(aline);
        if (aline.length() > 0) {
            aTreeline.append(aline);
            if (aTreeline.length() > 0 && aTreeline[aTreeline.length()-1]==';') {
                tree.buildTreeFrNewick(aTreeline);
                if (tree.isProblematic==0) {
                    orderedTopology = tree.orderedTopStr();
                    itr = allTrees.find(orderedTopology);
                    if (itr == allTrees.end()) {
                        // unique tree
                        allTrees.insert(pair<string,TreeInfo*>(orderedTopology,new TreeInfo(tree.topStr(), tree.weighting, "")));
                    } else {
                        // not unique
                        (itr->second)->weight += tree.weighting;
                    }
                } else {
                    // the tree is problematic, thus return it into the file "[input file].discard.trees"
                    if (fout6 == NULL) {
                        fout6 = new ofstream;
                        fout6->open(outFile5.c_str());
                    }
                    (*fout6) << aTreeline << endl;
                }
                aTreeline="";
            }
        }
    }
    fin.close();
	
    // sort all the unique trees according to the frequencies
    vector<pair<double,TreeInfo*> > freq_treeinfo;
    for (itr=allTrees.begin(); itr!=allTrees.end(); itr++) {
        freq_treeinfo.push_back(pair<double,TreeInfo*>((itr->second)->weight,itr->second));
    }
    sort(freq_treeinfo.begin(),freq_treeinfo.end(),myfunction);
    
    // print out the unique trees
    for (i=0; i<(int)freq_treeinfo.size(); i++) {
        fout1 << freq_treeinfo[i].second->topStr << ";" << endl;
    }
    fout1.close();
    
    // print out the unique trees with frequencies
    for (i=0; i<(int)freq_treeinfo.size(); i++) {
        fout2 << freq_treeinfo[i].second->topStr << "[" << freq_treeinfo[i].second->weight << "];" << endl;
    }
    fout2.close();

    if (fout6 != NULL) {
        fout6->close();
    }

	cout << "Number of unique trees : " << freq_treeinfo.size() << endl;
    fout0 << "Number of unique trees : " << freq_treeinfo.size() << endl;

    //===============================================================================//
    // For computation of MAST & KAST
    //===============================================================================//

    vector<string> uniqueTreeSet;
    for (i=0; i<(int)freq_treeinfo.size(); i++) {
        uniqueTreeSet.push_back(freq_treeinfo[i].second->topStr);
    }

    vector<string> MASTs;
    string KAST;
	int MAST_score;
	int KAST_score;
    int totalMAST;
	vector<int> MASTunionLeaves;
	vector<int> MASTinterLeaves;
	vector<string> leafNames;
	
	// to compute MAST and KAST of all the input trees
	// if "outputMAST" is 1, then output all the MAST trees.
	obtainMASTKAST(uniqueTreeSet, MASTs, KAST, MAST_score, KAST_score, totalMAST, outputMAST, MASTunionLeaves, MASTinterLeaves, leafNames);

    // show MAST results
    cout << endl;
    cout << "MAST result" << endl;
    cout << "-----------" << endl;
	cout << "# of maximum agreement subtree(s): " << totalMAST << endl;
    cout << "size of maximum agreement subtree: " << MAST_score << endl;

    fout0 << endl;
    fout0 << "MAST result" << endl;
    fout0 << "-----------" << endl;
	fout0 << "# of maximum agreement subtree(s): " << totalMAST << endl;
    fout0 << "size of maximum agreement subtree: " << MAST_score << endl;
    
	if (outputMAST) {
		for (i=0; i<(int)MASTs.size(); i++) {
			fout3 << MASTs[i] << endl;
		}
	}
        
    // show KAST results
    cout << endl;
    cout << "KAST result" << endl;
    cout << "-----------" << endl;
    if (KAST == "") {
        cout << "# of KAST: 0" << endl;
    } else {
        cout << "# of KAST(s): 1" << endl;
        cout << "size of KAST: " << KAST_score << endl;
        fout4 << KAST << endl;
    }

    fout0 << endl;
    fout0 << "KAST result" << endl;
    fout0 << "-----------" << endl;
    if (KAST == "") {
        fout0 << "# of KAST: 0" << endl;
    } else {
        fout0 << "size of KAST: " << KAST_score << endl;
        fout0 << "   # of KAST: 1" << endl;
    }
	
	// for each leaf, to show how many MAST/KAST containing the leaf
    // according to the descending order of the MAST #
    multimap<int,int> order;
    multimap<int,int>::reverse_iterator ritr;
	for (i=0; i<(int)leafNames.size(); i++) {
        order.insert(pair<int,int>(MASTunionLeaves[i],i));
	}
    fout5 << "# For each leaf, the number of MAST/KAST containing the leaf" << endl;
	fout5 << "All\tMAST\tKAST" << endl;
    for (ritr=order.rbegin(); ritr!=order.rend(); ritr++) {
        i = ritr->second;
        fout5 << leafNames[i] << "\t" << MASTunionLeaves[i] << "\t" << MASTinterLeaves[i] << endl;
    }
    
    /*
	// first show the leaves contained in KAST
	for (i=(int)leafNames.size()-1; i>=0; i--) {
		if (MASTinterLeaves[i]) {
			fout5 << leafNames[i] << "\t" << MASTunionLeaves[i] << "\t1" << endl;
		}
	}
	// then show the leaves contained in MAST but not KAST
	for (i=(int)leafNames.size()-1; i>=0; i--) {
		if ((MASTinterLeaves[i]==0) && (MASTunionLeaves[i] > 0)) {
			fout5 << leafNames[i] << "\t" << MASTunionLeaves[i] << "\t0" << endl;
		}
	}
	// then show the other leaves
	for (i=(int)leafNames.size()-1; i>=0; i--) {
		if ((MASTinterLeaves[i]==0) && (MASTunionLeaves[i]==0)) {
			fout5 << leafNames[i] << "\t0\t0" << endl;
		}
	}
    */

    
    //===============================================================================//
    // To show all the output files
    //===============================================================================//
    
    cout << endl;
    cout << "Output files" << endl;
    cout << "------------" << endl;

    int maxFileNameLen = (int) outFile0.length();
    if (maxFileNameLen < (int) outFile1.length())
        maxFileNameLen = (int) outFile1.length();
    if (maxFileNameLen < (int) outFile2.length())
        maxFileNameLen = (int) outFile2.length();
    if (outputMAST && maxFileNameLen < (int) outFile3.length())
        maxFileNameLen = (int) outFile3.length();
    if (maxFileNameLen < (int) outFile4.length())
        maxFileNameLen = (int) outFile4.length();
    if (maxFileNameLen < (int) outFile5.length())
        maxFileNameLen = (int) outFile5.length();
    if (fout6 != NULL && maxFileNameLen < (int)outFile6.length())
        maxFileNameLen = (int) outFile6.length();
        
    if ((int)outFile0.length() < maxFileNameLen)
        outFile0 = outFile0 + string(maxFileNameLen-outFile0.length(), ' ');
    cout << outFile0 << " : the summary of the result" << endl;
    
    if ((int)outFile1.length() < maxFileNameLen)
        outFile1 = outFile1 + string(maxFileNameLen-outFile1.length(), ' ');
    cout << outFile1 << " : the unique binary trees in newick format" << endl;

    if ((int)outFile2.length() < maxFileNameLen)
        outFile2 = outFile2 + string(maxFileNameLen-outFile2.length(), ' ');
    cout << outFile2 << " : the unique binary trees with frequencies" << endl;

	if (outputMAST) {
		if ((int)outFile3.length() < maxFileNameLen)
			outFile3 = outFile3 + string(maxFileNameLen-outFile3.length(), ' ');
		cout << outFile3 << " : the maximum agreement subtrees" << endl;
	}

    if ((int)outFile4.length() < maxFileNameLen)
        outFile4 = outFile4 + string(maxFileNameLen-outFile4.length(), ' ');
    cout << outFile4 << " : the KAST subtrees" << endl;

    if ((int)outFile5.length() < maxFileNameLen)
        outFile5 = outFile5 + string(maxFileNameLen-outFile5.length(), ' ');
    cout << outFile5 << " : for each sequence, show how many MAST/KAST containing the sequence" << endl;

    if (fout6 != NULL) {
        if ((int)outFile6.length() < maxFileNameLen)
            outFile6 = outFile6 + string(maxFileNameLen-outFile6.length(), ' ');
        cout << outFile6 << " : the discarded polytomic trees" << endl;
    }
}
