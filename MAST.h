/*
 *
 *    MAST.h
 *    MAST
 *
 *    Copyright (C) 2013, CSIRO
 *
 *    This program is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU General Public License
 *    as published by the Free Software Foundation; either version 2
 *    of the License, or (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 *
 */


#ifndef __MAST__MAST__
#define __MAST__MAST__

#include <cstring>
#include <math.h>
#include <vector>
#include <set>
#include <string>
#include <iostream>
#include <cstdlib>

using namespace std;


typedef vector<pair<int,int> > PairList;

class Node;
class Tree;

class tripleProfile {

public:
    
    long leafNum;
    int wordSize; // number of bits per word i.e. sizeof(unsigned long) * 8
    long profileWordNum; // i.e. ceiling ( length of profile / wordSize )
    unsigned long* profile; // the x-th bit is referring to i | j k, where x = ( i * leafNum + j )* leafNum + k, and i,j,k start from zero

	// union of all leaves in all MASTs
	int* unionLeaves;
	
	// intersection of all leaves in all MASTs
	int* interLeaves;
    
    // total number of MASTs
    int totalMAST;
    
    tripleProfile(long leafNum); // constructor
    ~tripleProfile(); // destructor
    
    // set ( i | j k ) to 1
    void setTriple(int i, int j, int k);

    // get ( i | j k )
    int getTriple(int i, int j, int k);
    
    // get the set of k's such that ( i | j k ) is 1
    void getTripleSet(int i, int j, vector<int>& k_set);
    
    // get the set of (i,j) such that ( i | j k ) is 1
    // void getPairSet(vector<pair<int,int> >& pair_set);
    
    // reset
    void reset();
    
    // show the triples
    void print();
	
    // intersect with the triples of the other profile
    void intersect(tripleProfile& p);

	//=============================================================================
	// For MAST
	//=============================================================================
    
    // compute MAST table (maximum agreement subtree)
    // return MAST value
    int computeMAST(vector<pair<int,int> >& pair_set, vector<int>& commonLeaves);
    
    // get all results of MAST
    void getAllResults(set<string>& results, int rootNode, vector<string>& leafNames);

    // get the number of MASTs
    int getMASTNum(int rootNode);
    
    // show the MAST table
    void showMAST(vector<pair<int,int> >& pair_set);
    
	// get union and intersection of leaves of all MASTs
	// Prerequisite: MAST table has been computed
	void getUnionInterMASTleaves();

private:

    unsigned long mask; // -1
    
    unsigned long one; // 1
    
    unsigned long* k_array; // the array of k's for the specific values of i and j (internal use only)
    
    int k_array_length; // the length of k array

    // get the array of k's for the specific values of i and j
    void getKArray(int i, int j);
    
    // MAST table ( size : leafNum * leafNum )
    int* MAST;
    
    int MASTopt;
	
    // get the MAST subtrees rooted at lca(a, b)
    void getSubTrees(int a, int b, set<string>* subTrees, int closedWithBracket, vector<string>& leafNames);

	// get the number of MAST subtrees rooted at lca(a, b)
	int getSubTreeNum(int a, int b);
    
	// get the union and the interaction of the leaves of MAST subtree rooted at lca (a,b)
	void getUnionInterLeaves(int a, int b, int* uniLeaves, int* intLeaves, int& totalMAST);

	
};


// to compute MAST and KAST of all the input trees
// if "outputMAST" is 1, then get all the MAST trees.
void obtainMASTKAST(vector<string>& intrees, vector<string>& MASTs, string& KAST, int& MAST_score, int& KAST_score, int& totalMAST,
    int outputMAST, vector<int>& MASTunionLeaves, vector<int>& MASTinterLeaves, vector<string>& leafNames);

// to obtain the MAST of all the input trees
// void obtainMASTKAST(vector<string>& intrees, vector<string>& MASTs, string& KAST, int& MAST_score, int& KAST_score)


#endif
