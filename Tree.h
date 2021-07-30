//
//  Tree.h
//  unique_tree
//
//  Created by kfwong on 18/4/13.
//  Copyright (c) 2013 CSIRO. All rights reserved.
//

#ifndef __unique_tree__Tree__
#define __unique_tree__Tree__

// uncomment the following line to activate the debug message
#define DEBUG_MODE

#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <map>

using namespace std;

class tripleProfile;

class Node {
    // for unrooted binary tree
public:
    int isLeaf; // 1 if it is a leaf, otherwise 0
    string leafName; // empty if not a leaf
    Node* node1; // not NULL
    Node* node2; // NULL for a leaf
    Node* node3; // NULL for a leaf
    char isValidDist1; // 0 if distance to node 1 is not defined
    char isValidDist2; // 0 if distance to node 2 is not defined
    char isValidDist3; // 0 if distance to node 3 is not defined
    double dist1;
    double dist2;
    double dist3;
    int leafID; // ID of the leaf
    
    Node(); // constructor
    
    void setValues(int isLeaf, string leafName, Node* node1, Node* node2, Node* node3,
              char isValidDist1, char isValidDist2, char isValidDist3,
              double dist1, double dist2, double dist3);
    
    void print(Node* parentNode, bool showEdgeLen);

    // return the topology string
    string topStr(Node* parentNode);
    
    // return the ordered-topology in newick format
    string orderedTopStr(Node* parentNode, bool isStartInter, string& smallestLeaf);

	// extract the subtree with specific set of leaves
	// leaveSelection[i] = 1 if the leaf with ID i is selected; 0 otherwise
	void getSubTree(Node* parentNode, vector<int>& leaveSelection, string& s1, string& s2);
	
    //================================================================================
    // functions added for computation of maximum agreement subtrees
    //================================================================================

    // From this node, build triple profile and the pair list for MUST
    // leafList : to store all the leaves from this node
    void buildTripleProfile(Node* parentNode, tripleProfile& profile, vector<int>* leafList, vector<pair<int,int> >* pairList, int compute_pairList, vector<int>& leavesConsidered);
    
};

class Tree {
    // unrooted binary tree
public:
	int isRooted; // whether this is a rooted tree
    vector<Node*> leafList; // list of all leaves
    vector<Node*> interNodeList; // list of all internal nodes
    
    Node* parseNewickSubStr(string& str, size_t endPos, size_t& startPos, Node* parent,
                                  double& distToParent, char& isValidDistToParent);
        // input: 1. string; 2. end position of the string; 3. parent node (if available)
        // output: 1. the root of the subtree generated ending at str[endPos];
        //         2. the startPos of the subtree on the string
        //         3. the edge length between this node and the parent node if available
        //
        // assume that there is no space in the input string
    
// public:

    double weighting; // weight of the tree
    int isProblematic; //  a flag to indicate whether the tree is problematic

    Tree(); // constructor
    ~Tree(); // destructor
    // build the tree according to the input string with Newick format
    void buildTreeFrNewick(string& str);
    // print the whole tree
    void print(bool showEdgeLen);
    // return the topology string
    string topStr();
    // return the ordered-topology in newick format
    string orderedTopStr();
    // clear the tree
    void clear();
	// extract the subtree with specific set of leaves
	// leaveSelection[i] = 1 if the leaf with ID i is selected; 0 otherwise
	string getSubTree(vector<int>& leaveSelection);

    //================================================================================
    // functions added for computation of maximum agreement subtrees
    //================================================================================
    
    // build triple profile and the pair list for MUST
    // remove the terminal edge of the leaf with ID = "leafID"
    // and consider the intersection of the removed edge as a new root
    // thus the triple profile is according to the tree without the leaf with ID = "leafID"
    void buildTripleProfile(int leafID, tripleProfile& profile, vector<pair<int,int> >& pairList, int compute_pairList, vector<int>& leavesConsidered);

    // get the order of the leaves
    void getLeafOrder(map<string,int>& leafOrder, vector<string>& leafNames);
    // assign the ID of the leaves
    void assignLeafID(map<string,int>& leafOrder);


};

#endif /* defined(__unique_tree__Tree__) */
