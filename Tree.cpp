//
//  Tree.cpp
//  unique_tree
//
//  Created by kfwong on 18/4/13.
//  Copyright (c) 2013 CSIRO. All rights reserved.
//

#include "Tree.h"
#include "MAST.h"

double getLastDouble(string& str, size_t endPos, size_t& startPos) {
    // get the last double before the position endPos
    // return the double and the starting position of the double
    
    size_t j = endPos;
    while (j>0 && (!isdigit(str[j])))
        j--;
    size_t i = j;
    while (i>0 && (isdigit(str[i-1]) || str[i-1]=='.' || str[i-1]=='-' || str[i-1]=='e'))
        i--;
    if (j>=i) {
        startPos = i;
        double lastDouble = atof(str.substr(i,j-i+1).c_str());
        return lastDouble;
    } else {
        cerr << "Error! Double does not exist! " << str.substr(0,endPos+1) << endl;
        exit(1);
    }
}

void removeSpaceInsideStr(string& str) {
    size_t new_size = 0;
    for (size_t i=0; i<str.size(); i++) {
        if (str[i]!=' ') {
            if (new_size < i)
                str[new_size] = str[i];
            new_size++;
        }
    }
    if (new_size < str.size())
        str.resize(new_size);
}

bool isEdgeLen(string& str, size_t endPos) {
    // check whether the last double is an edge length
    while (endPos>0) {
        // stop if the character is ':' ',' ')' or '('
        if (str[endPos]==':' || str[endPos]==',' || str[endPos]==')' || str[endPos]=='(')
            break;
        endPos--;
    }
    return (str[endPos]==':');
}

Node::Node() {
    // constructor
    node1 = NULL;
    node2 = NULL;
    node3 = NULL;
    leafID = -1;
}

void Node::setValues(int isLeaf, string leafName, Node* node1, Node* node2, Node* node3,
               char isValidDist1, char isValidDist2, char isValidDist3,
                     double dist1, double dist2, double dist3) {
    this->isLeaf = isLeaf;
    this->leafName = leafName;
    this->node1 = node1;
    this->node2 = node2;
    this->node3 = node3;
    this->isValidDist1 = isValidDist1;
    this->isValidDist2 = isValidDist2;
    this->isValidDist3 = isValidDist3;
    this->dist1 = dist1;
    this->dist2 = dist2;
    this->dist3 = dist3;
}

void Node::print(Node* parentNode, bool showEdgeLen) {
    if (isLeaf == 1) {
        // a leaf
        cout << leafName;
        if (showEdgeLen && isValidDist1==1)
            cout << ":" << dist1;
    } else if (parentNode != NULL) {
        if (node1 == parentNode) {
            cout << "(";node2->print(this, showEdgeLen);cout<<",";node3->print(this,showEdgeLen);cout << ")";
            if (showEdgeLen && isValidDist1==1)
            {cout << ":";cout<<dist1;}
        } else if (node2 == parentNode) {
            cout << "(";node1->print(this, showEdgeLen);cout<<",";node3->print(this,showEdgeLen);cout << ")";
            if (showEdgeLen && isValidDist2==1)
            {cout << ":";cout<<dist2;}
        } else if (node3 == parentNode) {
            cout << "(";node1->print(this, showEdgeLen);cout<<",";node2->print(this,showEdgeLen);cout << ")";
            if (showEdgeLen && isValidDist3==1)
            {cout << ":";cout<<dist3;}
        } else {
            cerr << "[Node::print()] Error! The parentNode does not match with node1, node2 or node3" << endl;
            exit(1);
        }
    } else {
        // starting node
        cout << "(";
        if (node1 != NULL) {
            node1->print(this, showEdgeLen);cout<<",";
        }
        if (node2 != NULL) {
            node2->print(this,showEdgeLen);cout<<",";
        }
        // node 3 cannot be NULL
        node3->print(this,showEdgeLen);cout << ")";
    }
}

string Node::topStr(Node* parentNode) {
    // return the topology string
    
    if (isLeaf == 1) {
        // a leaf
        return leafName;
    } else if (parentNode != NULL) {
        if (node1 == parentNode) {
            return "(" + node2->topStr(this) + "," + node3->topStr(this) + ")";
        } else if (node2 == parentNode) {
            return "(" + node1->topStr(this) + "," + node3->topStr(this) + ")";
        } else if (node3 == parentNode) {
            return "(" + node1->topStr(this) + "," + node2->topStr(this) + ")";
        } else {
            cerr << "[topStr] Error! The parentNode does not match with node1, node2 or node 3" << endl << flush;
            exit(1);
        }
    } else {
        // starting node
        string resultingStr = "(";
        if (node1 != NULL) {
            resultingStr.append(node1->topStr(this) + ",");
        }
        if (node2 != NULL) {
            resultingStr.append(node2->topStr(this) + ",");
        }
        // node 3 cannot be NULL
        return resultingStr + node3->topStr(this) + ")";
    }
}

string Node::orderedTopStr(Node* parentNode, bool isStartInter, string& smallestLeaf) {
    // assume the tree is a unrooted tree
    // return the ordered-topology in newick format
    // if want to output a rooted tree, then set outputRootedTree to 1
    if (isLeaf == 1) {
        if (parentNode == node1 || node1 == NULL) {
            // a leaf
            smallestLeaf = leafName;
            return leafName;
        } else {
            // first starting point
            return node1->orderedTopStr(this, 1, smallestLeaf);
        }
    } else {
        string leftStr, rightStr;
        string parentStr = "";
        string lSmallLeaf = ""; string rSmallLeaf=""; string pSmallLeaf="";
        if (node1 == parentNode) {
            leftStr = node2->orderedTopStr(this,0,lSmallLeaf);
            rightStr = node3->orderedTopStr(this,0,rSmallLeaf);
            if (isStartInter) parentStr = node1->orderedTopStr(this,0,pSmallLeaf);
        } else if (node2 == parentNode) {
            if (node1 != NULL)
                leftStr = node1->orderedTopStr(this,0,lSmallLeaf);
            rightStr = node3->orderedTopStr(this,0,rSmallLeaf);
            if (isStartInter) parentStr = node2->orderedTopStr(this,0,pSmallLeaf);
        } else if (node3 == parentNode) {
            leftStr = node1->orderedTopStr(this,0,lSmallLeaf); rightStr = node2->orderedTopStr(this,0,rSmallLeaf);
            if (isStartInter) parentStr = node3->orderedTopStr(this,0,pSmallLeaf);
        } else {
            cerr << "[orderedTopStr] Error! ParentNode is unknown!" << endl;
            exit(1);
        }
        if (lSmallLeaf=="") {
            return "("+parentStr+","+rightStr+")";
        } else {
            int cmp = lSmallLeaf.compare(rSmallLeaf);
            if (cmp == 0) {
                cerr << "[orderedTopStr] Error! The left smallest leaf equals to the right smallest leaf!" << endl;
                exit(1);
            } else if (cmp < 0) {
                smallestLeaf = lSmallLeaf;
                if (isStartInter)
                    return "("+parentStr+","+leftStr+","+rightStr+")";
                else
                    return "("+leftStr+","+rightStr+")";
            } else {
                smallestLeaf = rSmallLeaf;
                if (isStartInter)
                    return "("+parentStr+","+rightStr+","+leftStr+")";
                else
                    return "("+rightStr+","+leftStr+")";
            }
        }
    }
}

Node* Tree::parseNewickSubStr(string& str, size_t endPos, size_t& startPos, Node* parent,
                              double& distToParent, char& isValidDistToParent) {
    // input: 1. string; 2. end position of the string; 3. parent node (if available)
    // output: 1. the root of the subtree generated ending at str[endPos];
    //         2. the startPos of the subtree on the string
    //         3. the edge length between this node and the parent node if available
    //
    // assume that there is no space in the input string
    
    if (str.size() == 0)
        return NULL;
    
    // There are three types: starting node, internal node or leaf
    
    size_t pos;
    Node* node1 = NULL;
    Node* node2 = NULL;
    Node* node3 = NULL;
    double dist1=0.0;
    double dist2=0.0;
    double dist3=0.0;
    char isValidDist1=0;
    char isValidDist2=0;
    char isValidDist3=0;
    
    if (str[endPos]==';') {
        // TYPE I : starting node
		// if there are three subnodes, then it is a unrooted tree
        
        // first check whether the tree has weighting
        if (str[endPos-1]==']') {
            // the tree has weighting
            weighting = getLastDouble(str, endPos-2, pos); // get the weighting
            endPos = pos-2;
        } else {
            endPos--;
        }
		
        // get dist1 if available
        if (isEdgeLen(str, endPos)) {
            distToParent = getLastDouble(str, endPos, pos);
            isValidDistToParent = 1;
            endPos = pos-2;
        } else {
            distToParent = 0.0;
            isValidDistToParent = 0;
        }
		
        if (str[endPos]==')')
            endPos--;
        
        Node* startNode = new Node();
        interNodeList.push_back(startNode);

        node3 = parseNewickSubStr(str, endPos, pos, startNode, dist3, isValidDist3);
        if (node3 == NULL) {
            isProblematic = 1;
            return NULL;
        }
        
        if (pos>2) {
            // more than one node
            node2 = parseNewickSubStr(str, pos-2, pos, startNode, dist2, isValidDist2);
            if (node2 == NULL) {
                isProblematic = 1;
                return NULL;
            }

            if (pos>2) {
                // unrooted tree
				isRooted = 0;
                node1 = parseNewickSubStr(str, pos-2, pos, startNode, dist1, isValidDist1);
                if (node1 == NULL) {
                    isProblematic = 1;
                    return NULL;
                }
            } else {
                // rooted tree
            }
        } else {
            // only one node
            // reset the node 1 of the node to NULL
            leafList[0]->node1=NULL;
            // and remove the interal node, which is just inserted, from the interNodeList
            interNodeList.clear();
        }
        
        
        // for starting node, pos should be 1
        if (pos > 1) {
            // cerr << "Non-binary tree : " << str << endl << flush;
            isProblematic = 1;
            return NULL;
            // exit(1);
        }
        if (pos > 0)
            startPos = pos-1;
        else
            startPos = pos;

        // setting the vales for this starting node
        startNode->setValues(0,"",node1,node2,node3,isValidDist1,isValidDist2,isValidDist3,dist1,dist2,dist3);
        
        return startNode;
    } else {
        
        // get dist1 if available
        if (isEdgeLen(str, endPos)) {
            dist1 = getLastDouble(str, endPos, pos);
            isValidDist1 = 1;
            endPos = pos-2;
        } else {
            dist1 = 0.0;
            isValidDist1 = 0;
        }
        distToParent = dist1;
        isValidDistToParent = isValidDist1;
        
        // skip the bootstrap value if any
        pos = str.find_last_of("(,)", endPos);
        if (str[pos]==')')
            endPos = pos;
        
        if (str[endPos] == ')') {
            
            // TYPE II : internal node
            
            Node* interNode = new Node();
            interNodeList.push_back(interNode);
            
            node3 = parseNewickSubStr(str, endPos-1, pos, interNode, dist3, isValidDist3);
            if (node3 == NULL) {
                isProblematic = 1;
                return NULL;
            }

            node2 = parseNewickSubStr(str, pos-2, pos, interNode, dist2, isValidDist2);
            if (node2 == NULL) {
                isProblematic = 1;
                return NULL;
            }

            node1 = parent;
            
            startPos = pos-1;

            interNode->setValues(0,"",node1,node2,node3,isValidDist1,isValidDist2,isValidDist3,dist1,dist2,dist3);
            
            return interNode;
            
        } else {
            
            // TYPE III : leaf
            
            Node* leaf = new Node();
            leafList.push_back(leaf);
            
            // get the starting position of the leaf
            startPos = str.find_last_of("(,)", endPos);
            if (startPos == string::npos) {
                startPos = 0;
            } else {
                startPos = startPos+1;
            }
            
            // get the leaf name
            string leafName = "";
            if (endPos-startPos+1 > 0)
                leafName = str.substr(startPos, endPos-startPos+1);
            
            leaf->setValues(1,leafName,parent,NULL,NULL,isValidDist1,0,0,dist1,0.0,0.0);
            
            return leaf;
            
        }
    }
}

void Tree::buildTreeFrNewick(string& str) {
    // build the tree according to the input string with Newick format
    
    removeSpaceInsideStr(str);
    // add the ";" at the end if necessary
    if (str.length() > 0 && str[str.length()-1]!=';')
        str.append(";");
    // cout << str << endl;
    size_t startPos;
    double distToParent;
    char isValidDistToParent;
    clear();
    parseNewickSubStr(str, str.length()-1, startPos, NULL,
                      distToParent, isValidDistToParent);
    
    // cout << "Number of leaves : " << leafList.size() << endl << flush;
    // cout << "Number of internal nodes : " << interNodeList.size() << endl << flush;
    
}

void Tree::print(bool showEdgeLen) {
    // print the whole tree
    
    if (leafList.size() == 0) {
        cout << ";";
    } else if (interNodeList.size() > 0) {
        Node* startNode = interNodeList[0];
        startNode->print(NULL,showEdgeLen);cout<<";"<<endl;
    } else {
        // there should be only one leaf
        if (leafList.size()>1) {
            cerr << "[print] Error! No internal node but more than one leaf!" << endl;
            exit(1);
        }
        cout << leafList[0]->leafName;
        if (showEdgeLen && leafList[0]->isValidDist1)
            cout << ":" << leafList[0]->dist1;
        cout << ";" << endl;
    }
}

string Tree::topStr() {
    // return the topology string
    
    if (leafList.size() == 0) {
        return "";
    } else if (interNodeList.size() > 0) {
        Node* startNode = interNodeList[0];
        return startNode->topStr(NULL);
    } else {
        // there should be only one leaf
        if (leafList.size()!=1) {
            cerr << "[topStr] Error! No internal node but more than one leaf!" << endl;
            exit(1);
        }
        return leafList[0]->leafName;
    }
}

string Tree::orderedTopStr() {
    
    // return the ordered-topology in newick format
	if (leafList.size() == 0) {
		return ";";
	} else if (leafList.size() == 1) {
		return leafList[0]->leafName + ";";
	}
	
	string smallestLeaf;
	if (isRooted) {
		return interNodeList[0]->orderedTopStr(NULL, 0, smallestLeaf);
	} else {
		int smallestLeafID = 0;
		for (int i=1; i<(int)leafList.size(); i++) {
			if (leafList[smallestLeafID]->leafName > leafList[i]->leafName)
				smallestLeafID = i;
		}
		return leafList[smallestLeafID]->orderedTopStr(NULL, 0, smallestLeaf) + ";";
	}
}

void Tree::clear() {
    int i;
    for (i=0;i<(int)leafList.size();i++) {
        delete(leafList[i]);
    }
    for (i=0;i<(int)interNodeList.size();i++) {
        delete(interNodeList[i]);
    }
    leafList.clear();
    interNodeList.clear();
    // reset the parameters
    weighting = 1.0;
    isProblematic = 0;
	isRooted = 1;
}

Tree::Tree() {
    // constructor
    weighting = 1.0;
    isProblematic = 0;
}

Tree::~Tree() {
    // destructor
    clear();
}

// extract the subtree with specific set of leaves
// leaveSelection[i] = 1 if the leaf with ID i is selected; 0 otherwise
string Tree::getSubTree(vector<int>& leaveSelection) {

	if (leafList.size()==1 && leaveSelection[leafList[0]->leafID]==1)
		return leafList[0]->leafName + ";";
	else if (interNodeList.size() == 0)
		return ";";
	else {
		string s1,s2;
		interNodeList[0]->getSubTree(NULL, leaveSelection, s1, s2);
		if (isRooted) {
			if (s2 != "") {
				return "(" + s1 + "," + s2 + ");";
			} else {
				return s1 + ";";
			}
		} else {
			if (s2 != "") {
				if (s2[0]=='(') {
					return "(" + s1 + "," + s2.substr(1,s2.length()-2) + ");";
				} else if (s1[0]=='(') {
					return "(" + s1.substr(1,s1.length()-2) + "," + s2 + ");";
				} else {
					return "(" + s1 + "," + s2 + ");";
				}
			} else {
				return s1 + ";";
			}
		}
	}
}


// extract the subtree with specific set of leaves
// leaveSelection[i] = 1 if the leaf with ID i is selected; 0 otherwise
void Node::getSubTree(Node* parentNode, vector<int>& leaveSelection, string& s1, string& s2) {
	
	s1=""; s2="";
	if (isLeaf) {
		if (leaveSelection[leafID]) {
			s1 = leafName;
		}
		return;
	}
	
	vector<pair<string,int> > leafLabels;
	int i;
	string ss1,ss2;
	
	if (node1 != parentNode) {
		node1->getSubTree(this, leaveSelection, ss1, ss2);
		if (ss1 != "")
			leafLabels.push_back(pair<string,int>(ss1,0));
		if (ss2 != "")
			leafLabels.push_back(pair<string,int>(ss2,0));
	}
	if (node2 != parentNode) {
		node2->getSubTree(this, leaveSelection, ss1, ss2);
		if (ss1 != "")
			leafLabels.push_back(pair<string,int>(ss1,1));
		if (ss2 != "")
			leafLabels.push_back(pair<string,int>(ss2,1));
	}
	if (node3 != parentNode) {
		node3->getSubTree(this, leaveSelection, ss1, ss2);
		if (ss1 != "")
			leafLabels.push_back(pair<string,int>(ss1,2));
		if (ss2 != "")
			leafLabels.push_back(pair<string,int>(ss2,2));
	}
	
	if (leafLabels.size() <= 2) {
		if (leafLabels.size()>0)
			s1 = leafLabels[0].first;
		if (leafLabels.size()>1)
			s2 = leafLabels[1].first;
		return;
	}
	
	// grouping the labels in the same group
	vector<int> deleted;
	deleted.push_back(0);
	for (i=1; i<(int)leafLabels.size(); i++) {
		if (!deleted[i-1] && (leafLabels[i-1].second==leafLabels[i].second)) {
			// same group
			leafLabels[i-1].first = "(" + leafLabels[i-1].first + "," + leafLabels[i].first + ")";
			deleted.push_back(1);
		} else {
			deleted.push_back(0);
		}
	}
	int newSize = 0;
	for (i=0; i<(int)leafLabels.size(); i++) {
		if (!deleted[i]) {
			if (newSize < i)
				leafLabels[newSize] = leafLabels[i];
			newSize++;
		}
	}
	if (newSize < (int)leafLabels.size())
		leafLabels.resize(newSize);

	s1 = leafLabels[0].first;
	if (leafLabels.size() > 2)
		s2 = "(" + leafLabels[1].first + "," + leafLabels[2].first + ")";
	else
		s2 = leafLabels[1].first;
	
}



//================================================================================
// functions related to maximum agreement subtree
//================================================================================

// From this node, build triple profile and the pair list for MUST
// leafList : to store all the leaves from this node
void Node::buildTripleProfile(Node* parentNode, tripleProfile& profile, vector<int>* leafList, vector<pair<int,int> >* pairList, int compute_pairList, vector<int>& leavesConsidered) {

	//   A triple (a,b,c) is used to represent a rooted triple a|bc, where lca(a,b) is the root,
	//   lca(b,c) is a descendant of lca(a,b), and a > c > b
    
    // clear leafList;
    leafList->clear();
    
    if (isLeaf) {
        if (leavesConsidered[leafID])
            leafList->push_back(leafID);
        return;
    }
        
    Node* child1 = node1;
    Node* child2 = node2;
    int i,j,k;
    
    if (parentNode==node1)
        child1 = node3;
    else if (parentNode == node2)
        child2 = node3;
    
    vector<int>* child_leafList = new vector<int>;
    
    if (child1!=NULL) {
        child1->buildTripleProfile(this, profile, leafList, pairList, compute_pairList, leavesConsidered);
    }
    
    if (child2!=NULL) {
        child2->buildTripleProfile(this, profile, child_leafList, pairList, compute_pairList, leavesConsidered);
    }
    
    if ((int)leafList->size() > 0 && (int)child_leafList->size() > 0) {
        
        // update the profile
        int compare_ij, compare_ik, compare_jk;

        // every two leaves from child1's leafList + every leaf from child2's leafList
		if ((int)leafList->size() > 1) {
			for (i=0; i<(int)leafList->size()-1; i++) {
				for (j=i+1; j<(int)leafList->size(); j++) {
					compare_ij = (leafList->at(i) < leafList->at(j));
					for (k=0; k<(int)child_leafList->size(); k++) {
						compare_ik = (leafList->at(i) < child_leafList->at(k));
						compare_jk = (leafList->at(j) < child_leafList->at(k));
						if ((compare_ij && compare_jk) || (!(compare_ij || compare_jk))) {
							// k > j > i or k < j < i
							profile.setTriple(child_leafList->at(k), leafList->at(i), leafList->at(j));
						} else if (((!compare_ij) && compare_ik) || (compare_ij && (!compare_ik))) {
							// k > i > j or k < i < j
							profile.setTriple(child_leafList->at(k), leafList->at(j), leafList->at(i));
						}
					}
				}
			}
		}
		
        // every two leaves from child2's leafList + every leaf from child1's leafList
		if ((int)child_leafList->size() > 1) {
			for (i=0; i<(int)child_leafList->size()-1; i++) {
				for (j=i+1; j<(int)child_leafList->size(); j++) {
					compare_ij = (child_leafList->at(i) < child_leafList->at(j));
					for (k=0; k<(int)leafList->size(); k++) {
						compare_ik = (child_leafList->at(i) < leafList->at(k));
						compare_jk = (child_leafList->at(j) < leafList->at(k));
						if ((compare_ij && compare_jk) || (!(compare_ij || compare_jk))) {
							// k > j > i or k < j < i
							profile.setTriple(leafList->at(k), child_leafList->at(i), child_leafList->at(j));
						} else if (((!compare_ij) && compare_ik) || (compare_ij && (!compare_ik))) {
							// k > i > j or k < i < j
							profile.setTriple(leafList->at(k), child_leafList->at(j), child_leafList->at(i));
						}
					}
				}
			}
		}
        
        // update the pairList
		// not for both lists have only one leaf
        if (compute_pairList && ((int)leafList->size() > 1 || (int)child_leafList->size() > 1)) {
            for (i=0; i<(int)leafList->size(); i++) {
                for (j=0; j<(int)child_leafList->size(); j++) {
                    if (leafList->at(i) < child_leafList->at(j))
                        pairList->push_back(pair<int,int>(leafList->at(i),child_leafList->at(j)));
                    else
                        pairList->push_back(pair<int,int>(child_leafList->at(j),leafList->at(i)));
                }
            }
        }
    }

    for (i=0; i<(int)child_leafList->size(); i++)
        leafList->push_back(child_leafList->at(i));
    child_leafList->clear();

    delete child_leafList;
    
}


void Tree::getLeafOrder(map<string,int>& leafOrder, vector<string>& leafNames) {
    // get the order of the leaves
    int i;
    for (i=0; i<(int)leafList.size(); i++) {
        leafOrder.insert(pair<string,int>(leafList[i]->leafName, i));
        leafNames.push_back(leafList[i]->leafName);
    }
}

void Tree::assignLeafID(map<string,int>& leafOrder) {
    // assign the ID of the leaves
    int i;
    map<string,int>::iterator itr;
    for (i=0; i<(int)leafList.size(); i++) {
        itr = leafOrder.find(leafList[i]->leafName);
        if (itr != leafOrder.end()) {
            leafList[i]->leafID = itr->second;
			// cout << leafList[i]->leafName << " => " << leafList[i]->leafID << endl;
        }
    }
}

// build triple profile and the pair list for MUST
// remove the terminal edge of the leaf with ID = "leafID"
// and consider the intersection of the removed edge as a new root
// thus the triple profile is according to the tree without the leaf with ID = "leafID"
void Tree::buildTripleProfile(int leafID, tripleProfile& profile, vector<pair<int,int> >& pairList, int compute_pairList, vector<int>& leavesConsidered) {
    
    // reset the triple-profile
    profile.reset();
	
	// reset the pairList
	if (compute_pairList)
		pairList.clear();
    
    // find the corresponding leaf to start
    int i;
    Node* startLeaf = NULL;
    for (i=0; i<(int)leafList.size(); i++) {
        if (leafList[i]->leafID == leafID) {
            startLeaf = leafList[i];
            break;
        }
    }
    vector<int>* leafList = new vector<int>;
    if (startLeaf != NULL && startLeaf->node1 != NULL) {
        // go to its parent and start building the triple profile
        startLeaf->node1->buildTripleProfile(startLeaf, profile, leafList, &pairList, compute_pairList, leavesConsidered);
    }
    leafList->clear();
    delete leafList;
    
}
