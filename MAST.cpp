/*
 *
 *    tripleProfile.cpp
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
 
 /*
  *   The calculation of MAST is based on the algorithm mentioned in the paper:
  *   Swenson K.M., Chen E., Pattengale N.D., Sankoff D. (2012) The Kernel of Maximum Agreement Subtrees, IEEE/ACM TCBB, 9(4).
  *
  *   We have made a small modification on the algorithm to enhance the speed of the program practically
  *
  *   Definitions:
  *   A set of rooted binary trees T = {T_1, T_2, ..., T_k},
  *   A set of labels L such that each x \in L labels exactly one leaf of each T_i.
  *   The set of labels L is sorted according to the leaf position in one of the trees, say T_1.
  *   T(a,b) where a <= b \in L = the set of all agreement subtrees where lca(a,b) is the root of the tree
  *          with leaves between a and b.
  *   M(a,b) which is a subset of T(a,b) = the set of maximum agreement subtrees where lca(a,b) is the root
  *          of the tree with leaves between a and b.
  *   MAST(a,b) = the number of leaves in any member of M(a,b).
  *
  *   a|bc or (a,b,c) = a rooted triple, where lca(a,b) is the root, lca(b,c) is a descendant of lca(a,b), a > c > b or a < c < b
  *
  *   R = the set of rooted triples common to all trees in T
  *   X_{a,b} = {x: a|bx \in R} union with {b} such that lca(a,b) is the root.
  *
  *   MAST_z(y) = max{MAST(y,x): x \in X_{z,y} and y < x, MAST(x,y) : x \in X_{z,y} and x < y} 
  *               to be the MAST of the leaves between z and y in a subtree on y's side of the root
  *   
  *   MAST(a,b) = MAST_b(a) + MAST_a(b)
  */
  

#include "MAST.h"
#include "Tree.h"

// integer to string
// minLen : minimum length of the integer.
// If the number of digits < minimum length, then the integer will be displayed with leading zeros, unless the integer is zero
string intToStr(int i, int minLen) {
    char* d2c = (char*) "0123456789"; // the array for digit to char
    if (i<0) {
        return "-" + intToStr(-i, minLen);
    } else if (i<10) {
        if (minLen > 0)
            return string(minLen-1, '0') + string(1,d2c[i]);
        else
            return string(1,d2c[i]);
    } else {
        return intToStr(i/10, minLen-1)+string(1,d2c[i%10]);
    }
}
 // constructor
tripleProfile::tripleProfile(long leafNum) {
    long numBytes;
    this->leafNum = leafNum;
    wordSize = sizeof(unsigned long) * 8;
    profileWordNum = ( leafNum * leafNum * leafNum + (long) wordSize - 1 )/ (long) wordSize;
    numBytes = profileWordNum * sizeof(unsigned long);
    profile = new unsigned long [profileWordNum];
    memset(profile, 0, numBytes);

    k_array_length = (int) (leafNum + wordSize - 1) / wordSize;
    k_array = new unsigned long[k_array_length];
    
    mask = (unsigned long) -1;
    one = (unsigned long) 1;

    // MAST table ( size : leafNum * leafNum )
    MAST = new int[leafNum*leafNum];
	unionLeaves = NULL;
	interLeaves = NULL;
}

// destructor
tripleProfile::~tripleProfile() {
    delete[] profile;
    delete[] k_array;
    delete[] MAST;
	if (unionLeaves != NULL)
		delete[] unionLeaves;
	if (interLeaves != NULL)
		delete[] interLeaves;
}

// reset
void tripleProfile::reset() {
    memset(profile, 0, profileWordNum * sizeof(unsigned long));
}


// set ( i | j k ) to 1
void tripleProfile::setTriple(int i, int j, int k) {
    
    long pos = ( i * leafNum + j ) * leafNum + k;
    long interPos = pos / wordSize;
    int intraPos = pos % wordSize;
    
    if (interPos >= profileWordNum) {
        cerr << "Error! interPos (" << interPos << ") >= profileWordNum (" << profileWordNum << ")" << endl;
        exit(1);
    }
    
    profile[interPos] = profile[interPos] | ( one << (wordSize - intraPos - 1) );
}

// get ( i | j k )
int tripleProfile::getTriple(int i, int j, int k) {
    long pos = ( i * leafNum + j ) * leafNum + k;
    long interPos = pos / wordSize;
    int intraPos = pos % wordSize;
    return (int) (profile[interPos] & ( one << (wordSize - intraPos - 1) ));
}


// get the array of k's for the specific values of i and j
void tripleProfile::getKArray(int i, int j) {

    // initialize the k_array
    memset(k_array, 0, k_array_length * sizeof(unsigned long));
    
    long pos1 = ( i * leafNum + j ) * leafNum ;
    long interPos1 = pos1 / wordSize;
    int intraPos1 = pos1 % wordSize;

    long pos2 = ( i * leafNum + j + 1 ) * leafNum - 1;
    long interPos2 = pos2 / wordSize;
    int intraPos2 = pos2 % wordSize;

    int trimSize;
    long k=0;
    long interPos;
    unsigned long t;
    
    for (interPos=interPos1; interPos<=interPos2; interPos++) {
        
        if (profile[interPos]) {
            t = profile[interPos];
            if (interPos==interPos2 && intraPos2 < wordSize-1) {
                // trim the end part
                trimSize = wordSize - intraPos2 - 1;
                t = t & (mask << trimSize);
            }
            if (t) {
                if (k > 0 && intraPos1 > 0)
                    k_array[k-1] |= (t >> (wordSize - intraPos1));
                if (k < k_array_length)
                    k_array[k] = t << intraPos1;
            }
        }
        k++;
        
    }
}


// get the set of k's such that ( i | j k ) is 1
void tripleProfile::getTripleSet(int i, int j, vector<int>& k_set) {
    
    // clear the set
    k_set.clear();
    
    getKArray(i, j);
    
    int k,p;
    int leadZeroNum;
    unsigned long t;
    for (p=0; p<k_array_length; p++) {
        t = k_array[p];
        k = p*wordSize;
        while (t) {
            leadZeroNum = __builtin_clzl(t);
            k += leadZeroNum;
            k_set.push_back(k);
            if (leadZeroNum + 1 >= wordSize)
                t = 0;
            else
                t = t << (leadZeroNum + 1);
            k++;
        }
    }
    
}

// show the triples
void tripleProfile::print() {
    cout << "List of items inside tripleProfile:" << endl;
    int i,j;
    unsigned long t;
    int leadZeroNum;
    int position;
    int first,second,third;
    for (i=0; i<profileWordNum; i++) {
        t = profile[i];
        j=0;
        while (t) {
            leadZeroNum = __builtin_clzl(t);
            j += leadZeroNum;
            position = i * wordSize + j;
            first = position / (leafNum * leafNum);
            second = (position / leafNum) % leafNum;
            third = position % leafNum;
            cout << "(" << first << "," << second << "," << third << ")" << endl;
            if (leadZeroNum + 1 >= wordSize)
                t = 0;
            else
                t = t << (leadZeroNum + 1);
            j++;
        }
    }
}

// compute MAST table (maximum agreement subtree)
int tripleProfile::computeMAST(vector<pair<int,int> >& pair_set, vector<int>& commonLeaves) {
    
    // initialize the whole MAST table
    if (leafNum == 1) {
        MASTopt = 1;
        return 1;
    }
        
    MASTopt = 2;
    
    int i, j;
    int a, b, c;
    vector<int> c_set;
    int MASTa, MASTb;
    
    // initialize MAST
    memset(MAST, 0, leafNum*leafNum*sizeof(int));
    
    
    for (i=0; i<leafNum; i++) {
        if (commonLeaves[i]==0)
            continue;
		/*
		// according to our definition, MAST(i,j) = 0 if j<i
        for (j=0; j<i; j++) {
            if (commonLeaves[j]==0)
                continue;
            MAST[i*leafNum + j] = 2;
        }
		*/
        MAST[i*leafNum + i] = 1;
        for (j=i+1; j<leafNum; j++) {
            if (commonLeaves[j]==0)
                continue;
            MAST[i*leafNum + j] = 2;
        }
    }
    
    for (i=0; i<(int)pair_set.size(); i++) {
        a = pair_set[i].first;
        b = pair_set[i].second;
        MASTa = 0;
        MASTb = 0;

        // compute MASTa
        // get the set of c's such that ( a c | b ) is 1
        getTripleSet(b, a, c_set);
        c_set.push_back(a);
        for (j=0; j<(int)c_set.size(); j++) {
            c = c_set[j];
            if (a <= c && MAST[a*leafNum + c] > MASTa)
                MASTa = MAST[a*leafNum + c];
            else if (a > c && MAST[c*leafNum + a] > MASTa)
                MASTa = MAST[c*leafNum + a];
        }

        // compute MASTb
        // get the set of c's such that ( a | b c ) is 1
        getTripleSet(a, b, c_set);
        c_set.push_back(b);
        for (j=0; j<(int)c_set.size(); j++) {
            c = c_set[j];
            if (b <= c && MAST[b*leafNum + c] > MASTb)
                MASTb = MAST[b*leafNum + c];
            else if (b > c && MAST[c*leafNum + b] > MASTb)
                MASTb = MAST[c*leafNum + b];
        }
        
        MAST[a*leafNum + b] = MASTa + MASTb;
        
        if (MAST[a*leafNum + b] > MASTopt)
            MASTopt = MAST[a*leafNum + b];
    }
    
    return MASTopt;
}

// show the MAST table
void tripleProfile::showMAST(vector<pair<int,int> >& pair_set) {
    int i;
    int a,b;
    cout << "MAST values:" << endl;
    for (i=0; i<(int)pair_set.size(); i++) {
        a = pair_set[i].first;
        b = pair_set[i].second;
        cout << "(" << a << "," << b << "): " << MAST[a*leafNum + b] << endl;
    }
}


// intersect with the triples of the other profile
void tripleProfile::intersect(tripleProfile& p) {
    
    // assume the profile of p is the same size as this profile
    int i;
    for (i=0; i<profileWordNum; i++)
        profile[i] &= p.profile[i];
}

// integer to string
string intToStr(int i) {
    char* d2c = (char*) "0123456789"; // the array for digit to char
    if (i<0) {
        return "-" + intToStr(-i);
    } else if (i<10) {
        return string(1,d2c[i]);
    } else {
        return intToStr(i/10)+string(1,d2c[i%10]);
    }
}

// get the union and the interaction of the leaves of MAST subtree rooted at lca (a,b)
// and get the number of MASTs each leave appears (stored in the array uniLeaves)
// and get the total number of MASTs
void tripleProfile::getUnionInterLeaves(int a, int b, int* uniLeaves, int* intLeaves, int& totalMAST) {
	
	// assume that a <= b
    memset(uniLeaves, 0, leafNum * sizeof(int));
    memset(intLeaves, 0, leafNum * sizeof(int));

	if (a > b) {
		cerr << "[inside getSubKASTsize] Error! a (" << a << ") > b (" << b << ")" << endl;
		exit(1);
	}
	
    if (a==b) {
		uniLeaves[a] = 1;
		intLeaves[a] = 1;
        totalMAST = 1;
        return;
    }
	
	int* left_uniLeaves = new int[leafNum];
	int* left_intLeaves = new int[leafNum];
	int* right_uniLeaves = new int[leafNum];
	int* right_intLeaves = new int[leafNum];
    int left_totMAST = 0;
    int right_totMAST = 0;
	
    vector<int> c_set;
    int MASTa = 0;
    int MASTb = 0;
    int j,k;
    int c;
	
    // compute MASTa
    // get the set of c's such that ( b | a c) is 1
    getTripleSet(b, a, c_set);
    c_set.push_back(a);
    for (j=0; j<(int)c_set.size(); j++) {
        c = c_set[j];
		// since a < b, a <= c < b
        if (MAST[a*leafNum + c] > MASTa)
            MASTa = MAST[a*leafNum + c];
    }
	
	// for each optimal c,
	// get the interaction of the leaves of MASTs rooted at lca(a, c)
    // and get the number of MASTs each leave appears
	int* curr_uniLeaves = new int[leafNum];
	int* curr_intLeaves = new int[leafNum];
    int curr_totMAST;
	int first = 1;
    for (j=0; j<(int)c_set.size(); j++) {
        c = c_set[j];
        if (MAST[a*leafNum + c] == MASTa) {
			if (first) {
				getUnionInterLeaves(a, c, left_uniLeaves, left_intLeaves, left_totMAST);
				first = 0;
			} else {
				getUnionInterLeaves(a, c, curr_uniLeaves, curr_intLeaves, curr_totMAST);
				for (k=0; k<leafNum; k++) {
					left_uniLeaves[k] += curr_uniLeaves[k];
					left_intLeaves[k] = left_intLeaves[k] && curr_intLeaves[k];
				}
                left_totMAST += curr_totMAST;
			}
		}
    }

    // compute MASTb
    // get the set of c's such that ( a | b c ) is 1
    getTripleSet(a, b, c_set);
    c_set.push_back(b);
    for (j=0; j<(int)c_set.size(); j++) {
        c = c_set[j];
		// since a < b, a < c <= b
        if (MAST[c*leafNum + b] > MASTb)
            MASTb = MAST[c*leafNum + b];
    }

	// for each optimal c,
	// get the interaction of the leaves of MASTs rooted at lca(b, c)
    // get the number of MASTs each leave appears
	first = 1;
    for (j=0; j<(int)c_set.size(); j++) {
        c = c_set[j];
        if (MAST[c*leafNum + b] == MASTb) {
			if (first) {
				getUnionInterLeaves(c, b, right_uniLeaves, right_intLeaves, right_totMAST);
				first = 0;
			} else {
				getUnionInterLeaves(c, b, curr_uniLeaves, curr_intLeaves, curr_totMAST);
				for (k=0; k<leafNum; k++) {
					right_uniLeaves[k] += curr_uniLeaves[k];
					right_intLeaves[k] = right_intLeaves[k] && curr_intLeaves[k];
				}
                right_totMAST += curr_totMAST;
			}
		}
    }
	
	// compute the union and the interaction of the leaves of MASTs rooted at lca(a, b)
	for (k=0; k<leafNum; k++) {
		uniLeaves[k] = left_uniLeaves[k]*right_totMAST + right_uniLeaves[k]*left_totMAST;
		intLeaves[k] = left_intLeaves[k] || right_intLeaves[k];
	}
    totalMAST = left_totMAST * right_totMAST;

	/*
	// compute the size of KAST rooted at lca(a, b)
	int KASTsize = 0;
	for (k=0; k<leafNum; k++) {
		KASTsize += intLeaves[k];
	}
	cout << "a=" << a << " b=" << b << " KASTsize=" << KASTsize << endl;
	*/
	
	delete[] left_uniLeaves;
	delete[] left_intLeaves;
	delete[] right_uniLeaves;
	delete[] right_intLeaves;
	delete[] curr_uniLeaves;
	delete[] curr_intLeaves;
}


// get the MAST subtrees rooted at lca(a, b)
void tripleProfile::getSubTrees(int a, int b, set<string>* subTrees, int closedWithBracket, vector<string>& leafNames) {
    
	// assume that a <= b
	if (a > b) {
		cerr << "[inside tripleProfile] Error! a (" << a << ") > b (" << b << ")" << endl;
		exit(1);
	}
	
    if (a==b) {
        subTrees->insert(leafNames[a]);
		// subTrees->insert(intToStr(a));
        return;
    }
    
    vector<int> c_set;
    set<string>* left_subTrees = new set<string>;
    set<string>* right_subTrees = new set<string>;
    set<string>::iterator itr1, itr2;
    
    int MASTa = 0;
    int MASTb = 0;
    int j;
    int c;
    string str;
    
    // compute MASTa
    // get the set of c's such that ( b | a c) is 1
    getTripleSet(b, a, c_set);
    c_set.push_back(a);
    for (j=0; j<(int)c_set.size(); j++) {
        c = c_set[j];
		// since a < b, a <= c < b
        if (MAST[a*leafNum + c] > MASTa)
            MASTa = MAST[a*leafNum + c];
    }

    // get the optimal set of subtrees
    for (j=0; j<(int)c_set.size(); j++) {
        c = c_set[j];
        if (MAST[a*leafNum + c] == MASTa)
            getSubTrees(a, c, left_subTrees, 1, leafNames); // a <= c
    }

    // compute MASTb
    // get the set of c's such that ( a | b c ) is 1
    getTripleSet(a, b, c_set);
    c_set.push_back(b);
    for (j=0; j<(int)c_set.size(); j++) {
        c = c_set[j];
		// since a < b, a < c <= b
        if (MAST[c*leafNum + b] > MASTb)
            MASTb = MAST[c*leafNum + b];
    }

    // get the optimal set of subtrees
    for (j=0; j<(int)c_set.size(); j++) {
        c = c_set[j];
        if (MAST[c*leafNum + b] == MASTb)
            getSubTrees(c, b, right_subTrees, 1, leafNames);
    }
    
    // combine the left and the right subtrees
    for (itr1=left_subTrees->begin(); itr1!=left_subTrees->end(); itr1++) {
        for (itr2=right_subTrees->begin(); itr2!=right_subTrees->end(); itr2++) {
            str="";
            if (closedWithBracket)
                str.append("(");
            if ((*itr1) > (*itr2))
                str.append((*itr2) + "," + (*itr1));
            else
                str.append((*itr1) + "," + (*itr2));
            if (closedWithBracket)
                str.append(")");
            subTrees->insert(str);
        }
    }
    
	left_subTrees->clear();
    delete left_subTrees;
    
	right_subTrees->clear();
    delete right_subTrees;
    
}


// get the number of MAST subtrees rooted at lca(a, b)
int tripleProfile::getSubTreeNum(int a, int b) {
    
	// assume that a <= b
	if (a > b) {
		cerr << "[inside tripleProfile] Error! a (" << a << ") > b (" << b << ")" << endl;
		exit(1);
	}
	
    if (a==b) {
        return 1;
    }
    
    vector<int> c_set;
	int left_subTreeNum = 0;
	int right_subTreeNum = 0;
    
    int MASTa = 0;
    int MASTb = 0;
    int j;
    int c;
    string str;
    
    // compute MASTa
    // get the set of c's such that ( b | a c) is 1
    getTripleSet(b, a, c_set);
    c_set.push_back(a);
    for (j=0; j<(int)c_set.size(); j++) {
        c = c_set[j];
		// since a < b, a <= c < b
        if (MAST[a*leafNum + c] > MASTa)
            MASTa = MAST[a*leafNum + c];
    }

    // get the optimal set of subtrees
    for (j=0; j<(int)c_set.size(); j++) {
        c = c_set[j];
        if (MAST[a*leafNum + c] == MASTa)
			left_subTreeNum += getSubTreeNum(a, c); // a <= c
    }

    // compute MASTb
    // get the set of c's such that ( a | b c ) is 1
    getTripleSet(a, b, c_set);
    c_set.push_back(b);
    for (j=0; j<(int)c_set.size(); j++) {
        c = c_set[j];
		// since a < b, a < c <= b
        if (MAST[c*leafNum + b] > MASTb)
            MASTb = MAST[c*leafNum + b];
    }

    // get the optimal set of subtrees
    for (j=0; j<(int)c_set.size(); j++) {
        c = c_set[j];
        if (MAST[c*leafNum + b] == MASTb)
			right_subTreeNum += getSubTreeNum(c, b);
    }
    
    return left_subTreeNum * right_subTreeNum;
}

// get all results of MAST
void tripleProfile::getAllResults(set<string>& results, int rootNode, vector<string>& leafNames) {
    
    set<string>* subResults = new set<string>;
    set<string>::iterator itr;
    
    int i,j;
    for (i=0; i<leafNum; i++) {
        for (j=i; j<leafNum; j++) {
            if (MAST[i*leafNum+j] == MASTopt) {
                getSubTrees(i, j, subResults, 0, leafNames);
            }
        }
    }
    
    string str;
    Tree tree;
    
    for (itr=subResults->begin(); itr!=subResults->end(); itr++) {
        str = "(" + leafNames[rootNode] + "," + (*itr) + ")";
        tree.buildTreeFrNewick(str);
        results.insert(tree.orderedTopStr());
    }
        
    delete subResults;
                       
}

// get number of MAST
int tripleProfile::getMASTNum(int rootNode) {
    
	int num = 0;
    int i,j;
    for (i=0; i<leafNum; i++) {
        for (j=i; j<leafNum; j++) {
            if (MAST[i*leafNum+j] == MASTopt) {
				num += getSubTreeNum(i, j);
            }
        }
    }
	return num;
}


// get the intersection of leaves of all MASTs
// and get the number of MASTs each leave appear
// and get the total number of MASTs
// Prerequisite: MAST table has been computed
void tripleProfile::getUnionInterMASTleaves() {
	
	if (unionLeaves != NULL) {
		delete[] unionLeaves;
	}
	if (interLeaves != NULL) {
		delete[] interLeaves;
	}

	// allocate the memory for KAST, unionLeaves and interLeaves
	unionLeaves = new int[leafNum];
	interLeaves = new int[leafNum];

	int* currUniLeaves = new int[leafNum];
	int* currIntLeaves = new int[leafNum];
    int currTotMAST = 0;
	
	// initialize the arrays
	memset(unionLeaves, 0, leafNum*sizeof(int));
	memset(interLeaves, 0, leafNum*sizeof(int));
    totalMAST = 0;
	
    int i,j,k;
    int first = 1;
    for (i=0; i<leafNum; i++) {
        for (j=i; j<leafNum; j++) {
            if (MAST[i*leafNum+j] == MASTopt) {
                if (first) {
                    getUnionInterLeaves(i, j, unionLeaves, interLeaves, totalMAST);
                    first=0;
                } else {
                    getUnionInterLeaves(i, j, currUniLeaves, currIntLeaves, currTotMAST);
                    for (k=0; k<leafNum; k++) {
                        unionLeaves[k] += currUniLeaves[k];
                        interLeaves[k] = interLeaves[k] && currIntLeaves[k];
                    }
                    totalMAST += currTotMAST;
                }
            }
        }
    }
    
	delete[] currUniLeaves;
	delete[] currIntLeaves;
}



// to compute MAST and KAST of all the input trees
// if "outputMAST" is 1, then get all the MAST trees.
void obtainMASTKAST(vector<string>& intrees, vector<string>& MASTs, string& KAST, int& MAST_score, int& KAST_score, int& totalMAST,
    int outputMAST, vector<int>& MASTunionLeaves, vector<int>& MASTinterLeaves, vector<string>& leafNames) {
	
    MAST_score = 0;
    KAST_score = 0;
    totalMAST = 0;
    MASTs.clear();
    KAST = "";
	MASTunionLeaves.clear();
	MASTinterLeaves.clear();
	leafNames.clear();
	
    if (intrees.size() == 0) {
		return;
    }

    set<string> results;
    set<string>::iterator itr;
	
    vector<int> leafIDs;
    vector<int> commonLeaves;
    map<string,int> leafOrder;
    map<string,int>::iterator mapItr;
    string currLeaf;
    vector<Tree*> allTrees;
    Tree* tree;
	
	set<int> treeSizes;
	set<int>::iterator iitr;
	
	int isAllRootedTree = -1; // 0: all unrooted trees; 1: all rooted trees; 2: mixed; -1: initial value
    int i,j;

    
    // build all the trees
    for (i=(int)intrees.size()-1; i>=0; i--) {
        tree = new Tree();
        tree->buildTreeFrNewick(intrees[i]);
		if (isAllRootedTree == -1)
			isAllRootedTree = tree->isRooted;
		else if (isAllRootedTree != tree->isRooted)
			isAllRootedTree = 2;
        allTrees.push_back(tree);
    }
    
    // obtain the set of leaves appearing in any tree
    // and the set of leaves common in all trees
    for (i=0; i<(int)allTrees.size(); i++) {
        tree = allTrees[i];
		treeSizes.insert((int)tree->leafList.size());
        // examine all the leaves
        for (j=0; j<(int)tree->leafList.size(); j++) {
            currLeaf = tree->leafList[j]->leafName;
            mapItr = leafOrder.find(currLeaf);
            if (mapItr == leafOrder.end()) {
                leafOrder.insert(pair<string,int>(currLeaf,leafOrder.size()));
                leafNames.push_back(currLeaf);
                commonLeaves.push_back(1);
            } else {
                commonLeaves[mapItr->second]++;
            }
        }
    }
	
	cout << "Input tree size";
	if (treeSizes.size() > 1)
		cout << "s";
	cout << " : ";
	for (iitr=treeSizes.begin(); iitr!=treeSizes.end(); iitr++) {
		if (iitr!=treeSizes.begin())
			cout << ", ";
		cout << (*iitr);
	}
	
	// special case: # of trees = 1
	if (allTrees.size()==1) {
		if (intrees[0].length() > 0 && intrees[0].at(intrees[0].length()-1)!=';')
			intrees[0] += ";";
		MASTs.push_back(intrees[0]);
        KAST = intrees[0];
        MAST_score = (int)allTrees[0]->leafList.size();
        KAST_score = (int)allTrees[0]->leafList.size();
		for (i=0; i<(int)allTrees[0]->leafList.size(); i++) {
			MASTunionLeaves.push_back(1);
			MASTinterLeaves.push_back(1);
		}
		cout << endl;
        totalMAST = 1;
		return;
	}

	switch(isAllRootedTree) {
		case 0:
			cout << " (All trees are unrooted)" << endl;
			break;
		case 1:
			cerr << " (All trees are rooted)" << endl << endl;
			cerr << "Error! This module for computation of MAST and KAST is not designed for rooted trees." << endl;
			exit(1);
			break;
		case 2:
			cerr << " (Some trees are rooted and some are unrooted)" << endl << endl;
			cerr << "Error! This module for computation of MAST and KAST is not designed for rooted trees." << endl;
			exit(1);
			break;
	}


    for (i=0; i<(int)commonLeaves.size(); i++) {
        if (commonLeaves[i] == (int)intrees.size()) {
            commonLeaves[i] = 1;
        } else {
            commonLeaves[i] = 0;
        }
    }    

	// special case:
	// # of common leaves <= 3 if all trees are unrooted
	// # of common leaves <= 2 if all trees are rooted
	
    int numCommonLeaves = 0;
    for (i=0; i<(int)commonLeaves.size(); i++) {
        numCommonLeaves += commonLeaves[i];
    }
    if ((numCommonLeaves <= 3 && (isAllRootedTree==0)) || (numCommonLeaves <= 2 && (isAllRootedTree==1))) {
        if (numCommonLeaves > 1)
            KAST.append("(");
		int first = 1;
        for (i=0; i<(int)commonLeaves.size(); i++) {
            if (commonLeaves[i]) {
                if (first) {
                    KAST.append(leafNames[i]);
                    first = 0;
                } else {
                    KAST.append(","+leafNames[i]);
                }
            }
        }
        if (numCommonLeaves > 1)
            KAST.append(")");
		if (numCommonLeaves > 0) {
			KAST.append(";");
			MASTs.push_back(KAST);
		}
        KAST_score = numCommonLeaves;
        MAST_score = numCommonLeaves;
        totalMAST = 1;
		for (i=0; i<(int)commonLeaves.size(); i++) {
			MASTunionLeaves.push_back(commonLeaves[i]);
			MASTinterLeaves.push_back(commonLeaves[i]);
		}
		return;
    }
	
    // assign leaf ID
    for (i=0; i<(int)allTrees.size(); i++) {
        allTrees[i]->assignLeafID(leafOrder);
    }
    
    
	//================================================================================
	// Start computation of MAST and KAST
	//================================================================================
	
	// initialize MASTunionLeaves and MASTinterLeaves
	for (i=0; i<(int)commonLeaves.size(); i++) {
		MASTunionLeaves.push_back(0);
		MASTinterLeaves.push_back(commonLeaves[i]);
	}

    // the total number of leaves
    int numLeaves = (int) leafNames.size();
    
    int rootID;
    int treeID;
    string curr_unrooted_tree;
    int MAST_curr_score;
    tripleProfile currProfile(numLeaves);
	tripleProfile commonProfile(numLeaves);
	PairList pairList;
    
    // since MAST is only applied on a rooted tree
    // we need to consider the tree be rooted at each leaf one by one
    vector<int> leavesConsidered;
	leavesConsidered.assign(commonLeaves.begin(), commonLeaves.end());

    cout << "computing MAST and KAST";
    
    for (rootID=0; rootID<(int)leafOrder.size(); rootID++) {
        
		cout << "." << flush;
		
		// only consider those leaves common to all subtrees
        if (commonLeaves[rootID]==0)
            continue;
        
        // the leaf rootID should not be considered as a leaf inside the subtrees rooted at rootID
        leavesConsidered[rootID] = 0;
		
		// =============================================================
		// compute the profile common to all the trees rooted at rootID
		// =============================================================
        
		int compute_pairList = 1;
        allTrees[0]->buildTripleProfile(rootID, commonProfile, pairList, compute_pairList, leavesConsidered);
        
		compute_pairList = 0;
        for (treeID=1; treeID<(int)allTrees.size(); treeID++) {
            allTrees[treeID]->buildTripleProfile(rootID, currProfile, pairList, compute_pairList, leavesConsidered);
            commonProfile.intersect(currProfile);
        }
        
		// =============================================================
		// compute the MAST for all the trees rooted at rootID
		// =============================================================
		
        MAST_curr_score = commonProfile.computeMAST(pairList, leavesConsidered)+1;
        
		// =============================================================
		// if the size of MAST >= the size of the current optimal MAST,
		// then (1) get the intersection of the MASTs, 
        //          get the number of MASTs in which each leave appears,
        //          get the total number of MASTs
		//      (2) get the MAST if needed
		// =============================================================
		
        if (MAST_curr_score >= MAST_score) {
            if (MAST_curr_score > MAST_score) {
                MAST_score = MAST_curr_score;
                totalMAST = 0;
				// initialize MASTunionLeaves and MASTinterLeaves
				for (i=0; i<(int)commonLeaves.size(); i++) {
					MASTunionLeaves[i]=0;
					MASTinterLeaves[i]=commonLeaves[i];
				}
				if (outputMAST)
					results.clear();
            }
			// get the union and intersection of the MASTs
			commonProfile.getUnionInterMASTleaves();
			commonProfile.unionLeaves[rootID] = commonProfile.totalMAST;
			commonProfile.interLeaves[rootID] = 1;
			for (i=0; i<numLeaves; i++) {
				MASTunionLeaves[i] += commonProfile.unionLeaves[i];
				MASTinterLeaves[i] = (MASTinterLeaves[i] && commonProfile.interLeaves[i]);
			}
            totalMAST += commonProfile.totalMAST;
			if (outputMAST) {
				// get all result
				commonProfile.getAllResults(results, rootID, leafNames);
			}
            
            /*
            // for debugging
            // list out the content of unionLeaves
            cout << endl;
			for (i=0; i<numLeaves; i++) {
				cout << leafNames[i] << " " << MASTunionLeaves[i] << endl;
			}
            // list out all the MASTs
            for (itr=results.begin(); itr!=results.end(); itr++) {
                cout << (*itr) << endl;
            }
            exit(1);
            */
        }
        
        // reset leavesConsidered[rootID] to 1
        // leavesConsidered[rootID] = 1;
        
    }
    
	// get KAST
	for (i=0; i<numLeaves; i++)
		KAST_score += MASTinterLeaves[i];
	if (KAST_score > 0)
		KAST = allTrees[0]->getSubTree(MASTinterLeaves);
	
	// get MAST if necessary
	if (outputMAST) {
		for (itr=results.begin(); itr!=results.end(); itr++) {
			MASTs.push_back(*itr);
		}
	}
	
	cout << endl;
    
}
