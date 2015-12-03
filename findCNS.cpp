#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <math.h>
#include <sys/time.h>
#include <vector>
#include <sstream>
#include <cstdlib>
#include <string.h>
#include <algorithm>
#include <stdint.h>
//#include <forward_list>
#include <list>
#include <stack>
#include <limits.h>
#include <vector>
#define NINF INT_MIN

using namespace std;
void seqPos(int *s, int y, int *n, int *z);


// g++ -O3 -o splitMEM splitMEM.cc

//on yellowstone
// g++  -std=c++0x  -O3 -o splitMEM splitMEM.cc

// dot -Tpdf tree.dot -o tree.pdf

// Default length, can override at runtime
int Kmer_Len = 30;
string CDG_Filename = "cdg.dot"; //for output of compressed de Bruijn graph in dot format

int DEBUG = 0;
int VERBOSE = 0;
int VERIFY = 0;

typedef uint32_t treeint;
typedef uint64_t treeintLarge;

treeintLarge skippedbases = 0;
treeintLarge skippedextensions = 0;

treeint numKmerLens = 0;
treeint numNodesWithTable = 0;
treeint numEntriesInAuxTables = 0;

// Settings for linear time alg
bool FORCEROOT = false;
bool DOJUMP = true;
bool DOINTERNALSKIP = true;
bool DOPHASETRICK = true;

bool MEM = false;

const int basecount = 6;
int b2i(char base)
{
    switch (base)
    {
        case '$' : return 0;
        case 'A' : return 1;
        case 'C' : return 2;
        case 'G' : return 3;
        case 'N' : return 4;
        case 'T' : return 5;
            
        default:
            cerr << "Unknown base: " << base << endl;
            return b2i('N');
    };
}

class MerVertex_t;

class SuffixNode
{
public:
    static treeintLarge s_nodecount;
    
    SuffixNode(
               treeint  s,
               treeint  e,
               SuffixNode * p,
               SuffixNode * x)
    : m_start(s),
    m_end(e),
    m_parent(p)//,
    {
        s_nodecount++;
        m_isCopy = false;
        m_suffixTable = new SuffixNode*[1]; //for first suffix link, replaced with a larger table later
        
        m_suffixTable[0] = x;
        m_strdepth = 0;
        
        for (int i = 0; i < basecount; i++)
        {
            m_children[i] = NULL;
        }
    }
    
    //copy constructor
    SuffixNode(const SuffixNode & node)
    {
        m_start = node.m_start;
        m_end = node.m_end;
        m_SA_start = node.m_SA_start;
        m_SA_end = node.m_SA_end;
        
        m_parent = node.m_parent;
        for (int i = 0; i < basecount; i++)
        {
            m_children[i] = NULL;
        }
        m_isCopy = true;
        
        m_suffixTable = node.m_suffixTable;
        m_LMAproximityTable= node.m_LMAproximityTable;
        m_LMAprox_nodeTable = node.m_LMAprox_nodeTable;
        
        m_MEM = node.m_MEM;
        m_LMA = node.m_LMA;
        m_strdepth = node.m_strdepth;
        
    }
    
    void createSuffixLinkTable(treeint numSuffixLinks, bool otherTables = true)
    {
        if(numSuffixLinks==0)
            numSuffixLinks = 1; //min number of entries in table
        
        m_suffixTable = new SuffixNode*[numSuffixLinks];
        
        if(otherTables)
        {
            numNodesWithTable++;
            
            m_LMAproximityTable = new treeint[numSuffixLinks];
            m_LMAprox_nodeTable = new SuffixNode*[numSuffixLinks];
        }
        else
        {
            m_LMAproximityTable = NULL;
            m_LMAprox_nodeTable = NULL;
            
        }
    }
    
    void deleteSuffixLinkTable()
    {
        SuffixNode * suffixLink = m_suffixTable[0];
        
        if(m_suffixTable)
        {
            delete [] m_suffixTable;
        }
        if(m_LMAproximityTable)
        {
            delete [] m_LMAproximityTable;
            m_LMAproximityTable = NULL;
        }
        if(m_LMAprox_nodeTable)
        {
            delete [] m_LMAprox_nodeTable;
            m_LMAprox_nodeTable = NULL;
            
        }
        
        createSuffixLinkTable(1, false);
        m_suffixTable[0] = suffixLink;
        
    }
    
    ~SuffixNode()
    {
        if(!m_isCopy)
        {
            if(m_suffixTable)
                delete [] m_suffixTable;
            if(m_LMAproximityTable)
                delete [] m_LMAproximityTable;
            if(m_LMAprox_nodeTable)
                delete [] m_LMAprox_nodeTable;
        }
        
        for (int i = 0; i < basecount; i++)
        {
            if (m_children[i]) { delete m_children[i]; }
        }
    }
    
    string str(const string & s)
    {
        return s.substr(m_start, m_end-m_start+1);
    }
    
    treeint len(int i=-1)
    {
        if (i != -1)
        {
            if (i < m_end)
            {
                return i - m_start + 1;
            }
        }
        
        return m_end - m_start + 1;
    }
    
    ostream & printLabel(ostream & os, const string & str)
    {
        string seq = str.substr(m_start, m_end-m_start+1);
        
        if (m_start == m_end && m_start == 0)
        {
            os << "\"ROOT\"";
        }
        else
        {
            os << "\"" << seq;
            if(m_MEM)
            {
                os << " MEM ";
            }
            os << ":"
            << " [" << m_start
            << "," << m_end << "]"
            << "}\"";
        }
        
        return os;
    }
    
    
    ostream & printNodeLabel(ostream & os)
    {
        os << m_start << m_end << m_strdepth ;
        
        return os;
    }
    
    ostream & printEdgeLabel(ostream & os, const string & str)
    {
        string seq = str.substr(m_start, m_end-m_start+1);
        os << "\"" << seq << (MEM?(m_MEM ? " MEM " :""):"")
        << " SA [" << m_SA_start
        << ", " << m_SA_end <<"]"
        << " [" << m_start
        << "," << m_end << "]\"";
        
        return os;
    }
    
    //returns lowest marked ancestor (if this node is marked, it is the lowest)
    //if no marked ancestor, NULL
    SuffixNode* LMA()
    {
        return m_LMA;
    }
    
    //return num suffix links in table based on strdepth
    inline int nodeNumSuffixLinks()
    {
        if(m_strdepth == 0)
            return 0; //for root
        
        return (int)(ceil(log(m_strdepth) / log(2)));
    }
    
    treeint  m_start;
    treeint  m_end;
    treeint  m_strdepth;
    
    bool m_isCopy; //set to true for copy constructor. then suffix link table will not be deleted when the function returns
    bool m_MEM;
    bool m_prevChar[basecount];
    
    typedef SuffixNode * SuffixNodePtr;
    
    SuffixNode * m_parent;
    SuffixNode * m_children [basecount];
    SuffixNodePtr * m_suffixTable;
    
    treeint * m_LMAproximityTable;  //[i] stores smallest m_LMAproximity encountered in m_suffixTable[i]
    SuffixNodePtr * m_LMAprox_nodeTable; //stores LMAnode corresponding to m_LMAproximityTable entries
    
    SuffixNode * m_LMA;
    //interval in suffix array that corresponds to this node
    treeint m_SA_start;
    treeint m_SA_end;
};
treeintLarge SuffixNode::s_nodecount(1);


ostream & operator<< (ostream & os, SuffixNode * n)
{
    return n->printNodeLabel(os);
}


class SuffixTree
{
public:
    SuffixTree(const string & s)
    : m_nodecount(0), m_string(s), m_maxMEMstrdepth(0)
    {
        m_root = new SuffixNode(0,0,NULL,NULL);
        m_root->m_suffixTable[0] = m_root;
        m_suffixArray = new treeint[s.length()];
    }
    
    treeint * m_suffixArray;
    SuffixNode * m_root;
    treeintLarge m_nodecount;
    string m_string;
    treeint m_maxMEMstrdepth; //max strdepth at any MEM node; only nodes whose strdepth \leq maxMEMstrdepth have suffix link tables and skipped MEM tables, each one is the size of their strdepth
    //std::forward_list<SuffixNodeMark> nodesWithSuffixSkips;
    SuffixNode ** nodesWithSuffixSkips;
    treeint nextNodeWithSuffixSkips; //next position in array to fill
    treeint numNodesWithSuffixSkips;
    
    
    void resetMaxMemStrdepth()
    {
        m_maxMEMstrdepth = 0;
    }
    
    void dumpNode(SuffixNode * node)
    {
        int children = 0;
        for (int i = 0; i < basecount; i++)
        {
            SuffixNode * child = node->m_children[i];
            if (child)
            {
                children++;
                
                cout << " " << node << "->" << child;
                
                cout << " [minlen=" << child->len() << ", label=";
                child->printEdgeLabel(cout, m_string) << "]" << endl;
                
                dumpNode(child);
            }
        }
        
        if (node->m_suffixTable[0])
        {
            cout << " " << node << " -> " << node->m_suffixTable[0]
            << " [style=dotted, constraint=false]" << endl;
        }
        
        if (children == 0)
        {
            cout << " " << node << " [shape=box, label=";
            node->printLabel(cout, m_string) << "]" << endl;
        }
        else
        {
            cout << " " << node << " [label=";
            node->printLabel(cout, m_string) << "]" << endl;
        }
    }
    
    void dumpTree()
    {
        cerr << "Dumping tree" << endl;
        cout << "digraph G {" << endl;
        cout << " size=\"7.5,10\"" << endl;
        cout << " center=true" << endl;
        cout << " label=\"Suffix tree of \'" << m_string << "\' len:"
        << m_string.length()-1 << " nc:"
        << m_nodecount << "\"" << endl;
        dumpNode(m_root);
        cout << "}" << endl;
    }
    
    void dumpNodeText(ostream & out, SuffixNode * n, treeint depth)
    {
        for (int b = 0; b < basecount; b++)
        {
            if (n->m_children[b])
            {
                for (treeint i = 0; i < depth; i++)
                {
                    out << " ";
                }
                out << " ";
                out << n->m_children[b]->str(m_string) << endl;
                dumpNodeText(out, n->m_children[b], depth+1);
            }
        }
    }
    
    void dumpTreeText(ostream & out)
    {
        out << "Suffix Tree len=" << m_string.length()-1 << endl;
        out << "String: \"" << m_string << "\"" << endl;
        out << "+" << endl;
        dumpNodeText(out, m_root, 0);
    }
    
    
    void dumpTreeSorted(ostream & out, SuffixNode * node, const string & pathstring)
    {
        int c = 0;
        
        string mystring = node->str(m_string);
        string ps(pathstring);
        ps.append(mystring);
        
        for (int i = 0; i < basecount; i++)
        {
            if (node->m_children[i])
            {
                c++;
                dumpTreeSorted(out, node->m_children[i], ps);
            }
        }
        
        if (c == 0)
        {
            out << ps << endl;
        }
    }
    
    
    SuffixNode * createNode(treeint s, treeint e, SuffixNode * p, SuffixNode * x)
    {
        SuffixNode * retval = new SuffixNode(s, e, p, x);
        m_nodecount++;
        
        return retval;
    }
    
    //error to call this function with k<1.
    bool fillKhopSuffixLinkedList(SuffixNode * node, treeint k)
    {
        bool hasSuffixSkip = false;
        //calculate num suffix links at node: log base 2 of strdepth at node
        int nodeSuffixLinks =  node->nodeNumSuffixLinks();
        
        int children = 0;
        for (int b = 0; b < basecount; b++)
            if(node->m_children[b])
                children++;
        
        //SM added last clause 8/6/14
        if(children == 0 || node->m_strdepth > m_maxMEMstrdepth || node->m_strdepth < Kmer_Len)  //no tables to fill for this node
            return hasSuffixSkip;
        
        if(nodeSuffixLinks > k)
        {
            hasSuffixSkip = true;
            if(node->m_suffixTable[k-1])
            {
                //calculate suffix link that is 2^k hops away from node
                SuffixNode * kNeighbor = node->m_suffixTable[k-1]; //really the (k-1)Neighbor
                int kNeighborSuffixLinks =  kNeighbor->nodeNumSuffixLinks();
                
                if(kNeighborSuffixLinks > k-1)
                {
                    node->m_suffixTable[k] = kNeighbor->m_suffixTable[k-1];
                    
                    if(kNeighbor->m_LMAproximityTable[k-1] < node->m_LMAproximityTable[k-1])
                    {
                        node->m_LMAproximityTable[k] = kNeighbor->m_LMAproximityTable[k-1];
                        node->m_LMAprox_nodeTable[k] = kNeighbor->m_LMAprox_nodeTable[k-1];
                    }
                    else
                    {
                        node->m_LMAproximityTable[k] = node->m_LMAproximityTable[k-1];
                        node->m_LMAprox_nodeTable[k] = node->m_LMAprox_nodeTable[k-1];
                    }
                }
                else
                {
                    node->m_suffixTable[k] = NULL;
                    node->m_LMAproximityTable[k] = node->m_LMAproximityTable[0];
                    node->m_LMAprox_nodeTable[k] = node->m_LMAprox_nodeTable[0];
                    //return; //chidren won't be any deeper as far as strdepth //SM added 6/12/14

                }
            }
            else  //only applies to root
            {
                node->m_suffixTable[k] = NULL;
            }
        }
        
        return hasSuffixSkip;
     }
    
    //error to call this function with k<1.
    void fillKhopSuffixNode(SuffixNode * node, treeint k, bool fillLinkedList)
    {
        int nodeSuffixLinks =  node->nodeNumSuffixLinks();
        
        int children = 0;
        for (int b = 0; b < basecount; b++)
            if(node->m_children[b])
                children++;
        
        
        if(children == 0 || node->m_strdepth > m_maxMEMstrdepth )  //no tables to fill for leaf
            return;
        
        
        if(nodeSuffixLinks > k)
        {
            if(fillLinkedList)
            {
                //nodesWithSuffixSkips.push_front(node);  //SM added 6/16/14
                nodesWithSuffixSkips[nextNodeWithSuffixSkips++] = node;
            }
            
            numNodesWithSuffixSkips++;
            
            if(node->m_suffixTable[k-1])
            {
                //calculate suffix link that is 2^k hops away from node
                SuffixNode * kNeighbor = node->m_suffixTable[k-1]; //really the (k-1)Neighbor
                int kNeighborSuffixLinks =  kNeighbor->nodeNumSuffixLinks();
                
                if(kNeighborSuffixLinks > k-1)
                {
                    node->m_suffixTable[k] = kNeighbor->m_suffixTable[k-1];
                    
                    if(kNeighbor->m_LMAproximityTable[k-1] < node->m_LMAproximityTable[k-1])
                    {
                        node->m_LMAproximityTable[k] = kNeighbor->m_LMAproximityTable[k-1];
                        node->m_LMAprox_nodeTable[k] = kNeighbor->m_LMAprox_nodeTable[k-1];
                    }
                    else
                    {
                        node->m_LMAproximityTable[k] = node->m_LMAproximityTable[k-1];
                        node->m_LMAprox_nodeTable[k] = node->m_LMAprox_nodeTable[k-1];
                    }
                }
                else
                {
                    node->m_suffixTable[k] = NULL;
                    node->m_LMAproximityTable[k] = node->m_LMAproximityTable[0];
                    node->m_LMAprox_nodeTable[k] = node->m_LMAprox_nodeTable[0];
                }
            }
            else  //only applies to root
            {
                node->m_suffixTable[k] = NULL;
            }
        }
        
        
        for (int b = 0; b < basecount; b++)
        {
            if(node->m_children[b])
            {
                fillKhopSuffixNode(node->m_children[b], k, fillLinkedList);
            }
        }
    }
    
 
    
    //recursively perform DFS to calc k-hop suffix links for each node
    void fillKhopSuffix(treeint k, bool fillLinkedList)
    {
        if(fillLinkedList) //dynamically allocate array
        {
            nodesWithSuffixSkips = new SuffixNode*[m_nodecount];
            nextNodeWithSuffixSkips = 0; //array position to fill next
        }
            
        fillKhopSuffixNode(m_root, k, fillLinkedList);
    }
    
   
    
    //this function uses linked list instead of multiple recursive DFSs over suffix tree
    //this function populates the table m_SuffixLinkTable and associated skipped LMA tables
    //entry i is the suffix link that is 2^i hops away
    //traverse suffix tree floor(log n) = m_numSuffixLinks times to populate the table
    void fillSuffixTableNonRecursively()
    {
        //timeval starttime;
        //timeval endtime;
        bool hasSuffixSkip;
        treeint nextSNodeSetup; //where up to in array
        treeint nextSNodeReplace; //where to copy to in array
        treeint totalSNodesWithSuffixSkips; //where working array ends
        //forward_list<SuffixNodeMark>::iterator it;
        
        treeint seconds;
        treeint microseconds;
        double elapsed;

        
        
        //begin with 1 since m_suffixLInkTable[0] is calcuated as part of Ukkonen's construction algorithm
        unsigned int numSuffixLinks = (int)( ceil(log(m_maxMEMstrdepth) / log(2))); //floor(log(m_maxMEMstrdepth - Kmer_Len) / log(2));
        
        cerr<<" numSuffixLinks="<<numSuffixLinks<<endl;
        
        
        //fill first suffix skip entries by traversing entire tree and filling nodesWithSuffixSkips linked list with only nodes that have suffix skips
        int k = 1;
        
        numNodesWithSuffixSkips = 0;
        
        //gettimeofday(&starttime, NULL);

        fillKhopSuffix(k, true);
       
        
        totalSNodesWithSuffixSkips = nextNodeWithSuffixSkips;
        
        cerr<<" total nodes with suffix skips ="<<nextNodeWithSuffixSkips<<endl;
        cerr<<" numNodesWithSuffixSkips = "<<numNodesWithSuffixSkips<<endl;
        cerr<<" now using linked lists"<<endl;
        
        
        for(k = 2; k < numSuffixLinks; k++)
        {
            numNodesWithSuffixSkips = 0;
            
            //gettimeofday(&starttime, NULL);

            nextSNodeReplace = 0;
           
            for( nextSNodeSetup = 0; nextSNodeSetup<totalSNodesWithSuffixSkips; nextSNodeSetup++)
            {
                hasSuffixSkip = fillKhopSuffixLinkedList(nodesWithSuffixSkips[nextSNodeSetup], k);
                if(hasSuffixSkip)
                {
                    nodesWithSuffixSkips[nextSNodeReplace++] = nodesWithSuffixSkips[nextSNodeSetup];
                    numNodesWithSuffixSkips++;
                }
                //else
                    //nothing to do: loop advances nextSNodeSetup
            
            }
            
            totalSNodesWithSuffixSkips = numNodesWithSuffixSkips;
           
            
            
            cerr<<" numNodesWithSuffixSkips = "<<numNodesWithSuffixSkips<<endl;
            
          
        }
        
        delete [] nodesWithSuffixSkips;
        
    }
    
       
    void preprocessLMAnode(SuffixNode * node, SuffixNode * markedNode)
    {
        if(node->m_MEM && node->m_strdepth >= Kmer_Len)
        {
            markedNode = node;
        }
        node->m_LMA = markedNode;
        
        for (int b = 0; b < basecount; b++)
        {
            if(node->m_children[b])
            {
                preprocessLMAnode(node->m_children[b], markedNode);
            }
        }
        
        
    }
    
    //mark internal nodes that represent MEMs of length >=minMEM by setting m_MEM = true
    void preprocessLMA( )
    {
        //recursively perform DFS on suffix tree and set LMA to point to lowest ancestor that is marked as a MEM
        preprocessLMAnode(m_root, NULL);
    }
    
    
    //recursively unmark all MEM nodes
    //can stop when reach a node that is long enough to be a MEM but isn't
    void unmarkMEMnode(SuffixNode * node, int minMEM)
    {
        for (int b = 0; b < basecount; b++)
        {
            if(node->m_children[b])
                unmarkMEMnode(node->m_children[b], minMEM);
        }
        
        node->m_MEM = false;
        node->m_LMA = NULL;
        
        node->deleteSuffixLinkTable();
        
        
    }
    
    void unmarkMEMnodes(int minMEM)
    {
        unmarkMEMnode(m_root, minMEM);  //root cannot be marked node
    }
    
    
    //recursively mark node if it's a MEM, fill suffix array along the way with same DFS
    void markMEMnode(SuffixNode * node, treeint minMEM, treeint * nextSAentry, bool setupSA)
    {
        int children = 0;
        for (int b = 0; b < basecount; b++)
        {
            if(node->m_children[b])
            {
                children++;
            }
        }
        
        if(setupSA)
        {
            //when visit a node the first time: set string depth in m_strdepth to reflect the length of the labels on the path from the root to this node's end
            if(node == m_root)  //set string depth = length of labels on path from root to this node's end
                node->m_strdepth = 0;
            else
                node->m_strdepth = node->m_parent->m_strdepth + node->len();
            
            //fill next suffix array element for this leaf
            if(children == 0)
            {
                m_suffixArray[*nextSAentry] = node->m_start - node->m_parent->m_strdepth;
                
                node->m_SA_start = *nextSAentry;
                node->m_SA_end = *nextSAentry;
                
                (*nextSAentry)++;
            }
        }
        
        //perform DFS recursively
        bool firstChild = true;
        
        
        for (int b = 0; b < basecount; b++)
        {
            if(node->m_children[b])
            {
                markMEMnode(node->m_children[b], minMEM, nextSAentry, setupSA);
                for(int c = 0; c < basecount; c++)
                {
                    if(node->m_children[b]->m_prevChar[c])
                        node->m_prevChar[c] = true;
                    
                }
                if(setupSA)
                {
                    if(firstChild) //first child
                    {
                        node->m_SA_start = node->m_children[b]->m_SA_start;
                        node->m_SA_end = node->m_children[b]->m_SA_end;
                        firstChild = false;
                    }
                    else
                    {
                        if(node->m_children[b]->m_SA_start < node->m_SA_start)
                            node->m_SA_start = node->m_children[b]->m_SA_start;
                        
                        if(node->m_children[b]->m_SA_end > node->m_SA_end)
                            node->m_SA_end = node->m_children[b]->m_SA_end;
                    }
                }
            }
        }
        
        //when visit node the second time: set m_prevChar to true for its children's values (its value if leaf).  if more than one array element is set to true, set m_MEM to true if m_strdepth >= minMEM
        if(node == m_root)
        {
            node->m_MEM = false;
        }
        else if(children == 0) //leaf, cannot be MEM
        {
            node->m_MEM = false;
            treeint prevCharPos = node->m_end - node->m_strdepth;
            //if(prevCharPos >= 0) //ignore position 1, treat as $ before string so will be considered maximal (can't extend beyond beginning of string)
            {
                char prevChar;
                if(prevCharPos == 0)
                    prevChar = '$';
                else
                    prevChar = m_string[prevCharPos];
                
                int prevCharNum = b2i(prevChar);
                node->m_prevChar[prevCharNum] = true;
            }
        }
        else //internal node
        {
            if(node->m_strdepth >= minMEM) //see if this node is left maximal
            {
                
                int numPrevChars = 0;
                for(int b = 0; b < basecount; b++)
                {
                    if(node->m_prevChar[b])
                        numPrevChars++;
                }
                if(numPrevChars>=2) //can short-circuit OR if a child is a MEM
                {
                    node->m_MEM = true;
                    
                    numKmerLens++;
                    
                    if(node->m_strdepth > m_maxMEMstrdepth)
                        m_maxMEMstrdepth = node->m_strdepth;
                    
                    //when visit MEM node the second time, also propagate up the start positions from the children, adjusted to account for the length of this node
                    
                    
                }
                else
                    node->m_MEM = false;
            }
            
        }
        
    }
    
    //mark internal nodes that represent MEMs of length >=minMEM by setting m_MEM = true
    //fill suffix array with same DFS in suffix tree
    void markMEMnodes(int minMEM, bool setupSA)
    {
        cerr<<" marking MEM nodes in suffix tree of at least "<<minMEM<<" bp"<<endl;
        treeint nextSAentry = 0;
        //recursively perform DFS on suffix tree and mark internal nodes that are left maximal as MEMs
        markMEMnode(m_root, minMEM, &nextSAentry, setupSA);
    }
    
    
    
    void createAuxTablesAtNode(SuffixNode * node)
    {
        SuffixNode * suffixLink = node->m_suffixTable[0];
        if(node->m_suffixTable != NULL)
            delete [] node->m_suffixTable;
        
        int children = 0;
        for (int b = 0; b < basecount; b++)
        {
            if(node->m_children[b])
            {
                children++;
            }
        }
        
        if(node == m_root || children == 0 || node->m_strdepth > m_maxMEMstrdepth)  //leaves don't need tables for suffix links and skipped LMA
        {
            node->createSuffixLinkTable(1, false);
        }
        else
        {
            //don't go through all suffix links until root.  stop when strdepth of node is < Kmer_Len
            int numSuffixLinks = node->nodeNumSuffixLinks();
            if(numSuffixLinks==0)
                numSuffixLinks = 1; //min table size
            
            node->createSuffixLinkTable(numSuffixLinks, true);
        }
        
        node->m_suffixTable[0] = suffixLink;
        
        
        //from here comes from fillFirstLMAproxEntry_Node
    
        if(node != m_root && children > 0 && node->m_strdepth <= m_maxMEMstrdepth)
        {
            int minProx, thisProx, sProx;
            SuffixNode * snode = node->m_suffixTable[0];
            
            if(node->m_LMA)
                thisProx = node->m_strdepth - node->m_LMA->m_strdepth;  //this node's info is first entry
            else
                thisProx = node->m_strdepth;
            
            
            if(snode)
            {
                if(snode->m_LMA)
                {
                    sProx = snode->m_strdepth - snode->m_LMA->m_strdepth;  //this node's info is first entry
                }
                else
                {
                    sProx = snode->m_strdepth;
                }
                if(sProx < thisProx)
                    minProx = sProx;
                else
                    minProx = thisProx;
            }
            else
            {
                minProx = thisProx;
            }
            node->m_LMAproximityTable[0] = minProx;
            
            if(minProx == thisProx)
            {
                node->m_LMAprox_nodeTable[0] = node->LMA();
            }
            else
            {
                node->m_LMAprox_nodeTable[0] = snode->LMA();
            }
        } //to here
        
        for (int b = 0; b < basecount; b++)
        {
            if(node->m_children[b])
                createAuxTablesAtNode(node->m_children[b]);
        }
    }

    
    void createAuxTables()
    {
        numNodesWithTable = 0; //reset variable so can use in limited output as tables are created
        cerr<<" creating aux tables at nodes"<<endl;
        //recursively perform DFS on suffix tree and create suffix link and skipped LMA tables when appropriate
        createAuxTablesAtNode(m_root);
    }
    
     
  
    
    void buildUkkonen()
    {
        treeint len = m_string.length() - 1; // length of the string, not of the buffer
        char base = m_string[1];
        
        if (DEBUG)
        {
            cerr << "Building Ukkonen Tree for ";
            cerr << "string of len: " << len << endl;
        }
        
        // Construct T1
        SuffixNode * node = createNode(1, len, m_root, NULL);
        m_root->m_children[b2i(base)] = node;
        SuffixNode * firstleaf = node;
        SuffixNode * lastleaf = node;
        
        if (DEBUG)
        {
            cerr << "Phase 1 Child: ";
            node->printLabel(cerr, m_string) << endl;
        }
        
        treeint startj = 2;
        
        // phase i+1
        for (int i = 2; i <= len; i++)
        {
            DEBUG = 0;
            
            
            // Start at the last leaf created which will allow easy
            // access to the node for startj
            node = lastleaf;
            treeint nodewalk = 0;
            
            // Keep track of last internal nodes created in split so we can add suffix links
            SuffixNode * splitnode = NULL;
            
            if (!DOPHASETRICK)
            {
                startj = 2;
                node = firstleaf;
            }
            
            if (DEBUG)
            {
                char next = m_string[i];
                cerr << endl;
                cerr << i << ".0 " << "Phase " << i << " adding " << next << " starting with " << startj << endl;
                
                string beta = m_string.substr(1, i);
                cerr << i << ".1" << " Extension 1:  [implicit]" << endl;
            }
            
            for (treeint j = startj; j <= i; j++)
            {
                // Goal: Ensure S[j .. i] (beta) is in the suffix tree
                // Precondition: S[j-1 .. i] (alpha) is in the suffix tree "near" node
                //               All Internal nodes have a suffix link
                
                // Idea: 1) Remember where alpha is in the tree relative to node
                //       2) Walk up the tree w bases until we get to a node with a suffix link.
                //       3) Follow suffix link which shifts the path from S[j-1..i] to S[j..i]
                //       4) Walk down tree in new location ensuring S[i-w .. i] is in tree
                
                // Notes: 1) All internal nodes have a suffix link by next extension
                //        2) Any time we walk up to root, have to check S[j..i]
                //        3) Suffix [1..i] is always present so start extension j with 2
                
                treeint betapos = i; // The first position in string we need to check in tree
                
                if (DEBUG)
                {
                    cerr << endl;
                    string beta = m_string.substr(j, i-j+1);
                    cerr << i << "." << j << " Phase " << i << " Extension " << j << " bp:" << betapos << endl;
                    
                    cerr << i << "." << j << "  Walking up from n:";
                    cerr << " nw: " << nodewalk << endl;
                }
                
                if (node == m_root)
                {
                    // If we are at root, we have to check the full string s[j..i] anyways
                }
                else
                {
                    if (nodewalk)
                    {
                        // partially walked down node->child, but didn't switch to child
                        // Match at i=6 on left... nodewalk=2, at 5 after suffix link
                        // 5 = i-2+1
                        //                 o ----- o
                        //               5 A       A 5  <-
                        //            -> 6 T       T 6
                        
                        betapos -= nodewalk-1;
                        
                        if (DEBUG)
                        {
                            cerr << i << "." << j << "   Adjusted nw: " << nodewalk << endl;
                        }
                    }
                    else
                    {
                        // Exactly at a node or leaf.
                        // Walk up to parent, subtracting length of that edge
                        treeint len = node->len(i);
                        betapos -= len-1;
                        node = node->m_parent;
                        
                        if (DEBUG)
                        {
                            cerr << i << "." << j << "   Adjusted len: " << len << endl;
                        }
                    }
                    
                    if (DEBUG)
                    {
                        cerr << i << "." << j << "   parent bp: " << betapos <<  " n:";
                        cerr<<endl;
                    }
                    
                    if (node->m_suffixTable[0] == NULL)
                    {
                        // Subtract entire edge length
                        betapos -= node->len(i);
                        node = node->m_parent;
                        
                        if (DEBUG)
                        {
                            cerr << i << "." << j << "   grandparent bp: " << betapos << " n:";
                            cerr << endl;
                        }
                        
                        if (node->m_suffixTable[0] == NULL)
                        {
                            cerr << "Missing suffix link!!! ";
                            exit(1);
                        }
                    }
                }
                
                // jump across suffix link
                node = node->m_suffixTable[0];
                if (node == m_root) { betapos = j; } // have to check full string
                
                if (DEBUG)
                {
                    cerr << i << "." << j << "  Starting to walk down from bp: " << betapos << " to " << i << " n:";
                    cerr << endl;
                }
                
                if (FORCEROOT && node != m_root)
                {
                    node = m_root;
                    betapos = j;
                    
                    if (DEBUG)
                    {
                        cerr << i << "." << j << " AtRoot bp: " << betapos << endl;
                    }
                }
                
                bool done = false;
                startj = j+1; // assume this extension should be skipped in the next phase
                
                while ((betapos <= i) && !done)
                {
                    char base = m_string[betapos];
                    char b = b2i(base);
                    SuffixNode * child = node->m_children[b];
                    
                    if (DEBUG)
                    {
                        cerr << i << "." << j << "  node betapos: " << betapos << "[" << base << "] n:";
                        cerr << endl;
                    }
                    
                    if (!child)
                    {
                        if (splitnode && betapos == splitnode->m_start)
                        {
                            if (DEBUG)
                            {
                                cerr << i << "." << j << "   Add SL1: ";
                            }
                            
                            splitnode->m_parent->m_suffixTable[0] = node;
                            splitnode = NULL;
                        }
                        
                        SuffixNode * newnode = createNode(betapos, len, node, NULL);
                        node->m_children[b] = newnode;
                        lastleaf = newnode;
                        
                        if (DEBUG)
                        {
                            cerr << i << "." << j << "   New Node: ";
                        }
                        
                        node = newnode;
                        
                        // This is the first base that differs, but the edgelength to
                        // i may be longer. Therefore set nodewalk to 0, so the entire
                        // edge is subtracted.
                        nodewalk = 0;
                        done = true;
                        break;
                    }
                    else
                    {
                        treeint nodepos = child->m_start;
                        nodewalk = 0;
                        
                        char nodebase = m_string[nodepos];
                        
                        if (nodebase != base)
                        {
                            char nb = m_string[nodepos];
                            cerr << "ERROR: first base on edge doesn't match edge label" << endl;
                            cerr << "       nb: " << nb << " base: " << base << endl;
                            exit(1);
                        }
                        
                        // By construction, the string from j-1 to betapos to i-1
                        // must already by present in the suffix tree
                        // Therefore, we can skip checking every character, and zoom
                        // to exactly the right character, possibly skipping the entire edge
                        
                        if (DOJUMP)
                        {
                            treeint mustmatch = i-1 - betapos + 1;
                            treeint childlen = child->len(i);
                            
                            if (mustmatch >= childlen)
                            {
                                betapos += childlen;
                                nodepos += childlen;
                                
                                skippedbases += childlen;
                                
                                if (DEBUG)
                                {
                                    cerr << i << "." << j << "   Edge Jump by: " << childlen << " new bp: " << betapos << " np: " << nodepos << endl;
                                }
                                
                                if (nodepos != child->m_end+1)
                                {
                                    cerr << "ERROR: jump should have skipped entire edge, but didn't!" << endl;
                                    exit(1);
                                }
                            }
                            else if (mustmatch)
                            {
                                betapos += mustmatch;
                                nodepos += mustmatch;
                                nodewalk += mustmatch;
                                
                                skippedbases += mustmatch;
                                
                                if (DEBUG)
                                {
                                    cerr << i << "." << j << "   Partial Jump by: " << mustmatch << " new bp: " << betapos << " np: " << nodepos << endl;
                                }
                                
                                if (VERIFY)
                                {
                                    if (m_string[betapos-1] != m_string[nodepos-1])
                                    {
                                        cerr << "ERROR: jump should have matched at least the mustmatch-1 characters" << endl;
                                        cerr << "s[bp-1]: " << m_string[betapos-1] << " s[np-1]: " << m_string[nodepos-1] << endl;
                                        exit(1);
                                    }
                                }
                            }
                        }
                        
                        while (nodepos <= child->m_end && betapos <= i)
                        {
                            nodebase = m_string[nodepos];
                            
                            if (VERBOSE)
                            {
                                cerr << i << "." << j << "   child bp: " << betapos << "[" << m_string[betapos]
                                << "] nb [" << nodebase << "]" << endl;
                            }
                            
                            if (m_string[betapos] == nodebase)
                            {
                                if (splitnode && betapos == splitnode->m_start)
                                {
                                    if (DEBUG)
                                    {
                                        cerr << i << "." << j << "   Add SL2: ";
                                    }
                                    
                                    splitnode->m_parent->m_suffixTable[0] = node;
                                    splitnode = NULL;
                                }
                                
                                nodepos++; betapos++; nodewalk++;
                                
                                if (betapos == i+1)
                                {
                                    if (DEBUG)
                                    {
                                        cerr << i << "." << j << "    Internal edge match nw: " << nodewalk << endl;
                                    }
                                    
                                    if ((nodewalk == child->len(i)) && (child->m_end == len))
                                    {
                                        // we walked the whole edge to leaf, implicit rule I extension
                                        if (DEBUG)
                                        {
                                            cerr << i << "." << j << "    Leaf Node, Implicit Rule I Extension" << endl;
                                        }
                                    }
                                    else
                                    {
                                        // "Real" rule III implicit extension
                                        
                                        // The j-1 extension was the last explicit extension in this round
                                        // Start the next round at the last explicit extension
                                        if (DOPHASETRICK)
                                        {
                                            startj = j;
                                            
                                            treeint skip = startj - 2;
                                            
                                            if (DEBUG)
                                            {
                                                cerr << i << "." << j << "    Implicit Extension... start next phase at " << startj << ", saved " << skip << endl;
                                            }
                                            
                                            skippedextensions += skip;
                                        }
                                        
                                        if (DOINTERNALSKIP)
                                        {
                                            // Since we hit an internal match on a non-leaf, we know every other
                                            // extension in this phase will also hit an internal match.
                                            
                                            // Have to be careful since leafs get the full string immediately, but
                                            // they really have a Rule 1 extension
                                            
                                            treeint skip = i-j;
                                            
                                            if (DEBUG)
                                            {
                                                cerr << i << "." << j << "    Implicit Extension... skipping rest of phase, saved " << skip << endl;
                                            }
                                            
                                            skippedextensions += skip;
                                            j = i+1;
                                        }
                                    }
                                    
                                    done = true;
                                }
                            }
                            else
                            {
                                if (DEBUG)
                                {
                                    cerr << i << "." << j << "   Spliting ";
                                }
                                
                                // Split is a copy of the child with the end shifted
                                // Then adjust start of child
                                SuffixNode * split = createNode(child->m_start, nodepos-1, node, NULL);
                                
                                split->m_children[b2i(nodebase)] = child;
                                child->m_start = nodepos;
                                child->m_parent = split;
                                
                                if (DEBUG)
                                {
                                    cerr << " => ";
                                    //split->printLabel(cerr, m_string) << " + ";
                                    //child->printLabel(cerr, m_string) << endl;
                                }
                                
                                node->m_children[b] = split;
                                node = split;
                                
                                if (splitnode && betapos == splitnode->m_start)
                                {
                                    if (DEBUG)
                                    {
                                        cerr << i << "." << j << "   Add SL3: ";
                                        //splitnode->m_parent->printLabel(cerr, m_string) << " sl-> ";
                                        //node->printLabel(cerr, m_string) << endl;
                                    }
                                    
                                    splitnode->m_parent->m_suffixTable[0] = split;
                                    splitnode = NULL;
                                }
                                
                                // Now create the new node
                                base = m_string[betapos];
                                b = b2i(base);
                                SuffixNode * newnode = createNode(betapos, len, split, NULL);
                                lastleaf = newnode;
                                
                                split->m_children[b] = newnode;
                                splitnode = newnode;
                                
                                node = newnode;
                                
                                if (DEBUG)
                                {
                                    cerr << i << "." << j << "   Split New Node: ";
                                    //newnode->printLabel(cerr, m_string)
                                    cerr << endl;
                                }
                                
                                // This is the first base that differs, but the edgelength to
                                // i may be longer. Therefore set nodewalk to 0, so the entire
                                // edge is subtracted.
                                nodewalk = 0;
                                done = true;
                                break;
                            }
                        }
                    }
                    
                    if (!done) { node = child; }
                }
            }
            
            /*
             if (VERIFY)
             {
             verifySuffixLinks();
             }*/
        }
    }
};


SuffixTree * buildUkkonenSuffixTree(const string & s)
{
    SuffixTree * tree = new SuffixTree(s);
    tree->buildUkkonen();
    
    return tree;
}



void printHelp()
{
    /*cout << "splitMEM [options]"                                         << endl
     << "  -h                             help"                                           << endl
     //<< "  -dot                       Print tree in dot form"                         << endl
     //<< "  -sort                      Print tree in sorted text form"                 << endl
     << "  -file <file>                   Load sequence from file"                        << endl
     << "  -mem  <kmer length>            Indicates single kmer length to use; smallest MEM that is marked in the suffix tree"                                   << endl
     << "  -manyMEMs <kmerLenFile>        Create graph for several kmer legnths listed in file, each on separate line"<<endl
     << "  -cdg <output file name>        Name of file compressed de Bruijn graph will be written to.  Prefix in manyMEMs mode."<<endl
     << "  -multiFa                       Input file is a mult fasta file "<<endl
     << "  -printGenome <outputFaFile>    Prints the fasta file as it is represented in memory"<<endl;
     */
    cout << "splitMEM [options] -file <fasta.file>"                     << endl
    << "  -h               help"                                   << endl
    << "  -file <file>     Load sequence from file"                << endl
    << "  -mem <len>       Locate MEMs at least this long "        << endl
    << "  -manyMEMs <file> File of minimum MEM lengths"            << endl
    << "  -cdg <out>       Filename of compressed de Bruijn graph" << endl
    << "  -multiFa         Indicates the input file is a multifasta file for pan-genome analysis" << endl;
    
    exit(0);
}



int OPT_DisplaySeq = 0;
int OPT_DisplayStarts = 1;
int OPT_SeqToDisplay = 8;
int OPT_DisplayStats = 0;
int OPT_PrintGraphs = 0;


// class for storing a mer and its neighbors
class MerVertex_t
{
    static treeintLarge NODECOUNT; // give all the nodes a unique id
public:
    MerVertex_t(treeint len)
    : node_m(NODECOUNT++) {
        length_m = len;
    }
    
    MerVertex_t(treeint startpos, treeint len)
    : node_m(NODECOUNT++) {
        starts_m.push_back(startpos);
        length_m = len;
    }
    
    MerVertex_t(treeint startpos, treeint len, bool noNodecount)
    : node_m(-1) {
        starts_m.push_back(startpos);
        length_m = len;
    }
    
    long   node_m; //unique node id
    vector<MerVertex_t *> successor_m;
    vector<MerVertex_t *> predecessor_m;
    vector<treeint> starts_m;
    treeint length_m;
    
    treeint addStartPos(treeint startpos)
    {
        starts_m.push_back(startpos);
        return starts_m.size()-1;
    }
    
    static treeintLarge getNodeCount()
    {
        return NODECOUNT;
    }
    
    treeint getNumEdges()
    {
        
        treeint numEdges = 0;
        numEdges = successor_m.size();
        return numEdges;
    }
    
    
    // Return the subsequence stored in the node
    string str(const string & seq)
    {
        string r = seq.substr(starts_m[0]-1, length_m);
        return r;
    }
    
    static void resetNodeCount()
    {
        NODECOUNT = 0;
    }
};

treeintLarge MerVertex_t::NODECOUNT=0;

//data type to wrap together a node pointer and which startpos index it represents
//this way the set of NodePos_t objects is searchable by startpos.  a vector of MerVertex_t is not searchable by startpos
class NodePos_t
{
public:
    MerVertex_t* nodePtr_m;
    treeint startPosIdx_m;
    
    NodePos_t(MerVertex_t* nptr, treeintLarge idx = 0){
        nodePtr_m = nptr;
        startPosIdx_m = idx;
    }
    
    void changeNodePtr(MerVertex_t* nptr)
    {
        nodePtr_m = nptr;
    }
    
    void changeStartIdx(int newIdx)
    {
        startPosIdx_m = newIdx;
    }
    
    bool operator ==(const NodePos_t * aNode)
    {
        if(getNodeBegin() == aNode->getNodeBegin())
            return true;
        else
            return false;
    }
    
    bool operator >(const NodePos_t * aNode)
    {
        if(getNodeBegin() > aNode->getNodeBegin())
            return true;
        else
            return false;
    }
    
    //get endPos based on startPosIdx of this nodePos
    treeint getNodeEnd() const
    {
        return getNodeBegin()+nodePtr_m->length_m-1;
    }
    
    //get startPos based on startPosIdx of this nodePos
    treeint getNodeBegin() const
    {
        return nodePtr_m->starts_m[startPosIdx_m];
    }
    
};


bool SortLens(const treeint & a, const treeint & b)
{
    return b < a;
}

bool SortStarts(const treeint & a, const treeint & b)
{
    return a < b;
}

bool SortNodePos(NodePos_t* a, NodePos_t* b)
{
    return a->nodePtr_m->starts_m[a->startPosIdx_m] < b->nodePtr_m->starts_m[b->startPosIdx_m];
}

bool SortMerVertex(MerVertex_t* a, MerVertex_t* b)
{
    return a->node_m < b->node_m;
}

typedef vector<NodePos_t *> NodeTable_t;
typedef vector<MerVertex_t *> RepeatNodeTable_t;
typedef vector<MerVertex_t *> UniqueNodeTable_t;


// Class for storing the entire compressed de Bruijn graph
class deBruijnGraph_t
{
public:
    
    deBruijnGraph_t(const string & tag,
                    const string & seq)
    : tag_m(tag), seq_m(seq) {}
    
    const string & tag_m;
    const string & seq_m;
    
    RepeatNodeTable_t repeatNodes_m;
    UniqueNodeTable_t uniqueNodes_m;
    MerVertex_t * firstNode;  //for traversal
    
    ~deBruijnGraph_t()
    {
        for(size_t i=0; i<repeatNodes_m.size(); i++)
            delete repeatNodes_m[i];
        for(size_t i=0; i<uniqueNodes_m.size(); i++)
            delete uniqueNodes_m[i];
        
    }
    void createParentChildRelationship(MerVertex_t * parentNode, MerVertex_t * childNode)
    {
        parentNode->successor_m.push_back(childNode);
        childNode->predecessor_m.push_back(parentNode);
    }
    
    
    //create node_pos for each start in the repeat nodes, sort the set of node_pos objects in order of increasing startpos
    void sortRepeatNodeStartsToNodeTable(NodeTable_t& nodePosSorted)
    {
        MerVertex_t * newNode;
        NodePos_t * newNodePos;
        
        //loop through each start pos in each repeat node and create an entry in nodes_m for each
        for(treeint i=0; i<repeatNodes_m.size(); i++)
        {
            for(treeint j=0; j<repeatNodes_m[i]->starts_m.size(); j++)
            {
                newNodePos = new NodePos_t(repeatNodes_m[i], j);
                nodePosSorted.push_back(newNodePos);
            }
        }
        
        //sort starts in nodes_m
        sort(nodePosSorted.begin(), nodePosSorted.end(), SortNodePos);
    }
    
    //this fucntion creates edges and unique nodes along the way
    //can loop through nodes, search for successor for each start and make edges
    //but this function loops through sorted starts in nodes_m and make edges without searching.
    //either way, when the successor start is too far away, create a linking uniqueNode
    //when this function is called, it assumes nodes_m is already sorted by start position in increasing order
    void createEdgesAndUniqueNodes()
    {
        treeint thisStart, thisLength, nextStart;
        MerVertex_t* successor;
        MerVertex_t * lastNode = NULL;
        treeint i = 0;
        NodeTable_t nodePosSorted;  // each node can have several values for startpos so need to know which pos each link refers to
        
        cerr<<" there are "<<repeatNodes_m.size()<<" repeat nodes ";
        
        sortRepeatNodeStartsToNodeTable(nodePosSorted);
        cerr<<" rep. "<<nodePosSorted.size()<<" starts"<<endl;
        cerr<<" repeat nodes have been sorted "<<endl;
        
        //make sure there is a start node for the beginning of the sequence, if the first start pos is not 1
        if(nodePosSorted.size()>0 && nodePosSorted[0]->getNodeBegin()>1) // first position is considered 1
        {
            successor = new MerVertex_t(1, (nodePosSorted[0]->getNodeBegin())+Kmer_Len-1-1);
            //begin at 1 not 0 so that initial s is excluded
            lastNode = createUniqueNode(successor);
        }
        
        
        if(nodePosSorted.size()>0)
        {
            for(i = 0; i<nodePosSorted.size()-1; i++)
            {
                
                if(lastNode != NULL)
                {
                    createParentChildRelationship(lastNode, nodePosSorted[i]->nodePtr_m);
                }
                lastNode = NULL;
                
                //find successor for this start position
                thisStart = nodePosSorted[i]->getNodeBegin();
                thisLength = nodePosSorted[i]->nodePtr_m->length_m;
                nextStart = thisStart + thisLength - Kmer_Len + 1;
                
                
                //if this is the last nodePos, see if it goes until the end of the sequence or we need a linking node
                //if it is the next start in the sorted list, the next nodePos is the successor and create edge
                if(nodePosSorted[i+1]->getNodeBegin() == nextStart)
                {
                    successor = nodePosSorted[i+1]->nodePtr_m;
                }
                //if it is NOT the next start in the sorted list, create uniqueNode to link, since it must be in between
                else
                {
                    successor = new MerVertex_t(nextStart, (nodePosSorted[i+1]->getNodeBegin())-nextStart+Kmer_Len-1);
                    lastNode = createUniqueNode(successor);
                }
                
                //create edge to successor
                createParentChildRelationship(nodePosSorted[i]->nodePtr_m, successor);
            }
        }
        
        if(i==0) //no repeat nodes, degenerate case of 1 node in graph
        {
            successor = new MerVertex_t(1, seq_m.length()-2);
            createUniqueNode(successor);
            
        }
        else
        {
            if(lastNode != NULL)
            {
                createParentChildRelationship(lastNode, nodePosSorted[i]->nodePtr_m);
            }
            lastNode = NULL;
            
            if(nodePosSorted[i]->getNodeEnd() != seq_m.length())
            {
                thisStart = nodePosSorted[i]->getNodeBegin();
                thisLength = nodePosSorted[i]->nodePtr_m->length_m;
                nextStart = thisStart + thisLength - Kmer_Len + 1;
                
                //create node that goes from the end of this one until the end of the genome
                successor = new MerVertex_t(nextStart, seq_m.length() - nextStart);
                createUniqueNode(successor);
            }
        }
        
        if(nodePosSorted.size() > 1)
            //create edge to successor
            createParentChildRelationship(nodePosSorted[i]->nodePtr_m, successor);
        
        for(size_t i=0; i<nodePosSorted.size(); i++)
            delete nodePosSorted[i];
        
        
    }
    
    
    void createRepeatNodesFromSuffixTree(SuffixTree* stree)
    {
        //recursively perform DFS on suffix tree and break down MEMs and their subtrees
        createRepeatNodesFromMEM(stree->m_root, stree);
    }
    
    
    void createRepeatNode(SuffixNode* memNode, int offset, int length, SuffixTree * ST)
    {
        MerVertex_t* newNode;
        
        newNode = new MerVertex_t(ST->m_suffixArray[memNode->m_SA_start] + offset, length); //for first start, length
        for(int i = 1; i < memNode->m_SA_end - memNode->m_SA_start + 1; i++)
        {
            newNode->addStartPos(ST->m_suffixArray[memNode->m_SA_start+i] + offset);
        }
        repeatNodes_m.push_back(newNode);
        
    }
    
    MerVertex_t* createUniqueNode(MerVertex_t* newNode)
    {
        uniqueNodes_m.push_back(newNode);
        size_t pos = uniqueNodes_m.size() - 1;
        return uniqueNodes_m[pos];
    }
    
    //create repeat nodes from suffix tree traversal
    //break down MEMs to repeat nodes - find kmers shared among MEMs
    void createRepeatNodesFromMEM(SuffixNode * origMemNode, SuffixTree * ST)//, int extendLeft)
    {
        //speedup: if !MEM and long enough for Kmer_Len, none of descendandts are MEMs so don't bother
        if(!origMemNode->m_MEM && origMemNode->m_strdepth >= Kmer_Len )
            return;
        
        int children = 0;
        //go through marked nodes (DFS with recursive calls), can use BFS
        for (int b = 0; b < basecount; b++)
        {
            if(origMemNode->m_children[b])
            {
                children++;
                
                createRepeatNodesFromMEM(origMemNode->m_children[b], ST);
            }
        }
        
        
        //todo: second clause is redundant
        if(!origMemNode->m_MEM || origMemNode->m_strdepth < Kmer_Len)
            return;
        
        int origOffset = 0;
        int extendLeft = 0;
        
        //always extend MEM node to root so that it's the entire MEM.
        //if there is a prefix that is another MEM, it will be removed like any substring MEM (that is not a proper suffix)
        if(origMemNode->m_parent != ST->m_root)
        {
            extendLeft = origMemNode->m_strdepth - origMemNode->len();
        }
        
        
        //don't change the memNode. instead, navigate to appropriate location in suffix tree by following suffix links for extendLeft chars away from beginning of path from root
        //copy relevant information so as not to modify the original node in the suffix tree when extending left
        //this extended node represents the entire MEM from the root, with appropriate starts, but it maintains the original LMA
        //the parent is no longer the parent. might be more accurate to make the root the parent
        
        //use copy constructor
        SuffixNode memNodeExt(*origMemNode);
        memNodeExt.m_start -= extendLeft;
        
        SuffixNode * memNode = &memNodeExt;
        
        memNode->m_strdepth = memNode->len();
        
        bool needLastNode = true;
        SuffixNode* snode = memNode;
        SuffixNode* LMAnode;
        
        int skippedChars = 0; //count number of chars that get removed when follow suffix links - so that we an create a repeat node for them
        int offset = 0; //how far are from beginning of node the next new node should begin, useful for calculating starts for nodes that begin in middle of MEM that are not another MEM
        int rOffset = 0; //how far from right end of MEM the last repeatNode should end
        
        int slinkToTraverse = 0; //increments offset - takes care of prefix MEM (could be  a prefix to a suffix of the MEM)
        int slinkTraversing = 0; //how much of slinkToTraverse is being used with one suffix skip
        int slinkIndex = 0; //what index in suffix link table is being used for slinkTraversing
        bool endOuterLoop = false;
        
        //begin with modified node. then use suffix links to navigate to successively smaller suffixes of the entire MEM
        snode = memNode;
        
        while(snode->m_strdepth >= Kmer_Len)
        {
            //we don't  have to consider if the parent is marked - once a parent is marked, all internal node descendants are also marked.  we can simply see if the node itself is marked
            //if this node is marked - remove beginning of node that is MEM and continue with rest of it. at the same time, make a repeat node to represent whatever has been skipped with suffix links (need to know length and the starts)
            //start positions are already in this node's set of starts so nothing to change about it, just need to carve it out of longer node
            //if snode is not marked, check if parent is marked (LMA query)
            
            if(skippedChars == 0 && offset == 0)  //first LMA query for this MEM
            {
                LMAnode = snode->m_parent->LMA();
            }
            else
            {
                LMAnode = snode->LMA();
            }
            
            if(LMAnode != NULL)
            {
                
                if(LMAnode->m_strdepth + offset + skippedChars == memNode->m_strdepth)
                {
                    rOffset = LMAnode->m_strdepth - Kmer_Len + 1;
                    break; //end loop
                }
                else
                {
                    if(skippedChars >= 1) //link is long enough to be repeat node
                    {
                        createRepeatNode(memNode, offset+origOffset, skippedChars+Kmer_Len-1, ST);
                        
                    }
                    rOffset = 0;
                    offset += skippedChars; //skip newNode
                    slinkToTraverse = LMAnode->m_strdepth-Kmer_Len+1;
                    
                    
                    skippedChars = 0;
                }
            }
            
            //use m_suffixTable appropriate entry to traverse slinkToTraverse in floor(log length) time
            //follow suffix link(s)
            
            endOuterLoop = false;
            
            if(slinkToTraverse > 0)
            {
                //don't want to check for skipped LMA nodes for first snode so handle separately, before loop
                offset++;
                slinkToTraverse--;
                snode = snode->m_suffixTable[0];
                
                
                //already did extra slinkToTraverse-- at beginning of loop
                if(snode->m_MEM)
                {
                    needLastNode = false;
                    break; //end outer loop for MEM - will be processing another MEM for no reason (like in tandem repeat)
                }
                
                while(slinkToTraverse>0)
                {
                    slinkIndex = (int)(floor(log(slinkToTraverse) / log(2)));
                    slinkTraversing = (int)(pow(2.0, slinkIndex));
                    
                    
                    //todo: can remove this if statement and go to next (nested one)
                    if(snode->m_LMAprox_nodeTable[slinkIndex] != NULL)
                    {
                        if(snode->m_LMAproximityTable[slinkIndex] == 0) //at end of MEM that we are currently working with
                        {
                            
                            rOffset = snode->m_LMAprox_nodeTable[slinkIndex]->m_strdepth - Kmer_Len + 1;
                            endOuterLoop = true;
                            needLastNode = false;
                        }
                        else if(snode->m_LMAproximityTable[slinkIndex] < snode->m_strdepth-slinkToTraverse-Kmer_Len+1)
                        {
                            slinkToTraverse += (snode->m_strdepth-slinkToTraverse - snode->m_LMAproximityTable[slinkIndex] - Kmer_Len+1);
                            
                        }
                        
                    }
                    
                    
                    //adjust offset in bulk
                    offset += slinkTraversing;
                    snode = snode->m_suffixTable[slinkIndex]; //traverse suffix skips so that algorithm is O(n log n) time instead of O(n^2)
                    
                    slinkToTraverse -= slinkTraversing;
                    
                    if(endOuterLoop)
                        goto afterMemLoop;
                    
                }
                
            }
            else
            {
                snode = snode->m_suffixTable[0];             //follow suffix link
                skippedChars++;
                //cerr<<" inc skippedChars"<<endl;
            }
            
        }
        
    afterMemLoop:
        //finish processing node and make repeat node from wherever until the end (there are no more MEMs), unless it's a MEM
        if(needLastNode && memNode->len()-offset-rOffset>=Kmer_Len)
        {
            createRepeatNode(memNode, offset+origOffset, memNode->len()-offset-rOffset, ST);
        }
        
    }
    
    
    //takes suffix tree with MEMs marked and creates compressed de Bruijn graph
    void construct(SuffixTree* stree)
    {
        cerr<<" constructing CDG"<<endl;
        createRepeatNodesFromSuffixTree(stree);
        cerr<<" created repeat nodes"<<endl;
        
        createEdgesAndUniqueNodes();
        cerr<<" created unique nodes and all edges"<<endl;
    }
    
    // Returns number of nodes in the graph
    treeintLarge nodeCount()
    {
        return MerVertex_t::getNodeCount();
        
    }
    
    void printStats()
    {
        treeint numNodes = nodeCount();
        treeint numEdges = 0;
        treeint totalspan = 0;
        vector<treeint> lengths;
        treeint max = 0;
        treeint l;
        
        for(treeint i=0; i<repeatNodes_m.size(); i++)
        {
            numEdges += repeatNodes_m[i]->getNumEdges();
            l = repeatNodes_m[i]->length_m;
            
            totalspan += l;
            lengths.push_back(l);
            
            if (l > max) { max = l; }
        }
        
        for(treeint i=0; i<uniqueNodes_m.size(); i++)
        {
            numEdges += uniqueNodes_m[i]->getNumEdges();
            l = uniqueNodes_m[i]->length_m;
            
            totalspan += l;
            lengths.push_back(l);
            
            if (l > max) { max = l; }
        }
        
        sort(lengths.begin(), lengths.end(), SortLens);
        
        treeint target = totalspan/2;
        treeint sum = 0;
        treeint n50cnt = 0;
        
        for (treeint i = 0; i < lengths.size(); i++)
        {
            sum += lengths[i];
            
            if (sum >= target)
            {
                n50cnt = i;
                break;
            }
        }
        
        cerr << "n="          << numNodes
        //<< " or n="      << numNodesLong
        << " m="         << numEdges
        << " totalspan=" << totalspan
        << " max="       << max
        << " mean="      << ((double) totalspan)/ ((double) numNodes)
        << " n50="       << lengths[n50cnt]
        << " n50cnt="    << n50cnt
        << endl;
        
        
    }
    
    
    
    // Print graph in DOT format
    //starts are only sorted when lineNum=0, for final graph
    void print()
    {
        streambuf *coutbuf = cout.rdbuf(); //save old buf
        
        
        stringstream sstm;
        
        sstm << CDG_Filename;
        string dotFileName = sstm.str();
        ofstream out(dotFileName.c_str());
        //redirecting all ouptut since suffix tree goes to cout
        
        cout.rdbuf(out.rdbuf());
        
        cerr<<" redirected output to "<<dotFileName<<endl;
        
        int numNodes = 0;
        
        cout << "digraph G {" << endl;
        
        
        for(treeint i = 0; i<repeatNodes_m.size(); i++)
        {
            
            if(OPT_DisplayStarts && OPT_DisplaySeq)
            {
                string s = seq_m.substr(repeatNodes_m[i]->starts_m[0], repeatNodes_m[i]->length_m);
                int slen = repeatNodes_m[i]->length_m;
                
                if (slen > OPT_SeqToDisplay)
                {
                    string q = s.substr(0,OPT_SeqToDisplay/2);
                    q += "...";
                    q += s.substr(s.length()-OPT_SeqToDisplay/2, OPT_SeqToDisplay/2);
                    s = q;
                }
                
                cout << "  " << repeatNodes_m[i]->node_m << " [label=\"" << s << "\\n";
                //can sort but then invalidates nodePos_t objects, which is fine if the program will end after this function call (it won't prevent destructor from running)
                sort(repeatNodes_m[i]->starts_m.begin(), repeatNodes_m[i]->starts_m.end(), SortStarts);
                
                for (int j = 0; j < repeatNodes_m[i]->starts_m.size(); j++)
                {
                    if (j > 0) { cout << ","; }
                    cout << repeatNodes_m[i]->starts_m[j]; // uses 1-based coordinates
                }
                
                cout << ":" << repeatNodes_m[i]->length_m << "\"]" << endl;
                
            }
            else if (OPT_DisplayStarts)
            {
                cout << "  " << repeatNodes_m[i]->node_m<< " [label=\"";
                //can sort but then invalidates nodePos_t objects, which is fine if the program will end after this function call (it won't prevent destructor from running)
                sort(repeatNodes_m[i]->starts_m.begin(), repeatNodes_m[i]->starts_m.end(), SortStarts);
                
                for (int j = 0; j < repeatNodes_m[i]->starts_m.size(); j++)
                {
                    if (j > 0) { cout << ","; }
                    cout << repeatNodes_m[i]->starts_m[j]; // uses 1-based coordinates
                }
                
                cout << ":" << repeatNodes_m[i]->length_m << "\"]" << endl;
            }
            else if (OPT_DisplaySeq)
            {
                string s = seq_m.substr(repeatNodes_m[i]->starts_m[0]-1, repeatNodes_m[i]->length_m);
                int slen = repeatNodes_m[i]->length_m;
                
                if (slen > OPT_SeqToDisplay)
                {
                    string q = s.substr(0,OPT_SeqToDisplay/2);
                    q += "...";
                    q += s.substr(s.length()-OPT_SeqToDisplay/2, OPT_SeqToDisplay/2);
                    s = q;
                }
                cout << "  " << repeatNodes_m[i]->node_m << " [label=\"" << s <<  "\"]" << endl;
            }
            else
            {
                cout << "  " << repeatNodes_m[i]->node_m << " [label=\"" << repeatNodes_m[i]->length_m << "\"]" << endl;
            }
            
            numNodes++;
            
            for (int j = 0; j < repeatNodes_m[i]->successor_m.size(); j++)
            {
                cout << "  " << repeatNodes_m[i]->node_m << " -> " << repeatNodes_m[i]->successor_m[j]->node_m << endl;
            }
            
        }
        
        for(treeint i = 0; i<uniqueNodes_m.size(); i++)
        {
            
            if(OPT_DisplayStarts && OPT_DisplaySeq)
            {
                string s = seq_m.substr(uniqueNodes_m[i]->starts_m[0], uniqueNodes_m[i]->length_m);
                int slen = uniqueNodes_m[i]->length_m;
                if (slen > OPT_SeqToDisplay)
                {
                    string q = s.substr(0,OPT_SeqToDisplay/2);
                    q += "...";
                    q += s.substr(s.length()-OPT_SeqToDisplay/2, OPT_SeqToDisplay/2);
                    s = q;
                }
                
                cout << "  " << uniqueNodes_m[i]->node_m << " [label=\"" << s << "\\n";//<<  "\"]" << endl;
                //can sort but then invalidates nodePos_t objects, which is fine if the program will end after this function call (it won't prevent destructor from running)
                
                sort(uniqueNodes_m[i]->starts_m.begin(), uniqueNodes_m[i]->starts_m.end(), SortStarts);
                
                for (int j = 0; j < uniqueNodes_m[i]->starts_m.size(); j++)
                {
                    if (j > 0) { cout << ","; }
                    cout << uniqueNodes_m[i]->starts_m[j]; // uses 1-based coordinates
                }
                
                cout << ":" << uniqueNodes_m[i]->length_m << "\"]" << endl;
                
            }
            else if (OPT_DisplayStarts)
            {
                cout << "  " << uniqueNodes_m[i]->node_m<< " [label=\"";
                //can sort but then invalidates nodePos_t objects, which is fine if the program will end after this function call (it won't prevent destructor from running)
                
                sort(uniqueNodes_m[i]->starts_m.begin(), uniqueNodes_m[i]->starts_m.end(), SortStarts);
                
                for (int j = 0; j < uniqueNodes_m[i]->starts_m.size(); j++)
                {
                    if (j > 0) { cout << ","; }
                    cout << uniqueNodes_m[i]->starts_m[j]; // uses 1-based coordinates
                }
                
                cout << ":" << uniqueNodes_m[i]->length_m << "\"]" << endl;
            }
            else if (OPT_DisplaySeq)
            {
                string s = seq_m.substr(uniqueNodes_m[i]->starts_m[0]-1, uniqueNodes_m[i]->length_m);
                int slen = uniqueNodes_m[i]->length_m;
                
                if (slen > OPT_SeqToDisplay)
                {
                    string q = s.substr(0,OPT_SeqToDisplay/2);
                    q += "...";
                    q += s.substr(s.length()-OPT_SeqToDisplay/2, OPT_SeqToDisplay/2);
                    s = q;
                }
                cout << "  " << uniqueNodes_m[i]->node_m << " [label=\"" << s <<  "\"]" << endl;
            }
            else
            {
                cout << "  " << uniqueNodes_m[i]->node_m << " [label=\"" << uniqueNodes_m[i]->length_m << "\"]" << endl;
            }
            
            numNodes++;
            
            for (int j = 0; j < uniqueNodes_m[i]->successor_m.size(); j++)
            {
                cout << "  " << uniqueNodes_m[i]->node_m << " -> " << uniqueNodes_m[i]->successor_m[j]->node_m << endl;
            }
            
        }
        
        cout << "}" << endl;
        
        cout.rdbuf(coutbuf); //reset to standard output again
        
    }
    
};


// For now, just create the deBruijn Graph, and print it out
//void ComputeComplexity(const string & tag, const string & seq)
//todo: SM removed tag variable from parameter list, should get it from fasta file in main and pass in
void ComputeComplexity(const string & seq, SuffixTree * tree)
{
    treeint  i, j, n;
    
    n = seq . length ();
    
    if  (n < Kmer_Len)
    {
        cerr<<" string len = "<<"  Kmer_Len = "<<Kmer_Len<<endl;
        return;
    }
    
    //todo: add tag variable back into function call
    //deBruijnGraph_t graph(tag, seq);
    deBruijnGraph_t graph("", seq);
    graph.construct(tree);
    
    
    cerr << graph.nodeCount() << " nodes."<< endl;
    
    graph.printStats();
    
    if (!OPT_DisplayStats)
    {
        graph.print();
    }
    
    return;
}
bool string_has_all_of_the_same_chars(const std::string& s) {
    //cout <<"Inside string_has_all_of...\n";
    return s.find_first_not_of(s[0]) == std::string::npos;
}
bool containsN(string str){
   return str.find("N") != std::string::npos;
}   
bool containsX(string str){
   return str.find("X") != std::string::npos; 
}

/* New printMEMnode() function */
void printMEMnode(SuffixNode * node, SuffixTree * ST, int *sLength, string seq[], string MEM[], int *index, int **matrix)
{
        //speedup: if !MEM and long enough for Kmer_Len, none of descendants are MEMs so don't bother
        if(!node->m_MEM && node->m_strdepth >= Kmer_Len )
            return;

        for (int b = 0; b < basecount; b++)
        {
            if(node->m_children[b])
            {
                printMEMnode(node->m_children[b], ST, sLength, seq, MEM, index, matrix);
            }
        }
        
        int extendLeft = 0;

        if(node->m_MEM && node->m_strdepth >= Kmer_Len)
        {
            //todo:
            //this is a MEM, need to print it out
            if(node->m_parent != ST->m_root)
            {
                 extendLeft = node->m_strdepth - node->len();
            }


            SuffixNode memNodeExt(*node);
            memNodeExt.m_start -= extendLeft;
        
            SuffixNode * memNode = &memNodeExt;
            memNode->m_strdepth = memNode->len();
        
            //first start
            //ST->m_suffixArray[memNode->m_SA_start], length
            //for each start in MEM node
           // cout<<"Detecting the MEMS......\n";
            string s;
            int t = *index;	
			int d = t;
			 int count[3]={0, 0, 0};
            for(int i = 0; i < memNode->m_SA_end - memNode->m_SA_start + 1; i++)
            {
              int y = ST->m_suffixArray[memNode->m_SA_start+i];
		      int n=0, z=0;
					
		      seqPos(sLength, y, &n, &z);
			  if (i==0) 
			  {
			    if((z - n ) == 0) s = seq[n].substr(z-n,memNode->m_strdepth);
				else s = seq[n].substr(z-n-1,memNode->m_strdepth);
			  }  
			  if(!(string_has_all_of_the_same_chars(s) || containsN(s) || containsX(s))){
		      //cout<< seq[n].substr(z-n-1,memNode->m_strdepth)<< " Seq_"<<n+1<<"\t"<<z<<"\t"<< memNode->m_strdepth << endl;
              //if(!string_has_all_of_the_same_chars(s))
			    //matrix[t][n]=z;
				//printf("n=%d\n", n);
				if(n==0) count[0]+=1;
				if(n==1) count[1]+=1;
				if(n==2) count[2]+=1;
			  }
		    //cout << ST->m_suffixArray[memNode->m_SA_start+i] << "\t" << memNode->m_strdepth << endl;
            }	
           // if(!(count[0]==0 && count[1]==0 && count[2]==0)) 			
			 // cout<<"count1:"<<count[0]<<" "<<"count2:"<<count[1]<<" "<<"count3:"<<count[2]<<"\n";
			int *s0 = (int *) malloc(sizeof(int)*count[0]+1);
			int *s1 = (int *) malloc(sizeof(int)*count[1]+1);
			int *s2 = (int *) malloc(sizeof(int)*count[2]+1);
			int c1=0, c2=0, c3=0;
            for(int i = 0; i < memNode->m_SA_end - memNode->m_SA_start + 1; i++)
            {
              int y = ST->m_suffixArray[memNode->m_SA_start+i];
		      int n=0, z=0;
					
		      seqPos(sLength, y, &n, &z);
			  if (i==0) 
			  {
				//s = seq[n].substr(z-n-1,memNode->m_strdepth);
				if((z - n ) == 0) s = seq[n].substr(z-n,memNode->m_strdepth);
				else s = seq[n].substr(z-n-1,memNode->m_strdepth);
			  }  
			   if(!(string_has_all_of_the_same_chars(s) || containsN(s) || containsX(s))){
		        //cout<< seq[n].substr(z-n-1,memNode->m_strdepth)<< " Seq_"<<n+1<<"\t"<<z<<"\t"<< memNode->m_strdepth << endl;
              //if(!string_has_all_of_the_same_chars(s))
			    //matrix[t][n]=z;
				//printf("n=%d\n", n);
				if(n==0){s0[c1]=z; c1++;}
				if(n==1){s1[c2]=z; c2++;}
				if(n==2){s2[c3]=z; c3++;}
			  }
			 // if(!(c1==0 && c2==0 && c3==0)) 			
			//cout<<"c1:"<<c1<<" "<<"c2:"<<c2<<" "<<"c3:"<<c3<<"\n";
		    //cout << ST->m_suffixArray[memNode->m_SA_start+i] << "\t" << memNode->m_strdepth << endl;
            }
			//if(!(c1==0 && c2==0 && c3==0)) 			
			//cout<<"c1:"<<c1<<" "<<"c2:"<<c2<<" "<<"c3:"<<c3<<"\n";
			int p=0, q=0, r=0;
			for(p=0; p<=c1; p++){
			  if(c1>0 && p==c1) break;
			  for(q=0; q<=c2; q++){
			    if(c2>0 && q==c2) break;
			    for(r=0; r<=c3; r++){
				  if(c3>0 && r==c3) break;
				  d++;
				 // if (d !=1) printf("d= %d\n", d);
				  if(c1 != 0 && p!=c1) matrix[d][0]=s0[p];
				  else matrix[d][0]=0;	
				  if(c2 != 0 && q!=c2) matrix[d][1]=s1[q];
				  else  matrix[d][1]=0;
				  if(c3 != 0 && r !=c3) matrix[d][2]=s2[r];
				  else matrix[d][2]=0;
				  //d++;
				} 
              }            
			}
				int aa = d - t;
		//	int t = *index;
			//MEM[t]= (char *) malloc(sizeof(char)*strlen(s)+1);
			 if(!(string_has_all_of_the_same_chars(s) || containsN(s) || containsX(s))){
			  //cout<<"\n";
			  for(p=0; p<aa; p++){
			    MEM[t]= s;
			    t++;
			    *index=t;
			  }	
            } 
     //up to here
     //use code like the following to get starts of MEM.  
     //print MEMid (can be calculated on the spot),start, and length of each MEM and print in text file (3 columns, tab separated)
     
/*            
     newNode = new MerVertex_t(ST->m_suffixArray[memNode->m_SA_start] + offset, length); //for first start, length
     for(int i = 1; i < memNode->m_SA_end - memNode->m_SA_start + 1; i++)
     {
             newNode->addStartPos(ST->m_suffixArray[memNode->m_SA_start+i] + offset);
     }
  */      

        }
        
    }   
	

/* New countMEMnode() function */
	void countMEMnode(SuffixNode * node, SuffixTree * ST, int *sLength, string seq[],int *cnt)
{
        //speedup: if !MEM and long enough for Kmer_Len, none of descendants are MEMs so don't bother
        if(!node->m_MEM && node->m_strdepth >= Kmer_Len )
            return;

        for (int b = 0; b < basecount; b++)
        {
            if(node->m_children[b])
            {
			    //*cnt++;
				//int temp = *cnt;
			    //cout <<"cnt:"<<temp<<"\n"; 
                countMEMnode(node->m_children[b], ST, sLength, seq, cnt);
            }
        }
        int extendLeft = 0;
        if(node->m_MEM && node->m_strdepth >= Kmer_Len)
        {
            //todo:
            //this is a MEM, need to print it out
            if(node->m_parent != ST->m_root)
            {
                 extendLeft = node->m_strdepth - node->len();
            }
            SuffixNode memNodeExt(*node);
            memNodeExt.m_start -= extendLeft;
        
            SuffixNode * memNode = &memNodeExt;
            memNode->m_strdepth = memNode->len();
        
            //first start
            //ST->m_suffixArray[memNode->m_SA_start], length
            //for each start in MEM node
            //cout<<"Detecting the MEMS......\n";
            //cnt++;
			int t = *cnt;
			int count[3]={0, 0, 0};
			string s;
            for(int i = 0; i < memNode->m_SA_end - memNode->m_SA_start + 1; i++)
            {
                    int y = ST->m_suffixArray[memNode->m_SA_start+i];
		            int n=0, z=0;
					int xyz = 0;
					
		            seqPos(sLength, y, &n, &z);
			        if (i==0) 
			        {
					    /*if(xyz == 0 ) {
						cout<< seq[n]<<"\n";
						cout<<*sLength<<"\t"<<y<<"\t"<<z<<"\t"<<n<<"\t"<<(z-n-1)<<"\n";
						cout<<memNode->m_strdepth<<"\n";
						string tempstr = seq[n].substr(z-n-1,memNode->m_strdepth);
						xyz++;
						}*/
						if(z-n-1 < 0){
						  //printf("testing %d %d \n", z, n);
						  //printf("%d \n", y);
						  //cout<<memNode->m_strdepth<<"\n";
						  s = seq[n].substr(z-n,memNode->m_strdepth);
						}
			        	else s = seq[n].substr(z-n-1,memNode->m_strdepth);
			       }  
				   if(!(string_has_all_of_the_same_chars(s) || containsN(s) || containsX(s))){
					  
		             //cout<< seq[n].substr(z-n-1,memNode->m_strdepth)<< " Seq_"<<n+1<<"\t"<<z<<"\t"<< memNode->m_strdepth << endl;
                     //cout << "testing \n";
				     if(n==0) count[0]+=1;
				     if(n==1) count[1]+=1;
				     if(n==2) count[2]+=1;
					 if(count[0] == 0) count[0]=1;
					 if(count[1] == 0) count[1]=1;
					 if(count[2] == 0) count[2]=1;
					 //int t1=0;
					 t=t+(count[0]*count[1]*count[2]); *cnt = t;
			       }
			       
            }
             //if(!string_has_all_of_the_same_chars(s)){t++; *cnt = t;}
				
        }
        
    }    	
void seqPos(int *s, int y, int *n, int *z)
{
  int x=0, i=0;
  while(true)
  {
    *z = y - x;
    if((y-x) <= s[i]) break;
    x = x + s[i];
    i++;
  } 
   *n = i; 
}

void printMEMs(SuffixTree * tree, int *sLength, string seq[], string MEM[], int **matrix)
{
     //recursively perform DFS on suffix tree to find all MEMs
	 int index=0;
     printMEMnode(tree->m_root, tree, sLength, seq, MEM, &index, matrix);
	 /*int cnt=0;
	 countMEMnode(tree->m_root, tree, &cnt);
	 cout<<"The no. of MEM nodes is:"<<cnt<<"\n";*/
           
}
// find single source longest distances in a DAG

 
// Graph is represented using adjacency list. Every node of adjacency list
// contains vertex number of the vertex to which edge connects. It also
// contains weight of the edge
class AdjListNode
{
    int v;
    int weight;
public:
    AdjListNode(int _v, int _w)  { v = _v;  weight = _w;}
	
    int getV()       {  return v;  }
    int getWeight()  {  return weight; }
};

class revAdjListNode
{
    int v;
    int weight;
public:
    revAdjListNode(int _v, int _w)  { v = _v;  weight = _w;}
	
    int getV()       {  return v;  }
    int getWeight()  {  return weight; }
};
 
// Class to represent a graph using adjacency list representation
class Graph
{
    int V;    // No. of verticesí
    
    // Pointer to an array containing adjacency lists
    list<AdjListNode> *adj;
	list<revAdjListNode> *revAdj;
 
    // A function used by longestPath
    void topologicalSortUtil(int v, bool visited[], stack<int> &Stack);
public:
    Graph(int V);   // Constructor
 
    // function to add an edge to graph
    void addEdge(int u, int v, int weight);
 
    // Finds longest distances from given source vertex
    void longestPath(int s, int **matrix, int *ml, char *cID);
};
 
Graph::Graph(int V) // Constructor
{
    this->V = V;
	
    adj = new list<AdjListNode>[V];
	revAdj = new list<revAdjListNode>[V];
}
 
void Graph::addEdge(int u, int v, int weight)
{
    AdjListNode node(v, weight);
	revAdjListNode node1(u, weight);
    adj[u].push_back(node); // Add v to uís list
	revAdj[v].push_back(node1); 
	
}
 
// A recursive function used by longestPath.
void Graph::topologicalSortUtil(int v, bool visited[], stack<int> &Stack)
{
    // Mark the current node as visited
    visited[v] = true;
 
    // Recur for all the vertices adjacent to this vertex
    list<AdjListNode>::iterator i;
    for (i = adj[v].begin(); i != adj[v].end(); ++i)
    {
        AdjListNode node = *i;
        if (!visited[node.getV()])
            topologicalSortUtil(node.getV(), visited, Stack);
    }
 
    // Push current vertex to stack which stores topological sort
    Stack.push(v);
}
 
// The function to find longest distances from a given vertex. It uses
// recursive topologicalSortUtil() to get topological sorting.
void Graph::longestPath(int s, int **matrix, int *ml, char *cID)
{
    char fmem1[] = "LPMEMs";
	strcat(fmem1, cID);
	strcat(fmem1, ".txt");
	char fmem2[] = "mem";
	strcat(fmem2, cID);
	strcat(fmem2, ".txt");
	FILE *lpMEMs = fopen(fmem1, "w");
	char temp[100];
	int num = 0;
	FILE *f = fopen(fmem2, "r");
	
	//int dist[V];
	char Gobe[] = "Gobe";
	strcat(Gobe, cID);
	strcat(Gobe, ".csv");
	FILE *fg = fopen(Gobe, "a+");
	
	while(fgets(temp, 100, f)!= NULL) {
	  num++;
	}
	fclose(f);
	int ** mat = (int **) malloc(sizeof(int *) * num+1);
	int k=0;
	for(k=0; k <num+1; k++)
	  mat[k] = (int *) malloc(sizeof(int)*7);
	f = fopen(fmem2, "r");
	k=0;
	while(fgets(temp, 100, f)!= NULL) {
	  int ind = 0;
	  if(k >= 1){
	  
	    char* token = strtok(temp, "\t");
	    mat[k][ind] = atoi(token);
		//printf("token:\n");
        while (token) {
         // printf("%s\n", token);
          token = strtok(NULL, "\t");
		  ind++;
		  mat[k][ind] = atoi(token);
		  if(ind == 4) break;
        }
		//printf("\n");
	    ind++;
	    mat[k][ind] = fabs(mat[k][2] - mat[k][3]) + fabs(mat[k][2] - mat[k][4]) + fabs(mat[k][3] - mat[k][4]);
	  }
	  k++;
	}
	fclose(f);
	printf("Displaying the matrix inside Logest Path \n");
    for(k=0; k<num; k++){
	  printf("%d\t%d\t%d\t%d\t%d\t%d\n", mat[k][0], mat[k][1], mat[k][2], mat[k][3], mat[k][4], mat[k][5]);
    }
	
    stack<int> Stack;
    int dist[V];
	//cout<<"No of vertices "<<V<<"\n";
    //cout<<"for testing"<<V<<"\n"; 
    // Mark all the vertices as not visited
    bool *visited = new bool[V];
    for (int i = 0; i < V; i++)
        visited[i] = false;
 
    // Call the recursive helper function to store Topological Sort
    // starting from all vertices one by one
    for (int i = 0; i < V; i++)
        if (visited[i] == false)
            topologicalSortUtil(i, visited, Stack);
    //cout<<"for testing -- 2\n";
    // Initialize distances to all vertices as infinite and distance
    // to source as 0
    for (int i = 0; i < V; i++)
        dist[i] = NINF;
    dist[s] = 0;
 
    // Process vertices in topological order
    while (Stack.empty() == false)
    {
        // Get the next vertex from topological order
        int u = Stack.top();
        Stack.pop();
 
        // Update distances of all adjacent vertices
        list<AdjListNode>::iterator i;
        if (dist[u] != NINF)
        {
          for (i = adj[u].begin(); i != adj[u].end(); ++i)
             if (dist[i->getV()] < dist[u] + i->getWeight())
                dist[i->getV()] = dist[u] + i->getWeight();
        }
    }
 
    // Print the calculated longest distances
	int max = s;
	FILE *MEMList = NULL;
	//char temp[100];
    for (int i = 0; i < V; i++)
	{
        //(dist[i] == NINF)? cout << "INF ": cout << dist[i] << " ";
	    if(dist[i] >= dist[s]) max = i; 
	}	
	cout<<"The destination vertex is "<<max<<" with score "<<dist[max]<<"\n";
	
	int flag = 0;
	/*MEMList = fopen(fmem2, "r");
	int flag = 0;
	while(fgets(temp, 100, MEMList)!= NULL) {
	  int token = atoi(strtok(temp, " "));
	  if(token == max){
	    printf("token is %d\n", token);
		printf("MEM is %s\n", temp);
	    flag = 1;
		fprintf(lpMEMs, "%s", temp);
	  }
	  if (flag == 1) break;
    }
    fclose(MEMList);*/
	
	int start = max;
	int st = max;
	cout<<max<<"("<<dist[max]<<")";
	list<revAdjListNode>::iterator i1;
	int last = 0;
	while(start != s){
	    int d = -999;
	    for(i1 = revAdj[start].begin(); i1 != revAdj[start].end(); ++i1){
	      //if(dist[i1->getV()] > d) d = i1->getV();
		    if(dist[i1->getV()] > d){
		      d = dist[i1->getV()];
		      st = i1->getV() ;
		    }
	    }
	    cout<<"<-( ";
	  
	  //cout<<st<<" ";
	   int z=0; int dd = 0; int st1 = 0;
	    for(i1 = revAdj[start].begin(); i1 != revAdj[start].end(); ++i1){
		    cout<< " t_s: "<<i1->getV()<<"-"<<dist[i1->getV()]<<" t_e ";
	        if(dist[i1->getV()] == d) {
			    //cout<< "test_start:"<<i1->getV()<<"test_end ";
			  if(last < 1){	 
		        if(z == 0){
			      dd = mat[i1->getV()+1][5];
			      /*int d1 = mat[start+1][2] - mat[i1->getV()+1][2];
			      int d2 = mat[start+1][3] - mat[i1->getV()+1][3];
			      int d3 = mat[start+1][4] - mat[i1->getV()+1][4];
			      dd = fabs(d1-d2)+fabs(d2-d3)+fabs(d1-d3);*/
			      st1 = i1->getV();
		        }	
		        else {
		            /*int d1 = mat[start+1][2] - mat[i1->getV()+1][2];
			        int d2 = mat[start+1][3] - mat[i1->getV()+1][3];
			        int d3 = mat[start+1][4] - mat[i1->getV()+1][4];  
			        int ddd = fabs(d1-d2)+fabs(d2-d3)+fabs(d1-d3);
			        if(ddd< dd){
			         dd = ddd;
			         st1 = i1->getV();
			        }*/
		  	        if(mat[i1->getV()+1][5] < dd){
			          dd = mat[i1->getV()+1][5];
			          st1 = i1->getV();
			        }
			    }
			  }
			  
			  /* added 10_14 */
			  if(last >= 1){	 
		        if(z == 0){
			      //dd = mat[i1->getV()+1][5];
			      int d1 = mat[start+1][2] - mat[i1->getV()+1][2];
			      int d2 = mat[start+1][3] - mat[i1->getV()+1][3];
			      int d3 = mat[start+1][4] - mat[i1->getV()+1][4];
			      dd = fabs(d1-d2)+fabs(d2-d3)+fabs(d1-d3);
			      st1 = i1->getV();
		        }	
		        else {
		            int d1 = mat[start+1][2] - mat[i1->getV()+1][2];
			        int d2 = mat[start+1][3] - mat[i1->getV()+1][3];
			        int d3 = mat[start+1][4] - mat[i1->getV()+1][4];  
			        int ddd = fabs(d1-d2)+fabs(d2-d3)+fabs(d1-d3);
			        if(ddd< dd){
			         dd = ddd;
			         st1 = i1->getV();
			        }
		  	        
			    }
			  }
		      z++;
				
		     
		   }
			//cout<<"st1:"<<st1;
			
		
	    }
		//fprintf(lpMEMs, "\n");
		MEMList = fopen(fmem2, "r");
	    flag = 0; int aaa = 0;
	    while(fgets(temp, 100, MEMList)!= NULL) {
	        int token = atoi(strtok(temp, " "));
            if((token-1) == st1){
	          flag = 1;
		      fprintf(lpMEMs, "%s", temp);
	        }
			aaa++;
	        if (flag == 1) break;
	        
	    }
        fclose(MEMList);
	    cout<<" "<<st1<<" ";
		
		int i=0, j=0;
			for(i=0; i<3; i++){
			  for(j=i+1; j<3; j++){
			    fprintf(fg, "1,S_%d,%d,%d,HSP,+", i, matrix[st1][i],matrix[st1][i]+ml[st1]);
				fprintf(fg, "\n");
				fprintf(fg, "1,S_%d,%d,%d,HSP,+", j, matrix[st1][j],matrix[st1][j]+ml[st1]);
				fprintf(fg, "\n");
			 }	
			} 
	    cout<<")";
	    //start = d;
	    start = st1;
	    last++;
	    //cout<<"<-"<<st<<"("<<dist[st]<<")";
	}
	fclose(lpMEMs);
	fclose(fg);
	delete [] visited;
	cout<<"\n";
	
}
 void LP(char *fname, char *cID)
{
   
    FILE *instream = fopen(fname, "r");
    if(instream == NULL) {
      fprintf(stderr, "Unable to open file: %s\n", fname);
      exit(1);
    }
    
    printf("file name is %s\n", fname);
    char temp[100];
	int no_of_lines=0;
	while(fgets(temp, 100, instream)!= NULL) {
	  no_of_lines++;
	}
	fclose(instream);
	int count=0;
	int **matrix = (int **) malloc(sizeof(int *)*no_of_lines +1);
	int i = 0;
	for(i=0; i<no_of_lines; i++)
	  matrix[i] = (int *) malloc(sizeof(int)*4);
	//int matrix[100][4];
	int *ml = (int *) malloc(sizeof(int)*no_of_lines+1);
	//int ml[100];
	instream = fopen(fname, "r");
	while(fgets(temp, 100, instream)!= NULL) {
	  
	  temp[strlen(temp)-1]='\0';
	  //printf("%s\n", temp);
	  if(count!=0){
	    int m = atoi(strtok(temp, "\t"));
	//	printf("%d\n", m);
        int l = atoi(strtok(NULL, "\t"));
		//printf("%d\n", l);
        int x = atoi(strtok(NULL, "\t"));
        //printf("%d\n", x);
	    int y = atoi(strtok(NULL, "\t"));
            //printf("%d\n", y);
		//if(count==1) printf("%d %d %d %d", m,l,x,y);
	    int z = atoi(strtok(NULL, "\t"));
            //printf("%d\n", z);
		//if(count==1)printf("%d %d %d %d %d", m,l,x,y,z);
	    ml[count-1]=l;
		matrix[count-1][0]=x; 
		matrix[count-1][1]=y; 
		matrix[count-1][2]=z;
	  }
	  count++;
	}
    printf("count=%d\n", count);
	FILE *f;
	char fn[] = "adj.txt";
        printf("fn =%s\n", fn);
	f=fopen(fn, "w");
	count--;
	Graph g(count+2);
	int j=0;
	i = 0;
	int newCount = count+2;
	//int **adjMat = (int **) malloc(sizeof(int *)*count+1);
	int **adjMat = (int **) malloc(sizeof(int *)*newCount+1);
	for(i=0; i<newCount; i++){
	  adjMat[i] = (int *) malloc(sizeof(int)*newCount+1);
	  for(j=0; j<newCount; j++) adjMat[i][j]=0;
	}
	int **deg = (int **) malloc(sizeof(int *)*newCount+1);
	for(i=0; i<newCount; i++)
	  deg[i] = (int *) malloc(sizeof(int)*2+1);
	for(i=0; i<newCount; i++){
	  deg[i][0]=0;
	  deg[i][1]=0;
	}
	int p,q,r;
	cout << "cnt = "<<count<<"\n";
	int v1 = count, v2 = count+1;
	int c3=0, s=0;
	for(p=0; p<count; p++){
	    if(p<2) {printf("TESTING MATRIX %d %d %d\n", matrix[p][0], matrix[p][1], matrix[p][2]);}
	    for(q=p+1; q<count; q++){
		    int c2=0;
			int flag1 = 0, flag2=0, flag3=0; 
			for(r=0; r<3; r++){
			    if(matrix[p][r] > matrix[q][r] && matrix[p][r] > (matrix[q][r]+ml[q]) ) flag1++;
				if(matrix[p][r] < matrix[q][r] && (matrix[p][r] + ml[p]) < matrix[q][r]) flag2++;
			}
			c3++;
			if(flag1 == 3 ){
			    if(p==0) printf("TESTING p = 0, q = %d\n", q);
				// q ----> p Edge 
				g.addEdge(q, p, ml[p]);
                adjMat[q][p]=1;
				deg[q][1]+=1;  //out-degree of q
				deg[p][0]+=1;  // in-degree of p
				
				// v1 -----> q Edge 
				if(adjMat[v1][q]==0){
				  g.addEdge(v1, q, ml[q]);
				  deg[q][0]+=1;
				  deg[v1][1]+=1;
				  adjMat[v1][q]=1;
				}  
				
				// p -----> v2 Edge
				if(adjMat[p][v2]==0){
				  g.addEdge(p, v2, 1);
				  deg[v2][0]+=1;
				  deg[p][1]+=1;
				  adjMat[p][v2]=1;
                }				  
			}
            if(flag2 == 3){
			     if(p==0 ) printf("TESTING p = 0, q = %d\n", q);
			    // Edge  p ----> q  
			    g.addEdge(p, q, ml[q]);
				adjMat[p][q]=1;
				deg[p][1]+=1;
				deg[q][0]+=1;
				
				// Edge v1 -----> p 
				if(adjMat[v1][p] == 0){
				  g.addEdge(v1, p, ml[p]);
				  deg[p][0]+=1;
				  deg[v1][1]+=1;
				  adjMat[v1][p]=1;
				}  
				
				
				// Edge q -----> v2 
				if(adjMat[q][v2]==0){
				  g.addEdge(q, v2, 1);
      			  deg[q][1]+=1;
				  deg[v2][0]+=1;
				  adjMat[q][v2]=1;
				}  
				
			}				
			//if(c3==1) s = p;
		}
			  
	}
	
	for(i=0; i<newCount; i++){
	  for(j=0; j<newCount; j++)
	    fprintf(f, "%d ", adjMat[i][j]);
	  fprintf(f, "\n");
	}  
	
	 
    printf("Printing in-degree and out-degrees\n");
	for(i=0; i<newCount; i++){
      printf("%d\t%d\t%d\n", i, deg[i][0], deg[i][1]);
	  //printf("\n");
    }
	g.longestPath(count, matrix, ml, cID); 
	/*for(i=0; i<count; i++){
     
	    cout << "Following are longest distances from source vertex " << i <<" \n";
         g.longestPath(i);
		 printf("\n");
	
	 
    }*/
    //return 0;
 }
 
 
 
 
 
 /* Count number of genes in each gene */
  void findGeneCounts(int *gCount) {
    
    FILE *f1, *f2, *f3;
	  char fn1[]="Maize.anno";
	  char fn2[]="Sorghum.anno";
	  char fn3[]="Setaria.anno";
	  
	  f1 = fopen(fn1, "r");
	  f2 = fopen(fn2, "r");
	  f3 = fopen(fn3, "r");
	  int c1=0, c2=0, c3=0;
	  int size = 1000;
	  char *tempBuffer = (char *) malloc(sizeof(char) * size);
	  while(fgets(tempBuffer, size, f1)!=NULL){
		tempBuffer[strlen(tempBuffer)-1]='\0';
		if(tempBuffer[0] == '<' || tempBuffer[0] == '>') c1++;
	  }
	  fclose(f1);
	  while(fgets(tempBuffer, size, f2)!=NULL){
		tempBuffer[strlen(tempBuffer)-1]='\0';
		if(tempBuffer[0] == '<' || tempBuffer[0] == '>') c2++;
	  }
	  fclose(f2);
	  while(fgets(tempBuffer, size, f3)!=NULL){
		tempBuffer[strlen(tempBuffer)-1]='\0';
		if(tempBuffer[0] == '<' || tempBuffer[0] == '>') c3++;
	  }
	  fclose(f3);
	  printf("c1=%d c2=%d c3=%d\n", c1, c2, c3);
      gCount[0]=c1; gCount[1]=c2; gCount[2]=c3;
  
  }
  
  
 /* Find Gene Locations */
 void findGeneLocations(int **gene1, int **gene2, int **gene3, int *gCount){
  FILE *f1, *f2, *f3;
  char fn1[]="Maize.anno";
  char fn2[]="Sorghum.anno";
  char fn3[]="Setaria.anno";
  
  f1 = fopen(fn1, "r");
  f2 = fopen(fn2, "r");
  f3 = fopen(fn3, "r");
  int c1=0, c2=0, c3=0;
  int size = 1000;
  char *tempBuffer = (char *) malloc(sizeof(char) * size);
  c1=0;
  while(fgets(tempBuffer, size, f1)!=NULL){
    tempBuffer[strlen(tempBuffer)-1]='\0';
    if(tempBuffer[0] == '<' || tempBuffer[0] == '>'){
      char *a = strtok(tempBuffer, " ");
      int b = atoi(strtok(NULL, " "));
      int c = atoi(strtok(NULL, " "));
      if(b<0) b=0;
      if(c<0) c=0;
      gene1[c1][0] = b;
      gene1[c1][1] = c;
      c1++; 
    }
  }
  fclose(f1);
  c2=0;
  while(fgets(tempBuffer, size, f2)!=NULL){
    tempBuffer[strlen(tempBuffer)-1]='\0';
    if(tempBuffer[0] == '<' || tempBuffer[0] == '>'){
      char *a = strtok(tempBuffer, " ");
      int b = atoi(strtok(NULL, " ")); 
      int c = atoi(strtok(NULL, " "));
       if(b<0) b=0;
      if(c<0) c=0;

      gene2[c2][0] = b;
      gene2[c2][1] = c;
      c2++;
    }
  }
  fclose(f2);
  c3=0;
  while(fgets(tempBuffer, size, f3)!=NULL){
    tempBuffer[strlen(tempBuffer)-1]='\0';
    if(tempBuffer[0] == '<' || tempBuffer[0] == '>'){
      char *a = strtok(tempBuffer, " ");
      int b = atoi(strtok(NULL, " ")); 
      int c = atoi(strtok(NULL, " "));
       if(b<0) b=0;
      if(c<0) c=0;

      gene3[c3][0] = b;
      gene3[c3][1] = c;
      c3++;
    }
  }
  fclose(f3);
  
 }
 
 
 void getChrLocations(string seqName[], string chrLoc[][4]){
 
 cout << seqName[0] << " " << seqName[1] << " " << seqName[2] <<"\n"; 
 
 
string gff1 = "Oryza_sativa_v6.1_dna.gff3";
string gff2 = "Sbicolor_79_gene_exons.gff3";
string gff3 = "Sitalica_164_gene_exons.gff3";
//ifstream file1(gff1);
//ifstream file1(gff2);
//ifstream file1(gff3);
ifstream file1("Oryza_sativa_v6.1_dna.gff3");
ifstream file2("Sbicolor_79_gene_exons.gff3");
ifstream file3("Sitalica_164_gene_exons.gff3"); 
 
/* string gff2 = "Sbicolor_79_gene.gff3";
 string gff1 = "Sitalica_164_gene.gff3";
 string gff3 = "Zmays_284_6a.gene.gff3";
 ifstream file1("Sitalica_164_gene.gff3");
 ifstream file2("Sbicolor_79_gene.gff3");
 ifstream file3("Zmays_284_6a.gene.gff3");*/
 string line;
 size_t flag;
 //string chrLoc[3][3]; 
 char *l1, *l2, *l3;
 if(file1.is_open()){
   while(getline(file1, line)){
     if(line.find("gene") != string::npos && line.find("ID="+seqName[0]) != string::npos)
     {
       //cout << line << "\n";
       l1 = new char[line.length()+1];
       strcpy(l1, line.c_str());
       //cout<<l1 <<"\n";
       char * tok = strtok(l1, "\t");
	   chrLoc[0][3] = tok;
       int i = 1;
       while(tok){
         //if(i==4 || i==5 || i==7) cout<<tok<<" ";
         if(i==4) chrLoc[0][0] = tok;
         if(i==5) chrLoc[0][1] = tok;
         if(i==7) chrLoc[0][2] = tok;
         //cout<<i<< " - " << tok<<"\n ";
         tok = strtok(NULL, "\t");
         i++;
       }
     }
   }
   file1.close();
 }
 
 if(file2.is_open()){
   while(getline(file2, line)){
     if(line.find("gene") != string::npos && line.find("ID="+seqName[1]) != string::npos){
       //cout << line << "\n";
       l2 = new char[line.length()+1];
       strcpy(l2, line.c_str());
       //cout<<l2 <<"\n";
       char * tok = strtok(l2, "\t");
	   chrLoc[1][3] = tok;
       int i = 1;
       while(tok){
         //if(i==4 || i==5 || i==7) cout<<tok<<" ";
         if(i==4) chrLoc[1][0] = tok;
         if(i==5) chrLoc[1][1] = tok;
         if(i==7) chrLoc[1][2] = tok;
         //cout<<i<< " - " << tok<<"\n ";
         tok = strtok(NULL, "\t");
         i++;
       }
     }
   }
   file2.close();
 }
 //cout<< "\n";
 if(file3.is_open()){
   while(getline(file3, line)){
     if(line.find("gene") != string::npos && line.find("ID="+seqName[2]) != string::npos){
       //cout << line << "\n";
       l3 = new char[line.length()+1];
       strcpy(l3, line.c_str());
       //cout<<l3 <<"\n";
       char * tok = strtok(l3, "\t");
	   chrLoc[2][3] = tok;
       int i = 1;
       while(tok){
         //if(i==4 || i==5 || i==7) cout<<tok<<" ";
         if(i==4) chrLoc[2][0] = tok;
         if(i==5) chrLoc[2][1] = tok;
         if(i==7) chrLoc[2][2] = tok;
         //cout<<i<< " - " << tok<<"\n ";
         tok = strtok(NULL, "\t");
         i++;
       }
     }
   }
   file3.close();
 }
 
 
 }
 
 
class Node
{
public:
	void AddLink(int id)
	{
		next.push_back(id);
	}

public:
	vector <int> next;
};

void FindAllPathsAt(vector <Node> &all_nodes, int id, vector < vector<int> > &all_paths, vector <int> tmp)
{
    tmp.push_back(id);

    if(all_nodes[id].next.size() == 0) {
        all_paths.push_back(tmp);
        return;
    }

    for(size_t i=0; i < all_nodes[id].next.size(); i++) {
        vector <int> tmp2(tmp);
        FindAllPathsAt(all_nodes, all_nodes[id].next[i], all_paths, tmp2);
    }
}

void PrintPaths(const vector < vector<int> > &all_paths, int **matrix, int count) 
{    
    int i,j;
    /*for(i=0; i<count-1; i++){
      for(j=0; j<5; j++)
        printf("%d ", matrix[i][j]);
      printf("\n");  
    }*/
    
    int num_rows = sizeof(matrix) / sizeof(matrix[0]);
    int num_cols = sizeof(matrix[0]) / sizeof(matrix[0][0]);
    
    cout <<" rows: "<< num_rows << "cols:"<< num_cols; 
    ofstream myfile ("allPaths.txt");
    size_t max = 0;
    for(size_t i=0; i < all_paths.size(); i++) {
      if(max < all_paths[i].size()) max = all_paths[i].size();
    }
    cout << "max=" << max;
   if (myfile.is_open()){
    
	for(size_t i=0; i < all_paths.size(); i++) {
        // Don't print node if it points to nothing
        if(all_paths[i].size() == 1) {
            continue;
        }
        
        if(all_paths[i].size() == max)
	    {
	    	myfile<<all_paths[i][0];
	    	//cout << all_paths[i][0];
		
		  int score1 = 0; 
		  int score2 = 0;
		  int prev = 0;
		  for(size_t j=1; j < all_paths[i].size(); j++) {
			//if(all_paths[i][j].size() == max)
			  //cout << " -- > " << all_paths[i][j];
			  myfile<<"-"<<all_paths[i][j];
			  //if(j != all_paths[i].size()) score1 += matrix[all_paths[i][j]][0];
			  /*if(j != 1 || j != all_paths[i].size()){
			    int x = matrix[all_paths[i][j]][1] - matrix[prev][1];
			    int y = matrix[all_paths[i][j]][2] - matrix[prev][2];
			    int z = matrix[all_paths[i][j]][3] - matrix[prev][3];
			    score2 += fabs(x-y)+fabs(x-z)+fabs(y-z);          
			    
			  }*/
			  prev = all_paths[i][j];
		   }
		   myfile<< " ("<<score1<<")";
		   myfile<< " ("<<score2<<")";
		   myfile<<"\n";
		   //cout << endl;
		}

		//cout << endl;
	 }
	}
	myfile.close();
}

void allPath(char *fname, char *cID)
{
   
    FILE *instream = fopen(fname, "r");
    if(instream == NULL) {
      fprintf(stderr, "Unable to open file: %s\n", fname);
      exit(1);
    }
    
    printf("file name is %s\n", fname);
    char temp[100];
	int no_of_lines=0;
	while(fgets(temp, 100, instream)!= NULL) {
	  no_of_lines++;
	}
	fclose(instream);
	int count=0;
	int **matrix = (int **) malloc(sizeof(int *)*no_of_lines +1);
	int i = 0;
	for(i=0; i<no_of_lines; i++)
	  matrix[i] = (int *) malloc(sizeof(int)*5);
	//int matrix[100][4];
	int *ml = (int *) malloc(sizeof(int)*no_of_lines+1);
	//int ml[100];
	instream = fopen(fname, "r");
	while(fgets(temp, 100, instream)!= NULL) {
	  
	  temp[strlen(temp)-1]='\0';
	  //printf("%s\n", temp);
	  if(count!=0){
	    int m = atoi(strtok(temp, "\t"));
	//	printf("%d\n", m);
        int l = atoi(strtok(NULL, "\t"));
		//printf("%d\n", l);
        int x = atoi(strtok(NULL, "\t"));
        //printf("%d\n", x);
	    int y = atoi(strtok(NULL, "\t"));
            //printf("%d\n", y);
		//if(count==1) printf("%d %d %d %d", m,l,x,y);
	    int z = atoi(strtok(NULL, "\t"));
            //printf("%d\n", z);
		//if(count==1)printf("%d %d %d %d %d", m,l,x,y,z);
	    ml[count-1]=l;
	    matrix[count-1][0]=l;
		matrix[count-1][1]=x; 
		matrix[count-1][2]=y; 
		matrix[count-1][3]=z;
	  }
	  count++;
	}
    printf("count=%d\n", count);
	FILE *f;
	char fn[] = "adj.txt";
        printf("fn =%s\n", fn);
	f=fopen(fn, "w");
	count--;
	//Graph g(count+2);
	vector <Node> all_nodes(count+2);
	int j=0;
	i = 0;
	int newCount = count+2;
	//int **adjMat = (int **) malloc(sizeof(int *)*count+1);
	int **adjMat = (int **) malloc(sizeof(int *)*newCount+1);
	for(i=0; i<newCount; i++){
	  adjMat[i] = (int *) malloc(sizeof(int)*newCount+1);
	  for(j=0; j<newCount; j++) adjMat[i][j]=0;
	}
	int **deg = (int **) malloc(sizeof(int *)*newCount+1);
	for(i=0; i<newCount; i++)
	  deg[i] = (int *) malloc(sizeof(int)*2+1);
	for(i=0; i<newCount; i++){
	  deg[i][0]=0;
	  deg[i][1]=0;
	}
	int p,q,r;
	cout << "cnt = "<<count<<"\n";
	int v1 = count, v2 = count+1;
	int c3=0, s=0;
	for(p=0; p<count; p++){
	    if(p<2) {printf("TESTING MATRIX %d %d %d\n", matrix[p][0], matrix[p][1], matrix[p][2]);}
	    for(q=p+1; q<count; q++){
		    int c2=0;
			int flag1 = 0, flag2=0, flag3=0; 
			for(r=0; r<3; r++){
			    if(matrix[p][r] > matrix[q][r]) flag1++;
				if(matrix[p][r] < matrix[q][r]) flag2++;
			}
			c3++;
			if(flag1 == 3 ){
			    if(p==0) printf("TESTING p = 0, q = %d\n", q);
				// q ----> p Edge 
				//g.addEdge(q, p, ml[p]);
				all_nodes[q].AddLink(p);
                adjMat[q][p]=1;
				deg[q][1]+=1;  //out-degree of q
				deg[p][0]+=1;  // in-degree of p
				
				// v1 -----> q Edge 
				if(adjMat[v1][q]==0){
				  //g.addEdge(v1, q, ml[q]);
				  all_nodes[v1].AddLink(q);
				  deg[q][0]+=1;
				  deg[v1][1]+=1;
				  adjMat[v1][q]=1;
				}  
				
				// p -----> v2 Edge
				if(adjMat[p][v2]==0){
				  //g.addEdge(p, v2, 1);
				  all_nodes[p].AddLink(v2);
				  deg[v2][0]+=1;
				  deg[p][1]+=1;
				  adjMat[p][v2]=1;
                }				  
			}
            if(flag2 == 3){
			     if(p==0 ) printf("TESTING p = 0, q = %d\n", q);
			    // Edge  p ----> q  
			    //g.addEdge(p, q, ml[q]);
			     all_nodes[p].AddLink(q);
				adjMat[p][q]=1;
				deg[p][1]+=1;
				deg[q][0]+=1;
				
				// Edge v1 -----> p 
				if(adjMat[v1][p] == 0){
				  //g.addEdge(v1, p, ml[p]);
				   all_nodes[v1].AddLink(p);
				  deg[p][0]+=1;
				  deg[v1][1]+=1;
				  adjMat[v1][p]=1;
				}  
				
				
				// Edge q -----> v2 
				if(adjMat[q][v2]==0){
				  //g.addEdge(q, v2, 1);
				   all_nodes[q].AddLink(v2);
      			  deg[q][1]+=1;
				  deg[v2][0]+=1;
				  adjMat[q][v2]=1;
				}  
				
			}				
			//if(c3==1) s = p;
		}
			  
	}
	
	for(i=0; i<newCount; i++){
	  for(j=0; j<newCount; j++)
	    fprintf(f, "%d ", adjMat[i][j]);
	  fprintf(f, "\n");
	}  
	
	 
    printf("Printing in-degree and out-degrees\n");
	for(i=0; i<newCount; i++){
      printf("%d\t%d\t%d\n", i, deg[i][0], deg[i][1]);
	  //printf("\n");
    }
	//g.longestPath(count, matrix, ml, cID); 
	vector <int> tmp;
	vector <vector<int> > all_paths;
	FindAllPathsAt(all_nodes, v1, all_paths, tmp);
	cout<< "All paths at node "<< v1 << endl;
	PrintPaths(all_paths, matrix, count);
	/*for(i=0; i<count; i++){
     
	    cout << "Following are longest distances from source vertex " << i <<" \n";
         g.longestPath(i);
		 printf("\n");
	
	 
    }*/
    //return 0;
 }
 
 
 /* START OF MAIN() Function */

int main(int argc, char ** argv)
{
    string S = "sGCTAGCTAATATATATATATATATC$";
    bool txt = false;
    bool dot = false;
    bool sort = false;
    bool count = true;
    bool manyMEMs = false;
    bool multiFasta = false;
    bool printModGenome = false;
    int minMEM = 1;
    
    string filename, kmerLenFile;
    string modGenomeFilename;
    string CDGprefix;
    
    try{
        
        for (int c = 1; c < argc; c++)
        {
            
            if (!strcmp(argv[c], "-h"))       { printHelp(); }
            else if (!strcmp(argv[c], "-dot"))     { dot = true; }
            else if (!strcmp(argv[c], "-txt"))     { txt = true; }
            else if (!strcmp(argv[c], "-sort"))    { sort = true; }
            else if (!strcmp(argv[c], "-mem"))     { MEM = true; minMEM = atoi(argv[c+1]); c++; }
            else if (!strcmp(argv[c], "-manyMEMs")){ manyMEMs = true; kmerLenFile = argv[c+1]; c++;}
            else if (!strcmp(argv[c], "-cdg"))     { CDG_Filename = argv[c+1]; c++; }
            else if (!strcmp(argv[c], "-file"))    { filename = argv[c+1]; c++; }
            else if (!strcmp(argv[c], "-multiFa")) { multiFasta = true; }
            else if (!strcmp(argv[c], "-printGenome"))   {printModGenome=true; modGenomeFilename = argv[c+1]; c++;}
        }
        cout<<"filename: "<< filename << "\n";
        if (filename.empty())
        {
          printHelp();
        }

            SuffixTree * tree;

            cerr << "Loading " << filename << endl;
            
            S="s";
            
            ifstream file;
            file.open(filename.c_str());
            
            std::ofstream fastaStartOutfile;
            
            
            string buffer;
            
            if(multiFasta)
            {
                string fileName_str= CDG_Filename + "fastaPos.txt";
                fastaStartOutfile.open(fileName_str.c_str());
            }
            int *sLength;
            int no_seq=0;           
            while(getline(file, buffer))
            {
                
                //cerr<<"\n The length of "<<buffer<<" is :"<<buffer.length()<<"\n"; 

                if (buffer[0] == '>')
                {
                    if(!isspace(buffer[1])) no_seq++;
                    if(S.length()==1) //beginning of first line of file
                        cerr << buffer << endl;
                    //otherwise insert N
                    else
                        S += 'N';
                    
                    if(multiFasta)
                    {
                        fastaStartOutfile << S.length() <<endl;
                    }
                }
                else
                {
                    for (int i = 0; i < buffer.length(); i++)
                    {
                        char b = toupper(buffer[i]);
                        if (b == 'A' || b == 'C' || b == 'G' || b=='T' )
                        {
                            S += b;
                        }
                        else //catch ambiguity codes
                        {
                            S += 'A';
                        }
                    }
                }
                //cerr<<"length:"<<buffer.length()<<"\n";
            }
            sLength = (int *) malloc(sizeof(int) * no_seq);
			printf("The no. of sequences is %d\n", no_seq); 
            int max=0;
                        
            file.close();
			
            
            file.open(filename.c_str());

            std::ofstream fastaStartOutfile1;


           // string buffer;
				if(multiFasta)
            {
                string fileName_str= CDG_Filename + "fastaPos.txt";
                fastaStartOutfile1.open(fileName_str.c_str());
            }
            int i=0,j=0,k=0;
            //string[] seq = new string[no_seq];
            string seq[10];
			string seqName[10];
			int **geneLoc = (int **) malloc(sizeof(int *) * (no_seq+1));
			int c = 0;
			for(c=0; c<no_seq; c++) 
			  geneLoc[c] = (int *) malloc(sizeof(int)*3);
			
			int x1 = filename.find_last_of("_");
			//int y1 = filename.find_first_of(".");
			int y1 = filename.find_last_of(".");
			string  ID1 = filename.substr(x1+1, y1-x1-1);
			cout<<"ID = " << ID1 <<"\n";
			string  fileGobe= "Gobe"+ID1+".csv";
			
			char *fileG = new char[fileGobe.length()+1];
			strcpy(fileG, fileGobe.c_str());
			char Gobe[100];
            strcpy(Gobe, fileG);
			delete [] fileG;
			FILE *csv = fopen(Gobe, "w");
			
			
            while(getline(file, buffer))
            {

                //cerr<<"\n The length of "<<buffer<<" is :"<<buffer.length()<<"\n"; 

                if (buffer[0] != '>')
                {
                
                  sLength[i]=buffer.length();
                  seq[i]=buffer;
                  i++;
                }
				else if(!isspace(buffer[1])){
				  int len = buffer.length();
				  seqName[j] = buffer.substr(1, len);
				  cout<<seqName[j]<<"\n";
				  j++;
				}
				else {
				  char gene[buffer.length()+1];
				  strcpy(gene, buffer.c_str());
				  printf("%s\n", gene);
				  int aa = 0;
				  char* token = strtok(gene, " ");
				  
				  while (token) {
				    aa++;
				    printf("token: %d\n", atoi(token));
				    if(aa == 2) geneLoc[k][0] = atoi(token);
					if(aa == 3) geneLoc[k][1] = atoi(token);
					
					token = strtok(NULL, " ");
				  }
				  k++;
				
				}
				
				  
				
                //cerr<<"length:"<<buffer.length()<<"\n";
            }
            file.close();
			string chrLoc[3][4]; 
            getChrLocations(seqName, chrLoc);
            /*for(i=0; i<no_seq; i++){
              for(j=0; j<no_seq; j++)
                cout << chrLoc[i][j]<<" ";
              cout<< "\n";
            }*/
			
			i=0;
			// print the first few lines into Gobe file
			for(j=0; j<no_seq; j++){
			   fprintf(csv, "S_%d,S_%d,0,%d,track,", j, j,sLength[j]);
			   fprintf(csv,"\n");
            }
			for(i=0; i<no_seq; i++){
			  printf("%d %d\n", geneLoc[i][0], geneLoc[i][1]);
			  fprintf(csv, "%c%c_Gene,S_%d,%d,%d,gene,+",seqName[i][0],seqName[i][1],i,geneLoc[i][0],geneLoc[i][1]);
              fprintf(csv, "\n");			  
			 } 
             j=0;
            for(j=0; j<no_seq; j++) if(max > sLength[j]) max = sLength[j];
            for(j=0; j<no_seq; j++){
              cout<<"length of sequence"<<j+1<<" is: "<<sLength[j]<<"\n";
              //cout<<"S"<<j+1<<" :"<<seq[j]<<"\n";
            }
			fclose(csv);
			i=0;
			for(i=0; i<no_seq; i++)
			  printf("%d %d\n", geneLoc[i][0], geneLoc[i][1]);
			  
            j=0;
            for(j=0; j<i; j++) if(max > sLength[j]) max = sLength[j];
            for(j=0; j<i; j++){
              cout<<"length of "<<j+1<<" sequence is: "<<sLength[j]<<"\n";
              //cout<<"S"<<j+1<<" :"<<seq[j]<<"\n";
            }
            if(multiFasta)
                fastaStartOutfile.close();
            
            S += "$";
            
        if(printModGenome)
        {
            cerr<<"printing genome string to file"<<endl;
            //save S in file
            ofstream modGenomeFile;
            modGenomeFile.open(modGenomeFilename.c_str());
            modGenomeFile << S;
            modGenomeFile.close();
        }
        
        
        cerr << "Creating Suffix Tree for string of length " << S.length()-2 << endl;
        
        /*
        timeval starttime;
        timeval endtime;
        
        gettimeofday(&starttime, NULL);
        */
        tree = buildUkkonenSuffixTree(S);
        
        /*       
        gettimeofday(&endtime, NULL);
        
        treeint seconds = endtime.tv_sec - starttime.tv_sec;
        treeint microseconds = seconds*1000000 + endtime.tv_usec - starttime.tv_usec;
        double elapsed = microseconds / 1000000.0;
        
        cerr << "Suffix Tree Construction time, microseconds: " << microseconds  << endl;
        cerr << "Suffix Tree Construction time, seconds: " << elapsed << endl;
        */
        
        int nodesize = sizeof(SuffixNode);
        treeintLarge totalsize = nodesize * SuffixNode::s_nodecount;
        
        cerr << "Total nodes: " << SuffixNode::s_nodecount << endl;
        cerr << "Total space: " << totalsize << endl;
        cerr << "Node size: "   << nodesize << endl;
        
        double bytesperbase = totalsize / ((double) S.length()-2);
        cerr << "Bytes/Base: " << bytesperbase << endl;
        
        
        //construct compressed de Bruijn graph for sequence
        int numKmerLens = 0;
        int * kmerLens;
        
        if(MEM)
        {
            numKmerLens = 1;
            kmerLens = new int[numKmerLens];
            kmerLens[0] = minMEM;
        }
        else if(manyMEMs)
        {
            //store K values in text file, with different K on each line (blank line is ignored)
            // command line argument for manyMEMs that is the filename
            // scan the file twice: once to count number of lines and store in numKmerLens
            // then create array kmerLens = new int[numKmerLens];
            // then scan file a second time to poplate array of kmerLens with actual values of k
            
            
            if (!kmerLenFile.empty())
            {
                ifstream kfile;
                string buffer;
                
                kfile.open(kmerLenFile.c_str());
                
                while(kfile >> buffer)
                {
                    if(buffer != "")
                        numKmerLens++;
                }
                //return to beginning of file
                kfile.clear();
                kfile.seekg(0, ios::beg);
                
                kmerLens = new int[numKmerLens];
                
                int i = 0;
                for(int i = 0; i<numKmerLens; i++)
                {
                    kfile >> buffer;
                    if(buffer != "")
                    {
                        int nextK = atoi(buffer.c_str());
                        kmerLens[i] = nextK;
                    }
                    
                }
                
                kfile.close();
            }
            
            CDGprefix = CDG_Filename;
        }
        
        cerr<<"numKmerLens = "<<numKmerLens<<endl;
        
        for(int i = 0; i < numKmerLens; i++)
        {
            Kmer_Len = kmerLens[i];
            cerr<<endl<<" Kmer_Len = "<<Kmer_Len<<endl;
            cerr<<" maxMEMstrdepth="<<tree->m_maxMEMstrdepth<<endl;
 
           
            if(i==0)
                tree->markMEMnodes(Kmer_Len, true);
            else
                tree->markMEMnodes(Kmer_Len, false);
            
            
            tree->preprocessLMA();
            cerr<<" finished marking MEM nodes"<<endl;
            
           
            cerr<<" finished constructing suffix tree and marking MEMs"<<endl;
            
            
			int cnt=0;
			//cout<<"sLength="<<sLength<<"\n";
	        countMEMnode(tree->m_root, tree, sLength, seq, &cnt);
	        cout<<"The no. of MEM nodes is:"<<cnt<<"\n";
			
			string MEM[cnt];
			
		    int **matrix = (int **) malloc(sizeof(int *) * cnt);
            int p = 0;	
            for(p = 0; p < cnt; p++)
              matrix[p] = (int *) malloc(sizeof(int)*no_seq);
			printMEMs(tree, sLength, seq, MEM, matrix);
			
			// printing All MEMs 
			//cout<<"printing All MEMs "<<"\n";
			int countNonEmptyMEMs = 0;
            for(p = 0; p < cnt; p++)
			  if(!MEM[p].empty())
			  {
			    //cout<<MEM[p]<<"\n"; 
				countNonEmptyMEMs++;
			  }
			  int q = 0, r;  
			/*cout <<"Number of non empty MEMs "<< countNonEmptyMEMs << "\n";
			cout<<"---------------------------"<<"\n";  
			int q = 0, r;  
			for(p=0; p<cnt; p++){
			  for(q=0; q<no_seq; q++){
			    printf("%d ", matrix[p][q]);
			  }
			  printf("\n");
			}*/
			
			// removing MEMs belong to Gene (put 0) 
			
			
		    for(p=0; p<cnt; p++){
			  for(q=0; q<no_seq; q++){
			    if(q == 0){ 
				  //for(r=0; r<gCount[0]; r++)
				   if(matrix[p][q] >= geneLoc[q][0] && matrix[p][q] <= geneLoc[q][1])		
				     matrix[p][q] = 0;
                }				  
				if(q == 1){ 
				  //for(r=0; r<gCount[1]; r++)
				   if(matrix[p][q] >= geneLoc[q][0] && matrix[p][q] <= geneLoc[q][1])		
				     matrix[p][q] = 0;
                }	
				if(q == 2){ 
				  //for(r=0; r<gCount[2]; r++)
				   if(matrix[p][q] >= geneLoc[q][0] && matrix[p][q] <= geneLoc[q][1])		
				     matrix[p][q] = 0;
                }			
			  }
			}  
			
			// removing MEMS which are on both sides of gene 
			for(p=0; p<cnt; p++){
			  int cc1 = 0, cc2 = 0;
			  for(q=0; q<no_seq; q++){
			    if(matrix[p][q] < geneLoc[q][0]) cc1++;
				if(matrix[p][q] > geneLoc[q][1]) cc2++;
                 				
			  }
			  if(cc1 != no_seq && cc2 != no_seq){
			    for(q=0; q<no_seq; q++){
				  matrix[p][q] = 0;
				}
			  }
			}  
			
			// This is Specific for this CNS project as file names has a pattern  
			int x = filename.find_last_of("_");
			//int y = filename.find_first_of(".");
			int y = filename.find_last_of(".");
			string  ID = filename.substr(x+1, y-x-1);
			char *cID = new char[ID.length()+1];
			strcpy(cID, ID.c_str());
			
			cout << "cID = " << cID<<"\n";
			
			
			string  fileM1= "mem"+ID+".txt";
			string  fileM2= "mems"+ID+".txt";
			
			char *fileC1 = new char[fileM1.length()+1];
			strcpy(fileC1, fileM1.c_str());
			char mfname[100];
            strcpy(mfname, fileC1);
			delete [] fileC1;
			
			char *fileC2 = new char[fileM2.length()+1];
			strcpy(fileC2, fileM2.c_str());
			char mems[100];
            strcpy(mems, fileC2);
			delete [] fileC2;
            
            // q=0;
			//char mfname[] = "mem.txt";
			//char mems[] = "mems.txt";
			FILE *mem = fopen(mfname, "w");
			FILE *memALL = fopen(mems, "w");
			fprintf(mem, "MEM\tLength\tS_1\tS_2\tS_3\n");
            cout<<"MEM\t";
			cout<<"Length\t";
            for(p=0; p<no_seq; p++)
              cout<<"S_"<<p+1<<"\t";
            cout<<"\n";
			int l = 1;
            for(p=0; p<cnt; p++){
			  int u = 0, count1 = 0, count2=0;
			  for(q=0; q<no_seq; q++)
			  {
				if(matrix[p][q] !=0) u++;
                //if(matrix[p][q] >=0 && matrix[p][q] <=10000) count1++;
				//if(matrix[p][q] >14000) count2++;
			  }
			  //if(u >= 2){	//if at least two sequences share a MEM
			  //if(u==no_seq && (count1 == no_seq || count2 == no_seq)){ // if all sequences share the MEM
              if(u==no_seq ){			   
			    cout<<l;
			    //fprintf(mem, "%d\t", p); 
				fprintf(mem, "%d\t", l);
				cout<<"\t"<<MEM[p].length();
				char *cstr = new char[MEM[p].length() + 1];
                strcpy(cstr, MEM[p].c_str());
				fprintf(memALL, "%s", cstr);
				fprintf(memALL, "\n");
				delete [] cstr;
				fprintf(mem, "%d\t", MEM[p].length());
                for(q=0; q<no_seq; q++){
                  cout<<"\t"<<matrix[p][q];
				  fprintf(mem, "%d\t", matrix[p][q]);
				}  
				l++; 
                cout <<"\n";
				fprintf(mem, "\n");
			  }
            }
			fclose(mem);
			fclose(memALL);
			LP(mfname, cID);
			/*int **memMatrix = (int **) malloc(sizeof(int *)* l);
			for(i=0; i<l; i++)
			  memMatrix[i] = (int *) malloc(sizeof(int)*(no_seq+2));
			mem = fopen(mfname, "r");
			char temp[500];
			l=0;
			while(fgets(temp, 500, mem) != NULL){
			  temp[strlen(temp)-1] = '\0';
			  int p = 0;
			  int tok = atoi(strtok(temp, "\t"));
			  while(tok){
			    tok = atoi(strtok(NULL, "\t")); 
			    memMatrix[l][p] = tok;
			    p++;
			  }
			  l++;
			}  
			fclose(mem);*/
			//allPath(mfname, cID);
			delete [] sLength;
			//delete [] seqName;
			delete [] geneLoc;
			//delete [] MEM;
			delete [] matrix;
			
			
			/* Create CSV files for Visualization */
			char csv1[] = "MEM1";
			strcat(csv1, cID);
			strcat(csv1, ".csv");
			
			char csv2[] = "MEM2";
			strcat(csv2, cID);
			strcat(csv2, ".csv");
			
			
			char lpmemTxt[] = "LPMEMs";
			strcat(lpmemTxt, cID);
			strcat(lpmemTxt, ".txt");
			
			FILE *cs1 = fopen(csv1, "w");
			FILE *cs2 = fopen(csv2, "w");
			i=0;
			fprintf(cs2, "%s,", "Length");
			for(i=0; i<no_seq; i++){
			  if(i < (no_seq-1)){
			    fprintf(cs1, "%c%c,", seqName[i][0], seqName[i][1]); 
			    fprintf(cs2, "%c%c,", seqName[i][0], seqName[i][1]); 
			  }
			  else{
			    fprintf(cs1, "%c%c\n", seqName[i][0], seqName[i][1]); 
			    fprintf(cs2, "%c%c\n", seqName[i][0], seqName[i][1]);  
			  }
			    
			}
			
			fprintf(cs2, "%d,", 0);
			for(i=0; i<no_seq; i++){
			  if(i<no_seq-1){
			    fprintf(cs1, "%d,", 0);
			    fprintf(cs2, "%d,", 0);
			  }
			  else{
			    fprintf(cs1, "%d\n", 0);
			    fprintf(cs2, "%d\n", 0);
			  }
			}
			FILE *cs3 = fopen(lpmemTxt, "r");
			printf("%s\n", lpmemTxt);
			//char * tempb = (char *) malloc(sizeof(char)*1000);
			char tempb[1000];
			while(fgets(tempb, 1000, cs3) != NULL){
			  tempb[strlen(tempb)-1] = '\0';
			  int x1 = atoi(strtok(tempb, "\t"));
			  int x2 = atoi(strtok(NULL, "\t"));
			  fprintf(cs2, "%d,", x2);
			  int x3 = atoi(strtok(NULL, "\t"));
			  fprintf(cs2, "%d,", x3);
			  fprintf(cs1, "%d,", x3);
			  int x4 = atoi(strtok(NULL, "\t"));
			  fprintf(cs2, "%d,", x4);
			  fprintf(cs1, "%d,", x4);
			  int x5 = atoi(strtok(NULL, "\t"));
			  fprintf(cs2, "%d\n", x5);
			  fprintf(cs1, "%d\n", x5);
			}
			
			fprintf(cs2, "%d,", 0);
			for(i=0; i<no_seq; i++){
			  if(i<no_seq-1){
			    fprintf(cs1, "%d,", 30000);
			    fprintf(cs2, "%d,", 30000);
			  }
			  else{
			    fprintf(cs1, "%d", 30000);
			    fprintf(cs2, "%d", 30000);
			  }
			}
			
			
			fclose(cs1);
			fclose(cs2);
			fclose(cs3);
			cout << "max = "<<max <<"\n";
			
			char csv3[] = "CHRMEM_";
			strcat(csv3, cID);
			strcat(csv3, ".csv");
			
			char csv4[] = "MEM2";
			strcat(csv4, cID);
			strcat(csv4, ".csv");
			
			FILE *cs4 = fopen(csv4, "r");
			FILE *cs5 = fopen(csv3, "w");
            i = 0;
            cout<< "before while \n";
            //cout << seqName[0] << " " << seqName[1] << " " << seqName[2] <<"\n"; 
            fprintf(cs5, "%s\n", seqName[0].c_str());
            fprintf(cs5, "%s\n", seqName[1].c_str());
            fprintf(cs5, "%s\n", seqName[2].c_str());
            while(fgets(tempb, 1000, cs4) != NULL){
              tempb[strlen(tempb)-1] = '\0';
              if(i!=0){
                int x1 = atoi(strtok(tempb, ","));
                int x2 = atoi(strtok(NULL, ","));
                int x3 = atoi(strtok(NULL, ","));
                int x4 = atoi(strtok(NULL, ","));
				cout << " x2 = "<< x2 << " x3 = "<< x3 << " x4 = "<< x4 <<"\n";
                //if(!(x1 == 0 && x2 == 0 && x3 == 0 && x4 == 0) && !(x1 == 0 && x2 == 30000 && x3 ==30000  && x4 == 30000)){
                if(x1 != 0){  
                  if(chrLoc[0][no_seq-1][0] == '+'){
                    x2 = x2 + atoi(chrLoc[0][0].c_str()) - 10000;                    
                  }
                  else{
                    x2 = atoi(chrLoc[0][1].c_str()) + 10000 - x2 - x1;
                  }
                  
                  if(chrLoc[1][no_seq-1][0] == '+'){
                    x3 = x3 + atoi(chrLoc[1][0].c_str()) - 10000;                    
                  }
                  else{
                    x3 = atoi(chrLoc[1][1].c_str()) + 10000 - x3 - x1;
                  }
                  
                  if(chrLoc[2][no_seq-1][0] == '+'){
                    x4 = x4 + atoi(chrLoc[2][0].c_str()) - 10000;                    
                  }
                  else{
                    x4 = atoi(chrLoc[2][1].c_str()) + 10000 - x4 - x1;
                  }
				  cout << "x2 = "<< x2 << "x3 = "<< x3 << "x4 = "<< x4 <<"\n";
                  fprintf(cs5, "%d,%d,%d,%d,%s,%s,%s\n", x1, x2, x3, x4, chrLoc[0][3].c_str(), chrLoc[1][3].c_str(), chrLoc[2][3].c_str());
                  
                }
                //else fprintf(cs5, "%s\n", tempb);
              }
              else fprintf(cs5, "%s,Chr1,Chr2,Chr3\n", tempb);
			  i++;
            }			
			
			fclose(cs4);
			fclose(cs5);
			delete [] cID;
			
			//delete [] cID;
			
           
        }
        
        delete [] kmerLens;
        
        
        if (txt)  { tree->dumpTreeText(cout); }
        else if (dot)  { tree->dumpTree(); cerr<<" dumping ST to dot file";}
        else if (sort) { tree->dumpTreeSorted(cout, tree->m_root, ""); }
        
        
    }
    catch(std::bad_alloc& ba){
        cerr<<"bad_alloc caught "<<ba.what()<<endl;
        
    }
    catch(exception& ex){
        cerr<<" some other exception caught"<<ex.what() <<endl;
    }
    
}
