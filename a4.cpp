#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <map>
#include <cmath>
#include <algorithm>


// Format checker just assumes you have Alarm.bif and Solved_Alarm.bif (your file) in current directory
using namespace std;

// Our graph consists of a list of nodes where each node is represented as follows:
class Graph_Node{

private:
	string Node_Name;  // Variable name 
    int index;  //index in network
	vector<int> Children; // Children of a particular node - these are index of nodes in graph.
	vector<string> Parents; // Parents of a particular node- note these are names of parents
	int nvalues;  // Number of categories a variable represented by this node can take
	vector<string> values; // Categories of possible values
	vector<float> CPT; // conditional probability table as a 1-d array . Look for BIF format to understand its meaning

public:
	// Constructor- a node is initialised with its name and its categories
    Graph_Node(){
        Node_Name = "";
        nvalues = 0;
    }
    Graph_Node(string name,int n,vector<string> vals)
	{
		Node_Name=name;
	
		nvalues=n;
		values=vals;
		

	}
	string get_name()
	{
		return Node_Name;
	}
    int get_index()
	{
		return index;
	}
	vector<int> get_children()
	{
		return Children;
	}
	vector<string> get_Parents()
	{
		return Parents;
	}
	vector<float> get_CPT()
	{
		return CPT;
	}
	int get_nvalues()
	{
		return nvalues;
	}
	vector<string> get_values()
	{
		return values;
	}
    void set_index(int val){
        index = val;
    }
	void set_CPT(vector<float> new_CPT)
	{
		CPT.clear();
		CPT=new_CPT;
	}
    void set_Parents(vector<string> Parent_Nodes)
    {
        Parents.clear();
        Parents=Parent_Nodes;
    }
    // add another node in a graph as a child of this node
    int add_child(int new_child_index )
    {
        for(int i=0;i<Children.size();i++)
        {
            if(Children[i]==new_child_index)
                return 0;
        }
        Children.push_back(new_child_index);
        return 1;
    }



};


 // The whole network represted as a list of nodes
class network{

public:
	map<string, Graph_Node> Pres_Graph;
    map<int, string> indexToName;

	int addNode(Graph_Node node)
	{
		Pres_Graph[node.get_name()] = node;
		return 0;
	}
    
    
	int netSize()
	{
		return Pres_Graph.size();
	}
    // get the index of node with a given name
    int get_index(string val_name)
    {
        if (Pres_Graph.find(val_name) != Pres_Graph.end())
            return Pres_Graph[val_name].get_index();
        else
            return -1;
        // list<Graph_Node>::iterator listIt;
        // int count=0;
        // for(listIt=Pres_Graph.begin();listIt!=Pres_Graph.end();listIt++)
        // {
        //     if(listIt->get_name().compare(val_name)==0)
        //         return count;
        //     count++;
        // }
        // return -1;
    }
// get the node at nth index
    // list<Graph_Node>::iterator get_nth_node(int n)
    // {
    //    list<Graph_Node>::iterator listIt;
    //     int count=0;
    //     for(listIt=Pres_Graph.begin();listIt!=Pres_Graph.end();listIt++)
    //     {
    //         if(count==n)
    //             return listIt;
    //         count++;
    //     }
    //     return listIt; 
    // }

    //get the iterator of a node with a given name
    // list<Graph_Node>::iterator search_node(string val_name)
    // {
    //     list<Graph_Node>::iterator listIt;
    //     for(listIt=Pres_Graph.begin();listIt!=Pres_Graph.end();listIt++)
    //     {
    //         if(listIt->get_name().compare(val_name)==0)
    //             return listIt;
    //     }
    
    //         cout<<"node not found\n";
    //     return listIt;
    // }
	

};

network read_network()
{
	network Alarm;
	string line;
	int find=0;
  	ifstream myfile("alarm.bif"); 
  	string temp;
  	string name;
  	vector<string> values;
  	
    if (myfile.is_open())
    {
    	while (! myfile.eof() )
    	{
    		stringstream ss;
      		getline (myfile,line);
      		
      		
      		ss.str(line);
     		ss>>temp;
     		
     		
     		if(temp.compare("variable")==0)
     		{
                    
     				ss>>name;
     				getline (myfile,line);
                   
     				stringstream ss2;
     				ss2.str(line);
     				for(int i=0;i<4;i++)
     				{
     					
     					ss2>>temp;
     					
     					
     				}
     				values.clear();
     				while(temp.compare("};")!=0)
     				{
     					values.push_back(temp);
     					
     					ss2>>temp;
    				}
     				Graph_Node new_node(name,values.size(),values);
                    new_node.set_index(Alarm.netSize());
                    Alarm.indexToName[Alarm.netSize()] = name;
     				int pos=Alarm.addNode(new_node);

     				
     		}
     		else if(temp.compare("probability")==0)
     		{
                    
     				ss>>temp;
     				ss>>temp;
     				
                    // list<Graph_Node>::iterator listIt;
                    // list<Graph_Node>::iterator listIt1;
     				// listIt=Alarm.search_node(temp);
                    int index=Alarm.Pres_Graph[temp].get_index();
                    string initName = temp;
                    ss>>temp;
                    values.clear();
     				while(temp.compare(")")!=0)
     				{
                        // listIt1=Alarm.search_node(temp);
                        // listIt1->add_child(index);
                        Alarm.Pres_Graph[temp].add_child(index);
     					values.push_back(temp);
     					
     					ss>>temp;

    				}
                    Alarm.Pres_Graph[initName].set_Parents(values);
    				getline (myfile,line);
     				stringstream ss2;
                    
     				ss2.str(line);
     				ss2>> temp;
                    
     				ss2>> temp;
                    
     				vector<float> curr_CPT;
                    string::size_type sz;
     				while(temp.compare(";")!=0)
     				{
                        
     					curr_CPT.push_back(atof(temp.c_str()));
     					
     					ss2>>temp;
                       
                        

    				}
                    
                    Alarm.Pres_Graph[initName].set_CPT(curr_CPT);


     		}
            else
            {
                
            }
     		
     		

    		
    		
    	}
    	
    	if(find==1)
    	myfile.close();
  	}
  	
  	return Alarm;
}


void expectation(network& Alarm, vector<pair<vector<string>,float>>& data, vector<int> unknownIndex, vector<vector<double>>& unknownDistribution){
    int dataRows = data.size();
    for(int i=0;i<dataRows;i++){
        int index = unknownIndex[i];
        string name = Alarm.indexToName[index];
        int nvalues = Alarm.Pres_Graph[name].get_nvalues();
        double sum = 0;
        for(int j=0;j<nvalues;j++){
            double prob = 1;
            for(auto it = Alarm.Pres_Graph.begin(); it!=Alarm.Pres_Graph.end(); it++){
                vector<string> Parents = it->second.get_Parents();
                int nuValues = it->second.get_nvalues();
                vector<int> cptSizePrefix(Parents.size(),nuValues);
                // vector<int> cptSizePrefix(Parents.size(),nvalues);
                for(int k=Parents.size()-2;k>=0;k--){
                    int numValues = Alarm.Pres_Graph[Parents[k+1]].get_nvalues();
                    cptSizePrefix[k] = cptSizePrefix[k+1]*numValues;
                }
                int location;
                if(name == it->second.get_name())
                    location = j;
                else{
                    vector<string> values = it->second.get_values();
                    for(int k=0;k<values.size();k++){
                        if(values[k] == data[i].first[it->second.get_index()]){
                            location = k;
                            break;
                        }
                    }
                }
                for(int k=0;k<Parents.size();k++){
                    if(Parents[k] == name){
                        location+=j*cptSizePrefix[k];
                        continue;
                    }
                    else{
                        vector<string> values = Alarm.Pres_Graph[Parents[k]].get_values();
                        int index = Alarm.Pres_Graph[Parents[k]].get_index();
                        for(int l=0;l<values.size();l++){
                            if(values[l] == data[i].first[index]){
                                location+=l*cptSizePrefix[k];
                                break;
                            }
                        }
                    }
                }
                vector<float> cpt = it->second.get_CPT();
                double val = cpt[location];
                prob*=val;
                // cout<<prob<<" ";
            }
            // cout<<endl;
            unknownDistribution[i][j] = prob;
            sum+=prob;
        }
        for(int j=0;j<nvalues;j++)
            unknownDistribution[i][j]/=sum;
    }
}

void maximization(network& Alarm, vector<pair<vector<string>,float>> data, vector<int> unknownIndex, vector<vector<double>> unknownDistribution){
    int variables = Alarm.netSize();
    cout<<variables<<endl;
    for(int i=0;i<variables;i++){
        cout<<i<<endl;
        string name = Alarm.indexToName[i];
        int nvalues = Alarm.Pres_Graph[name].get_nvalues();
        vector<string> Parents = Alarm.Pres_Graph[name].get_Parents();
        vector<int> ParentIndices;
        int numParents = Parents.size();
        for(int j=0;j<numParents;j++){
            ParentIndices.push_back(Alarm.Pres_Graph[Parents[j]].get_index());
        }
        cout<<"B"<<endl;
        sort(ParentIndices.begin(),ParentIndices.end());
        vector<int> cptSizePrefix(numParents,nvalues);
        for(int j=numParents-2;j>=0;j--){
            int numValues = Alarm.Pres_Graph[Parents[j+1]].get_nvalues();
            cptSizePrefix[j] = cptSizePrefix[j+1]*numValues;
        }
        cout<<"C"<<endl;
        long size;
        if(numParents>0)
            size = Alarm.Pres_Graph[Parents[0]].get_nvalues()*cptSizePrefix[0];
        else
            size = nvalues;
        vector<double> numEvents(size,0);
        cout<<data.size()<<endl;
        for(int j=0;j<data.size();j++){
            cout<<j<<endl;
            long location = 0;
            bool flag = false;
            for(int k=0;k<numParents;k++){
                if (unknownIndex[j] != ParentIndices[k]){
                    vector<string> values = Alarm.Pres_Graph[Alarm.indexToName[ParentIndices[k]]].get_values();
                    for(int l=0;l<values.size();l++){
                        if(values[l] == data[j].first[ParentIndices[k]]){
                            location += l*cptSizePrefix[k];
                            break;
                        }
                    }
                }
                else
                    flag = true;
            }
            cout<<"D1"<<endl;
            if(i != unknownIndex[j]){
                vector<string> values = Alarm.Pres_Graph[name].get_values();
                cout<<"E1"<<endl;
                for(int l=0;l<values.size();l++){
                    if(values[l] == data[j].first[i]){
                        location += l;
                        break;
                    }
                }
                cout<<"E2"<<endl;
            }
            else
                flag = true;
            if(flag == false)
                numEvents[location]+=1;
            else{
                cout<<"F1"<<endl;
                for(int k=0;k<numParents;k++){
                    if (unknownIndex[j] == ParentIndices[k]){
                        int nvalues = Alarm.Pres_Graph[Alarm.indexToName[ParentIndices[k]]].get_nvalues();
                        cout<<"H1"<<endl;
                        for(int l=0;l<nvalues;l++){
                            long loc = location+l*cptSizePrefix[ParentIndices[k]];
                            numEvents[loc]+=unknownDistribution[j][l];
                        }
                        cout<<"H2"<<endl;
                        break;
                    }
                }
                cout<<"F2"<<endl;
                if(i == unknownIndex[j]){
                    int nvalues = Alarm.Pres_Graph[name].get_nvalues();
                    cout<<"G1"<<endl;
                    for(int l=0;l<nvalues;l++){
                        long loc = location+l;
                        numEvents[loc]+=unknownDistribution[j][l];
                    }
                    cout<<"G2"<<endl;
                }
            }
            cout<<"D2"<<endl;
        }
        cout<<"E"<<endl;
        double sum = 0;
        vector<float> CPT(numEvents.size(),0);
        for(int j=0;j<numEvents.size();j++){
            sum+=numEvents[j];
            if (j%nvalues==nvalues-1){
                for(int k=j-nvalues+1;k<=j;k++){
                    CPT[k] = numEvents[k]/sum;
                }
            }
        }
        cout<<"F"<<endl;
        Alarm.Pres_Graph[name].set_CPT(CPT);
    }
}


int main()
{
	network Alarm;
	Alarm=read_network();
    int variables = Alarm.netSize();
    vector<pair<vector<string>,float>> data;
    vector<int> unknownIndex;
    ifstream dataFile("records.dat", ios::in);

    string line;
    bool mode = false;
    while(getline(dataFile, line)){
        int len = line.length();
        string val = "";
        vector<string> vals;
        int count = 0;
        for(int i=0;i<len;i++){
            if (line[i] == '\"'){
                mode = !mode;
                if (mode == false){
                    vals.push_back(val);
                    if (val=="?")
                        unknownIndex.push_back(count);
                    val = "";
                    count++;
                }
            }
            if (mode == true){
                if(line[i]!='\"')
                    val+=line[i];
            }
        }
        data.push_back(make_pair(vals,1));
    }

    vector<vector<double>> unknownDistribution;
    for(int i=0;i<data.size();i++){
        int nvalues = Alarm.Pres_Graph[Alarm.indexToName[unknownIndex[i]]].get_nvalues();
        vector<double> temp(nvalues,0);
        unknownDistribution.push_back(temp);
    }

    for(auto it = Alarm.Pres_Graph.begin(); it != Alarm.Pres_Graph.end(); it++){
        int parents = it->second.get_Parents().size();
        int tableSize = pow(2,parents+1);
        int values = it->second.get_nvalues();
        vector<float> initCPT(tableSize, 1/values);
        it->second.set_CPT(initCPT);
    }

    for(int i=0;i<1;i++){
        expectation(Alarm, data, unknownIndex, unknownDistribution);
        cout<<"A"<<endl;
        maximization(Alarm, data, unknownIndex, unknownDistribution);
    }
    // for(int i=0;i<unknownDistribution.size();i++){
    //     for(int j=0;j<unknownDistribution[i].size();j++)
    //         cout<<unknownDistribution[i][j]<<" ";
    //     cout<<endl;
    // }

    
// Example: to do something
	cout<<"Perfect! Hurrah! \n";
	
}