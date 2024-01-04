//
// File: transposeBranchLengthsByPercentile.cpp
// Created by: Bastien Boussau
// Created on: Oct Mon 24 18:50 2005
//
//COMPILATION: g++ -lbpp-core -lbpp-phyl transposeBranchLengthsByPercentile.cpp -o transposeBranchLengthsByPercentile 


/*
Copyright or Copr. CNRS

This software is a computer program whose purpose is to simulate sequence
data according to a phylogenetic tree and an evolutionary model.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/

// From the STL:
#include <iostream>
#include <fstream>
#include <iomanip>
#include <set>
#include <map>


using namespace std;

// From PhylLib:
#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Phyl/TreeTemplateTools.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Io/Newick.h>


// From Utils:
#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Numeric/Random/RandomTools.h>

using namespace bpp;

/**
 * @brief Read trees from an input file, without segment annotations.
 */
void readTrees(ifstream& file, vector<Tree*>& trees, vector<double>& pos) throw (Exception)
{
  string line = "";
  double begin, end;
  string::size_type index1, index2, index3;
  double previousPos = 0;
  pos.push_back(0);
  string newickStr;
  while(!file.eof())
  {
    string tmp = TextTools::removeSurroundingWhiteSpaces(FileTools::getNextLine(file));
    if(tmp.size() == 0 || tmp.substr(0, 1) == "#") continue;
    line += tmp;
        
    index1 = line.find_first_of(" \t");
    index2 = line.find_first_of(" \t", index1 + 1);
    begin  = TextTools::toDouble(line.substr(0, index1));
    end    = TextTools::toDouble(line.substr(index1 + 1, index2 - index1 - 1));
    index3 = line.find_first_of(";", index2 + 1);
    try {
    while(index3 == string::npos)
    {
      if(file.eof()) throw Exception("Error when parsing tree file: incomplete tree.");
      line += FileTools::getNextLine(file);
      index3 = line.find_first_of(";", index3);
    }
    newickStr = line.substr(index2 + 1, index3 - index2);
    }
    catch (exception &e) {
      std::cout<<e.what()<<std::endl;
    }
    try {
        //std::cout << newickStr <<std::endl;
    TreeTemplate<Node>* t = TreeTemplateTools::parenthesisToTree(newickStr, true);
        //std::cout << TreeTemplateTools::treeToParenthesis(*t, true) <<std::endl;

       trees.push_back(t);
    }
    catch(exception &e) {
      std::cout<<e.what()<<std::endl;
    }
    line = line.substr(index3 + 1);
  }
}

void help()
{
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
  (*ApplicationTools::message << "transposeBranchLengthsByPercentile parameter1_name=parameter1_value"    ).endLine();
  (*ApplicationTools::message << "      parameter2_name=parameter2_value ... param=option_file").endLine();
  (*ApplicationTools::message).endLine();
  (*ApplicationTools::message << "  Refer to the Bio++ Program Suite Manual for a list of available options.").endLine();
  (*ApplicationTools::message << "ex: ./transposeBranchLengthsByPercentile input.tree.file=treefile output.file=outname").endLine();
  (*ApplicationTools::message << "__________________________________________________________________________" ).endLine();
}

int main(int args, char ** argv)
{   
  cout << "******************************************************************" << endl;
  cout << "*            Bio++ Branch lengths transposer, version 1.1.0      *" << endl;
  cout << "* Author: B. Boussau                                             *" << endl;
  cout << "*                                           Last Modif. 08/08/09 *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;
 

  if(args == 1)
  {
    help();
    return 0;
  }
  
  try {
      vector <double> positions;
  	BppApplication transposeBranchLengthsByPercentile(args, argv, "transposeBranchLengthsByPercentile");
  	transposeBranchLengthsByPercentile.startTimer();

	vector<TreeTemplate<Node>*> trees;

    //TreeTemplate<Node>* tree = new TreeTemplate<Node> (*PhylogeneticsApplicationTools::getTree(transposeBranchLengthsByPercentile.getParams()));
    
    string treesPath = ApplicationTools::getAFilePath("input.trees.file", transposeBranchLengthsByPercentile.getParams(), false, true);

    ApplicationTools::displayResult("Trees file", treesPath);
    ifstream treesFile(treesPath.c_str(), ios::in);
    vector<Tree*> trees2;
    readTrees(treesFile, trees2, positions);
    for (int i = 0 ; i< trees2.size() ; i++) 
    {
      trees.push_back(new TreeTemplate<Node> (*trees2[i]));
      delete trees2[i];
    }

      std::map<std::string, std::vector < double> > nameToDistances;
      std::set< double > internalBranches ; 

      std::vector< std::string > names = trees[0]->getLeavesNames();
      for (int i = 0 ; i< names.size() ; i++) 
      {
          nameToDistances[names[i]] = std::vector<double>();
      }
      
    for (int i = 0 ; i< trees.size() ; i++) 
    {
      vector <int> nodeIds = trees[i]->getNodesId();
      for (int j = 0 ; j < nodeIds.size() ; j++) 
        {
        if (trees[i]->getNode(nodeIds[j])->hasFather()) 
          {
          //out<< trees[i]->getNode(nodeIds[j])->getDistanceToFather();
          if (trees[i]->getNode(nodeIds[j])->isLeaf()) 
            {
                nameToDistances[trees[i]->getNode(nodeIds[j])->getName()].push_back(trees[i]->getNode(nodeIds[j])->getDistanceToFather());
	      //out << "\tno\t" << trees[i]->getNode(nodeIds[j])->getName() <<endl; 
          }
          else 
            {
                internalBranches.insert(trees[i]->getNode(nodeIds[j])->getDistanceToFather());
           // out << "\tyes\tNA" <<endl; 
          }
        }

        }
      
    }
      
//////////////
      std::vector< double > internalBranchesVec ; 
      set<double>::iterator it = internalBranches.begin();
      set<double>::iterator end = internalBranches.end();
      for ( ; it != end; ++it ) {
          internalBranchesVec.push_back(*it); //the vector should be ordered...
      }
      
//////////////      
      vector<TreeTemplate<Node>*> treesToChange;

      string treesPath2 = ApplicationTools::getAFilePath("input.trees.to.change.file", transposeBranchLengthsByPercentile.getParams(), false, true);
      
      ApplicationTools::displayResult("Trees to change file", treesPath2);
      ifstream treesFile2(treesPath2.c_str(), ios::in);
      vector<Tree*> trees3;
      readTrees(treesFile2, trees3, positions);
      for (int i = 0 ; i< trees3.size() ; i++) 
      {
          treesToChange.push_back(new TreeTemplate<Node> (*trees3[i]));
          delete trees3[i];
      }
          
      string outName = ApplicationTools::getAFilePath("output.file", transposeBranchLengthsByPercentile.getParams(), true, false);
      
      std::ofstream out (outName.c_str(), std::ios::out); 
  //    out<< "Length\tInternal\tName"<<endl;

      std::set< double > internalBranchesToChange ;
      
      for (int i = 0 ; i< treesToChange.size() ; i++) 
      {
          vector <int> nodeIds = treesToChange[i]->getNodesId();
          for (int j = 0 ; j < nodeIds.size() ; j++) 
          {
              if (treesToChange[i]->getNode(nodeIds[j])->hasFather()) 
              {
                  if (treesToChange[i]->getNode(nodeIds[j])->isLeaf()) 
                  {
                  }
                  else 
                  {
                      internalBranchesToChange.insert(treesToChange[i]->getNode(nodeIds[j])->getDistanceToFather());
                  }
              }
          }
      }
      //////////////
      std::map< double, double > internalBranchesToChangeMap ; // map between the branch length and its position (from smaller to larger) 
      it = internalBranchesToChange.begin();
      end = internalBranchesToChange.end();
      int i=0;
      double siz = (double)internalBranchesToChange.size();
      for ( ; it != end; ++it, ++i ) {
          internalBranchesToChangeMap[*it] = (double) i / siz;
      }      
      //////////////      

      
      
      
      for (int i = 0 ; i< treesToChange.size() ; i++) 
      {
          vector <int> nodeIds = treesToChange[i]->getNodesId();
          for (int j = 0 ; j < nodeIds.size() ; j++) 
          {
              if (treesToChange[i]->getNode(nodeIds[j])->hasFather()) 
              {
                  if (treesToChange[i]->getNode(nodeIds[j])->isLeaf()) 
                  {
                      int randInt = RandomTools::giveIntRandomNumberBetweenZeroAndEntry (nameToDistances[treesToChange[i]->getNode(nodeIds[j])->getName()].size() );
                      treesToChange[i]->getNode(nodeIds[j])->setDistanceToFather( nameToDistances[treesToChange[i]->getNode(nodeIds[j])->getName()][randInt] );
                  }
                  else 
                  {
                      double dist = treesToChange[i]->getNode(nodeIds[j])->getDistanceToFather() ; 
                      double percentile = internalBranchesToChangeMap[dist] ; 
                      int index = (int) (percentile * internalBranchesVec.size());
                      treesToChange[i]->getNode(nodeIds[j])->setDistanceToFather( internalBranchesVec[index] );
                  }
              }
              
          }
          out << TreeTemplateTools::treeToParenthesis(*treesToChange[i], false) ;
      }

      
      
      

      out.close();
     /* 
      vector <double> listOfBranchLengths;

  vector <double> maxBranchLengths;
      for (int i = 0 ; i< trees.size() ; i++) {
        cout <<"Tree "<<i<<endl;
        VectorTools::append(listOfBranchLengths, trees[i]->getBranchLengths());
        double maxi = VectorTools::max(trees[i]->getBranchLengths());
        maxBranchLengths.push_back(maxi);
      }
    vector <string> internal;
           for(unsigned int i = 0; i < tree->getRootNode()->getNumberOfSons(); i++)
             {
                 Vdouble sonBrLen = TreeTemplateTools::getBranchLengths(* tree->getRootNode()->getSon(i));
                 for(unsigned int j = 0; j < sonBrLen.size(); j++) brLen.push_back(sonBrLen[j]);
               }
           return brLen;
    
    
    
  
   string outName = "branchLengths";
              out<< "Length"<<endl;
            for (int j = 0 ; j < listOfBranchLengths.size(); j++) {
              out<< listOfBranchLengths[j]<<endl;
            }
            out.close();
    
    outName = "maxBranchLengths";
    out.open (outName.c_str(), std::ios::out); 
    out<< "maxLength"<<endl;
    for (int j = 0 ; j < maxBranchLengths.size(); j++) {
      out<< maxBranchLengths[j]<<endl;
    }
    out.close();
*/
    
      
      
      
      
  for (unsigned int i = 0; i < trees.size(); i++)
    delete trees[i];

  //  delete tree;
  transposeBranchLengthsByPercentile.done();

  }
  catch (exception& e)
  {
    cout << e.what() << endl;
    return 1;
  }

  return 0;
}

