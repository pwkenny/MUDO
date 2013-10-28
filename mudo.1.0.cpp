//  MUDO (MolecUle eDitOr)
//  Peter W. Kenny, IQSC-USP
//  
//  MUDO is the successor to the Leatherface molecule editor [Kenny & Sadowski (2005) Structure modification in
//  in chemical databases. Methods and principles in medicinal chemistry. In: Oprea T (ed) Chemoinformatics in 
//  drug discovery, vol 23, pp 271–285. doi: http://dx.doi.org/10.1002/3527603743.ch11].  A manuscript 
//  (Automated molecule edting in molecular design) has been submitted to JCAMD that describes MUDO and how
//   it can be applied in molecular design (structure standardization; setting tautomer and protonation states
//   for virtual screening; identification of matched molecular pairs). Structural transforms are applied using 
//  SMIRKS notation and the software has been created using OpenEye's OEChem toolkit.


#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <sstream>
#include <string>
#include <algorithm>
#include <iterator>
#include "openeye.h"
#include "oechem.h"


using namespace OEChem;
using namespace OESystem;
using namespace std;



const char *InterfaceData = 
"!BRIEF UsingOEInterfaceHelp [-o] <output> [-i] <input> [-s] <smirks> [-v] <vbind> [-m] <mode> [-e] <exhaustive> [-c] <canonical>  \n"
"!PARAMETER -input\n"
"  !ALIAS -i\n"
"  !TYPE string\n"
"  !REQUIRED true\n"
"  !BRIEF Input file\n"
"  !DETAIL\n"
"Input file of molecules (suffix will determine file type)\n"
"!END\n"
"!PARAMETER -output\n"
"  !ALIAS -o\n"
"  !TYPE string\n"
"  !REQUIRED true\n"
"  !BRIEF The output file\n"
"  !DETAIL\n"
"Output file of molecules (suffix will determine file type)\n"
"!END\n"
"!PARAMETER -smirks\n"
"  !ALIAS -s\n"
"  !TYPE string\n"
"  !REQUIRED true\n"
"  !BRIEF SMIRKS\n"
"  !DETAIL\n"
"SMIRKS definitions \n"
"!END\n"
"!PARAMETER -vectorbind\n"
"  !ALIAS -v\n"
"  !TYPE string\n"
"  !REQUIRED false\n"
"  !BRIEF Vector bindings\n"
"  !DETAIL\n"
"Optional file vector bindings\n"
"!END\n"
"!PARAMETER -mode\n"
"  !ALIAS -m\n"
"  !TYPE string\n"
"  !REQUIRED true\n"
"  !BRIEF Mode\n"
"  !DETAIL\n"
"Mode of operation (options: normal, link or enum ) \n"
"!END\n"
"!PARAMETER -exhaustive\n"
"  !ALIAS -e\n"
"  !TYPE bool\n"
"  !DEFAULT false\n"
"  !REQUIRED false\n"
"  !BRIEF Exhaustive\n"
"  !DETAIL\n"
"Determines whether enumeration is exhaustive (true or false)\n"
"!END\n"
"!PARAMETER -canonical\n"
"  !ALIAS -c\n"
"  !DEFAULT false\n"
"  !TYPE bool\n"
"  !REQUIRED false\n"
"  !BRIEF Canonical\n"
"  !DETAIL\n"
"Only output canonical form when enumerating (i.e. for duplicate recognition)\n"
"!END\n"
"!PARAMETER -debug\n"
"  !ALIAS -d\n"
"  !TYPE bool\n"
"  !DEFAULT false\n"
"  !REQUIRED false\n"
"  !BRIEF Debug\n"
"  !DETAIL\n"
"Debug flag\n"
"!END";


int main(int argc , char** argv)
{

 OEInterface itf(InterfaceData, argc, argv);

 oemolistream infil;
 oemolostream outfil;
 ifstream smkfil;
 ifstream vbfil;

 OEGraphMol mol, refmol ;
 OEIter<OEMolBase> product;


 OEUniMolecularRxn *lf_edit ;
 OELibraryGen *lf_enum ;

 vector< pair<string, string> > definitions;
 list<string> curr_forms , new_forms , add_lists ;   
 string  initial_smiles , line , mode , smiles , smirksin , title, vbname , vbdef ;
 unsigned int  i , i_mol ,i_smrk , j, mode_norm , mode_enum , n_edit , n_form , n_pass , vbuse ;
 bool canon  , debug , exhaustive;
 
 //unsigned int smiflag = OESMILESFlag::ISOMERIC;



 //  Command line stuff

 mode = itf.Get<std::string>("-m");

 if (itf.Has<bool>("-c")) canon = itf.Get<bool>("-c");
 if (itf.Has<bool>("-e")) exhaustive = itf.Get<bool>("-e");
 if (itf.Has<bool>("-d")) debug = itf.Get<bool>("-d");


 infil.open(itf.Get<std::string>("-i"));
 outfil.open(itf.Get<std::string>("-o"));
 smkfil.open(itf.Get<std::string>("-s").c_str());
  if (itf.Has<string>("-v"))
   {
     vbfil.open(itf.Get<std::string>("-v").c_str());
     if (!vbfil) OEThrow.Error("Unable to open vector binding file\n",smirksin.c_str());
   }
  if (!infil) OEThrow.Error("Unable to open input file\n",smirksin.c_str());
  if (!outfil) OEThrow.Error("Unable to open output file\n",smirksin.c_str());
  if (!smkfil) OEThrow.Error("Unable to open SMIRKS file\n",smirksin.c_str());

  if ( mode == "enum" ) cout << "Mode: " << mode << "\n" ;


 
// Read vector bindings (Thanks, Krysztina, for showing me how to do this)

 if (vbfil)
   {
     vbuse = 1;
     while(getline(vbfil,line))
        {
          if (!vbfil.eof())
              {  
                if ( line.find("#") != 0) 
                    {
                      istringstream iss(line);
                      iss >> vbname >> vbdef; 
                      definitions.push_back(pair<string,string>( vbname, vbdef ));
                    }
               }
         }
   }

// Read SMIRKS file to determine number of SMIRKS patterns

i_smrk = 0;
while(getline(smkfil,line))
   {
     if (!smkfil.eof())
        {  
           if ( line.find("#") != 0) 
              {
                 istringstream iss(line);
                 iss >> smirksin;
                 i_smrk++;
              }
        }
    }


 smkfil.close();
 smkfil.open(itf.Get<std::string>("-s").c_str());
 if ( mode == "enum" || mode == "link" ) 
   { 
     lf_enum = new OELibraryGen [i_smrk];
   }
 else
   {
     lf_edit = new OEUniMolecularRxn [i_smrk];
   }
 
    

// Read SMIRKS file again to get the SMIRKS

i_smrk = 0;
while(getline(smkfil,line))
   {
     if (!smkfil.eof())
        {  
           if ( line.find("#") != 0) 
              {
                 istringstream iss(line);
                 if (vbuse) 
                   {
		    OESmartsLexReplace(line, definitions);
                    smirksin = line;
                   } 
                 smirksin = line;
                 if ( mode == "enum" || mode == "link" )
		   {
                      if (!lf_enum[i_smrk].Init(smirksin.c_str()))
                        OEThrow.Error("Unable to parse %s",smirksin.c_str());
                   }
                 else
		   {
                     if (!lf_edit[i_smrk].Init(smirksin.c_str()))
                        OEThrow.Error("Unable to parse %s",smirksin.c_str());
		   }
                 i_smrk++;
              }
        }
    }

// Edit molecules in normal mode


 if ( mode == "normal" )
   {
       while (OEReadMolecule(infil, mol) )
         {
           n_pass = 0; n_edit = 0;
           while ( ( n_pass == 0 ) || ( ( n_edit > 0 ) && exhaustive ) )
             {
               n_edit = 0; n_pass++;
               for ( i = 0 ; i < i_smrk ; i++ )
                {
                  if (lf_edit[i](mol)) n_edit++ ;
                }
              }
           OEWriteMolecule(outfil, mol);
	 }
   }
 
 // Edit molecules in link mode

 if  ( mode == "link" )
   {
    i_mol = 0;
    while (OEReadMolecule(infil, mol) )
         {
           if ( i_mol == 0 )
	     {
               refmol = mol;
               mol.Clear();
               i_mol++;
	     }
           else
             {
                for ( i = 0 ; i < i_smrk ; i++ )
                    {
		      //lf_enum[i].SetValenceCorrection(false);
                     lf_enum[i].SetStartingMaterial(refmol,0);
                     lf_enum[i].SetStartingMaterial(mol,1);
                     if (! OEHasResidues(mol))
		         {
                            OEPerceiveResidues(mol, OEPreserveResInfo::All);
                         }
                     for (product = lf_enum[i].GetProducts();product;++product)
 		         {
                          if ( debug )
		             cout << "In link loop " <<  "\n";
			  OEWriteMolecule(outfil,product);
		         }
                }
             }       
 	 }
   }

// Edit molecules in enumeration mode

 while (OEReadMolecule(infil, mol) && (mode == "enum"))
    {
      curr_forms.clear();
      n_form = 0 ; n_pass = 0; 
      title = mol.GetTitle();
      OECreateIsoSmiString(smiles,mol);
      initial_smiles = smiles;
      curr_forms.push_front(smiles);
      if ( debug )
          cout << "chek smiles 0: " << smiles << "\n";

      while ( ( n_form != curr_forms.size() && ( ( n_pass < 1 ) || exhaustive ) ) ) 
        {
          if (debug )
              cout << "chek smiles 1: " << smiles << "\n";
	  n_form = curr_forms.size();
          n_pass++;
          for ( i = 0 ; i < i_smrk ; i++ )
            {
              new_forms.clear();
              list<string>::iterator p = curr_forms.begin() ;
              while ( p != curr_forms.end() )
      	      {
                  mol.Clear();
      	          OEParseSmiles(mol,*p);
                  lf_enum[i].SetStartingMaterial(mol,0,false);
                  for (product = lf_enum[i].GetProducts();product;++product)
                   {
                     OECreateIsoSmiString(smiles,product);
                     new_forms.push_back(smiles);
                     if (debug )                 
      		         cout << "product smiles = " << i << " " << smiles << "\n";
                   }
                  if ( debug )
                     cout << "chek " << *p << " " << new_forms.size() << " " << curr_forms.size() <<  "\n" ;
                  p++;
                }
              curr_forms.sort();
              new_forms.sort();
      	      curr_forms.merge(new_forms);
              curr_forms.unique();
              if ( debug )
                   cout << "List length: " << curr_forms.size() << "\n";
              list<string>::iterator q = curr_forms.begin() ;
              while ( q != curr_forms.end() )
      	         {
                    q++;
                  }
              if ( debug )
                 cout << "smirks chek " << i_smrk << " " << i << "\n";
            }
	}

      //  Convert SMILES in enumeration list to molecules and write to output file and only process
      //  the first member of the list if 'canonical' flag is set

     j = 0;
     list<string>::iterator q = curr_forms.begin() ;
     while ( (q != curr_forms.end()) && ( j == 0 ) )
	{
	  if ( canon )
	     j++;
           mol.Clear();
           if ( initial_smiles != *q || exhaustive )
             {
               OEParseSmiles(mol,*q);
               mol.SetTitle(title);
               OEWriteMolecule(outfil, mol);
             }
           q++;
              
        }
    }
  return 0;
}
