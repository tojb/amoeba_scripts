#!/usr/bin/env python3

import numpy as np
import sys

"""
This is a script which I have used to convert a pdb file to a tinker .txyz
file, while keeping atoms in the same order and with the same names as in the
pdb.  I wrote this because the canonical pdbxyz program distributed with tinker
was giving me problems, and I could not debug it.  Probably the current pdbxyz
at the time you are reading this is actually fine for whatever your purposes
are: you should try that first, and then this.  

It has been tested only for proteins, and only for some residues, and only
for the amoebabio18.prm amoeba forcefield.

My intention with this is that if I decide I like using amoeba,
I can keep adding aliases for pdb atom names to amoeba atom types as needed,
 and eventually it should be complete and robust.

Apart from assigning pdb atom names to atom types from the amoeba forcefield file,
the other computational work is to define bonding.

Bonding is done in a two-stage process: we have some expected bonding information 
stored in the hash "preferred nebs", and once the preferred nebs have been 
searched for, remaining neighbours up to the expected valence of the amoeba
atom type are assigned by distance.  If your bonding comes out wrong, then
add more information for the affected atom types to "preferred_nebs", below.

Check the output txyz carefully, and if you need box info you have to add that 
by hand.

Joshua T. Berryman, 2020.

"""


##defaults
infile      = "ilqins_4x2x12_ipqVac.pdb"
param_file  = "amoebabio18.prm"
outFileName = "ilqins_4x2x12_ipqVac_generated.xyz"

##allow replacement of defaults
if len(sys.argv) > 1:
    infile      = sys.argv[1]
if len(sys.argv) > 2:
    param_file  = sys.argv[2]
if len(sys.argv) > 3:
    outFileName = sys.argv[3]


shortNames = {\
"Glycine":       "GLY",\
"Alanine":       "ALA",\
"Valine":        "VAL",\
"Leucine":       "LEU",\
"Isoleucine":    "ILE",\
"Serine":        "SER",\
"Threonine":     "THR",\
"Cysteine":      "CYS",\
"Proline":       "PRO",\
"Phenylalanine": "PHE",\
"Tyrosine":      "TYR",\
"Tryptophan":    "TRP",\
"Histidine":     "HIS",\
"Aspartic Acid": "ASP",\
"Asparagine":    "ASN",\
"Glutamine":     "GLN",\
"Glutamatic Acid": "GLU",\
"Methionine":    "MET",\
"Lysine":        "LYS",\
"Arginine":      "ARG",\
"Ornithine":     "ORN"}
longNames = {}
for l in shortNames.keys():
    longNames[shortNames[l]] = l
  
##put some bonding information here
##so we are not relying only on distances
preferred_nebs = {\
"N"     : {"CA", "H"},
"NH3+"  : {"CA", "H3N+"},
"C"  : {"O", "CA"},
"O"  : {"C"},
"CA" : {"N", "C", "CB", "HA", "HA1", "HA2", "NH3+"},
"CB" : {"CA", "CG", "CG1", "CG2", "HB", "HB1", "HB2"},
"CG" : {"CB", "CD", "CD1", "CD2", "HG", "HG1", "HG2"},
"CG1" : {"CB", "CD", "CD1", "CD2", "HG", "HG1"},
"CG2" : {"CB", "CD", "CD2", "HG", "HG2"},
"HG1" : {"CG1"},
"HG2" : {"CG2"},
"OXT" : {"COO-"},
}




 
##This is my attempt to generate an amoeba xyz file
##based on an (AMBER) pdb file and the amoebabio18.prm amoeba forcefield.

class AmoebaAtRecord:
    def __init__(self, typeId, entryId, mass, valence, Z, atName, at_longName, res_longName):
        self.typeId  = typeId
        self.entryId = entryId
        self.mass    = mass
        self.valence = valence
        self.Z       = Z
        self.atName  = atName
        self.at_longName  = at_longName
        self.res_longName = res_longName
        return
    def __repr__(self):
        s = ""
        for k in self.__dict__:
            s += "%s: %s, " % (k, self.__dict__[k])
        return s

def read_amoebaParm( param_file ):
    f = open(param_file, "r")

    res_by_longName = {}
    alias_biotypes  = {}

    alias_lut = {}

    for line in f:
       L = line.split()
       if len(L) < 1: continue
       ##atom type definition
       if L[0] == "atom": 
          ll = line.split('"')
          l  = len(L)
          entryId = int(L[1])
          typeId  = int(L[2])
          atName  = L[3]
          valence = int(L[l-1])
          mass    = float(L[l-2])
          Z       = int(L[l-3])
          name_string   = ll[1]
          res_longName  = name_string.split()
          res_longName  = " ".join(res_longName[:len(res_longName)-1])
          at_longName   = name_string.split()[-1]
          
          at = AmoebaAtRecord(typeId, entryId,\
                              mass, valence, Z, atName, at_longName, res_longName)
          if res_longName not in res_by_longName:
              res_by_longName[res_longName] = [at]
          else:
              res_by_longName[res_longName].append(at)
          alias_lut[entryId] = at 


       ##aliases of atom dype definitions
       if L[0] == "biotype": 
          ll = line.split('"')
          res_longName  = ll[1]
          at_longName   = L[2]
          assign_type   = int(L[-1])
          try:
              alias_biotypes[res_longName+at_longName] = alias_lut[assign_type]
          except:
              #print("no alias for: ", line)
              alias_biotypes[res_longName+at_longName] = None
 
    return res_by_longName, alias_biotypes

def readPdbFile( pdb_file ):
   f       = open( pdb_file, "r" )
   isTer   = True
   resList = []
   prevRes = -1
   for line in f:
       if "ATOM" in line:

          ##list of candiate atom names, starting with the obvious one 
          atName  = [line[11:16].strip()]
          resName =  line[16:21].strip()

          resId   = int(line[21:30])
          x       = float(line[30:38])        
          y       = float(line[38:46])        
          z       = float(line[46:54])

          if resId != prevRes:
              resList.append([])
              prevRes = resId
          resList[-1].append({"atNames": atName,\
                           "baseAtName": atName[0],\
                              "resName": resName,\
                                "resId": resId,
                                  "crd": np.array([x,y,z])})
   return resList

def match_atom_types(res_by_longName, alias_biotypes, resList_pdb):

   xyzAtTypes = [] 
   for res in resList_pdb:
       if res[0]["resName"] not in longNames:
           print("residue name %s not recognised, add it to the dict at the top of this file." %\
               res[0]["resName"])

       ##build a list of valid residue names that might be in the amoeba parm.
       is_nTer      = False
       is_cTer      = False

       ##This is where we alias amber/pdb atom names for amoeba atom names
       ##mostly just a case of hydrogen naming conventions.
       ##
       ##Amber atom names need to be unique, amoeba not so much.

       ##this can change, so refresh from default at start of each atom.
       res_longName = [longNames[ res[0]["resName"] ]]

       for at in res:

           ##important one: backbone amide hydrogen
           if "H" in at["atNames"]:
               at["atNames"].append("HN")

           ##N-terminus hydrogens
           if "H3" in at["atNames"]: 
               res_longName.append("N-Terminal "+res[0]["resName"])
               res_longName.append("N-Terminal")
               is_nTer = True

           ##C-terminus
           if "OXT" in at["atNames"]: 
               res_longName.append("C-Terminal "+res[0]["resName"])
               is_cTer = True

           if "HB2" in at["atNames"]:
               at["atNames"].append("HB")
           if "HB3" in at["atNames"]:
               at["atNames"].append("HB")
           if "HG21" in at["atNames"]:
               at["atNames"].append("HG2")
           if "HG22" in at["atNames"]:
               at["atNames"].append("HG2")
           if "HG23" in at["atNames"]:
               at["atNames"].append("HG2")
           if "HG12" in at["atNames"]:
               at["atNames"].append("HG1")
           if "HG13" in at["atNames"]:
               at["atNames"].append("HG1")
           if "HD11" in at["atNames"] or "HD12" in at["atNames"] or "HD13" in at["atNames"]:
               at["atNames"].append("HD1")
           if "HD21" in at["atNames"] or "HD22" in at["atNames"] or "HD23" in at["atNames"]:
               at["atNames"].append("HD2")
           if "HE11" in at["atNames"] or "HE12" in at["atNames"] or "HE13" in at["atNames"]:
               at["atNames"].append("HE1")
           if "HE21" in at["atNames"] or "HE22" in at["atNames"] or "HE23" in at["atNames"]:
               at["atNames"].append("HE2")
           if "CD1" in at["atNames"] and "Isoleucine" in res_longName:
               at["atNames"].append("CD")
           if ("HG2" in at["atNames"] or "HG3" in at["atNames"]) and "Glutamine" in res_longName:
               at["atNames"].append("HG")
           if "HD11" in at["atNames"] and "Isoleucine" in res_longName:
               at["atNames"].append("HD")
           if "HD12" in at["atNames"] and "Isoleucine" in res_longName:
               at["atNames"].append("HD")
           if "HD13" in at["atNames"] and "Isoleucine" in res_longName:
               at["atNames"].append("HD")

       if is_nTer:
           for i_at in range(len(res)):
               if "N" in res[i_at]["atNames"]:
                   res[i_at]["atNames"].append("NH3+")
                   res[i_at]["atNames"].remove("N")
               if "H1" in res[i_at]["atNames"] or "H2" in res[i_at]["atNames"] \
                                               or "H3" in res[i_at]["atNames"]:
                   res[i_at]["atNames"].append("H3N+")
       if is_cTer:
           for i_at in range(len(res)):
               if "O" in res[i_at]["atNames"]:
                   res[i_at]["atNames"] = ["OXT"] + res[i_at]["atNames"]
               if "C" in res[i_at]["atNames"]:
                   res_longName.append("C-Terminal")
                   res[i_at]["atNames"] = ["COO-"]

       print("mapping pdb: %s to amoeba:" % res[0]["resName"], res_longName)

       for at in res:
           amoebaType = None
           
           for pdbName in at["atNames"]:

               if amoebaType is not None: break

               ##test for aliases
               if amoebaType is None:
                   for try_res in res_longName:
                       #print("didn't match res, at:", try_res, pdbName, " so trying aliases.")
                       if try_res+pdbName in alias_biotypes:
                           amoebaType = alias_biotypes[try_res+pdbName]
                           break

               ##test for normal atom definitions
               if amoebaType is None:
                   for try_res in res_longName:
                       if try_res in res_by_longName:
                           for at_type in res_by_longName[try_res]:
                               if at_type.at_longName == pdbName:
                                   amoebaType = at_type
                                   break           

           if amoebaType is None:
               print("No match for: %s %s. At this point, need to compare the script and the amoeba parm file for near-misses." % (res_longName, pdbName))
               raise ValueError          
           else:
               print("Match for: %s %s is: " % (res_longName, pdbName), amoebaType)
               xyzAtTypes.append(amoebaType)                      
   return xyzAtTypes

def build_nebKeys(key):
   i, j, k  = key.split("_")
   i = int(i)
   j = int(j)
   k = int(k)
   nebs     = [] 
   for ii in [-1,0,1]:
     for jj in [-1,0,1]:
       for kk in [-1,0,1]:
          nebs.append("%i_%i_%i" % (i+ii,j+jj,k+kk))
   return nebs 

def buildNebLists(xyzAtTypes, resList_pdb):
   i_at  = -1
   rcut  =  8.  ##maximum covalent bond
   rcut2 = rcut*rcut

   ##build a cell list, geometric LUT
   lut      = {}  
   atList   = []
   nebLists = []
   for res in resList_pdb:
      for at in res:
         atList.append(at)        
         i_at += 1
         nebLists.append([])
         atType = xyzAtTypes[i_at]
         crds   = at["crd"]
         key    = "%i_%i_%i" % (int(crds[0]/rcut), int(crds[1]/rcut), int(crds[2]/rcut))
         if key in lut:
             lut[key].append(i_at)
         else:
             lut[key] = [i_at]

   ##make list of nearest neighbours inside rcut
   ##and fill valence from that
   i_at  = -1
   for at in atList: 
      i_at   += 1 
      atType  = xyzAtTypes[i_at]
      atNames = at["atNames"]
      resName = at["resName"]
      resId   = at["resId"]
      crds    = at["crd"]
      valence = atType.valence
      key     = "%i_%i_%i" %\
           (int(crds[0]/rcut), int(crds[1]/rcut), int(crds[2]/rcut))
      nebKeys = build_nebKeys(key)
      nebAts  = []
      for k in nebKeys:
          if k not in lut: continue
          for ii_at in lut[k]:
              if ii_at != i_at:
                  dx  = atList[ii_at]["crd"] - crds
                  dr2 = np.dot(dx, dx)
                  nebAts.append( (ii_at, dr2) )
      nebAts.sort(key = lambda tup:tup[1])
      c_neb = 0

      ##two-pass: check for preferred nebs first
      for i_neb_r in nebAts:
          i_neb      = i_neb_r[0]
          r2_neb     = i_neb_r[1]
          if r2_neb > 8*8:
              break
          if atList[i_neb]["resId"] != resId:
              continue
          nebAtNames = atList[i_neb]["atNames"]
          gotThisOne = False
          for n in nebAtNames:
              for atName_1 in atNames:
                 if atName_1 not in preferred_nebs: continue
                 if n in preferred_nebs[atName_1]:
                    nebLists[i_at].append(i_neb)
                    c_neb += 1
                    gotThisOne = True
                    break
              if gotThisOne:
                 break
          if c_neb >= valence:
              break

      ##second pass, just use up non-assigned near atoms
      if c_neb < valence:
          for i_neb_r in nebAts:
              i_neb      = i_neb_r[0]
              r2_neb     = i_neb_r[1]
              if i_neb not in nebLists[i_at]:
                 nebLists[i_at].append(i_neb)
                 c_neb += 1
                 if c_neb >= valence:
                    break
      nebLists[i_at].sort()

   ##check that every neighbour relation is bi-directional
   print("consistency-checking neb lists")
   for i in range(len(atList)):
      print(atList[i]["resName"], atList[i]["atNames"])
      for i_neb in nebLists[i]:
         print("   -->", atList[i_neb]["atNames"])
         if i not in nebLists[i_neb]:
             print("WARNING atom %i -> %i but %i not -> %i" %\
                    (i, i_neb, i_neb, i))
             print(atList[i])
             for ii in nebLists[i]:
                 print("  -", atList[ii])
             print()
             print(atList[i_neb])
             for ii in nebLists[i_neb]:
                 print("  -", atList[ii])
             raise ValueError
   return nebLists
    
def write_xyzFile(outFileName, resList_pdb, xyzAtTypes, nebLists):
    f = open(outFileName, "w")
    f.write("%i\n" % len(xyzAtTypes))
    i_at = -1
    for res in resList_pdb:
        for at in res:
            i_at += 1
            line  = "%i" % (i_at+1)
            line  = line.rjust(6)
            line  = line + "  "

            ##sticking with the original amber atom name? Good idea?
            line  = line + at["baseAtName"] 
            
            ##xyz coordinates
            line  = line.ljust(13)
            for x in at["crd"][:]:
                crdstr = ("%10.6f" % x).rjust(12)
                line   = line + crdstr

            ##careful: not the overall type, but the individual 
            ##atom name entry, which seems to be over-redundant.            
            typeStr = "%i" % xyzAtTypes[i_at].entryId

            line    = line + typeStr.rjust(6)
            for neb in nebLists[i_at]:
                line = line + (str(neb+1)).rjust(8)
            f.write("%s\n" % line)

res_by_longName, alias_biotypes = read_amoebaParm(param_file)
#for k in res_by_longName:
#   print("======read res %s" % k)
#   for at in res_by_longName[k]:
#        print(at)

#for k in alias_biotypes:
#   print("======read alias %s <-> " % k, alias_biotypes[k])
resList_pdb                     = readPdbFile( infile )
print("read %i residues from %s" % (len(resList_pdb), infile))

##define the atom types, valence etc for the atoms in the pdb
xyzAtTypes = match_atom_types(res_by_longName, alias_biotypes, resList_pdb)

##build the neghbour lists based on valence of the xyz types and 
##coords from the pdb.
print("building neighbour lists to define bonds...")
nebLists   = buildNebLists( xyzAtTypes, resList_pdb ) 

###now ready to write out the xyz file:
write_xyzFile(outFileName, resList_pdb, xyzAtTypes, nebLists)

 
           
           
          
