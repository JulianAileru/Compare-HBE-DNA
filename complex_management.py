#1st function used to process pdb file by removing heteroatom,selecting top occupancy residues, and deleting all heteroatoms : writes processed pdb file 
def anatomize(exp="XRAY",pdb="pdb_file",directory="/Users/julianaileru/GuoLab/complexes_05_21_24"):
    import os 
    import subprocess
    try:
        remove_insertion = ["pdb_fixinsert", f"{os.path.join(directory,pdb)}"] # delete all insertion codes 
        process2 = subprocess.Popen(
        remove_insertion, stdout=subprocess.PIPE)
        delhetatms = ["pdb_delhetatm"] #removes all HETATM
        process3 = subprocess.Popen(delhetatms, stdin=process2.stdout, stdout=subprocess.PIPE)
        top_occupancy = ["pdb_selaltloc"] #selects the label with the highest occupancy value for each atom 
        process4 = subprocess.Popen(top_occupancy, stdin=process3.stdout, stdout=subprocess.PIPE)
        with open(f"{pdb[:-4]}_processed.pdb", "w") as f1:
            for line in process4.stdout:
                f1.write(line.decode("utf-8"))
        return f'{f1.name}'
    except:
        raise ValueError("pdb-tools by haddocking required")
#2nd : Optional Function to map residue id of processed_pdb file back to original file. Works by chain only, this is okay since were are only interested in 1 chain at a time. 
def mapping_insertion_code(pdb_file, processed_pdb_file, chain_id="A"):
    from Bio.PDB import PDBParser
    import pandas as pd
    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', None)
    pd.set_option('display.max_colwidth', 1)

    original_pdb = {}
    modified_pdb = {}
    original_pdb_2 = {}
    modified_pdb_2 = {}
    letter = "ABCDEGHIJKLMNOPQRSTUWXYZ"
    p = PDBParser()
    structure = p.get_structure("original_file", f"{pdb_file}")
    for model in structure:
        for chain in model:
            if chain.get_id() == f"{chain_id}":
                for residue in chain:
                    #print(residue.resname)
                    for atom in residue:
                        original_pdb[atom.get_serial_number()] = (residue.id[1:],residue.resname,atom.name)
    structure2 = p.get_structure("processed_file",f"{processed_pdb_file}")
    for model in structure2:
            for chain_iter in model:  # Changed variable name to avoid conflict
                if chain_iter.get_id() == f"{chain_id}":  # Changed variable name
                    for residue in chain_iter:
                        for atom in residue:
                            modified_pdb[atom.get_serial_number()] = (residue.id[1:])

    for key, value in modified_pdb.items():
            value = value[0]
            modified_pdb_2[key] = value

    for key,value in original_pdb.items():
        res_id,res_name,atom_name = value
        res_num,icode=res_id
        x = str(atom_name)+" "+str(res_name)+" "+str(res_num)+str(icode)
        original_pdb_2[key] = x
    for key, value in modified_pdb.items():
            value = value[0]
            modified_pdb_2[key] = value

    original_to_mod = {(x, "".join(original_pdb_2[x])): modified_pdb_2[x] for x in modified_pdb_2}
    df = pd.DataFrame(original_to_mod.values(),index=original_to_mod.keys())
    df.reset_index(inplace=True)
    df.drop(columns=["level_0"], inplace=True)
    df.rename(columns={0:"new_resi"},inplace=True)
    df[["ATOM","res","resi",""]] = df.level_1.str.split(" ",expand=True)
    df.drop(columns = ["level_1"])
    df2 = df[["res","resi","new_resi"]]
    df2 = df2.drop_duplicates()
    df2.rename(columns={"res":"RESIDUE","resi":"RESI","new_resi":"NEW_RESI"},inplace=True)
    return df2

#3rd function, runs reduced on processed file, add hydrogens to XRAY file followed by renumbering of atoms. 
def reduce(processed_pdb_file,path_to_reduce="/projects/cindel/jaileru/FIRST-6.2.1-bin-32-gcc3.4.3/reduce.3.14.080821.linuxi386",renumber=True):
    import subprocess
    import os
    if renumber:
        process1 = subprocess.run([path_to_reduce, "-BUILD",processed_pdb_file],stdout=subprocess.PIPE,text=True)
        with open(f"{processed_pdb_file[:-14]}_reduced.pdb","w") as f1:
            print(f"writing reduced file : {f1.name}")
            f1.write(process1.stdout)
        process2 = subprocess.run(["pdb_reatom","-1",f"{processed_pdb_file[:-14]}_reduced.pdb"],stdout=subprocess.PIPE,text=True)
        with open(f"{processed_pdb_file[:-14]}_Rrenumbered.pdb","w") as f2:
            print(f"writing renumbered file : {f2.name}")
            f2.write(process2.stdout)
    else: #to be built if you did not want to renumber atoms. 
        pass
    return f"{processed_pdb_file[:-14]}_Rrenumbered.pdb"

#4th function: run FIRST on reduce/renumbered file : generate hydrogen bond enegeries between different atoms (inidcated by atom number)
def FIRST(Rrenumber_file,HBE_file_name,first_path = "/projects/cindel/jaileru/FIRST-6.2.1-bin-32-gcc3.4.3/FIRST",lib_path="/projects/cindel/jaileru/FIRST-6.2.1-bin-32-gcc3.4.3/"):
    import subprocess
    import os
    FIRST_cmd = [first_path,"-L",lib_path,"-hbout","-E","-0.5","-non",Rrenumber_file]
    p = subprocess.run(FIRST_cmd)
    subprocess.run(["mv","hbonds.out",f'{HBE_file_name}.out'])
    print(f"hydrogen bond energies written to file: {HBE_file_name}.out")
    return f'{HBE_file_name}.out'

#5th function: parse output of first to including information about hydrogen bond from pdb file 
def parse_HBE(hb_outfile,Rrenumbered_file):
    import pandas as pd
    import os 
    import warnings
    warnings.filterwarnings("ignore",category=FutureWarning)
    pdb_id = str(Rrenumbered_file[:-16])
    ATOM_lines = []
    hb_lines = []
    with open(Rrenumbered_file,"r") as f1:
        for line in f1:
            if line.startswith("ATOM"):
                line = line.split()
                ATOM_lines.append([f'{pdb_id}',line[4],line[3],line[1],line[2]])
    with open(hb_outfile,"r") as f2:
        for line in f2:
            line = line.split()
            hb_lines.append(line)
    hb_df = pd.DataFrame(hb_lines)
    hb_df = hb_df.rename(columns={0:"hb_atom1", 1:"hb_atom2",2:"HBE",3:"unknown"})
    ATOM_df = pd.DataFrame(ATOM_lines)
    ATOM_df = ATOM_df.rename(columns = {0:"pdb_id",1:"chain_id",2:"residue_id",3:"atom_number",4:"atom_id"})
    atom_dict = {}
    hb_dict = {}
    atom_number_atom_properties = {}
    for index,row in ATOM_df.iterrows():
    #print(row)
        pdb_id = row[0]
        chain_id = row[1]
        residue_id = row[2]
        atom_number = row[3]
        atom_id = row[4]
        atom_dict[atom_number] = [atom_number,pdb_id,chain_id,residue_id,atom_id]
        atom_number_atom_properties[atom_number] = [atom_id,atom_number]
    #print(atom_id)
    with open(f"HB_{pdb_id}_parsed.out","w") as f1:
        for index1,row1 in hb_df.iterrows():
            hb_atom1 = row1[0]
            hb_atom2 = row1[1]
            HBE = row1[2]
            hb_dict[hb_atom1,hb_atom2] = [hb_atom1,hb_atom2,HBE]
        for key,value in hb_dict.items():
            for atom_number,value1 in atom_dict.items():
                if atom_number == key[0]:
                    f1.write(" ".join(atom_dict[atom_number]) + " " + " ".join(atom_dict[hb_dict[key][1]]) + " " + hb_dict[key][2] + "\n")
    return f"HB_{pdb_id}_parsed.out"
    
#6th function: subdivides HBE information into different categories of interest 
def filter_file(file="/Users/julianaileru/GuoLab/parsed_HB_1b8i_processed_reduced_renumbered.out",protein_chain="A"):
    import pandas as pd 
    #Inital Information
    chains = []
    dna_chains = []
    dna_residues = ["DT","DA","DC","DG"]
    protein_residues = ["ALA","ARG","ASN","ASP","CYS","GLU","GLN","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"]
    dna_bases = ['N1', 'N2', 'N3',  'N4', 'N6', 'N7', 'N9','C2', 'C4', 'C5', 'C6', 'C7', 'C8','O2', 'O3', 'O4', 'O5','H1', 'H2', 'H3', 'H5', 'H8', 'H21', 'H22', 'H41', 'H42','H6', 'H61', 'H62', 'H71', 'H72', 'H73']
    protein_atoms = []
    dna_atoms = [] 
    with open(f"{file}","r") as f1:
        for line in f1:
            line = line.strip()
            line = line.split(" ")
            atom1,pdb1,chain1,residue1,atomid1,atom2,pdb2,chain2,residue2,atomid2,HBE = line
            
            #Get Chain IDs 
            chains.append(chain1)
            chains.append(chain2)
            #Get DNA Chain IDs
            if residue1 in dna_residues:
                dna_chains.append(chain1)
            if residue2 in dna_residues:
                dna_chains.append(chain2)
            #Get DNA atoms
            for i in set(dna_chains):
                if (chain1 == i and residue1 in dna_residues):
                    dna_atoms.append(atomid1)
                if (chain2 == i and residue2 in dna_residues):
                    dna_atoms.append(atomid2)
            #Get Protein atoms 
            if (chain1 == protein_chain and residue1 in protein_residues):
                protein_atoms.append(atomid1)
            if (chain2 == protein_chain and residue2 in protein_residues):
                protein_atoms.append(atomid2)
    #Unique Chains
    chains = set(chains) 
    #Unique Atoms
    dna_atoms = set(dna_atoms)
    protein_atoms = set(protein_atoms)
    protein_backbone = ['N', 'CA', 'C', 'O', 'H', 'H1', 'H2', 'H3', 'HA']
    dna_backbone = [x.strip() for x in dna_atoms if x not in dna_bases]
    protein_sidechain = [x for x in protein_atoms if x not in protein_backbone]
    
    #Filter for Interface
    Interface = []
    BB_BB = []
    SC_base = []
    Mixed = []

    with open(f"{file}","r") as f1:
        for line in f1:
            line = line.strip()
            line = line.split(" ")
            atom1,pdb1,chain1,residue1,atomid1,atom2,pdb2,chain2,residue2,atomid2,HBE = line
            #Filter for Interface ('protein-DNA')
            if (chain1 == protein_chain and chain2 in dna_chains): 
                Interface.append(line)
                #Filter for BB-BB ('protein-DNA')
                if (atomid1 in protein_backbone) and (atomid2 in dna_backbone):
                    BB_BB.append(line)
                    #print(line)
                #Filter for SC-base ('protein-DNA')
                if (atomid1 in protein_sidechain) and (atomid2 in dna_bases):
                    #print(line)
                    SC_base.append(line)
                #Filter for Mixed ('protein-DNA')
                if (atomid1 in protein_sidechain) and (atomid2 in dna_backbone): #SC-BB
                    Mixed.append(line)
                    #print(line)
                if (atomid1 in protein_backbone) and (atomid2 in dna_bases): #BB-Base
                    Mixed.append(line)
                    #print(line)
                
                
            #Filter for Interface ('DNA-protein')
            if (chain1 in dna_chains and chain2 == protein_chain): 
                Interface.append(line)
                #Filter for BB-BB ('DNA-protein')
                if (atomid1 in dna_backbone) and (atomid2 in protein_backbone):
                    BB_BB.append(line)
                    #print(line)
                #Filter for SC-base ('DNA-protein')
                if (atomid1 in dna_bases) and (atomid2 in protein_sidechain):
                    #print(line)
                    SC_base.append(line)
                #Filter for Mixed ('DNA-protein') -- BB-SC base-BB
                if (atomid1 in dna_backbone) and (atomid2 in protein_sidechain): #BB-SC
                    Mixed.append(line)
                if (atomid1 in dna_bases) and (atomid2 in protein_backbone): #base-BB
                    Mixed.append(line)
    BB_BB_df = pd.DataFrame(BB_BB,columns=["atom1","pdb1","chain1","residue1","atomid1","atom2","pdb2","chain2","residue2","atomid2","HBE"])
    SC_base_df = pd.DataFrame(SC_base,columns=["atom1","pdb1","chain1","residue1","atomid1","atom2","pdb2","chain2","residue2","atomid2","HBE"])
    Mixed_df = pd.DataFrame(Mixed,columns=["atom1","pdb1","chain1","residue1","atomid1","atom2","pdb2","chain2","residue2","atomid2","HBE"])
    Interface_df = pd.DataFrame(Interface,columns=["atom1","pdb1","chain1","residue1","atomid1","atom2","pdb2","chain2","residue2","atomid2","HBE"])
    #return f'DNA Chains: {",".join(list(dna_chains))} | Protein Chain selected: {protein_chain}'
    print("returned BB_BB_df","SC_base_df","Mixed_df","Interface_df")
    x = [BB_BB_df,SC_base_df,Mixed_df,Interface_df]
    return x

def graph_HBE(directory):
    import os
    import subprocess
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    import warnings
    import os 
    warnings.simplefilter(action='ignore', category=FutureWarning)
    if not os.path.isdir("HBE_figures"):
        subprocess.run(["mkdir","HBE_figures"],check=True)
    else:
        pass
        #print("files written to ...")
    BB_BB = []
    Mixed = []
    SC_base = []
    Interface = []
    directory = os.listdir()
    directory = [x for x in directory if x.endswith("_HBE_sort")]

    for subdirectory in directory:
        lst = os.listdir(subdirectory)
        for csv in lst:
            if csv == "BB_BB.csv":
                temp1 = pd.read_csv(f'{subdirectory}/{csv}')
                BB_BB.append(temp1)
            if csv == "Mixed.csv":
                temp2 = pd.read_csv(f'{subdirectory}/{csv}')
                Mixed.append(temp2)
            if csv == "SC_base.csv":
                temp3 = pd.read_csv(f'{subdirectory}/{csv}')
                SC_base.append(temp3)
            if csv == "Interface.csv":
                temp4 = pd.read_csv(f'{subdirectory}/{csv}')
                Interface.append(temp4)
    BB_BB_df = pd.concat(BB_BB,axis=0,ignore_index=True)
    SC_base_df = pd.concat(SC_base,axis=0,ignore_index=True)
    Mixed_df = pd.concat(Mixed,axis=0,ignore_index=True)
    Interface_df = pd.concat(Interface,axis=0,ignore_index=True)
    #HBE categories (Malik and Guo 2022)
    dfs = [BB_BB_df, SC_base_df, Mixed_df,Interface_df]
    df_name = ["BB_BB_df","SC_base_df","Mixed_df","Interface_df"]
    for i,j in zip(dfs,df_name):
        i.attrs = {"name":j}
    for df in dfs[:-1]:
        #print(df.attrs.get("name"))
        df["Category"] = np.nan
        df.loc[(df["HBE"].astype(float) < -0.1) & (df["HBE"].astype(float) > -0.6), "Category"] = "I"
        df.loc[(df["HBE"].astype(float) < -0.6) & (df["HBE"].astype(float) > -1.0), "Category"] = "II"
        df.loc[(df["HBE"].astype(float) < -1.0) & (df["HBE"].astype(float) > -1.5), "Category"] = "III"
        df.loc[(df["HBE"].astype(float) < -1.5), "Category"] = "IV"
        df.dropna(subset=["Category"], inplace=True)
    for i in range(len(dfs[:-1])):
        df = dfs[i]
        df_name = df.attrs.get("name")
        I = sum(df[df["Category"] == "I"].value_counts().values)
        II = sum(df[df["Category"] == "II"].value_counts().values)
        III = sum(df[df["Category"] == "III"].value_counts().values)
        IV = sum(df[df["Category"] == "IV"].value_counts().values)
        category_total = sum([I,II,III,IV])
        fig, ax = plt.subplots()
        ax.bar(x=["I", "II", "III", "IV"], height=[(I / category_total)*100, (II / category_total)*100, (III / category_total)*100, (IV / category_total)*100])
        ax.set_title(f'{df_name[:-3]}')
        ax.set_xlabel("HBE Categories")
        ax.set_ylabel("Percentage of Hydrogen Bonds")
        plt.savefig(f'HBE_figures/{df_name[:-3]}')
   
def organize(directory):
    import os
    import subprocess
    print("moving files to correct directories...")
    directory_path = os.path.abspath(directory)
    directory = os.listdir(directory)
    subdirectory = [x for x in directory if os.path.isdir(x) and x.endswith("_HBE_sort")]
    for subdir in subdirectory:
        subdir_path = os.path.join(directory_path,subdir)
        items = [x for x in directory if x.startswith(subdir[:4]) and os.path.isfile(x)]
        items.append("".join([x for x in directory if x.startswith(f'HB_{subdir[:4]}')]))
        for i in items:
            subprocess.run(["mv",f'{os.path.join(directory_path,i)}',subdir_path])
    print("done")
