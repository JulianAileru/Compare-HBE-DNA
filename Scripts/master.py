from scripts.complex_management import anatomize,mapping_insertion_code,reduce,FIRST,parse_HBE,filter_file,graph_HBE,organize
import os
import subprocess
directory = os.listdir()
complexes = {}
with open("pdb_list.txt","r") as f1:
    for line in f1:
        line = line.split(" ")
        complexes[f'{line[0]}.pdb'] = line[1].strip()

for i,j in complexes.items():
    #path = '/Users/julianaileru/GuoLab/complexes_05_21_24/'
    #full_path = os.path.join(path,i)
    process1 = anatomize(pdb=i,directory=f'{os.getcwd()}')
    process2 = reduce(process1)
    process3 = FIRST(process2,HBE_file_name=i[:-4])
    process4 = parse_HBE(hb_outfile=process3,Rrenumbered_file=process2)
    process5 = filter_file(file=process4,protein_chain=j) #specify chain_id default = "A"
    try:
        subprocess.run(["mkdir",f"{i[:-4]}_HBE_sort"],check=True)
        for q,j in zip(process5,["BB_BB","SC_base","Mixed","Interface"]):
            q.to_csv(f"{i[:-4]}_HBE_sort/{j}.csv")
    except:
        print(f"directory: {i[:-4]}_finals already exists, refusing to overwrite files in directory")

graph_HBE(os.getcwd())

organize(os.getcwd())
