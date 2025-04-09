import subprocess
with open("pdb_list.txt","r") as f1:
    for line in f1:
        line = line.split(" ")
        with open (f"data/{line[0]}.pdb","w") as f2:
            subprocess.run(["pdb_fetch",line[0]],stdout=f2)


