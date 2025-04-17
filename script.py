
import pymol2
import nglview as nv
import py3Dmol
import MDAnalysis as mda
import os
import sys
import shutil
import subprocess
import warnings
import pandas as pd
from openbabel import openbabel

from rdkit import Chem
from rdkit.Chem import AllChem, Draw, rdMolTransforms
from rdkit.Chem.Draw import MolsToGridImage, IPythonConsole

aa_dict = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"
}

# den här kommer inte behövas om vi tar bort at man ska generera alla mutationer. 
aa_list = [
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL"
]

# Kanske vi borde skapa en metod för att göra smiles till .mol format

class KEX():

    def __init__(self, filename):
        self.raw_filename = filename
        self.pdb_filenames = []
        self.pdbqt_filenames = []
        self.ligand_filenames = []

        self.create_directories()
        #self.clean_up_dir("pdb")
        #self.clean_up_dir("pdbqt")
        #self.clean_up_dir("docking_results")
        self.clean_up_dir("ligands")
        
        new_filename = self.clean_up()
        self.pdb_filenames.append(new_filename) 
        self.starting_enzyme = self.pdb_filenames[0]

        self.copy_starting_enzyme_into_pdb_dir()


    def create_directories(self):
        self.pdb_dir = os.path.join(os.getcwd(), "pdb")
        os.makedirs(self.pdb_dir, exist_ok = True) 
        
        self.pdbqt_dir = os.path.join(os.getcwd(), "pdbqt")
        os.makedirs(self.pdbqt_dir, exist_ok = True)

        self.ligands_dir = os.path.join(os.getcwd(), "ligands")
        os.makedirs(self.ligands_dir, exist_ok = True)

        self.docking_results_dir = os.path.join(os.getcwd(), "docking_results")
        os.makedirs(self.docking_results_dir, exist_ok = True)

        


    def copy_starting_enzyme_into_pdb_dir(self):        
        source = os.path.join(os.getcwd(), self.starting_enzyme)
        destination = os.path.join(self.pdb_dir, self.starting_enzyme)
        shutil.copy2(source, destination)
        
    
    
    
    def clean_up_dir(self, directory): # kanske run den här i __init__
        
        if directory == "pdb":
            directory = self.pdb_dir
            self.pdb_filenames = []
            
        elif directory == "docking_results":
            directory = self.docking_results_dir

        elif directory == "ligands":
            directory = self.ligands_dir
            self.ligand_filenames = []
            
        elif directory == "pdbqt":
            directory = self.pdbqt_dir
            self.pdbqt_filenames = []
        
        for filename in os.listdir(directory):           
            file_path = os.path.join(directory, filename)
            if os.path.isfile(file_path):
                os.remove(file_path)
        #self.copy_starting_enzyme_into_pdb_dir()
    
    def viz(self): # Saga
        pass

    
    
    def find_molecule_coordinates(self, molecule_name, chain): # Ebba

        with pymol2.PyMOL() as pm:
            pm.cmd.load(self.raw_filename, "protein")

            inhibitor_selection = f"resn {molecule_name} and chain {chain}"
            center = pm.cmd.centerofmass(inhibitor_selection)
            self.docking_center = center # vi får fundera på om vi vill ha kvar den här, nu används return i 
            return center


    
    
    def clean_up(self): # denna innehåller massa extra grejer som är specifika för vår pdb
        clean_filename = self.raw_filename.replace(".pdb", "_clean.pdb")
        with open(self.raw_filename, "r") as infile, open(clean_filename, "w") as outfile:
            for line in infile:
                if line.startswith(("ATOM", "TER")):
                    outfile.write(line) 
                elif line.startswith("HETATM") and "N10" in line:
                    line = "ATOM  " + line[6:]  # Gör ingen skillnad någon skillnad, hela N10 tas bort av pdb2pqr när den inte känner igen N10 
                    #line = line.replace("N10", "SER") # Denna rad gör att pqr filen inte kan skapas eftersom SER inte ser ut så som vi ger den
                    outfile.write(line)
        return clean_filename
        
    def extract_N10_atoms(self, filename):
        N10_lines = []
        with open(filename, 'r') as infile:
            for line in infile:
                if "N10" in line:
                    N10_lines.append(line)
        self.N10_lines_list = N10_lines

    def get_last_atom_line(self, filepath):
        with open(filepath, "r") as f:
            lines = f.readlines()
            for line in reversed(lines):
                if line.startswith(("ATOM", "HETATM")):
                    return line
        return None 

    def append_N10_atoms_to_pdbqt(self, filename):
        with open(f"{self.pdbqt_dir}/{filename}", 'a') as infile:
            last_atom_line = self.get_last_atom_line(f"{self.pdbqt_dir}/{filename}")
            last_id = int(last_atom_line[6:11]) 
            atom_id = last_id+1
            for line in self.N10_lines_list:
                
                atom_name = line[12:16].strip()
                atom_name_formatted = atom_name.ljust(4)
                res_name = "N10"
                chain_id = line[21]
                res_num = line[22:26].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                element = line[76:78].strip() if len(line) >= 78 else atom_name[0]
                infile.write(
    f"ATOM  {atom_id:5d} {atom_name_formatted}{res_name:>4} {chain_id}{int(res_num):4d}    "
    f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00     0.000{element:>2}\n"
)
                atom_id += 1

                
    
    def create_molecule(self, smiles, name, show_structure = False):
        """creates a molecule in .mol from smiles"""
        molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
        AllChem.EmbedMolecule(molecule)
        AllChem.MMFFOptimizeMolecule(molecule)
        # Convert to MOL format
        mol = Chem.MolToMolBlock(molecule)
        # Visualize with py3Dmol
        if show_structure:
            viewer = py3Dmol.view(width=400, height=400)
            viewer.addModel(mol, "mol")
            viewer.setStyle({"stick": {}})
            viewer.zoomTo()
            viewer.show()
        # Save as MOL file
        mol_filename = f"{name}.mol"
        Chem.MolToMolFile(molecule, mol_filename)
        print(f"Molecule saved as {mol_filename}")
    
    
        
    def mutations(self, subunits = 'All', positions = None, mutations = None, label = None): 
        # Fixa så man kan generera alla kombinationer för två residue positioner
        # Snygga till koden
        # Kolla på hur man kan utcekla subinits delen så man kan göra ex en mutation på en subunit och en helt annan mutation på en annan subunit
        
        def apply_mutation():
            pm.cmd.get_wizard().set_mode(mutation)
            pm.cmd.get_wizard().do_select(f"/MutProt//{chain}/{pos}/")
            pm.cmd.get_wizard().apply()              
        
        
        # Kolla upp något snyggt sätt att hantera felaktiga inputs
        if positions == None or mutations == None: 
            return "Some error message"
        
        if type(positions) == int:
            positions = [positions]
        
        if subunits != 'All' and type(subunits) != list:
            return "Some error message"

        
        original_stdout = sys.stdout
        log_file = open("mutations_log.txt", "w")
        sys.stdout = log_file
        
        with pymol2.PyMOL() as pm:
            pm.cmd.load(self.starting_enzyme, "MutProt")
            pm.cmd.wizard("mutagenesis")
            pm.cmd.refresh_wizard()
        
            if subunits == 'All':
                chains = pm.cmd.get_chains("MutProt")
            else:
                chains = subunits 

            if mutations == "All":
                mutations = aa_list
            
            new_filename = [] 

            if len(positions) == len(mutations):
                for i, pos in enumerate(positions):                                        
                    mutation = mutations[i]
                    for chain in chains:
                        model = pm.cmd.get_model(f"/MutProt//{chain}/{pos}")
                        starting_aa = model.atom[0].resn
                        apply_mutation()
                    mutation_str = f"{aa_dict[starting_aa]}{pos}{aa_dict[mutation]}"
                    new_filename.append(mutation_str)
                
                if label is not None:
                    new_filename.append(label)
                new_filename = f"{'_'.join(new_filename)}.pdb"
                self.pdb_filenames.append(new_filename)


                output_path = os.path.join(self.pdb_dir, new_filename)
                pm.cmd.save(output_path, "MutProt")

            else:
                for pos in positions:
                    for mutation in mutations:
                        new_filename = []

                        pm.cmd.delete("MutProt")
                        pm.cmd.load(self.starting_enzyme, "MutProt")
                        pm.cmd.wizard("mutagenesis")
                        pm.cmd.refresh_wizard()
                    
                        for chain in chains:
                            model = pm.cmd.get_model(f"/MutProt//{chain}/{pos}")
                            starting_aa = model.atom[0].resn
                            apply_mutation()
                    
                        
                        if starting_aa == mutation:    
                            continue
                        pm.cmd.set_wizard()
                    
                        mutation_str = f"{aa_dict[starting_aa]}{pos}{aa_dict[mutation]}"
                        new_filename.append(mutation_str)               
                    
                        new_filename = f"{'_'.join(new_filename)}.pdb"
                        self.pdb_filenames.append(new_filename)

                        output_path = os.path.join(self.pdb_dir, new_filename)
                        pm.cmd.save(output_path, "MutProt")
            
        


        sys.stdout = original_stdout  
        log_file.close()
        
           
             

    
    def create_pdbqt(self):
       
        for filename in self.pdb_filenames:
            filename = filename[:-4]
            subprocess.run(f"pdb2pqr30 --keep-chain --with-ph 7.4 --ff=PARSE {self.pdb_dir}/{filename}.pdb {self.pdbqt_dir}/{filename}.pqr -q --log-level INFO") # CRITICAL
            u = mda.Universe(f"{self.pdbqt_dir}/{filename}.pqr")
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", UserWarning)
                na

            # Removes the first two lines, makes the file work in Vina
            with open(f"{self.pdbqt_dir}/{filename}_temp.pdbqt") as rf, open(f"{self.pdbqt_dir}/{filename}.pdbqt", 'w') as wf:
                for i, line in enumerate(rf):
                    if i >= 2:
                        wf.write(line)
       
            os.remove(f"{self.pdbqt_dir}/{filename}_temp.pdbqt")
            os.remove(f"{self.pdbqt_dir}/{filename}.log")
            os.remove(f"{self.pdbqt_dir}/{filename}.pqr")
            self.pdbqt_filenames.append(f"{filename}.pdbqt")

    def create_pdbqt_no_charges(self):  
        for filename in self.pdb_filenames:
            filename = filename[:-4]
            
            with pymol2.PyMOL() as pm:
                pm.cmd.load(f"{self.pdb_dir}/{filename}.pdb", "mol")
                pm.cmd.h_add("mol")
                pm.cmd.save(f"{self.pdb_dir}/{filename}.pdb", "mol")
            
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", UserWarning)
                u = mda.Universe(f"{self.pdb_dir}/{filename}.pdb")
                u.atoms.write(f'{self.pdbqt_dir}/{filename}_temp.pdbqt')
            
            with open(f"{self.pdbqt_dir}/{filename}_temp.pdbqt") as rf, open(f"{self.pdbqt_dir}/{filename}.pdbqt", 'w') as wf:
                for i, line in enumerate(rf):
                    if i >= 2:
                        wf.write(line)
            
            os.remove(f"{self.pdbqt_dir}/{filename}_temp.pdbqt")
            self.pdbqt_filenames.append(f"{filename}.pdbqt")
    
    
    def mol_to_pdbqt(self, filename):
        filename = filename[:-4]
        subprocess.run(f"obabel {filename}.mol -O {filename}.pdbqt -h --partialcharge gasteiger")
        self.ligand_filenames.append(f"{filename}.pdbqt")
        
        source = os.path.join(os.getcwd(), f"{filename}.pdbqt")
        destination = os.path.join(self.ligands_dir, f"{filename}.pdbqt")
        shutil.copy2(source, destination)
        os.remove(f"{filename}.pdbqt")

    def mol_to_pdbqt_new(self, filename, add_ligand = False):
        filename = filename[:-4]
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("mol", "pdbqt")
        
        mol = openbabel.OBMol()
        obConversion.ReadFile(mol, f"{filename}.mol")
        mol.AddHydrogens()

        # Lägg till laddningar
        charge_model = openbabel.OBChargeModel.FindType("gasteiger")
        charge_model.ComputeCharges(mol)
        obConversion.AddOption("h")
        obConversion.WriteFile(mol, f'{filename}.pdbqt')
        obConversion.CloseOutFile()

        if add_ligand:
            self.ligand_filenames.append(f"{filename}.pdbqt")
            source = os.path.join(os.getcwd(), f"{filename}.pdbqt")
            destination = os.path.join(self.ligands_dir, f"{filename}.pdbqt")
            shutil.copy2(source, destination)
            os.remove(f"{filename}.pdbqt")
    
    def add_ligand(self, filename): #lägg in den i funktionen över
        source = os.path.join(os.getcwd(), filename)
        destination = os.path.join(self.ligands_dir, filename)
        shutil.copy2(source, destination)
        # os.remove(filename)
        self.ligand_filenames.append(filename)
    
    def windows_docking(self, center, boxsize = 20): 
        cx = center[0]
        cy = center[1]
        cz = center[2]
        bx = boxsize
        by = boxsize
        bz = boxsize
        
        d = {}
        
        for ligand in self.ligand_filenames:
            results = []
            enzymes = []
            for enzyme in self.pdbqt_filenames:
                config = open('config.txt', mode='w')
                config.write(f"receptor={self.pdbqt_dir}/{enzyme}\n")
                config.write(f"ligand={self.ligands_dir}/{ligand}\n")
                config.write('center_x=')
                config.write(str(cx))
                config.write('\n')
                config.write('center_y=')
                config.write(str(cy))
                config.write('\n')
                config.write('center_z=')
                config.write(str(cz))
                config.write('\n')
                config.write('size_x=')
                config.write(str(bx))
                config.write('\n')
                config.write('size_y=')
                config.write(str(by))
                config.write('\n')
                config.write('size_z=')
                config.write(str(bz))
                config.write('\n')
                config.close()
                new_filename = f"docked_{ligand[:-6]}_in_{enzyme[:-6]}.pdbqt"
                r = subprocess.run(f'"vina.exe" --config config.txt --log log.txt --out {self.docking_results_dir}/{new_filename} --exhaustiveness 20 --num_modes 20 --energy_range 6', capture_output=True, text = True)
                
                #Lös detta bättre?
                os.remove("config.txt")
                os.remove("log.txt")
                
                # spara resultaten från docked filen
                
                enzymes.append(enzyme[:-6])

                with open(f"{self.docking_results_dir}/{new_filename}") as rf:
                    next(rf)
                    res = float(rf.readline().split()[3])
                
                results.append(res)
            
            d[f"{ligand[:-6]} (kcal/mol)"] = pd.Series(results, index = enzymes)
            df = pd.DataFrame(d)
        return df
    
    def os_docking(self):
        pass
        
    