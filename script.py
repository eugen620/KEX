# imorts
import pymol2
import nglview as nv
import MDAnalysis as mda
import os
import sys
import shutil
import subprocess
from dict_module import aa_dict, aa_list



class KEX():

    def __init__(self, filename):
        self.raw_filename = filename
        self.pdb_filenames = []
        self.pdbqt_filenames = []
        self.ligand_filenames = []

        # kör clean up och här appenda den till self.pdb_filenames

        self.pdb_filenames.append(filename) # *** Tillfällig kod *** vi vill egentligen inte lägga in den råa filen här, vi ska fixa en clean up method och sedan append.
        self.starting_enzyme = self.pdb_filenames[0]
        
        self.create_directories()
        self.copy_starting_enzyme_into_pdb_dir() # Skriv om den här sen så att den tar clean versionen av enzymet som skapas i klassen
        

    def create_directories(self):
        self.pdb_dir = os.path.join(os.getcwd(), "pdb")
        os.makedirs(self.pdb_dir, exist_ok = True) 
        
        self.pdbqt_dir = os.path.join(os.getcwd(), "pdbqt")
        os.makedirs(self.pdbqt_dir, exist_ok = True)

        


    def copy_starting_enzyme_into_pdb_dir(self):        
        source = os.path.join(os.getcwd(), self.starting_enzyme)
        destination = os.path.join(self.pdb_dir, self.starting_enzyme)
        shutil.copy2(source, destination)
    
    
    
    def clean_up_dir(self, directory): # kanske run den här i __init__
        
        if directory == "pdb":
            directory = self.pdb_dir
            self.pdb_filenames = []
            self.pdb_filenames.append(self.raw_filename)
            
            
        elif directory == "pdbqt":
            directory = self.pdbqt_dir
            self.pdbqt_filenames = []
        
        for filename in os.listdir(directory):           
            file_path = os.path.join(directory, filename)
            if os.path.isfile(file_path):
                os.remove(file_path)
        self.copy_starting_enzyme_into_pdb_dir()
    
    def viz(self): # Saga
        pass

    
    
    def find_molecule_coordinates(self): # Ebba
        # tar molekylens namn och chain som input
        # använd self.raw_filename
        # returnera koordinaterna för center, alternativt spara i attibut self.docking_center, kan användas senare då
        pass

    
    
    def clean_up(self): # Ebba
        # använd self.raw_filename
        # använd rensa alla molekyler (kan göras som i notebooken eller med ex mdanalysis)
        # returnera clean filen, kör funktionen i __init__
        pass

    
    
    def add_functional_group(self): # Saga
        # testa kolla pymol
        pass


    
    def mutations(self, subunits = 'All', positions = None, mutations = None): 
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
        
           
             

    
    def create_pdbqt(self): # Eugen
        
        for filename in self.pdb_filenames:
            filename = filename[:-4]
            subprocess.run(f"pdb2pqr30 --keep-chain --with-ph 7.4 --ff=PARSE {self.pdb_dir}/{filename}.pdb {self.pdbqt_dir}/{filename}.pqr -q --log-level CRITICAL")
            u = mda.Universe(f"{self.pdbqt_dir}/{filename}.pqr")
            u.atoms.write(f'{self.pdbqt_dir}/{filename}_temp.pdbqt')

            # Removes the first two lines, makes the file work in Vina
            with open(f"{self.pdbqt_dir}/{filename}_temp.pdbqt") as rf, open(f"{self.pdbqt_dir}/{filename}.pdbqt", 'w') as wf:
                for i, line in enumerate(rf):
                    if i >= 2:
                        wf.write(line)
       
            os.remove(f"{self.pdbqt_dir}/{filename}_temp.pdbqt")
            os.remove(f"{self.pdbqt_dir}/{filename}.log")
            os.remove(f"{self.pdbqt_dir}/{filename}.pqr")
            self.pdbqt_filenames.append(f"{filename}.pdbqt")
    
    
    def mol_to_pdbqt(self, filename):
        filename = filename[:-4]
        subprocess.run(f"obabel {filename}.mol -O {filename}.pdbqt --partialcharge gasteiger")
        self.ligand_filenames.append(f"{filename}.pdbqt")
    

    
    def windows_docking(self, center, boxsize = 20): # Eugen
        cx = center[0] # byta ut mot class attr när metoden find_molecule_coordinates är gjord
        cy = center[1]
        cz = center[2]
        bx = boxsize
        by = boxsize
        bz = boxsize

        
        for ligand in self.ligand_filenames:
            for enzyme in self.pdbqt_filenames:
                print(enzyme)
                config = open('config.txt', mode='w') # kanaske skapa olika filer och spara de i en egen mapp
                config.write(f"receptor={self.pdbqt_dir}/{enzyme}\n")
                config.write(f"ligand={ligand}\n")
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
                res = subprocess.run(f'"vina.exe" --config config.txt --log log.txt --out docked_{ligand[:-6]}_in_{enzyme[:-6]}.pdbqt --exhaustiveness 20 --num_modes 20 --energy_range 6', capture_output=True, text = True)
                print(res.stdout)
    
    def os_docking(self): # Ebba
        pass
        
    