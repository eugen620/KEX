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

        
        # kör clean up och här appenda den till self.pdb_filenames

        
        # *** Tillfällig kod ****
        self.pdb_filenames.append(filename) # vi vill egentligen inte lägga in den råa filen här, vi ska fixa en clean up method och sedan appen den clean versionen i append metoden.
        self.starting_enzyme = self.pdb_filenames[0]
        

        
        # Sjapar pdb och pdbqt mappar och attribut med filepath så man kan använda de i koden sen
        # Kanske ha denna del i mutations metoden, nu skapar den mappen oavsätt
        self.pdb_dir = os.path.join(os.getcwd(), "pdb")
        os.makedirs(self.pdb_dir, exist_ok = True) 

        self.pdbqt_dir = os.path.join(os.getcwd(), "pdbqt")
        os.makedirs(self.pdbqt_dir, exist_ok = True)

        # Den här delen tar just nu clean versionen från samma directory som notebooken är i och flyttar in den i pdb mappen.
        # Vi kan skriva om senare så att clean up metoden direkt sätter in clean cersionen i den mappen.
        source = os.path.join(os.getcwd(), self.starting_enzyme)
        destination = os.path.join(self.pdb_dir, self.starting_enzyme)
        shutil.copy2(source, destination)
    
    
    
    def clean_up_pdb_dir(self):
        for filename in os.listdir(self.pdb_dir):           
            file_path = os.path.join(self.pdb_dir, filename)
            if os.path.isfile(file_path):
                os.remove(file_path)
        
        # Förbättra den här delen sen
        self.pdb_filenames = []
        self.pdb_filenames.append(self.raw_filename)
        
    
    
    
    def viz(self): # Saga
        pass

    
    
    def find_molecule_coordinates(self): # Ebba
        # tar molekylens namn och chain som input
        # använd self.raw_filename
        # returnera koordinaterna för center
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
    
    def windows_docking(self): # Eugen
        pass

    
    def os_docking(self): # Ebba
        pass
        
    