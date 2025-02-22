# imorts
import pymol2
import nglview as nv
import MDAnalysis as mda
import os
import subprocess

from dict_module import aa_dict, aa_list

class KEX():

    def __init__(self, filename):
        self.raw_filename = filename
        self.pdb_filenames = []
        self.pdbqt_filenames = []

        # kör clean up och appenda den till self.pdb_filenames

        
        # *** Tillfällig kod ****
        self.pdb_filenames.append(filename)
        self.starting_enzyme = self.pdb_filenames[0]
        


        self.pdb_dir = os.path.join(os.getcwd(), "pdb")
        os.makedirs(self.pdb_dir, exist_ok = True) 


    def clean_up_pdb_dir(self):
        for filename in os.listdir(self.pdb_dir):           
            file_path = os.path.join(self.pdb_dir, filename)
            if os.path.isfile(file_path):
                os.remove(file_path)
        
    
    
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


      

    
    def mutations(self, subunits = 'All', positions = None, mutations = None): # Eugen
        
        def apply_mutation():
            pm.cmd.get_wizard().set_mode(mutation)
            pm.cmd.get_wizard().do_select(f"/MutProt//{chain}/{pos}/")
            pm.cmd.get_wizard().apply()              
        
        if positions == None or mutations == None:
            return "You have to enter information"
        
        if type(positions) == int:
            positions = [positions]
        
        if subunits != 'All' and type(subunits) != list:
            return "Enter the chains in a list"

                
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
                
                # spara mutationen i klassen
                new_filename = f"{'_'.join(new_filename)}.pdb"
                self.pdb_filenames.append(new_filename)

                # spara mutationen i pdb mappen
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
                            print(f"Här är det en odödig mutation, {starting_aa} till {mutation}")    
                            continue
                        pm.cmd.set_wizard()
                    
                        mutation_str = f"{aa_dict[starting_aa]}{pos}{aa_dict[mutation]}"
                        new_filename.append(mutation_str)               
                    
                        new_filename = f"{'_'.join(new_filename)}.pdb"
                        self.pdb_filenames.append(new_filename)

                        output_path = os.path.join(self.pdb_dir, new_filename)
                        pm.cmd.save(output_path, "MutProt")
            
                        
             

    
    def create_pdbqt(self): # Eugen
        # file handeling
        pass

    
    def windows_docking(self): # Eugen
        pass

    
    def os_docking(self): # Ebba
        pass
        
    