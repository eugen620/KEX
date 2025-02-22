# imorts
import pymol2
import nglview as nv
import MDAnalysis as mda
import os
import subprocess

from dict_module import aa_dict

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
                
            for chain in chains:
                for i, pos in enumerate(positions):
                    mutation = mutations[i]
                    
                    model = pm.cmd.get_model(f"/MutProt//{chain}/{pos}")
                    starting_aa = model.atom[0].resn

                    apply_mutation()
                

            pm.cmd.set_wizard()

            # spara mutationen i klassen
            new_filename = f"{aa_dict[starting_aa]}{pos}{aa_dict[mutation]}.pdb"
            self.pdb_filenames.append(new_filename)

            # spara mutationen i pdb mappen
            output_path = os.path.join(self.pdb_dir, new_filename)
            pm.cmd.save(output_path, "MutProt")
            
                        
             
            
                

    
    def create_pdbqt(self): # Eugen
        # file handeling
        pass

    
    def windows_docking(self): # Eugen
        pass

    
    def os_docking(self): # Ebba
        pass
        
    