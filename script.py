# imorts


class KEX():

    def __init__(self, filename):
        self.raw_filename = filename
        self.pdb_filenames = []
        self.pdbqt_filenames = []

        # kör clean up och appenda den till self.pdb_filenames


    
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

    
    def mutations(self): # Eugen
        pass

    
    def create_pdbqt(self): # Eugen
        # file handeling
        pass

    
    def windows_docking(self): # Eugen
        pass

    
    def os_docking(self): # Ebba
        pass
        
    