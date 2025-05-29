
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

from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout

from rdkit import Chem
from rdkit.Chem import AllChem, Draw, rdMolTransforms
from rdkit.Chem.Draw import MolsToGridImage, IPythonConsole

aa_dict = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"
}



class KEX():

    def __init__(self, filename):
        """
        Initialize a KEX object and prepare the working environment.

        Parameters
        ----------
        filename : str
            The name of the input PDB file.

        Attributes
        ----------
        raw_filename : str
            The original filename provided by the user.
        pdb_filenames : list
            List storing the cleaned PDB filenames.
        pdbqt_filenames : list
            List storing converted PDBQT filenames.
        ligand_filenames : list
            List of filenames for prepared ligand structures.
        starting_enzyme : str
            The filename of the initial enzyme structure.
        """
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
        """
        Create necessary directories for storing input and output files.
    
        This method creates the following folders in the current working directory:
        - 'pdb' for cleaned PDB files
        - 'pdbqt' for converted PDBQT files
        - 'ligands' for prepared ligand structures
        - 'docking_results' for output from docking simulations
    
        Existing directories are preserved.
        """
        self.pdb_dir = os.path.join(os.getcwd(), "pdb")
        os.makedirs(self.pdb_dir, exist_ok = True) 
        
        self.pdbqt_dir = os.path.join(os.getcwd(), "pdbqt")
        os.makedirs(self.pdbqt_dir, exist_ok = True)

        self.ligands_dir = os.path.join(os.getcwd(), "ligands")
        os.makedirs(self.ligands_dir, exist_ok = True)

        self.docking_results_dir = os.path.join(os.getcwd(), "docking_results")
        os.makedirs(self.docking_results_dir, exist_ok = True)

        


    def copy_starting_enzyme_into_pdb_dir(self):
        """
        Copy the starting enzyme structure file into the 'pdb' directory.
    
        This method copies the cleaned enzyme PDB file from the current working directory
        to the designated 'pdb' folder for further processing.
        """
        source = os.path.join(os.getcwd(), self.starting_enzyme)
        destination = os.path.join(self.pdb_dir, self.starting_enzyme)
        shutil.copy2(source, destination)
        
    
    
    
    def clean_up_dir(self, directory):
        """
        Remove all files from a specified working directory and reset related filename lists.
    
        Parameters
        ----------
        directory : str
            The name of the directory to clean. Valid options are:
            - "pdb": deletes all files in the 'pdb' folder and resets self.pdb_filenames
            - "pdbqt": deletes all files in the 'pdbqt' folder and resets self.pdbqt_filenames
            - "ligands": deletes all files in the 'ligands' folder and resets self.ligand_filenames
            - "docking_results": deletes all files in the 'docking_results' folder
    
        Notes
        -----
        This method only removes regular files (not subdirectories).
        """
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
    
    
    def find_molecule_coordinates(self, molecule_name, chain): # Ebba
        """
        Calculate the center of mass coordinates for a specific molecule within a given chain.
    
        Parameters
        ----------
        molecule_name : str
            The residue name of the molecule to locate (e.g., an inhibitor or ligand).
        chain : str
            The chain ID where the molecule is located.
    
        Returns
        -------
        tuple of float
            A 3D coordinate tuple representing the center of mass of the selected molecule.
    
        Side Effects
        ------------
        Sets the attribute `self.docking_center` to the computed coordinates.
        """
        with pymol2.PyMOL() as pm:
            pm.cmd.load(self.raw_filename, "protein")

            molecule_selection = f"resn {molecule_name} and chain {chain}"
            center = pm.cmd.centerofmass(molecule_selection)
            self.docking_center = center 
            return center


    
    
    def clean_up(self):
        """
        Create a cleaned version of the original PDB file containing only relevant structural lines.
    
        This method reads the raw PDB file and writes a new file that includes only lines
        starting with "ATOM" or "TER", effectively removing any unnecessary lines.
    
        Returns
        -------
        str
            The filename of the cleaned PDB file (with suffix '_clean.pdb').
        """
        clean_filename = self.raw_filename.replace(".pdb", "_clean.pdb")
        with open(self.raw_filename, "r") as infile, open(clean_filename, "w") as outfile:
            for line in infile:
                if line.startswith(("ATOM", "TER")):
                    outfile.write(line) 
        return clean_filename
        

                
    
    def create_molecule(self, smiles, name, show_structure = False):
        """
        Generate a 3D molecule from a SMILES string and save it as a .mol file.
    
        This method creates a hydrogenated molecule from a SMILES string, embeds it in 3D space,
        performs MMFF geometry optimization, and saves it as a .mol file. Optionally, it displays 
        the structure in an 3D viewer using py3Dmol.
    
        Parameters
        ----------
        smiles : str
            The SMILES string representing the molecular structure.
        name : str
            The base name used to save the .mol file (without extension).
        show_structure : bool, optional
            If True, displays the optimized molecule in a 3D viewer (default is False).
    
        Side Effects
        ------------
        - Saves a .mol file with the optimized molecular structure.
        - Optionally displays the structure in an interactive window.
        - Prints the file save confirmation to the console.
        """
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
            viewer.setBackgroundColor("black")
            viewer.zoomTo()
            viewer.show()
        # Save as MOL file
        mol_filename = f"{name}.mol"
        Chem.MolToMolFile(molecule, mol_filename)
        print(f"Molecule saved as {mol_filename}")
    
    
        
    def mutations(self, subunits='All', positions=None, mutations=None, label=None):
        """
        Apply one or more point mutations to specified residue positions using PyMOL.
    
        If subunits is set to 'All', mutations are applied to all chains. Otherwise, the user can
        specify which chains to mutate. The resulting structure is saved as a PDB file with a
        descriptive filename, and all PyMOL output is logged.
    
        Parameters
        ----------
        subunits : str or list of str
            'All' to mutate all chains, or a list of specific chain IDs.
        positions : int or list of int
            Residue positions to mutate.
        mutations : str or list of str
            Target residue names (3-letter codes) to mutate to.
        label : str, optional
            Optional label to include in the output filename.
        """
        
    
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
    
        sys.stdout = original_stdout  
        log_file.close()

        
           
    def molecular_dynamics(self):
        """
        Run energy minimization and short NVT molecular dynamics simulations for each PDB structure.
    
        For each structure in self.pdb_filenames, this method sets up an OpenMM simulation using
        the Amber14 force field, performs energy minimization, runs 10,000 steps of NVT dynamics,
        and saves the final structure as a new PDB file with '_md' appended to the name.
    
        Side Effects
        ------------
        - Modifies each entry in self.pdb_filenames to point to the MD-relaxed structure.
        - Writes simulation output to 'output.pdb' and logs to 'md_log.txt'.
        - Prints simulation progress to stdout.
        """
        for i, filename in enumerate(self.pdb_filenames):
            filename = filename[:-4]
            pdb = PDBFile(f"{self.pdb_dir}/{filename}.pdb")
            # Specify the forcefield
            forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
            modeller = Modeller(pdb.topology, pdb.positions)
            modeller.deleteWater()
            residues=modeller.addHydrogens(forcefield)

            # System setup
            system = forcefield.createSystem(modeller.topology, nonbondedMethod=NoCutoff, constraints=HBonds)
            integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
            simulation = Simulation(modeller.topology, system, integrator)
            simulation.context.setPositions(modeller.positions)

            print("Minimizing energy")
            simulation.minimizeEnergy()

            # Setup simulation
            simulation.reporters.append(PDBReporter('output.pdb', 1000))
            simulation.reporters.append(StateDataReporter(stdout, 100, step=True, potentialEnergy=True, temperature=True, volume=True))
            simulation.reporters.append(StateDataReporter("md_log.txt", 100, step=True, potentialEnergy=True, temperature=True, volume=True))

            print("Running NVT")
            simulation.step(10000)
            
            print("Minimizing energy")
            simulation.minimizeEnergy()
            
            positions = simulation.context.getState(getPositions=True).getPositions()
            PDBFile.writeFile(simulation.topology, positions, open(f'{self.pdb_dir}/{filename}_md.pdb', 'w'))
            self.pdb_filenames[i] = f'{filename}_md.pdb'
                
    def create_pdbqt(self):
        """
        Convert PDB structures to PDBQT format using pdb2pqr30 and MDAnalysis.
    
        For each file in self.pdb_filenames:
        - Runs pdb2pqr30 to assign charges and add hydrogens.
        - Converts the resulting .pqr file to .pdbqt using MDAnalysis.
        - Cleans up intermediate files.
        - Appends the final .pdbqt filename to self.pdbqt_filenames.
    
        Note
        ----
        Removes the first two lines of the .pdbqt to ensure compatibility with AutoDock Vina.
        """
       
        for filename in self.pdb_filenames:
            filename = filename[:-4]
            subprocess.run(f"pdb2pqr30 --keep-chain --with-ph 7.4 --ff=PARSE {self.pdb_dir}/{filename}.pdb {self.pdbqt_dir}/{filename}.pqr -q --log-level INFO") # CRITICAL
            u = mda.Universe(f"{self.pdbqt_dir}/{filename}.pqr")
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", UserWarning)
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

    
    


    def mol_to_pdbqt_new(self, filename, add_ligand = False):
        """
        Convert a .mol file to .pdbqt format using Open Babel and optionally move it to the ligand directory.
    
        Parameters
        ----------
        filename : str
            The name of the .mol file to convert (with extension).
        add_ligand : bool, optional
            If True, adds the resulting .pdbqt file to self.ligand_filenames and moves it to the 'ligands' folder.
    
        Side Effects
        ------------
        - Creates a .pdbqt file from the input .mol file using Gasteiger charges.
        - Optionally updates ligand_filenames and moves the file to the ligand directory.
        """
        
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
        """
        Perform molecular docking on Windows using AutoDock Vina for all enzyme-ligand combinations.
    
        For each ligand and enzyme pair, this method:
        - Creates a Vina configuration file centered at the specified coordinates.
        - Runs docking using 'vina.exe'.
        - Extracts the binding affinity from the resulting .pdbqt file.
        - Stores the results in a DataFrame with enzymes as rows and ligands as columns.
    
        Parameters
        ----------
        center : tuple of float
            The (x, y, z) coordinates of the center of the docking box.
        boxsize : int, optional
            The size of the docking box in Ångström along each axis (default is 20).
    
        Returns
        -------
        pd.DataFrame
            A DataFrame containing binding affinities (in kcal/mol) for each ligand-enzyme pair.
        """
        
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
    

        
    