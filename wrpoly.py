from pymatgen.core import Structure
from scipy.stats import skew
from pymatgen.core import Composition
from scipy.spatial import ConvexHull
import numpy as np
import pandas as pd
import math
from polyhedral_analysis.configuration import Configuration
from polyhedral_analysis.polyhedra_recipe import PolyhedraRecipe

"""
Cleaned version of wrpoly.py containing only functions needed for VizBATT.py
This module analyzes crystal structures to extract polyhedra information.
Main entry points:
- get_structure(filename): Load and simplify a crystal structure from a CIF file
- get_average_df(structure): Extract comprehensive polyhedra analysis data
"""

def get_transition_metal_sites(structure):
    """Get site indices of all transition metals, metalloids, and metals in the structure."""
    tms = [i for i, site in enumerate(structure) if site.species.elements[0].is_transition_metal \
        or site.species.elements[0].is_metalloid \
                or site.species.elements[0].is_metal]
    return tms

def get_alakali_sites(structure):
    """Get site indices of all alkali elements in the structure."""
    alkali = [i for i, site in enumerate(structure) if site.species.elements[0].is_alkali]
    return alkali

def get_alkaline_sites(structure):
    """Get site indices of all alkaline earth elements in the structure."""
    alkaline = [i for i, site in enumerate(structure) if site.species.elements[0].is_alkaline]
    return alkaline

def get_lanthanide_sites(structure):
    """Get site indices of all lanthanide elements in the structure."""
    landanide = [i for i, site in enumerate(structure) if site.species.elements[0].is_lanthanoid]
    return landanide

def get_actinide_sites(structure):
    """Get site indices of all actinide elements in the structure."""
    actinide = [i for i, site in enumerate(structure) if site.species.elements[0].is_actinoid]
    return actinide

def simplify_structure(struct):
    """
    Simplify a pymatgen Structure by selecting the most prevalent species for each site.
    If there's a tie, the first species is selected.
    """
    for i, site in enumerate(struct.sites):
        # Sort the species on the site by occupancy, and pick the most prevalent one.
        most_prevalent_species = sorted(site.species.items(), key=lambda x: x[1], reverse=True)[0][0]
        
        # Replace the site species with the most prevalent one.
        struct.replace(i, most_prevalent_species, properties=site.properties)
    
    return struct

def get_structure(filename):
    """
    Load a crystal structure from a file and simplify it.
    Returns None if loading fails.
    """
    try:
        structure = Structure.from_file(filename)
        return simplify_structure(structure)
    except:
        return None

def bj_get_vertex_atoms(structure):
    """Get list of non-central atom species (vertex atoms for polyhedra)."""
    try:
        sites = [i.species.elements[0].symbol for i in structure]
        tms_sites = bj_get_central_atoms(structure)
        species = list(set(sites))
        non_tms_species = list(set(species) - set(tms_sites))
    except:
        sites = [i.species.elements[0].element.value for i in structure]
        tms_sites = bj_get_central_atoms(structure)
        species = list(set(sites))
        non_tms_species = list(set(species) - set(tms_sites))
    return non_tms_species

def bj_get_central_atoms(structure):
    """Get list of central atom species for polyhedra analysis."""
    try:
        tms = get_transition_metal_sites(structure)
        tms_sites = np.unique([structure[i].specie.value for i in tms]).tolist()
        alakali_sites = get_alakali_sites(structure)
        alkaline_sites = get_alkaline_sites(structure)
        landanide_sites = get_lanthanide_sites(structure)
        actinide_sites = get_actinide_sites(structure)
        phos_sites = [i for i in range(len(structure)) if structure[i].specie.value == 'P']
        tms_sites = tms_sites + alakali_sites + alkaline_sites + landanide_sites + actinide_sites
    except:
        tms = get_transition_metal_sites(structure)
        tms_sites = np.unique([structure[i].species.elements[0].element.value for i in tms]).tolist()
        alakali_sites = get_alakali_sites(structure)
        alkaline_sites = get_alkaline_sites(structure)
        landanide_sites = get_lanthanide_sites(structure)
        actinide_sites = get_actinide_sites(structure)
        phos_sites = [i for i in range(len(structure)) if structure[i].species.elements[0].element.value  == 'P']
        tms_sites = tms_sites + alakali_sites + alkaline_sites + landanide_sites + actinide_sites
                
    return tms_sites

def get_bj_config(structure, radii=2.5):
    """Create a polyhedral analysis configuration for the structure."""
    central_atoms = bj_get_central_atoms(structure)
    vertex_atoms = bj_get_vertex_atoms(structure)
    print('radius cutoff: ', radii)
    try:
        recipe = PolyhedraRecipe(method='distance cutoff', 
                        coordination_cutoff=radii,
                        central_atoms=central_atoms,
                        vertex_atoms=vertex_atoms, 
                        n_neighbours=6)
        config = Configuration(structure=structure, recipes=[recipe])
        polyhedra = [i.best_fit_geometry['geometry'] for i in config.polyhedra]
        print('polyhedra: ', polyhedra)
        if len(polyhedra) > 0:
            return config
    except:
        print('failed to get config using radius cutoff: ', radii)
        print('switching to nearest neighbours method')
        recipe = PolyhedraRecipe(method='nearest neighbours', 
                        coordination_cutoff=6,
                        central_atoms=central_atoms,
                        vertex_atoms=vertex_atoms, 
                        n_neighbours=6)
        config = Configuration(structure=structure, recipes=[recipe])
        polyhedra = [i.best_fit_geometry['geometry'] for i in config.polyhedra]
        print('polyhedra: ', polyhedra)
        if len(polyhedra) > 0:
            return config
    return config

def scale_structure_for_analysis(structure, poly_radii=2.5, poly_factor=3):
    """
    Scale the structure to ensure it's large enough for polyhedra analysis.
    Creates supercell if any lattice dimension is smaller than poly_radii*poly_factor.
    """
    lattice = structure.lattice
    scale_axes = [i for i, j in enumerate(lattice.abc) if j < poly_radii*poly_factor]
    scale_factors = [poly_radii*poly_factor/np.linalg.norm(i) for i in lattice.abc if i < poly_radii*poly_factor]
    if len(scale_factors) == 0:
        scale_factors = [1]
    scale_factor = np.ceil(min(scale_factors))
    scales = [1, 1, 1]
    for i in scale_axes:
        scales[i] = scale_factor
    print(scales)
    structure.make_supercell(scales)
    return structure

def bj_structure(structure, poly_radii=2.5, poly_factor=3):
    """Create a polyhedral configuration from a structure."""
    structure = scale_structure_for_analysis(structure, poly_radii, poly_factor)
    config = get_bj_config(structure, radii=poly_radii)
    return config

def bj_get_poly_distances(poly):
    """Get distances from central atom to all vertex atoms in a polyhedron."""
    return [poly.central_atom.distance(poly.vertices[i]) for i in range(len(poly.vertices))]

def bj_get_average_poly_distance(poly):
    """Get average distance from central atom to vertex atoms."""
    return sum(bj_get_poly_distances(poly))/len(poly.vertices)

def bj_get_distortion_index(poly):
    """Calculate distortion index for a polyhedron."""
    return np.sum(np.abs(np.array(bj_get_poly_distances(poly)) - bj_get_average_poly_distance(poly)))/len(poly.vertices)/bj_get_average_poly_distance(poly)

def bj_get_l0(poly):
    """Calculate ideal bond length l0 for a polyhedron based on its volume and geometry."""
    vol = poly.volume
    struc_type = poly.best_fit_geometry['geometry'].lower()
    
    if 'octahedron' in struc_type or 'trigonal bipyramid' in struc_type:
        s = (vol * 6 * math.sqrt(2))**(1/3)
        l0 = math.sqrt(2/3) * s / 2

    elif 'tetrahedron' in struc_type:
        s = (vol * 6 * math.sqrt(2))**(1/3)
        l0 = math.sqrt(2/3) * s / 2

    elif 'trigonal prism' in struc_type:
        s = (vol * 8 / (math.sqrt(3)**3))**(1/3)
        l0 = (math.sqrt(3) / 4) * s

    elif 'square pyramid' in struc_type:
        s = (vol * 3 * math.sqrt(2))**(1/3)
        l0 = s / math.sqrt(2)
    
    return l0

def bj_quadratic_elongation(poly):
    """Calculate quadratic elongation for a polyhedron."""
    l0 = np.array(bj_get_l0(poly))
    l = np.array(bj_get_average_poly_distance(poly))
    quad_elongation = np.mean((l/l0)**2)
    return quad_elongation

def bj_get_bond_angles(poly):
    """Get all bond angles in a polyhedron, excluding angles > 160 degrees."""
    return [i for i in poly.angles() if i < 160]

def bond_angle_table():
    """Return ideal bond angles for common polyhedra geometries."""
    return {'Octahedron': 90,
            'Cube': 90,
            'Dodecahedron': 116.56505117707799,
            'Icosahedron': 138.18968510438353,
            'Tetrahedron': 109.47122063449069}

def bj_bond_angle_variance(poly):
    """Calculate bond angle variance for a polyhedron."""
    try:
        ideal_bond_angle = bond_angle_table()[poly.best_fit_geometry['geometry']]
        bond_angles = bj_get_bond_angles(poly)
        return np.sum((np.array(bond_angles) - ideal_bond_angle)**2)/(len(bond_angles)-1)
    except:
        print('failed to get ideal bond angle')
        return None

def get_poly_volume(poly):
    """Get volume of a polyhedron."""
    return poly.volume

def get_bj_df(structure):
    """
    Create a comprehensive DataFrame with polyhedra analysis data.
    Returns DataFrame with columns for polyhedra types, bond lengths, distortion indices,
    connectivity (corner/edge/face sharing), volumes, bond angles, and more.
    """
    bj_struc = bj_structure(structure)
    polys = [i for i in bj_struc.polyhedra]

    # Extracting information into separate lists
    poly_types = [i.best_fit_geometry['geometry'] for i in polys]
    bond_lengths = [bj_get_poly_distances(i) for i in polys]
    avg_bond_length = [bj_get_average_poly_distance(i) for i in polys]
    std_bond_length = [np.std(i) for i in bond_lengths]
    skew_bond_length = [skew(i) for i in bond_lengths]
    distortion_index = [bj_get_distortion_index(i) for i in polys]
    quadratic_elongation = [bj_quadratic_elongation(i) for i in polys]
    central_atoms = [i.central_atom.site.specie.value for i in polys]
    polyhedra_atoms = [[j.site.specie.value for j in i.vertices] for i in polys]
    polyhedra_atoms = ["".join(map(str, i)) for i in polyhedra_atoms]
    polyhedra_atoms = [Composition(i).formula for i in polyhedra_atoms]
    poly_formulas = ["".join(i) for i in polyhedra_atoms]
    poly_formulas = [Composition(i).formula for i in poly_formulas]

    # Creating a dictionary for DataFrame construction
    data = {
        'polyhedra_atoms': polyhedra_atoms,
        'polyhedra_formula': poly_formulas,
        'poly_types': poly_types,
        'bond_lengths': bond_lengths,
        'avg_bond_length': avg_bond_length,
        'std_bond_length': std_bond_length,
        'skew_bond_length': skew_bond_length,
        'distortion_index': distortion_index,
        'quadratic_elongation': quadratic_elongation,
        'central_atoms': central_atoms,
        'polyhedra': polys,
        'n_corner_pairs': [len(i.corner_sharing_neighbour_list()) for i in polys],
        'n_edge_pairs': [len(i.edge_sharing_neighbour_list()) for i in polys],
        'n_face_pairs': [len(i.face_sharing_neighbour_list()) for i in polys],
        'n_atoms_in_cell': [len(structure) for i in polys],
        'Formula': [structure.composition.reduced_formula for i in polys],
        'avg_bond_angles': [np.mean(bj_get_bond_angles(i)) for i in polys],
        'std_bond_length_angles': [np.std(bj_get_bond_angles(i)) for i in polys],
        'skew_bond_length_angles': [skew(bj_get_bond_angles(i)) for i in polys],
        'volumes': [get_poly_volume(i) for i in polys],
        'bond_angle_variance': [bj_bond_angle_variance(i) for i in polys],
        'space_group': [structure.get_space_group_info()[1] for i in polys]
    }

    # Creating DataFrame
    df = pd.DataFrame(data)
    return df

def get_average_df(structure):
    """
    Create aggregated DataFrame with average polyhedra properties grouped by central atom,
    formula, polyhedra type, and polyhedra formula.
    This is the main analysis function called by VizBATT.py.
    """
    try:
        df = get_bj_df(structure)

        # Convert relevant columns to numeric
        numeric_columns = ['avg_bond_length', 'std_bond_length', 'skew_bond_length', 'distortion_index',
                            'quadratic_elongation', 'n_corner_pairs', 'n_edge_pairs', 'n_face_pairs',
                            'n_atoms_in_cell', 'avg_bond_angles', 'std_bond_length_angles', 'skew_bond_length_angles',
                            'volumes', 'bond_angle_variance']

        # Print data types before conversion
        print("Data types before conversion:")
        print(df.dtypes)

        df[numeric_columns] = df[numeric_columns].apply(pd.to_numeric, errors='coerce')

        # Print data types after conversion
        print("\nData types after conversion:")
        print(df.dtypes)

        # Group by and calculate mean
        new_df = df.groupby(['central_atoms', 'Formula', 'poly_types', 'polyhedra_formula', 'avg_bond_length'])[numeric_columns].mean()
        sum_df = df.groupby(['central_atoms', 'Formula', 'poly_types', 'polyhedra_formula', 'avg_bond_length'])[numeric_columns].sum()
        new_df['n_corner_pairs'] = sum_df['n_corner_pairs'].values
        new_df['n_edge_pairs'] = sum_df['n_edge_pairs'].values
        new_df['n_face_pairs'] = sum_df['n_face_pairs'].values
        new_df['central_atoms'] = new_df.index.get_level_values('central_atoms') 
        new_df['polyhedra_formula'] = new_df.index.get_level_values('polyhedra_formula')
        new_df['poly_types'] = new_df.index.get_level_values('poly_types')
        new_df['n_atoms_per_cell'] = [structure.num_sites for i in range(len(new_df))]
        # Divide all n_*_pairs by n_atoms_per_cell
        new_df['n_corner_pairs'] = new_df['n_corner_pairs']/new_df['n_atoms_per_cell']
        new_df['n_edge_pairs'] = new_df['n_edge_pairs']/new_df['n_atoms_per_cell']
        new_df['n_face_pairs'] = new_df['n_face_pairs']/new_df['n_atoms_per_cell']
        print(df['Formula'].astype(str).unique())
        new_df['Formula'] = [df['Formula'].astype(str).unique()[0] for i in range(len(new_df))]
        new_df = new_df.reset_index(drop=True)
    except:
        print('error creating dataframe for ', structure.composition.reduced_formula)
        new_df = pd.DataFrame()
    return new_df
