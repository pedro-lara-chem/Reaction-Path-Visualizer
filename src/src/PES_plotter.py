# =============================================================================
#
# Potential Energy Surface (PES) Visualization Tool (Advanced)
#
# Description:
#   This script generates 1D potential energy surface diagrams for chemical
#   reactions. It creates interactive 3D plots, animated GIFs of deactivation
#   pathways, and static 2D plots for publication. This advanced version
#   can handle complex deactivation pathways involving multiple conical
#   intersections (CIs) and intersystem crossings (ISCs).
#
# Features:
#   - Interactive GUI (Tkinter) for entering PES data (minima, transition
#     states, CIs, ISCs, etc.).
#   - A test mode to bypass the GUI and use predefined data for rapid testing.
#   - A recursive pathfinding algorithm to trace multi-step deactivation
#     cascades (e.g., S2 -> S1 -> T1 -> S0).
#   - Choice between flat 'ribbon', 'parabolic', or 'gaussian' 3D PES meshes.
#   - 3D visualization of PES ribbons using PyVista.
#   - Generation of separate, animated GIFs for each deactivation pathway,
#     including a ball that traverses the path and changes color at each
#     crossing point (CI or ISC).
#   - Static 3D plots with a continuous series of curved, on-surface arrows
#     indicating the full deactivation route.
#   - Publication-quality 2D plots using Matplotlib with labeled points.
#
# How to Run:
#   - For help and a list of options:
#     python script_name --help
#
#   - For default flat ribbons in GUI mode:
#     python script_name
#
#   - To choose a specific mesh type (e.g., gaussian):
#     python script_name --mesh-type gaussian
#
#   - To choose a specific file for reading the data.
#     python script_name -f name_of_csv file
#
#   - For non-interactive test mode with predefined data:
#     python script_name --test
#
# =============================================================================

import numpy as np
import os
import re
import datetime
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq
import pyvista as pv 
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import argparse
import tkinter as tk
from tkinter import simpledialog, messagebox
import itertools
import pandas as pd
from typing import List, Tuple, Dict, Any, Optional
from adjustText import adjust_text

# Define type aliases for clarity
Point = Tuple[float, float]
Crossing = Dict[str, Any]


# --- StateBuilder Class ---
class StateBuilder:
    """
    Constructs a single electronic state by accumulating its critical points.

    This class acts as a "builder" for a potential energy state, collecting
    minima, transition states, and crossings before finalizing them into a
    structured dictionary ready for interpolation.
    """
    def __init__(self, state_name: str):
        """
        Initializes the StateBuilder.

        :param state_name: The name of the electronic state (e.g., 'S0', 'T1').
        """
        self.state_name: str = state_name
        self.minima: List[Point] = []
        self.transition_states: List[Point] = []
        self.crossings: List[Crossing] = []
        self.franck_condon: Optional[Dict[str, float]] = None
        # self.raw_points stores all points to be used for spline interpolation
        self.raw_points: List[Point] = []

    def add_minimum(self, rc: float, energy: float):
        """Adds a minimum energy point to the state."""
        self.minima.append((rc, energy))
        self.raw_points.append((rc, energy))

    def add_transition_state(self, rc: float, energy: float):
        """Adds a transition state (saddle point) to the state."""
        self.transition_states.append((rc, energy))
        self.raw_points.append((rc, energy))

    def add_guide_point(self, rc: float, energy: float):
        """Adds a guide point used only for spline interpolation."""
        self.raw_points.append((rc, energy))

    def add_crossing(self, rc: float, energy: float, coupled_state: str, id: str, type_: str):
        """Adds a crossing point (CI or ISC) to the state."""
        self.crossings.append({
            'rc': rc, 'energy': energy,
            'coupled_state': coupled_state,
            'id': id,
            'type': type_
        })
        self.raw_points.append((rc, energy))

    def set_franck_condon(self, rc: float, energy: float):
        """Sets the Franck-Condon point for the state."""
        self.franck_condon = {'rc': rc, 'energy': energy}
        self.raw_points.append((rc, energy))

    def finalize(self) -> Dict[str, Any]:
        """
        Finalizes the state data for interpolation.

        This method sorts all collected points by their reaction coordinate (RC).

        :return: A dictionary containing all processed data for the state.
        """
 
        # Sort all points by the reaction coordinate (RC) before interpolation
        self.raw_points.sort(key=lambda p: p[0])
        points = np.array(self.raw_points)

        return {
            'rc': points[:, 0] if len(points) > 0 else np.array([]),
            'energy': points[:, 1] if len(points) > 0 else np.array([]),
            'minima': self.minima,
            'TS': self.transition_states,
            'crossings': self.crossings,
            'franck_condon': self.franck_condon
        }


# --- GUI for Input ---
class PESInputGUI:
    """
    A class to create a simple graphical user interface (GUI) for collecting
    data about potential energy surface (PES) states.

    This class handles the user input for the number of singlet and triplet states,
    and for each state, it collects information about its minima, saddle points,
    conical intersections (CIs), intersystem crossings (ISCs), and Franck-Condon (FC) points.
    """
    def __init__(self, root):
        """
        Initializes the GUI, setting up data structures and color palettes.
        
        Args:
            root (tk.Tk): The root Tkinter window.
        """
        print("Initializing GUI...")
        self.root = root
        # Dictionary to store all the data for each state (e.g., 'S0', 'S1', 'T1').
        self.state_data = {}
        # Dictionary to assign a unique color to each state's plot.
        self.state_colors = {}
        # A lookup for crossing points by their unique ID to avoid duplication.
        self.ci_points_by_id = {}
        # A predefined colorblind-friendly color palette (based on Paul Tol's schemes) suitable for projectors.
        self.color_palette = [
            '#332288', '#117733', '#44AA99', '#88CCEE', 
            '#DDCC77', '#CC6677', '#AA4499', '#882255', 
            '#999933', '#661100'
        ]
        self.color_index = 0

    def ask_states(self):
        """Asks the user for the number of singlet and triplet states via dialog boxes."""
        self.n_singlets = int(simpledialog.askstring("Input", "How many singlet states?"))
        self.n_triplets = int(simpledialog.askstring("Input", "How many triplet states?"))

    def get_next_color(self):
        """Provides the next color from the palette in a cyclical manner."""
        color = self.color_palette[self.color_index % len(self.color_palette)]
        self.color_index += 1
        return color

    def collect_state_data(self, state):
        """
        Collects detailed data for a single electronic state via a series of dialog boxes.
        A summary is shown to the user for confirmation before proceeding.

        Args:
            state (str): The name of the state (e.g., "S0", "T1").

        Returns:
            dict: A dictionary containing all the data points and details for the state.
        """
        while True:
            builder = StateBuilder(state)
            summary_text = f"--- Summary for {state} ---\n"

            num_minima = int(simpledialog.askstring("Input", f"Enter number of minima for {state}"))
            num_saddles = int(simpledialog.askstring("Input", f"Enter number of saddle points (TS) for {state}"))
            num_guides = int(simpledialog.askstring("Input", f"Enter number of guide points for {state} (for interpolation shape)"))
            ci_present = simpledialog.askstring("Input", f"Does {state} have conical intersections (CIs)? (yes/no)").lower() == 'yes'
            isc_present = simpledialog.askstring("Input", f"Does {state} have intersystem crossings (ISCs)? (yes/no)").lower() == 'yes'

            if state != "S0" and simpledialog.askstring("Input", f"Do you want to define a Franck-Condon (FC) point for {state}? (yes or no)").lower() == 'yes':
                fc_rc = float(simpledialog.askstring("Input", "RC of FC point: "))
                fc_energy = float(simpledialog.askstring("Input", "Energy of FC point: "))
                builder.set_franck_condon(fc_rc, fc_energy)
                summary_text += f"FC Point: (RC={fc_rc}, E={fc_energy:.2f})\n"

            if num_minima > 0:
                summary_text += "\nMinima:\n"
                for i in range(num_minima):
                    rc = float(simpledialog.askstring("Input", f"RC of minimum {i+1}: "))
                    energy = float(simpledialog.askstring("Input", f"Energy of minimum {i+1}: "))
                    builder.add_minimum(rc, energy)
                    summary_text += f"  - Min {i+1}: (RC={rc}, E={energy:.2f})\n"

            if num_saddles > 0:
                summary_text += "\nSaddle Points (TS):\n"
                for i in range(num_saddles):
                    rc = float(simpledialog.askstring("Input", f"RC of saddle {i+1}: "))
                    energy = float(simpledialog.askstring("Input", f"Energy of saddle {i+1}: "))
                    builder.add_transition_state(rc, energy)
                    summary_text += f"  - TS {i+1}: (RC={rc}, E={energy:.2f})\n"

            if num_guides > 0:
                summary_text += "\nGuide Points:\n"
                for i in range(num_guides):
                    rc = float(simpledialog.askstring("Input", f"RC of guide point {i+1}: "))
                    energy = float(simpledialog.askstring("Input", f"Energy of guide point {i+1}: "))
                    builder.add_guide_point(rc, energy)
                    summary_text += f"  - Guide {i+1}: (RC={rc}, E={energy:.2f})\n"
            summary_text += "\nConical Intersections (CIs):\n"
            if ci_present:
                num_ci = int(simpledialog.askstring("Input", f"Enter number of CIs in {state}: "))
                if num_ci > 0:
                    for i in range(num_ci):
                        rc = float(simpledialog.askstring("Input", f"RC of CI {i+1}: "))
                        energy = float(simpledialog.askstring("Input", f"Energy of CI {i+1}: "))
                        coupled_state = simpledialog.askstring("Input", f"CI {i+1} connects to state: ")
                        ci_id = simpledialog.askstring("Input", f"Unique ID for CI {i+1}: ")
                        builder.add_crossing(rc, energy, coupled_state, ci_id, 'CI')
                        summary_text += f"  - CI {i+1} (ID: {ci_id}): (RC={rc}, E={energy:.2f}, connects to {coupled_state})\n"

            summary_text += "\nIntersystem Crossings (ISCs):\n"
            if isc_present:
                num_isc = int(simpledialog.askstring("Input", f"Enter number of ISCs in {state}: "))
                if num_isc > 0:
                    for i in range(num_isc):
                        rc = float(simpledialog.askstring("Input", f"RC of ISC {i+1}: "))
                        energy = float(simpledialog.askstring("Input", f"Energy of ISC {i+1}: "))
                        coupled_state = simpledialog.askstring("Input", f"ISC {i+1} connects to state: ")
                        isc_id = simpledialog.askstring("Input", f"Unique ID for ISC {i+1}: ")
                        builder.add_crossing(rc, energy, coupled_state, isc_id, 'ISC')
                        summary_text += f"  - ISC {i+1} (ID: {isc_id}): (RC={rc}, E={energy:.2f}, connects to {coupled_state})\n"

            if messagebox.askyesno("Confirm Data", f"{summary_text}\nIs this data correct?"):
                break

        return builder.finalize()

    
    def run(self):
        """Executes the main GUI workflow: ask for states, then collect data for each."""
        self.ask_states()
        # Loop through Singlet states
        for i in range(self.n_singlets):
            state = f"S{i}"; self.state_colors[state] = self.get_next_color(); self.state_data[state] = self.collect_state_data(state)
        # Loop through Triplet states
        for i in range(self.n_triplets):
            state = f"T{i+1}"; self.state_colors[state] = self.get_next_color(); self.state_data[state] = self.collect_state_data(state)
        print("\n--- All data collected. Closing GUI. ---")
        self.root.quit() # Close the (hidden) Tkinter root window to allow the script to proceed.

def get_data_from_file(filepath):
    """
    Reads PES data from a CSV file using pandas and formats it for the script.
    """
    print(f"--- FILE MODE ---\nLoading data from {filepath}...")
    try:
        df = pd.read_csv(filepath, sep=None, engine="python")
        df = df.where(pd.notnull(df), None)
        df.columns = df.columns.str.replace('[^A-Za-z0-9_]+', '', regex=True)
    except FileNotFoundError:
        print(f"Error: The file '{filepath}' was not found.")
        exit()
    except Exception as e:
        print(f"Error reading or parsing the file: {e}")
        exit()

    states = df['state'].unique()
    builders = {state: StateBuilder(state) for state in states}

    for _, row in df.iterrows():
        state = row['state']
        builder = builders[state]
        ptype, rc, energy = row['point_type'], row['rc'], row['energy']

        if ptype == 'minima':
            builder.add_minimum(rc, energy)
        elif ptype == 'TS':
            builder.add_transition_state(rc, energy)
        elif ptype == 'FC':
            builder.set_franck_condon(rc, energy)
        elif ptype == 'guide':
            builder.add_guide_point(rc, energy)
        elif ptype in ['CI', 'ISC']:
            builder.add_crossing(rc, energy, row['coupled_state'], row['id'], ptype)

    processed_state_data = {state: builder.finalize() for state, builder in builders.items()}

    color_palette = ['#D60270', 
                    '#9B4F96',  
                    '#0038A8', 
                    '#F58231', 
                    '#FFE119', 
                    '#469990',  
                    '#9A6324',  
                    '#808000',  
                    '#FFD8B1',  
                    '#A9A9A9',  
                    '#BFEF45',  
                    '#AAFFC3',  
                    '#DAA520',  
                    '#FFFAC8',  
                    '#7FFFD4',  
                    '#F0E68C',  
                    '#D2B48C',  
                    '#FF7F50',  
                    '#7FFF00', 
                    '#D2691E'   
                    ]
    state_colors = {state: color_palette[i % len(color_palette)] for i, state in enumerate(states)}

    class MockGUI:
        def __init__(self):
            self.state_data = processed_state_data
            self.state_colors = state_colors
    return MockGUI()

def get_test_data():
    """
    Provides a hardcoded set of data for testing, bypassing the GUI.
    This new data is designed to create an UNWANTED collision between S1 and T1
    to test the collision resolution function.
    """
    print("--- TEST MODE --- \nLoading predefined test data designed for collision detection...")
    raw_data = {
        'S0': {'minima': [(0.0, 0.0)], 'TS': [], 'crossings': [], 'fc': None},
        'S1': {
            'minima': [(2.0, 3.0)], 
            'TS': [], 
            'crossings': [
                {'rc': 1.0, 'energy': 4.5, 'coupled_state': 'S2', 'id': 'CI-S2S1', 'type': 'CI'}
            ], 
            'fc': None
        },
        'T1': {
            'minima': [(0.5, 4.0)], 
            'TS': [], 
            'crossings': [], 
            'fc': None
        },
        'S2': {
            'minima': [], 
            'TS': [], 
            'crossings': [
                {'rc': 1.0, 'energy': 4.5, 'coupled_state': 'S1', 'id': 'CI-S2S1', 'type': 'CI'}
            ], 
            'fc': {'rc': 0.0, 'energy': 5.5}
        }
    }
    
    processed_state_data = {}
    for state, data in raw_data.items():
        builder = StateBuilder(state)
        if data.get('fc'):
            builder.set_franck_condon(data['fc']['rc'], data['fc']['energy'])
        for rc, e in data.get('minima', []):
            builder.add_minimum(rc, e)
        for rc, e in data.get('TS', []):
            builder.add_transition_state(rc, e)
        for cross in data.get('crossings', []):
            builder.add_crossing(cross['rc'], cross['energy'], cross['coupled_state'], cross['id'], cross['type'])
        processed_state_data[state] = builder.finalize()
    
    state_colors = {'S0': '#332288', 'S1': '#117733', 'S2': '#44AA99', 'T1': '#88CCEE'}
    
    class MockGUI:
        def __init__(self):
            self.state_data = processed_state_data
            self.state_colors = state_colors
    return MockGUI()

def interpolation_1d(rc, energy, num_points=500, rc_range=None):
    """
    Performs 1D cubic spline interpolation on the given data points.
    
    Args:
        rc (np.array): Array of reaction coordinates.
        energy (np.array): Array of corresponding energies.
        num_points (int): The number of points to generate for the smooth curve.
        rc_range (tuple, optional): A tuple (min_rc, max_rc) to define the interpolation range.
        
    Returns:
        tuple: A tuple containing the new RC array, the new energy array, and the spline object.
    """
    # Handle cases with too few points for spline interpolation.
    if len(rc) < 2: 
        return np.array([]), np.array([]), None
    # Create the cubic spline function.
    cs = CubicSpline(rc, energy)

    # Define the range for the new, dense set of points.
    min_rc, max_rc = (rc_range if rc_range else (np.min(rc), np.max(rc)))
    if min_rc >= max_rc:
        return np.array([]), np.array([]), cs
        
    rc_new = np.linspace(min_rc, max_rc, num_points)
    return rc_new, cs(rc_new), cs


def create_ribbon_mesh(rc_path, energy_path, width=0.5, depth=0.0):
    """Creates a 3D flat ribbon mesh from a 1D path."""
    if len(rc_path) < 2: return None # Cannot create a mesh with less than 2 points
    y_center = 0.0
    # Define vertices for the front and back edges of the ribbon.
    v_front = np.column_stack((rc_path, np.full_like(rc_path, y_center - width / 2), energy_path))
    v_back = np.column_stack((rc_path, np.full_like(rc_path, y_center + width / 2), energy_path))
    mesh = pv.StructuredGrid()
    # Combine front and back vertices.
    mesh.points = np.vstack((v_front, v_back))
    # Define the grid dimensions (length of path, 2 points for width, 1 for depth).
    mesh.dimensions = [len(rc_path), 2, 1]
    return mesh

def create_parabolic_well_mesh(rc_path, energy_path, width=1.0, depth=0.5):
    """
    Creates a parabolic (harmonic oscillator) well-shaped 3D mesh.
    The energy increases quadratically in the direction perpendicular to the RC.
    """
    if len(rc_path) < 2: return None
    num_points_rc, num_points_perp = len(rc_path), 51 # Number of points along RC and perpendicular to it.
    y_coords = np.linspace(-width / 2, width / 2, num_points_perp)
    # Z-offset follows a parabolic curve: z = k*y^2
    z_offset = depth * (2 * y_coords / width)**2
    points = np.zeros((num_points_rc * num_points_perp, 3))
    # Create the grid of points.
    for i in range(num_points_rc):
        for j in range(num_points_perp):
            idx = i * num_points_perp + j
            points[idx] = [rc_path[i], y_coords[j], energy_path[i] + z_offset[j]]
    mesh = pv.StructuredGrid()
    mesh.points = points
    mesh.dimensions = [num_points_perp, num_points_rc, 1]
    return mesh

def create_gaussian_well_mesh(rc_path, energy_path, width=1.0, depth=0.5):
    """
    Creates a Gaussian well-shaped 3D mesh.
    The energy increases following a Gaussian curve perpendicular to the RC.
    """
    if len(rc_path) < 2: return None
    num_points_rc, num_points_perp = len(rc_path), 51
    y_coords = np.linspace(-width / 2, width / 2, num_points_perp)
    # Standard deviation (sigma) controls the "spread" of the well.
    sigma = width / 4  
    # Z-offset follows an inverted Gaussian function.
    z_offset = depth * (1 - np.exp(-0.5 * (y_coords / sigma)**2))
    points = np.zeros((num_points_rc * num_points_perp, 3))
    # Create the grid of points.
    for i in range(num_points_rc):
        for j in range(num_points_perp):
            idx = i * num_points_perp + j
            points[idx] = [rc_path[i], y_coords[j], energy_path[i] + z_offset[j]]
    mesh = pv.StructuredGrid()
    mesh.points = points
    mesh.dimensions = [num_points_perp, num_points_rc, 1]
    return mesh

def project_point_to_well(rc, z):
    """Projects a point onto the center (y=0) of any well-shaped mesh."""
    return (rc, 0.0, z)

def add_curved_arrow_on_surface(plotter, spline, start_rc, end_rc, mesh_type='ribbon', rc_range_span=None, force_plot=False):
    """
    Draws a curved arrow following the PES spline with a 3D cone arrowhead.

    Args:
        plotter (pv.Plotter): The PyVista plotter.
        spline (CubicSpline): Energy spline for the PES.
        start_rc (float): Start of arrow.
        end_rc (float): End of arrow (used to compute direction).
        mesh_type (str): 'ribbon', 'parabolic', or 'gaussian'.
        rc_range_span (float): Total RC range of the PES (used for dynamic cone size).
        force_plot (bool): If False, skip arrow if energy increases (used for CI/ISC).
    """
    if spline is None:
        return

    # Compute shortened arrow end (30%)
    visible_end_rc = start_rc + (end_rc - start_rc) * 0.3
    rc_path = np.linspace(start_rc, visible_end_rc, 20)

    try:
        energy_path = spline(rc_path)
    except ValueError:
        return

    # Get energies at start and end of full path (not shortened)
    try:
        z_start = float(spline(start_rc))
        z_end = float(spline(end_rc))
    except Exception:
        return

    # === CI/ISC logic: skip arrow if energy increases ===
    if not force_plot and z_end >= z_start:
        return  # Do not draw uphill arrows for CI/ISC segments

    # Project points to surface
    if mesh_type != 'ribbon':
        points = np.array([(rc, 0.0, z) for rc, z in zip(rc_path, energy_path)])
    else:
        points = np.column_stack((rc_path, np.zeros_like(rc_path), energy_path))

    # Draw arrow body
    plotter.add_mesh(pv.Spline(points), color='black', line_width=5)

    # Draw 3D cone arrowhead
    if len(points) >= 2:
        tip_point = points[-1]
        pre_tip = points[-2]
        direction = tip_point - pre_tip
        norm = np.linalg.norm(direction)
        if norm == 0:
            return
        direction /= norm

        # --- Dynamic cone scaling based on RC span ---
        if rc_range_span is None:
            rc_range_span = 5.0  # fallback
        rc_range_span = max(rc_range_span, 0.01)  # avoid divide-by-zero

        # Updated transition logic: between 0.5 and 2.0
        if rc_range_span <= 0.5:
            t = 0.0
        elif rc_range_span >= 2.0:
            t = 1.0
        else:
            t = (rc_range_span - 0.5) / (2.0 - 0.5)

        # Interpolate cone size
        height = 0.8 * (1 - t) + 0.15 * t
        radius = 0.4 * (1 - t) + 0.05 * t

        cone = pv.Cone(
            center=tip_point,
            direction=direction,
            height=height,
            radius=radius,
            resolution=20
        )
        plotter.add_mesh(cone, color='black')


def find_deactivation_path_recursive(start_rc: float, start_state: str, gui_data: dict, interpolated_paths: dict, visited_obstacles: set = None) -> list:
    """
    Finds a realistic deactivation path. It considers all critical points as
    potential barriers and correctly calculates the true barrier height when
    evaluating escape routes from a minimum.
    """
    if visited_obstacles is None:
        visited_obstacles = set()

    # --- Base Case: Stop if the current state is invalid or has no spline ---
    if start_state not in gui_data.state_data or start_state not in interpolated_paths:
        return []
    
    current_path_info = interpolated_paths[start_state]
    state_data = gui_data.state_data[start_state]
    spline = current_path_info.get('spline')
    if spline is None: return []

    path_segments = []
    start_energy = float(spline(start_rc))

    # --- Step 1: Build barrier list and find accessible minima ---
    potential_barriers = []
    potential_barriers.extend(state_data.get('minima', []))
    potential_barriers.extend(state_data.get('TS', []))
    potential_barriers.extend([(c['rc'], c['energy']) for c in state_data.get('crossings', [])])
    if fc_point := state_data.get('franck_condon'):
        potential_barriers.append((fc_point['rc'], fc_point['energy']))

    all_minima = state_data.get('minima', [])
    downhill_minima = [m for m in all_minima if m[1] < start_energy]
    
    accessible_minima = []
    for min_rc, min_e in downhill_minima:
        search_range = sorted([start_rc, min_rc])
        intervening_barriers = [b for b in potential_barriers if search_range[0] < b[0] < search_range[1]]
        is_accessible = not any(b[1] > start_energy for b in intervening_barriers)
        
        if is_accessible:
            accessible_minima.append((min_rc, min_e))

    target_min_rc = None
    if accessible_minima:
        nearest_min = min(accessible_minima, key=lambda m: abs(m[0] - start_rc))
        target_min_rc = nearest_min[0]

    # --- Step 2: Check for intercepting crossings on the path to a minimum ---
    available_crossings = [
        c for c in state_data.get('crossings', [])
        if c['id'] not in visited_obstacles and c['energy'] < start_energy
    ]
    
    if target_min_rc is not None:
       search_range_rc = sorted([start_rc, target_min_rc])
       crossings_on_path = [
           c for c in available_crossings 
           if search_range_rc[0] <= c['rc'] <= search_range_rc[1]
       ]

       if crossings_on_path:
           next_crossing = min(crossings_on_path, key=lambda c: abs(c['rc'] - start_rc))
           path_segments.append({
               "rc_start": start_rc, "rc_end": next_crossing['rc'], "spline": spline,
               "start_state": start_state, "end_state": next_crossing['coupled_state']
           })
           visited_obstacles.add(next_crossing['id'])
           further_path = find_deactivation_path_recursive(
               next_crossing['rc'], next_crossing['coupled_state'], gui_data, interpolated_paths, visited_obstacles
           )
           path_segments.extend(further_path)
           return path_segments

    # --- Step 3: Decay to the destination (Minimum OR Crossing) ---
    if target_min_rc is not None:
        path_segments.append({
            "rc_start": start_rc, "rc_end": target_min_rc, "spline": spline,
            "start_state": start_state, "end_state": start_state
        })
        start_rc = target_min_rc
        start_energy = float(spline(start_rc))
    else:
        accessible_crossings = []
        for cross in available_crossings:
            search_range = sorted([start_rc, cross['rc']])
            intervening_barriers = [b for b in potential_barriers if search_range[0] < b[0] < search_range[1]]
            is_accessible = not any(b[1] > start_energy for b in intervening_barriers)
            if is_accessible:
                accessible_crossings.append(cross)
        
        if accessible_crossings:
            next_crossing = min(accessible_crossings, key=lambda c: abs(c['rc'] - start_rc))
            path_segments.append({
                "rc_start": start_rc, "rc_end": next_crossing['rc'], "spline": spline,
                "start_state": start_state, "end_state": next_crossing['coupled_state']
            })
            visited_obstacles.add(next_crossing['id'])
            further_path = find_deactivation_path_recursive(
                next_crossing['rc'], next_crossing['coupled_state'], gui_data, interpolated_paths, visited_obstacles
            )
            path_segments.extend(further_path)
            return path_segments
            
    # --- Step 4: From a minimum, find the lowest TRUE escape barrier < 0.15 eV ---
    all_ts = state_data.get('TS', [])
    potential_obstacles_from_min = []
    potential_obstacles_from_min.extend(
        [{'rc': ts[0], 'energy': ts[1], 'type': 'TS', 'coupled_state': start_state, 'id': f"TS-{ts[0]}"}
         for ts in all_ts if f"TS-{ts[0]}" not in visited_obstacles]
    )
    potential_obstacles_from_min.extend([
        c for c in state_data.get('crossings', [])
        if c['id'] not in visited_obstacles
    ])

    if not potential_obstacles_from_min:
        return path_segments

    lowest_true_barrier = float('inf')
    best_escape_route = None

    for obstacle in potential_obstacles_from_min:
        # Find the highest point on the path from the minimum to the obstacle.
        search_range = sorted([start_rc, obstacle['rc']])
        path_barriers = [b for b in potential_barriers if search_range[0] < b[0] < search_range[1]]
        
        # The true barrier height is the max energy on the path, or the obstacle itself if the path is clear.
        max_energy_on_path = max([b[1] for b in path_barriers]) if path_barriers else -1
        true_barrier_energy = max(max_energy_on_path, obstacle['energy'])
        
        true_barrier = true_barrier_energy - start_energy
        
        if 0 < true_barrier < lowest_true_barrier:
            lowest_true_barrier = true_barrier
            best_escape_route = obstacle
    
    if best_escape_route and lowest_true_barrier < 0.15:
        path_segments.append({
            "rc_start": start_rc, "rc_end": best_escape_route['rc'], "spline": spline,
            "start_state": start_state, "end_state": best_escape_route['coupled_state']
        })
        visited_obstacles.add(best_escape_route['id'])
        further_path = find_deactivation_path_recursive(
            best_escape_route['rc'], best_escape_route['coupled_state'], gui_data, interpolated_paths, visited_obstacles
        )
        path_segments.extend(further_path)

    return path_segments

def _create_single_2d_plot(ax, interpolated_paths, gui_data, ylim=None):
    """
    Helper function to generate a single 2D plot on a given Matplotlib axis.
    Features:
      - Simulates bi-colored text by aligning two text objects at the center.
      - Energy label is centered relative to the marker, not the text parts.
      - Increased offsets to prevent labels from touching markers.
      - Deterministic placement (no adjust_text) for stability.
      - High z-order (10) for text to ensure it appears ON TOP of markers (zorder 5).
    """
    plotted_crossing_ids = set()
    plot_lines = []
    MARKER_SIZE = 22 
    
    # Text background to ensure readability over lines
    text_bbox = dict(boxstyle='round,pad=0.1', fc='white', alpha=0.7, ec='none')

    # Helper: Get state label with math formatting (e.g., S_1)
    def format_state_label_content(state_name):
        match = re.match(r'([ST])(\d+)', state_name)
        if match:
            return f'${match.group(1)}_{{{match.group(2)}}}$'
        return state_name

    # --- Step 1: Plot Curves ---
    for state, path in interpolated_paths.items():
        state_color = gui_data.state_colors[state]
        raw_data = gui_data.state_data[state]
        
        # Legend label
        label_for_legend = f"{format_state_label_content(state)}"

        if len(path['rc']) > 0:
            line, = ax.plot(path['rc'], path['energy'], label=label_for_legend, color=state_color, linewidth=2.5)
            plot_lines.append(line)

        # --- Step 2: Plot Points & Add Labels ---
        
        # 1. Minima -> Label BELOW (Offset -40)
        if calc_minima := path.get('calc_minima', []):
            for i, (m_rc, m_e, m_e_user) in enumerate(sorted(calc_minima, key=lambda m: m[0])):
                if ylim is None or (ylim[0] <= m_e <= ylim[1]):
                    suffix = f"-{chr(65 + i)}" if len(calc_minima) > 1 else ""
                    label = f"{format_state_label_content(state)}{suffix}\n{m_e_user:.2f}"
                    
                    # Marker (zorder=5)
                    ax.plot(m_rc, m_e, marker='o', color=state_color, markersize=MARKER_SIZE, zorder=5)
                    
                    # Label (zorder=10) -> Draws ON TOP of marker
                    ax.annotate(label, xy=(m_rc, m_e), xytext=(0, -40), 
                                textcoords='offset points', ha='center', va='top',
                                fontsize=18, color=state_color, fontweight='bold',
                                bbox=text_bbox, zorder=10)

        # 2. TS -> Label ABOVE (Offset +30)
        if calc_saddles := path.get('calc_saddles', []):
            for i, (ts_rc, ts_e, ts_e_user) in enumerate(sorted(calc_saddles, key=lambda t: t[0])):
                if ylim is None or (ylim[0] <= ts_e <= ylim[1]):
                    suffix = f"-{chr(65 + i)}" if len(calc_saddles) > 1 else ""
                    label = f"TS{suffix}\n{ts_e_user:.2f}"
                    
                    ax.plot(ts_rc, ts_e, marker='D', color=state_color, markersize=MARKER_SIZE, zorder=5)
                    
                    ax.annotate(label, xy=(ts_rc, ts_e), xytext=(0, 30), 
                                textcoords='offset points', ha='center', va='bottom',
                                fontsize=18, color=state_color, fontweight='bold',
                                bbox=text_bbox, zorder=10)

        # 3. FC -> Label ABOVE (Offset +30)
        if fc := raw_data.get('franck_condon'):
            if ylim is None or (ylim[0] <= fc['energy'] <= ylim[1]):
                label = f"FC {format_state_label_content(state)}\n{fc['energy']:.2f}"
                
                ax.plot(fc['rc'], fc['energy'], marker='s', color=state_color, markersize=MARKER_SIZE, zorder=5)
                
                ax.annotate(label, xy=(fc['rc'], fc['energy']), xytext=(0, 30), 
                            textcoords='offset points', ha='center', va='bottom',
                            fontsize=18, color=state_color, fontweight='bold',
                            bbox=text_bbox, zorder=10)

        # 4. Crossings -> Bi-Colored Label Trick
        if crossings := raw_data.get('crossings', []):
            for cross in crossings:
                if cross['id'] in plotted_crossing_ids: continue
                
                if ylim is None or (ylim[0] <= cross['energy'] <= ylim[1]):
                    s1, s2 = state, cross['coupled_state']
                    
                    if s1 not in gui_data.state_colors or s2 not in gui_data.state_colors:
                        continue 

                    # Determine hierarchy using global helper function
                    s1_order = get_state_ordering(s1)
                    s2_order = get_state_ordering(s2)
                    if s1_order > s2_order:
                        higher, lower = s1, s2
                    else:
                        higher, lower = s2, s1
                    
                    ch = gui_data.state_colors[higher] 
                    cl = gui_data.state_colors[lower]
                    lbl_h = format_state_label_content(higher)
                    lbl_l = format_state_label_content(lower)

                    # --- Markers (zorder=5) ---
                    ax.plot(cross['rc'], cross['energy'], marker='^', color=ch, markersize=MARKER_SIZE, 
                            markeredgecolor='black', markeredgewidth=0.5, zorder=5)
                    ax.plot(cross['rc'], cross['energy'], marker='v', color=cl, markersize=MARKER_SIZE, 
                            markeredgecolor='black', markeredgewidth=0.5, zorder=5)

                    # --- Bi-Color Text Logic (zorder=10) ---
                    # We position everything relative to the MARKER (cross['rc'], cross['energy'])
                    
                    y_offset_text = 45   # Top row (States)
                    y_offset_energy = 25 # Bottom row (Energy)

                    # 1. Left Text: "CI High" (Aligned Right to center 0)
                    # We nudge it left by -1 point to create a tiny visual gap
                    ax.annotate(f"{cross['type']} {lbl_h}", 
                                xy=(cross['rc'], cross['energy']), xytext=(-1, y_offset_text), 
                                textcoords='offset points', ha='right', va='bottom',
                                fontsize=18, color=ch, fontweight='bold',
                                bbox=text_bbox, zorder=10)

                    # 2. Right Text: "/Low" (Aligned Left to center 0)
                    # We nudge it right by +1 point
                    ax.annotate(f"/{lbl_l}", 
                                xy=(cross['rc'], cross['energy']), xytext=(1, y_offset_text), 
                                textcoords='offset points', ha='left', va='bottom',
                                fontsize=18, color=cl, fontweight='bold',
                                bbox=text_bbox, zorder=10)

                    # 3. Energy (Centered below text, relative to MARKER)
                    # This ensures the energy is centered on the point, regardless of text length above
                    ax.annotate(f"{cross['energy']:.2f}", 
                                xy=(cross['rc'], cross['energy']), xytext=(0, y_offset_energy), 
                                textcoords='offset points', ha='center', va='bottom',
                                fontsize=18, color=cl, fontweight='bold',
                                bbox=text_bbox, zorder=10)
                    
                    plotted_crossing_ids.add(cross['id'])

    # --- Step 3: Cleanup ---
    ax.axis('off')
    if ylim:
        ax.set_ylim(*ylim)

    ax.legend(handles=plot_lines, fontsize=20)
    
def plot_2d_matplotlib(interpolated_paths, gui_data, output_dir):
    """
    Creates two 2D plots:
    1. A full plot showing all states.
    2. A "zoomed-in" plot focusing on the excited states if they are
       energetically well-separated from the S0 ground state.
    """
    print("Generating 2D plots...")

    # --- Plot 1: Full Frame ---
    fig_full, ax_full = plt.subplots(figsize=(20, 16))
    _create_single_2d_plot(ax_full, interpolated_paths, gui_data)
    save_path_full = os.path.join(output_dir, 'PES_2D_Matplotlib_Full.png')
    plt.savefig(save_path_full, bbox_inches='tight', pad_inches=0.1, dpi=1000)
    plt.close(fig_full)
    print(f"Full 2D plot saved to {save_path_full}")

    # --- Plot 2: Conditional Zoomed Frame ---
    ZOOM_THRESHOLD = 2.0  # Min energy gap (in eV) to trigger the zoom plot
    ENERGY_BUFFER = 0.5   # Buffer (in eV) above and below the excited states

    s0_data = gui_data.state_data.get('S0')
    if not s0_data or not s0_data.get('minima'):
        print("  - Skipping zoomed plot: S0 minimum not defined.")
        return

    s0_min_energy = min(m[1] for m in s0_data['minima'])
    
    excited_state_energies = []
    for state, data in gui_data.state_data.items():
        if state != 'S0':
            excited_state_energies.extend(p[1] for p in data.get('minima', []))
            excited_state_energies.extend(p[1] for p in data.get('TS', []))
            excited_state_energies.extend(c['energy'] for c in data.get('crossings', []))
            if fc := data.get('franck_condon'):
                excited_state_energies.append(fc['energy'])

    if not excited_state_energies:
        print("  - Skipping zoomed plot: No excited state points found.")
        return

    min_excited_energy = min(excited_state_energies)
    max_excited_energy = max(excited_state_energies)

    if min_excited_energy - s0_min_energy > ZOOM_THRESHOLD:
        print(f"  - Energy gap ({min_excited_energy - s0_min_energy:.2f} eV) > threshold ({ZOOM_THRESHOLD:.2f} eV). Generating zoomed plot.")
        
        fig_zoom, ax_zoom = plt.subplots(figsize=(18, 18))
        
        # Define the y-axis limits for the zoom
        ylim = (min_excited_energy - ENERGY_BUFFER, max_excited_energy + ENERGY_BUFFER)
        
        _create_single_2d_plot(ax_zoom, interpolated_paths, gui_data, ylim=ylim)
        
        save_path_zoom = os.path.join(output_dir, 'PES_2D_Matplotlib_Zoom.png')
        plt.savefig(save_path_zoom, bbox_inches='tight', pad_inches=0.1, dpi=1000)
        plt.close(fig_zoom)
        print(f"Zoomed 2D plot saved to {save_path_zoom}")
    else:
        print("  - Skipping zoomed plot: Energy gap is too small.")

def create_deactivation_animations(output_dir, interpolated_paths, gui_data, mesh_type='ribbon'):
    """
    Creates a GIF animation for each deactivation path found, starting from a
    Franck-Condon point. ONLY PLOTS THE STATES INVOLVED IN THE PATH.
    """
    print("\n--- Creating Deactivation Animations ---")
    
    mesh_creators = {
        'ribbon': create_ribbon_mesh,
        'parabolic': create_parabolic_well_mesh,
        'gaussian': create_gaussian_well_mesh
    }
    mesh_func = mesh_creators.get(mesh_type, create_ribbon_mesh)

    for state_name, data in gui_data.state_data.items():
        if not (fc := data.get("franck_condon")): continue
        
        print(f"Processing animation for deactivation from {state_name}...")
        path_segments = find_deactivation_path_recursive(fc['rc'], state_name, gui_data, interpolated_paths)
        if not path_segments:
            print(f"  Warning: No valid deactivation path found starting from {state_name}. Skipping animation.")
            continue
            
        involved_states = set(seg['start_state'] for seg in path_segments)
        involved_states.add(path_segments[-1]['end_state'])
        print(f"  -> Path involves states: {', '.join(sorted(list(involved_states)))}")

        plotter = pv.Plotter(notebook=False, window_size=[1200, 800], off_screen=True)
        plotter.enable_lightkit()

        print(f"  -> Building static scene for {state_name} animation...")
        # --- Flags to track if point types are present in this specific animation ---
        fc_present, ci_present, isc_present, ts_present = False, False, False, False

        for s_name in involved_states:
            s_data = gui_data.state_data[s_name]
            if s_name in interpolated_paths and len(interpolated_paths[s_name]['rc']) > 0:
                path = interpolated_paths[s_name]
                mesh = mesh_func(path['rc'], path['energy'])
                if mesh:
                    plotter.add_mesh(mesh, color=gui_data.state_colors[s_name], label=s_name, emissive=True, opacity=0.75)
                
                is_well = mesh_type != 'ribbon'
                if s_data.get('franck_condon'):
                    fc_present = True
                    s_fc = s_data.get('franck_condon')
                    center = project_point_to_well(s_fc['rc'], s_fc['energy']) if is_well else (s_fc['rc'], 0.0, s_fc['energy'])
                    plotter.add_mesh(pv.Sphere(radius=0.08, center=center), color='magenta')
                
                if s_data.get('TS'):
                    ts_present = True
                    for ts_rc, ts_e in s_data.get('TS'):
                        center = project_point_to_well(ts_rc, ts_e) if is_well else (ts_rc, 0.0, ts_e)
                        plotter.add_mesh(pv.Sphere(radius=0.08, center=center), color='cyan')

                if s_cross := s_data.get('crossings'):
                    for cross in s_cross:
                        if cross['coupled_state'] in involved_states:
                            if cross['type'] == 'CI':
                                ci_present = True
                                color = 'red'
                            else: # ISC
                                isc_present = True
                                color = 'blue'
                            center = project_point_to_well(cross['rc'], cross['energy']) if is_well else (cross['rc'], 0.0, cross['energy'])
                            plotter.add_mesh(pv.Sphere(radius=0.1, center=center), color=color)
        
        # --- Add proxy actors for the legend, placing them near the action ---
        # Find a reference point from the involved states to place the proxies
        ref_center = (0, 0, 0) # Default fallback
        for s_name in involved_states:
            if s_name in interpolated_paths and len(interpolated_paths[s_name]['rc']) > 0:
                path = interpolated_paths[s_name]
                ref_center = (path['rc'][0], 0.0, path['energy'][0])
                break # Found a valid point, no need to search more

        if fc_present: plotter.add_mesh(pv.Sphere(radius=0.001, center=ref_center), color='magenta', label='FC Point')
        if ts_present: plotter.add_mesh(pv.Sphere(radius=0.001, center=ref_center), color='cyan', label='TS Point')
        if ci_present: plotter.add_mesh(pv.Sphere(radius=0.001, center=ref_center), color='red', label='CI Point')
        if isc_present: plotter.add_mesh(pv.Sphere(radius=0.001, center=ref_center), color='blue', label='ISC Point')

        plotter.add_title(f'Deactivation from {state_name}'); plotter.set_background('white'); plotter.add_legend()
        plotter.view_xz()
        plotter.camera.Azimuth(30) 
        plotter.camera.Elevation(30) 

        gif_path = os.path.join(output_dir, f"PES_3D_animation_{state_name}.gif")
        plotter.open_gif(gif_path)
        print(f"  -> Opened GIF for writing at {gif_path}")
        
        start_rc = path_segments[0]['rc_start']
        start_z = path_segments[0]['spline'](start_rc)
        initial_center = project_point_to_well(start_rc, start_z) if mesh_type != 'ribbon' else (start_rc, 0.0, start_z)

        smooth_sphere = pv.Sphere(radius=0.12, center=initial_center)
        ball_actor = plotter.add_mesh(smooth_sphere, color=gui_data.state_colors[state_name], emissive=True)

        for segment in path_segments:
            ball_actor.prop.color = gui_data.state_colors.get(segment['start_state'], 'white')
            end_point_desc = f"to {segment['end_state']}" if segment['end_state'] != segment['start_state'] else 'to minimum'
            print(f"    ...Animating segment from {segment['start_state']} {end_point_desc}.")
            
            for f in range(20):
                progress = f / 19
                rc_pos = segment['rc_start'] + progress * (segment['rc_end'] - segment['rc_start'])
                try:
                    z_pos = segment['spline'](rc_pos)
                    new_center = project_point_to_well(rc_pos, z_pos) if mesh_type != 'ribbon' else (rc_pos, 0.0, z_pos)
                    smooth_sphere.translate(np.array(new_center) - np.array(smooth_sphere.center), inplace=True)
                    plotter.render(); plotter.write_frame()
                except (ValueError, IndexError):
                    # This segment goes outside the new truncated range, stop animating it.
                    break
        plotter.close(); print(f"  -> Animation saved to {gif_path}")

def get_state_ordering(state_name):
    """Helper function to determine the energy ordering of states."""
    match = re.match(r'([ST])(\d+)', state_name)
    if not match:
        return 0 # S0
    type_val = 0 if match.group(1) == 'S' else 100 # Arbitrary large number to separate S and T
    return type_val + int(match.group(2))

def resolve_collisions_by_truncation(interpolated_paths, gui_data):
    """
    Detects and resolves unwanted spline crossings by truncating one of the splines.

    This function iteratively checks for intersections. If an unwanted crossing is found,
    it decides which state to truncate based on proximity to important features like
    transition states. It then recalculates the truncated spline and repeats the process.
    This version includes a safety check to prevent truncating segments that contain
    user-defined critical points (minima, TS, CIs, ISCs).
    """
    print("\n--- Checking for and resolving spline collisions by truncation ---")
    max_iterations = 10
    
    original_points = {
        state: np.column_stack((data['rc'], data['energy'])) 
        for state, data in gui_data.state_data.items() if len(data['rc']) > 0
    }

    for iteration in range(max_iterations):
        collision_found_in_iter = False
        state_names = sorted(list(interpolated_paths.keys()), key=get_state_ordering)
        
        for state1_name, state2_name in itertools.combinations(state_names, 2):
            if collision_found_in_iter: break
            path1, path2 = interpolated_paths[state1_name], interpolated_paths[state2_name]
            spline1, spline2 = path1.get('spline'), path2.get('spline')

            if spline1 is None or spline2 is None: continue

            min_rc = max(np.min(path1['rc']), np.min(path2['rc'])) if len(path1['rc']) > 0 and len(path2['rc']) > 0 else 0
            max_rc = min(np.max(path1['rc']), np.max(path2['rc'])) if len(path1['rc']) > 0 and len(path2['rc']) > 0 else 0
            if min_rc >= max_rc: continue

            spline_diff = lambda rc: spline1(rc) - spline2(rc)
            test_points = np.linspace(min_rc, max_rc, 500)
            signs = np.sign(spline_diff(test_points))
            
            for i in range(len(signs) - 1):
                if signs[i] != signs[i+1]:
                    try:
                        cross_rc = brentq(spline_diff, test_points[i], test_points[i+1])
                        is_known_crossing = False
                        cross_tolerance = 0.1
                        
                        all_crossings = gui_data.state_data[state1_name].get('crossings', []) + gui_data.state_data[state2_name].get('crossings', [])
                        for cross_info in all_crossings:
                            if cross_info['coupled_state'] in (state1_name, state2_name) and abs(cross_info['rc'] - cross_rc) < cross_tolerance:
                                is_known_crossing = True
                                break
                        
                        if not is_known_crossing:
                            ts1 = gui_data.state_data[state1_name].get('TS', [])
                            ts2 = gui_data.state_data[state2_name].get('TS', [])
                            dist1_to_ts = min([abs(ts[0] - cross_rc) for ts in ts1]) if ts1 else float('inf')
                            dist2_to_ts = min([abs(ts[0] - cross_rc) for ts in ts2]) if ts2 else float('inf')

                            state_to_truncate = state1_name if dist1_to_ts > dist2_to_ts else state2_name
                            if dist1_to_ts == dist2_to_ts:
                                state_to_truncate = state1_name if get_state_ordering(state1_name) > get_state_ordering(state2_name) else state2_name

                            print(f"  -> Unwanted crossing at RC={cross_rc:.2f} between {state1_name} and {state2_name}. Proposing to truncate {state_to_truncate}.")
                            
                            minima = gui_data.state_data[state_to_truncate].get('minima', [])
                            if not minima:
                                print(f"     ...Cannot truncate {state_to_truncate}, has no defined minimum. Skipping.")
                                continue

                            closest_min_rc = min(minima, key=lambda m: abs(m[0] - cross_rc))[0]
                            safety_factor = 0.8
                            new_boundary = closest_min_rc + (cross_rc - closest_min_rc) * safety_factor
                            current_min_rc, current_max_rc = np.min(interpolated_paths[state_to_truncate]['rc']), np.max(interpolated_paths[state_to_truncate]['rc'])

                            # Safety check to prevent removing critical points
                            critical_points_rc = [p[0] for p in gui_data.state_data[state_to_truncate].get('minima', [])] + \
                                                 [p[0] for p in gui_data.state_data[state_to_truncate].get('TS', [])] + \
                                                 [c['rc'] for c in gui_data.state_data[state_to_truncate].get('crossings', [])]
                            if fc := gui_data.state_data[state_to_truncate].get('franck_condon'):
                                critical_points_rc.append(fc['rc'])

                            removed_rc_range = (new_boundary, current_max_rc) if new_boundary > closest_min_rc else (current_min_rc, new_boundary)
                            
                            is_safe_to_truncate = True
                            for crit_rc in critical_points_rc:
                                if removed_rc_range[0] <= crit_rc <= removed_rc_range[1]:
                                    print(f"     ...Aborting truncation: would remove critical point at RC={crit_rc:.2f}.")
                                    is_safe_to_truncate = False
                                    break
                            
                            if not is_safe_to_truncate:
                                continue # Skip this truncation and move on

                            # Perform truncation
                            new_rc_range = (current_min_rc, new_boundary) if new_boundary > closest_min_rc else (new_boundary, current_max_rc)
                            rc_pts, energy_pts = original_points[state_to_truncate][:, 0], original_points[state_to_truncate][:, 1]
                            rc_path, energy_path, spline = interpolation_1d(rc_pts, energy_pts, rc_range=new_rc_range)
                            interpolated_paths[state_to_truncate].update({'rc': rc_path, 'energy': energy_path, 'spline': spline})
                            
                            collision_found_in_iter = True
                            break
                    
                    except ValueError:
                        continue
                if collision_found_in_iter: break
        
        if not collision_found_in_iter:
            print("--- All spline collisions resolved successfully. ---")
            return interpolated_paths
            
    print("--- WARNING: Max iterations reached. Some splines may still cross unintentionally. ---")
    return interpolated_paths

def main():
    """Main function to run the PES visualization script."""
    # --- 1. Parse Command-Line Arguments ---
    parser = argparse.ArgumentParser(description="Generate Potential Energy Surface (PES) plots and animations.", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-f", "--file", type=str, help="Path to a CSV file with PES data.")
    parser.add_argument("-t", "--test", action="store_true", help="Run in TEST_MODE with predefined data.")
    parser.add_argument("-m", "--mesh-type", type=str, default="ribbon", choices=['ribbon', 'parabolic', 'gaussian'], help="Type of 3D mesh for the PES.")
    args = parser.parse_args()

    # --- 2. Get Input Data
    if args.file:
        gui = get_data_from_file(args.file)
    elif args.test:
        gui = get_test_data()
    else:
        print("--- Running in GUI MODE ---")
        root = tk.Tk(); root.withdraw()
        gui = PESInputGUI(root); gui.run()

    # --- 3. Setup Output Directory ---
    if args.file:
        # Use the input file name for the directory name
        base_name = os.path.basename(args.file)
        # Remove the file extension (e.g., .csv)
        file_name_without_ext = os.path.splitext(base_name)[0]
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        output_dir = f"PES_Plot_{file_name_without_ext}_{timestamp}"
    else:
        # Fallback to timestamp for GUI or test mode
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        output_dir = f"PES_Plot_{timestamp}"
    os.makedirs(output_dir, exist_ok=True)
    print(f"\nOutput will be saved to directory: {output_dir}")

    # --- 4. Process and Interpolate Data ---
    print("\n--- Processing and Interpolating State Data ---")
    interpolated_paths = {}
    all_points_for_spline = {} # Store all raw points for collision correction
    for state_name, data in gui.state_data.items():
        if len(data['rc']) == 0: 
            print(f"  - Skipping state {state_name} (no data points).")
            continue
        print(f"  - Processing state {state_name}...")
        
        all_points = list(zip(data['rc'], data['energy']))

        # Add padding points to guide spline at endpoint minima
        if data.get('minima'):
            # Get the reaction coordinates of all significant user-defined points for this state
            significant_rcs = [p[0] for p in all_points]

            for min_rc, min_e in data.get('minima'):
                # Check if there are any OTHER significant points to the left or right

                has_points_left = any(rc < min_rc for rc in significant_rcs)
                has_points_right = any(rc > min_rc for rc in significant_rcs)

                # If a minimum is an endpoint among the significant points, add padding
                if not has_points_left:
                    all_points.append((min_rc - 0.1, min_e + 0.1))
                    print(f"  - Adding left padding point for minimum at RC={min_rc:.2f} on state {state_name}")

                if not has_points_right:
                    all_points.append((min_rc + 0.1, min_e + 0.1))
                    print(f"  - Adding right padding point for minimum at RC={min_rc:.2f} on state {state_name}")
        
        for min_rc, min_e in data.get('minima', []):
            all_points.append((min_rc - 0.01, min_e + 0.0001))
            all_points.append((min_rc + 0.01, min_e + 0.0001))

        for ts_rc, ts_e in data.get('TS', []):
            all_points.append((ts_rc - 0.01, ts_e - 0.0001))
            all_points.append((ts_rc + 0.01, ts_e - 0.0001))
            
        all_points_for_spline[state_name] = all_points
        
        unique_points = sorted(list(set(all_points)), key=lambda p: p[0])
        
        final_rc, final_energy = zip(*unique_points)
        final_rc, final_energy = np.array(final_rc), np.array(final_energy)
        
        rc_path, energy_path, spline = interpolation_1d(final_rc, final_energy)
        
        calculated_minima = [(m[0], m[1], m[1]) for m in data['minima']]
        calculated_saddles = [(ts[0], ts[1], ts[1]) for ts in data['TS']]

        interpolated_paths[state_name] = {
            'rc': rc_path, 
            'energy': energy_path, 
            'y': 0.0, 
            'calc_minima': calculated_minima, 
            'calc_saddles': calculated_saddles,
            'spline': spline
        }
    
    # --- Resolve Collisions ---
    interpolated_paths = resolve_collisions_by_truncation(interpolated_paths, gui)

    # --- 5. Generate Static 3D Plot ---
    print("\n--- Generating Static 3D Plot ---")
    mesh_creators = {
        'ribbon': create_ribbon_mesh,
        'parabolic': create_parabolic_well_mesh,
        'gaussian': create_gaussian_well_mesh
    }
    mesh_func = mesh_creators.get(args.mesh_type, create_ribbon_mesh)
    is_well = args.mesh_type != 'ribbon'

    static_plotter = pv.Plotter(notebook=False, window_size=[1200, 800])
    
    # --- Flags to track point types for the legend ---
    fc_present, ci_present, isc_present, ts_present = False, False, False, False

    for state_name, data in gui.state_data.items():
        if state_name in interpolated_paths and len(interpolated_paths[state_name]['rc']) > 0:
            print(f"  -> Adding mesh and markers for state {state_name}...")
            path_info = interpolated_paths[state_name]
            
            state_color = gui.state_colors[state_name]
            
            # 1. Create the top surface using the original ribbon/parabolic/gaussian function
            #    We don't need the version that bakes vertex colors anymore.
            top_surface = mesh_func(path_info['rc'], path_info['energy'], depth=2.0)

            if top_surface:
                # 2. Extract the surface and extrude to create a solid object
                poly_surface = top_surface.extract_surface()
                extrusion_thickness = 0.1
                solid_mesh = poly_surface.extrude((0, 0, -extrusion_thickness), capping=True)

                # 3. Add the mesh with a SINGLE color and make it FULLY OPAQUE
                static_plotter.enable_lightkit()
                static_plotter.add_mesh(
                    solid_mesh,
                    color=state_color,    # Use one uniform color
                    label=state_name,
                    opacity=1.0,          # Make it fully opaque
                    specular=0.3,         # Add a little shine
                    specular_power=10,
                    smooth_shading=True
                )

            
            if data.get('franck_condon'):
                fc_present = True
                fc_point = data.get('franck_condon')
                center = project_point_to_well(fc_point['rc'], fc_point['energy']) if is_well else (fc_point['rc'], 0.0, fc_point['energy'])
                static_plotter.add_mesh(pv.Sphere(radius=0.08, center=center), color='magenta')
            
            if data.get('TS'):
                ts_present = True
                for ts_rc, ts_e in data.get('TS'):
                    center = project_point_to_well(ts_rc, ts_e) if is_well else (ts_rc, 0.0, ts_e)
                    static_plotter.add_mesh(pv.Sphere(radius=0.08, center=center), color='cyan')
            
            if crossings := data.get('crossings'):
                for cross in crossings:
                    if cross['type'] == 'CI': ci_present = True; color = 'red'
                    else: isc_present = True; color = 'blue'
                    center = project_point_to_well(cross['rc'], cross['energy']) if is_well else (cross['rc'], 0.0, cross['energy'])
                    static_plotter.add_mesh(pv.Sphere(radius=0.1, center=center), color=color)

    print("  -> Adding curved deactivation arrows to static plot...")
    for state_name, data in gui.state_data.items():
        # Draw arrows for deactivation paths starting from Franck-Condon points
        if fc := data.get('franck_condon'):
            path_segments = find_deactivation_path_recursive(fc['rc'], state_name, gui, interpolated_paths)
            for segment in path_segments:
                start_rc = segment['rc_start']
                end_rc = segment['rc_end']
                # Shorten the arrow to 50% of its original length to make it indicative
                shortened_end_rc = start_rc + (end_rc - start_rc) * 0.5
                add_curved_arrow_on_surface(static_plotter, segment['spline'], start_rc, shortened_end_rc, mesh_type=args.mesh_type, force_plot=True)
        
        # Draw arrows from CIs and ISCs as starting points
        for cross in data.get('crossings', []):
            if cross['type'] in ('CI', 'ISC'):
                start_rc = cross['rc']
                start_state = state_name

                path_segments = find_deactivation_path_recursive(
                    start_rc, start_state, gui, interpolated_paths
                )

                for segment in path_segments:
                    start_rc = segment['rc_start']
                    end_rc = segment['rc_end']
                    shortened_end_rc = start_rc + (end_rc - start_rc) * 0.5  # optional shortening
                    add_curved_arrow_on_surface(static_plotter, segment['spline'], start_rc, shortened_end_rc, mesh_type=args.mesh_type, force_plot=False)
                
        # Draw arrows from each TS to its adjacent minima on the same state.
        if ts_points := data.get('TS'):
            if state_name in interpolated_paths and (all_minima := data.get('minima')):
                spline = interpolated_paths[state_name]['spline']
                if spline is None: continue

                for ts_rc, _ in ts_points:
                    # Find minima to the left and right of the current TS
                    minima_left = [m for m in all_minima if m[0] < ts_rc]
                    minima_right = [m for m in all_minima if m[0] > ts_rc]

                    # Draw arrow to the closest minimum on the left
                    if minima_left:
                        closest_left = max(minima_left, key=lambda m: m[0])
                        add_curved_arrow_on_surface(static_plotter, spline, ts_rc, closest_left[0], mesh_type=args.mesh_type, force_plot=True)

                    # Draw arrow to the closest minimum on the right
                    if minima_right:
                        closest_right = min(minima_right, key=lambda m: m[0])
                        add_curved_arrow_on_surface(static_plotter, spline, ts_rc, closest_right[0], mesh_type=args.mesh_type, force_plot= True)

    # --- Add proxy actors for the legend ---
    if fc_present: static_plotter.add_mesh(pv.Sphere(radius=0.001), color='magenta', label='FC Point')
    if ts_present: static_plotter.add_mesh(pv.Sphere(radius=0.001), color='cyan', label='TS Point')
    if ci_present: static_plotter.add_mesh(pv.Sphere(radius=0.001), color='red', label='CI Point')
    if isc_present: static_plotter.add_mesh(pv.Sphere(radius=0.001), color='blue', label='ISC Point')

    static_plotter.add_title('Potential Energy Surface Pathways'); static_plotter.set_background('white'); static_plotter.add_legend()
    static_plotter.view_xz(); static_plotter.camera.Elevation(15)
    gltf_file = os.path.join(output_dir, "PES_3D_plot.gltf")
    print(f"  -> Saving static 3D model to {gltf_file}...")
    static_plotter.export_gltf(gltf_file)
    static_plotter.close()

    # --- 6. Generate Animations and 2D Plot ---
    create_deactivation_animations(output_dir, interpolated_paths, gui, mesh_type=args.mesh_type)
    plot_2d_matplotlib(interpolated_paths, gui, output_dir)
    print("\n--- Script Finished Successfully ---")

# This ensures the main() function is called only when the script is executed directly.
if __name__ == "__main__":
    main()
