# Reaction-Path-Visualizer
A Python tool for visualizing chemical Potential Energy Surfaces (PES), calculating deactivation pathways, and generating 3D animations of deactivation paths, 3D PES and 2D PES.

# PES-Plotter: 3D Potential Energy Surface Visualization

![License](https://img.shields.io/badge/license-MIT-blue.svg) ![Python](https://img.shields.io/badge/python-3.8+-blue.svg)

**PES-Plotter** is an advanced visualization tool designed for computational chemists. It generates interactive 3D visualizations, 2D publication-quality plots, and animated deactivation pathways for complex chemical reactions involving Conical Intersections (CIs) and Intersystem Crossings (ISCs).

## 🌟 Features

* **Interactive 3D Meshes:** Visualize surfaces as Ribbons, Parabolic Wells, or Gaussian Wells using PyVista.
* **Pathway Finding:** Recursive algorithm to trace multi-step deactivation cascades (e.g., S2 → S1 → T1 → S0).
* **Collision Detection:** Automatically detects and resolves spline collisions between states using truncation logic.
* **Publication Ready:** Generates high-resolution 2D plots using Matplotlib with `adjustText` for clean labeling.
* **Animation:** Generates specific GIFs showing the ball-rolling deactivation mechanism.

## 🛠 Installation

1.  **Clone the repository:**
    ```bash
    git clone [https://github.com/YOUR_USERNAME/PES-Plotter.git](https://github.com/YOUR_USERNAME/PES-Plotter.git)
    cd PES-Plotter
    ```

2.  **Install dependencies:**
    ```bash
    pip install -r requirements.txt
    ```

## 🚀 Usage

You can run the script in several modes depending on your needs.

### 1. GUI Mode (Interactive)
This is the default mode. It opens a graphical interface (Tkinter) where you can manually enter the number of states, minima, transition states, and crossing points.
```bash
python src/Final_PES_plotter.py
```
### 2. File Mode (Batch Processing)
Use the -f flag to load reaction data directly from a CSV file. This is ideal for reproducing plots without re-entering data manually.
```bash
python src/Final_PES_plotter.py -f data/example.csv
```
An example CSV file is provided in the data/ folder to show the required format.

### 3. Test Mode
Use the --test flag to run the script with built-in dummy data. This is useful for quickly verifying that the installation is correct and all libraries are working.

```bash

python src/Final_PES_plotter.py --test
```
### 4. Advanced Options
You can customize the visualization style using the command line arguments.

Change Mesh Type: The default mesh is a flat ribbon. You can change this to 3D wells using the --mesh-type argument:

ribbon (Default)

parabolic (Harmonic oscillator style)

gaussian (Gaussian shaped wells)

Example:

```bash

python src/Final_PES_plotter.py --mesh-type gaussian
```
View Help: To see a full list of available arguments and options:

```bash

python src/Final_PES_plotter.py --help
```
## 📂 Output
The script automatically creates a new folder for every run, timestamped to prevent overwriting previous results (e.g., PES_Plot_2023-10-27_14-30-00).

Inside this folder, you will find:

PES_2D_Matplotlib_Full.png: High-resolution static 2D plot of the PES.

PES_2D_Matplotlib_Zoom.png: (Optional) Zoomed-in view of excited states if energetic separation allows.

PES_3D_plot.gltf: An interactive 3D model file. This can be opened in Windows 3D Viewer, Blender, or online GLTF viewers.

PES_3D_animation_*.gif: Individual animated GIFs for every deactivation pathway found (e.g., S2 to S0).

## 📄 License
This project is distributed under the MIT License.

You are free to use, modify, and distribute this software, provided that the original copyright notice and permission notice are included in all copies or substantial portions of the software. See the LICENSE file for more details.
