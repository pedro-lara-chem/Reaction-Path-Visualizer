
# Advanced Potential Energy Surface (PES) Visualization Tool

An advanced Python tool for generating, analyzing, and animating 1D Potential Energy Surface (PES) diagrams for complex chemical reactions. 

Whether you need publication-quality 2D Matplotlib plots, interactive 3D PyVista surfaces, or animated GIFs of multi-step deactivation cascades, this tool provides a comprehensive suite for photochemical and ground-state reaction visualization.

## ✨ Key Features

* **Multi-State Deactivation Cascades**: Features a recursive pathfinding algorithm to automatically trace and visualize complex cascades (e.g., S2 -> S1 -> T1 -> S0).
* **Multi-State Crossings (CIs & ISCs)**: Fully supports Conical Intersections (CIs) and Intersystem Crossings (ISCs) involving *more than two states* simultaneously.
* **Guide Points for Splines**: Add non-critical "guide points" to perfectly shape your cubic spline interpolations without them being flagged as Minima or Transition States.
* **Publication-Quality 2D Plots**: Generates high-DPI 2D Matplotlib graphs. Automatically produces a "zoomed-in" plot of excited states if they are energetically well-separated (>2.0 eV) from the S0 minimum.
* **Multicolored Labels**: Smart 2D plotting automatically generates beautifully aligned, multi-colored labels and split markers for state crossings.
* **Barriers Mode**: Optional flag to dynamically rename all "TS" (Transition State) labels to "Barrier" across all plots and legends.
* **3D Interactive Topologies**: Render PES surfaces as flat `ribbon`, `parabolic`, or `gaussian` meshes in 3D.
* **Automated Collision Resolution**: Intelligently detects and resolves unphysical accidental spline crossings by safely truncating curves before they collide.
* **Deactivation Animations**: Automatically exports `.gif` animations showing a ball rolling down the PES, changing color as it crosses states.


## 🛠 Installation

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/pedro-lara-chem/Reaction-Path-Visualizer.git
    cd Reaction-Path-Visualizer
    ```

2.  **Install dependencies:**
    ```bash
    pip install -r requirements.txt
    ```

## 🚀 Usage

You can run the script via an interactive GUI, or provide data rapidly using a CSV file.

### Basic Commands

**1. Interactive GUI Mode (Default)**
Run without arguments to open the Tkinter GUI and input states, minima, CIs, and ISCs manually.
`python PES_plotter.py`

**2. CSV File Mode**
Bypass the GUI and generate plots instantly from a prepared CSV file.
`python PES_plotter.py -f my_data.csv`

**3. Change 3D Mesh Type**
Choose between `ribbon` (default), `parabolic`, or `gaussian` energy wells for the 3D plots.
`python PES_plotter.py -m gaussian`

**4. Barriers Mode**
Replace the default "TS" nomenclature with "Barrier" in all visual outputs.
`python PES_plotter.py -b`

**5. Test Mode**
Run a hardcoded test scenario (useful for debugging collision resolution).
`python PES_plotter.py -t`

---

## 🛠️ Deep Dive: Advanced Functionalities

### 1. Guide Points
Standard spline interpolation can sometimes behave unpredictably between widely spaced critical points. You can define **Guide Points** (`point_type: guide`). These points strictly anchor the curve's shape without being mathematically treated as a minimum or transition state by the pathfinding algorithm.

### 2. Multi-State ISCs and CIs
Photochemical pathways frequently feature regions where three or more states become degenerate. The tool supports multi-state crossings. In your CSV, simply list the coupled states separated by commas (e.g., `S1, T1, S0`). The recursive pathfinder and 3D animator will correctly route the deactivation cascade through these complex nodes.

### 3. Smart Multicolored Labels (2D Plots)
When generating static 2D plots, overlapping text can ruin a figure. This tool utilizes `adjustText` to prevent label overlaps. Furthermore, for crossings involving multiple states, the tool creates **custom multi-colored text labels and alternating marker shapes** (e.g., a "CI S1/S0" label where "S1" is blue and "S0" is orange, perfectly stacked above the intersection).

### 4. Barriers Mode (`-b`)
Depending on your subfield (e.g., general kinetics vs. strict transition state theory), "Barrier" might be the preferred terminology over "TS" (Transition State). Passing the `-b` or `--barrier-labels` flag dynamically updates the Matplotlib text annotations, PyVista 3D legends, and exported files to use "Barrier".

### 5. Automated Collision Truncation
If the spline of S2 accidentally dips below S1 in a region where you haven't defined a crossing, the script will automatically run a brentq root-finding check. It identifies the unphysical intersection and truncates the least important spline (safeguarding your defined critical points) to ensure state order integrity.

---

## 📄 CSV Formatting Guide

To use the `-f` flag, format your CSV file like the example below. The tool uses Pandas, so it is flexible, but requires specific column headers: `state`, `point_type`, `rc`, `energy`, `coupled_state`, and `id`.

| state | point_type | rc  | energy | coupled_state | id        |
|-------|------------|-----|--------|---------------|-----------|
| S0    | minima     | 0.0 | 0.0    |               |           |
| S1    | FC         | 0.0 | 4.5    |               |           |
| S1    | minima     | 1.5 | 3.2    |               |           |
| S1    | TS         | 0.8 | 3.8    |               |           |
| S1    | guide      | 2.5 | 3.9    |               |           |
| S1    | CI         | 2.0 | 3.5    | "S2, S0"      | CI-S1S0   |
| T1    | ISC        | 1.8 | 3.0    | S1            | ISC-S1-T1 |

* **point_type options:** `minima`, `TS`, `FC`, `guide`, `CI`, `ISC`.
* **coupled_state:** Only required for CIs and ISCs. Can be a single state (`S0`) or comma-separated for multi-state with quotes ("S1, T1").
* **id:** A unique identifier for crossings to prevent duplicating logic when checking both interacting states.

## 📂 Output
The script automatically creates a new folder for every run, timestamped to prevent overwriting previous results (e.g., PES_Plot_2023-10-27_14-30-00).

Inside this folder, you will find:

* **PES_2D_Matplotlib_Full.png**: High-resolution static 2D plot of the PES.

* **PES_2D_Matplotlib_Zoom.png**: (Optional) Zoomed-in view of excited states if energetic separation allows.

* **PES_3D_plot.gltf**: An interactive 3D model file. This can be opened in Windows 3D Viewer, Blender, or online GLTF viewers.

* **PES_3D_animation_*.gif**: Individual animated GIFs for every deactivation pathway found (e.g., S2 to S0).

## 🤝 Contributing
Contributions, issues, and feature requests are welcome! Feel free to check the issues page.

## 📄 License
This project is distributed under the MIT License.

You are free to use, modify, and distribute this software, provided that the original copyright notice and permission notice are included in all copies or substantial portions of the software. See the LICENSE file for more details.
