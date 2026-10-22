```markdown
# Behavioural and Evolutionary Adaptations in Plant–Pollinator Networks

This repository contains the R simulation framework and analytical pipeline for modelling the eco-evolutionary dynamics of plant–pollinator mutualistic networks. It simulates population dynamics, behavioural adaptation, and evolutionary adaptation, followed by ecological network extraction and topological metric analysis.

---

## Repository Structure

| File | Description |
| :--- | :--- |
| `ModelSimulation.R` | Core simulation script implementing population dynamic models, behavioural and evolutionary adaptation mechanisms. |
| `NetworksAndMetricsGeneration.R` | Constructs bipartite plant–pollinator networks from simulation output and calculates ecological network metrics (e.g., persistence, connectance, nestedness, modularity). |
| `PostModelAnalysis.R` | Downstream statistical analysis, summary statistics, hypothesis testing, and data visualizations. |
| `sim_networks.rds` | Pre-computed 500 simulation outputs and network structures for quick analysis and reproduction. |

---

## Workflow Overview

The pipeline is organised sequentially from raw dynamical modelling to network synthesis and final visualization:

```text
       ┌──────────────────────┐
       │  ModelSimulation.R   │  ◄── Simulate population & eco-evolutionary dynamics
       └──────────┬───────────┘
                  │
                  ▼
┌───────────────────────────────────┐
│  NetworksAndMetricsGeneration.R   │  ◄── Generate networks & compute metrics
└─────────────────┬─────────────────┘     
                  │
                  ▼
       ┌──────────────────────┐
       │  PostModelAnalysis.R │  ◄── Run post-hoc statistics & generate figures
       └──────────────────────┘

```

1. **Simulate Dynamics (`ModelSimulation.R`)**: Runs the equations governing plant and pollinator abundances alongside behavioural and evolutionary over time.
2. **Network Generation & Metric Extraction (`NetworksAndMetricsGeneration.R`)**: Converts equilibrium or time-series interaction weights into bipartite network matrices and computes structural metrics.
3. **Downstream Analysis (`PostModelAnalysis.R`)**: Imports computed metrics (or loads `sim_networks.rds`) to produce summary tables and publication-ready plots.

---

## Prerequisites & Installation

Ensure you have R (>= 4.0.0) installed alongside the required dependency packages:

```r
install.packages(c(
  "bipartite",
  "igraph",
  "deSolve",
  "tidyverse",
  "ggplot2"
))

```

---

## Getting Started in RStudio

### 1. Clone the Repository via RStudio

Instead of using the command line, you can clone the repository directly inside RStudio:

1. Open **RStudio** and click on **File** > **New Project...**
2. Select **Version Control** ("Create a new project from a version control repository").
3. Select **Git**.
4. Paste the repository URL into the **Repository URL** field:
`https://github.com/mmyahaya/Behavioural-evolutionary-adaptations.git`
5. Choose your local directory for the project and click **Create Project**.

### 2. Run the Simulation Pipeline

Once the project is open in RStudio, you can run the scripts sequentially in your R console or source them:

```r
# Step 1: Run the simulation
source("ModelSimulation.R")

# Step 2: Build networks and extract metrics
source("NetworksAndMetricsGeneration.R")

# Step 3: Run downstream post-modelling analysis and plotting
source("PostModelAnalysis.R")

```

### 3. Using Pre-computed Networks

If you want to skip the intensive simulation steps and jump straight into the analysis, load the pre-computed outputs directly:

```r
sim_networks <- readRDS("sim_networks.rds")

```

---

## Author

* **Mukhtar Muhammed Yahaya** ([@mmyahaya](https://github.com/mmyahaya))

---

## License

This project is open-source. Please check the repository settings for specific usage permissions.

```

```
