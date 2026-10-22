# Behavioural and Evolutionary Adaptations in Plant–Pollinator Networks
This repository contains the R simulation framework and analytical pipeline for modelling the eco-evolutionary dynamics of plant–pollinator mutualistic networks. It simulates population dynamics, behavioural adaptation, and evolutionary adaptation, followed by ecological network extraction and topological metric analysis.

```markdown
---

## Repository Structure

| File | Description |
| :--- | :--- |
| `ModelSimulation.R` | Core simulation script implementing population dynamic models, behavioural and evolutionary adaptation mechanisms. |
| `NetworksAndMetricsGeneration.R` | Constructs bipartite plant–pollinator networks from simulation output and calculates ecological network metrics (e.g., persistence, connectance, nestedness, modularity). |
| `PostModelAnalysis.R` | Downstream statistical analysis, summary statistics, hypothesis testing, and data visualizations. |
| `sim_networks.rds` | Pre-computed 500 simulation outputs and network structures for demonstrtation. |

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
2. **Network Generation & Metric Extraction (`NetworksAndMetricsGeneration.R`)**: Converts equilibrium into bipartite network matrices and computes structural metrics.
3. **Postmodel Analysis (`PostModelAnalysis.R`)**: Perform statistical analyses on simulated networks and produce plots.


---

## Clone the Repository via RStudio

You can clone the repository directly inside RStudio:

1. Open **RStudio** and click on **File** > **New Project...**
2. Select **Version Control** ("Create a new project from a version control repository").
3. Select **Git**.
4. Paste the repository URL into the **Repository URL** field:
`https://github.com/mmyahaya/Behavioural-evolutionary-adaptations.git`
5. Choose your local directory for the project and click **Create Project**.


---

## Author

* **Mukhtar Muhammed Yahaya** ([@mmyahaya](https://github.com/mmyahaya))

---

## License

```
MIT License

Copyright (c) 2026 Mukhtar Muhammed Yahaya
```

