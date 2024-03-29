This repository contains all material related to the paper:

> [1] Ranciati, S., Roverato, A. Luati, A., (2021), *Fused graphical lasso for brain networks with symmetries.*, Journal of Royal Statistical Society: Series C, 70 (5), 1299- 1322;

### Main files
In the main folder of the repository, users can find the scripts that together implement procedures described in [1].
More specifically:

- **ADMM-symmetric-graphical-lasso.R**: main script, contains the Alternating Descending Method of Multipliers adapted for symmetric graphical lasso (*sgl*);
- **colored-graphical-models-generators.R**: collateral script, it provides functions to identify colored graphs from the concentration matrix estimated via *sgl*
- **misc-functions.R**: miscellanea script, where utilities function are stored.

---

### Folder: "example"
This folder contains the reproducible code for results reported in **Table 1** of [1].

By running the script **fit_sgl_example.R**, the input simulated data (**input_example.RData**) is automatically loaded, and the scripts in the main folder sourced. The results can be visualized by typing in R the name of the object *summary_perf* in the console.

For reproduciblity, users can also find the a saved version of the output if one wants to avoid running the code again (**output_example.RData**).
