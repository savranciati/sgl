# *sgl*: symmetric graphical lasso
R code repository for scripts implementing the algorithm discussed in:

> [1] Ranciati, S., Roverato, A. and Luati, A. (2020). *Fused graphical lasso for brain networks with symmetries.*

## Repository content
### Folder: main
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

---
# Known issues and Bugs

* 
