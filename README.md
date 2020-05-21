# Symmetric graphical lasso: *sgl*
R code repository for *sgl* scripts implementing the algorithm discussed in:

> [1] Ranciati, S., Roverato, A. and Luati, A. (2020). *Fused graphical lasso for brain networks with symmetries.*

## Content
### Folder: main
In the main folder of the repository, users can find three main scripts that together implement the procedures described in [1].
More specifically:

- **ADMM-symmetric-graphical-lasso.R**: main script, contains the Alternating Descending Method of Multipliers adapted for symmetric graphical lasso (*sgl*);
- **colored-graphical-models-generators.R**: collateral script, it provides functions to identify colored graphs from the concentration matrix estimated via *sgl*
- **misc_functions.R**: miscellanea script, where utilities function are stored.

### Folder: "example"
Simulated dataset used to produce results of **Table 1** reported in [1]; saved in .RData workspace format with name "**Table1_results.RData**".

For reproduciblity, users can also find the input data (**input_tab1.RData**) and a lite-version of the script (**reproduce_tab1.R**) used in the simulation study.
