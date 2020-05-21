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
Results of **Table 1** reported in [1], saved in .RData workspace format with name "**Table1_results.RData**".
Contents of the Table are stored inside the object "Tab1_res".

For reproduciblity, users can also find the input data (**input_tab1.RData**) and a lite-version of the script (**reproduce_tab1.R**) used in the simulation study.
