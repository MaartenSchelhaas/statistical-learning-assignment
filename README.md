# Statistical Learning Assignment

Group assignment implementing a hybrid kNN classifier in R. The assignment text is in `Assignment statistical learning.pdf` and the skeleton to implement is in `hybrid_kNN.R`.

## Requirements

Install these once:

- R 4.4.1 or newer (https://cran.r-project.org). Everyone should use the same minor version (4.4.x) to avoid package differences.
- RStudio (https://posit.co/download/rstudio-desktop).
- Windows only: Rtools matching your R version (https://cran.r-project.org/bin/windows/Rtools). Only needed if a package has to be compiled from source.

## First-time setup

1. Clone the repository.
2. Open `statistical-learning-assignment.Rproj` (double-click it, or in RStudio use File > Open Project). Always work from the project so the right package library is used.
3. The first time the project opens, renv installs itself automatically. Watch the console for messages.
4. In the R console, run:

   ```r
   renv::restore()
   ```

   This installs the exact package versions listed in `renv.lock` into a library local to this project. Answer `y` if it asks for confirmation. It is the R equivalent of `uv sync`.

5. Check that everything is in sync:

   ```r
   renv::status()
   ```

   It should report "No issues found".

## Daily workflow

After pulling changes from GitHub, run `renv::restore()` if `renv.lock` changed. Otherwise, just work.

### Adding a package

```r
renv::install("packagename")   # install it into the project library
library(packagename)           # use it in your code
renv::snapshot()               # record it in renv.lock
```

Commit the updated `renv.lock` together with the code that uses the package. `install.packages()` also works inside the project, but `renv::install()` is the safer habit.

### Running the code

Open `hybrid_kNN.R` in RStudio and run it with Source, or from the console:

```r
source("hybrid_kNN.R")
```

From a terminal in the project folder (Rscript may need to be on your PATH):

```
Rscript hybrid_kNN.R
```

## Git workflow

- `main` always holds working code. Do not commit to it directly.
- Create a branch per piece of work, for example `yourname/knn-distance`, and open a pull request into `main` when it is ready.
- Pull `main` into your branch regularly to avoid large merge conflicts.
- If two people change `renv.lock` at the same time, resolve the conflict by merging both sets of packages, then run `renv::restore()` followed by `renv::snapshot()`.

## Project files

| File | Purpose |
| --- | --- |
| `hybrid_kNN.R` | kNN skeleton to implement |
| `simulation.R` | Simulation study skeleton (second file to submit) |
| `report/` | Report source and PDF (max 8 pages, include an AI use statement) |
| `results/` | Simulation output (contents are git-ignored) |
| `Assignment statistical learning.pdf` | The assignment |
| `statistical-learning-assignment.Rproj` | RStudio project file, open this |
| `renv.lock` | Locked package versions (commit this) |
| `renv/` | renv activation files (the local library inside is ignored by git) |
| `.Rprofile` | Activates renv when R starts in this folder |

## Troubleshooting

- If `renv::restore()` fails on a package that needs compiling, install Rtools (Windows) and try again.
- If R says the project is out of sync after pulling, run `renv::restore()`.
- If you opened the folder without the `.Rproj` file and renv did not activate, close it and reopen via the `.Rproj` file.
