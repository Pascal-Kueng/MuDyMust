# MuDyMust

For project goals and background information, please refer to this [web-page](https://irhcollaborative.com/?p=6566).

---

## Requirements

- **Git** – Install from [git-scm.com](https://git-scm.com/).
- A Git client for repository management (e.g., [GitHub Desktop](https://desktop.github.com/)).
- **R** version 4.4 or newer.
- **R Studio** 

---

## Setup Instructions

1. **Clone the Repository:**  
   Use [GitHub Desktop](https://desktop.github.com/) or another Git client to clone the repository to your local machine.

2. **Open the Project:**  
   Open the project using the provided R project file, `00_MuDyMust.Rproj`.

3. **Restore the R Environment:**  
   In your R console, run the following commands:
   ```r
   install.packages("renv")
   renv::restore()
   ```

4. **Add Your Datasets:**  
   Place your datasets in the `/Datasets` folder. Ensure that you store them only in this folder to prevent accidental pushes to GitHub.

---

## When Working on the Repository

- *Branching*: Create a personal branch for your work.
- *Pull Updates*: Before you begin working, always pull the latest changes from the main branch.
- *Commits*: Commit your changes regularly and push them to GitHub.
- *Pull Requests*: When your changes are ready, create a pull request and request a code review.

---

