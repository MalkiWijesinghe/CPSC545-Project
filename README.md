## Dimension Reduction and Clustering Program

### Flow of the Program:

#### Setup:
- Imports necessary libraries for data manipulation, dimensionality reduction, clustering, and visualization.
- Defines functions for preprocessing data, performing various dimensionality reduction techniques, clustering evaluation, generating tSNE plots, and handling output.

#### Argument Parsing:
- Utilizes the argparse library to parse command-line arguments for `base_data_path` and `output_path`.
- Defaults are provided for these paths if no arguments are passed.

#### Function Definitions:
- Several functions are defined for preprocessing, performing dimensionality reduction (e.g., ICA, IPCA, NMF), clustering evaluation, tSNE plotting, and handling dimension reduction outputs.

#### Execution Section:
- Loads scRNA datasets, preprocesses them, retrieves singleR cell annotations, and visualizes the data with tSNE plots.

#### Dimensionality Reduction & Evaluation:
- Iterates through different dimensionality reduction techniques, performs the reduction, evaluates clustering using 10xGenomics and singleR annotations, and generates tSNE plots for each method.

#### Output Handling:
- Produces PDFs containing aggregated plots before and after dimensionality reduction for both 10xGenomics and singleR annotations.


### Execution Instructions:
1. **Environment Setup:**
   - Ensure R is installed along with required libraries (specified in the script).

2. **Running the Program:**
   - Open a terminal in the directory containing the script and execute the following command:
     ```bash
     Rscript PCA_UMAP_tSNE.R
     Rscript dim_reduction_and_clustering.R --base_data_path <path_to_input_data> --output_path <path_to_output_folder>
     ```
     - Replace `<path_to_input_data>` with the directory containing scRNA datasets.
     - Replace `<path_to_output_folder>` with the directory to store the program outputs.
     - If no arguments are passed, default paths are used.

#### File Descriptions:

- **dim_reduction_and_clustering.R:**
  - Main R script performing dimensionality reduction and clustering on scRNA datasets.
  - Contains functions for preprocessing, dimension reduction, evaluation, and visualization.

- **Aggregated Plots:**
  - Output directory storing aggregated PDF plots before and after dimensionality reduction for both 10xGenomics and singleR annotations.

#### Notes:
- Ensure the required packages mentioned in the script are installed to run the program successfully.
- Input scRNA datasets must be in a compatible format for the script to process.
