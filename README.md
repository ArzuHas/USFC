Overview
This MATLAB script performs calculations related to Unified Structural and Functional Connectivity (USFC) modeling for brain connectivity analysis. It processes subject-specific functional (FC) and structural (SC) connectivity matrices to compute various metrics and save the results for further analysis.

Prerequisites
MATLAB (R2021b or later recommended)
Required MATLAB toolboxes and functions:
niftiread (for NIfTI file reading)
Custom functions: importSubjIDs, saveAsEdge
Script Description
Constants Definition:

num_aal: Number of brain regions based on the atlas being used.
listpath, datapath, outputpath: File paths for subject lists, FC matrices, and output directory, respectively.

Data Loading:
Loads AAL atlas indices and functional and structural connectivity matrices.

Center of Mass and Distance Calculation:
Computes and saves the center of mass and distance matrices for the AAL atlas regions.

Main Computation Loop:

For each subject:
Loads FC and SC matrices.
Calculates minimum cost routes and saves matrices related to connectivity and route-specific data.
Saves the final matrices in edge file format for further analysis.
Usage

Set File Paths:
Update listpath, datapath, outputpath, and other file paths in the script with the actual paths to your data directories and files.

Run the Script:
Execute the script in MATLAB. The script will process data for each subject listed in the provided file and save the results in the specified output directories.
Detailed Steps

Loading Data:

AAL_ind_116to90.mat: Atlas with AAL regions compressed to 90 regions.
Functional and structural connectivity matrices for each subject.

Center of Mass Calculation:

Computes the center of mass for each brain region based on the AAL atlas.
Saves results to CMs_AAL_adult_90AAL.mat.

Distance Matrix Calculation:
Computes distance matrices based on center of mass data and saves them to Distance_CMs_AAL_adult_90AAL.mat.

Processing Each Subject:

Loads FC and SC matrices for each subject.
Computes minimum cost routes and route-specific matrices.
Saves the computed matrices in .mat files and .edge format for further analysis.

Output Files
.mat Files:

FC_Matrices_FDR_corrected_90AAL: Functional connectivity matrices.
Structural_Matrices_individual_nothr_90AAL: Structural connectivity matrices.
CMs_AAL_adult_90AAL.mat: Center of mass data.
Distance_CMs_AAL_adult_90AAL.mat: Distance matrix.
Cost_Route_all_Matrix.mat: All route costs and routes.
Route*_Matrix.mat: Route-specific matrices.
USFC_Matrix.mat, USFC_Matrix_abs.mat, USFC_Route_count_Matrix.mat: USFC matrices and route counts.
.edge Files:

SC_M*.edge, FC_M*.edge, USFC.edge, USFC_abs.edge, Route_counted.edge: Connectivity matrices in edge format.
Custom Functions
importSubjIDs: Imports a list of subject IDs.
saveAsEdge: Saves matrices in edge format for further analysis and plotting.
Notes
Ensure that all file paths and file names are updated to match your directory structure and filenames.
If additional custom functions are used, ensure they are in the MATLAB path.
For any issues or further assistance, please contact arzu.hassilemek@cshs.org.

