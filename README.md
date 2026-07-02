# COAD_ac4C
Machine learning identifies ac4C-related prognostic signature and TUBA1C as therapeutic target in COAD
Risk score calculation
The acRGBS risk score for each patient was calculated based on the expression levels of four hub genes and their corresponding coefficients derived from the final prognostic model. The formula was defined as follows:

Risk score = -0.3150266 × ExpSARAF -0.2220952 × ExpCDC42SE2 -0.25836728 × ExpTSPYL2 -0.05806783 × ExpTUBA1C

where Exp represents the normalized expression value of each gene and β represents the corresponding model coefficient. Patients were classified into high- and low-acRGBS groups according to the cutoff value determined in the training cohort.
