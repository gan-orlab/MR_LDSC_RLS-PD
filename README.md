# MR_LDSC_RLS-PD


**MR_RLS(exposure)-PD(outcome)**
We performed multi-variant two-sample Mendelian randomization to evaluate the effect of RLS on PD. We used GWAS significant SNPs that are associated with RLS (exposure) to create instrumental variables (IVs) in a second GWAS of PD (outcome). 
1. For the construction of our IVs, we selected the summary statistics from the GWAS meta-analysis of RLS using the R package “MRInstruments” which is a part of TwoSampleMR R package.
2.	We calculated F-statistics and R2 in exposure.
3.	For the outcome, we used PD GWAS full summary statistics.
4.	We performed clumping for the exposure and harmonization of datasets.
5.	Using the output from the 'mr' function, we generated a report containing tables and graphs summarising the results. 
6.	Steiger filtering was performed to exclude SNPs that explain more variance in the outcome than in the exposure.
7.	Report included an inverse-variance weighted (IVW) method, MR-Egger method, and Weighted median (WM).
8.	Heterogeneity was tested using Cochran’s Q test in the IVW and MR Egger methods.
9.	We performed MR-PRESSO test to detect horizontal pleiotropy and detect possible outliers
10.	For each method, we constructed funnel plots to be able to manually inspect results and notice pleiotropic outliers.
The MR_RLS(exposure)-PD(outcome).R is a script containing all the code we utilized in our study for MR.


**MR_PD(exposure)-RLS(outcome)**
We performed multi-variant two-sample Mendelian randomization to evaluate the effect of PD on RLS. We used GWAS significant SNPs that are associated with PD (exposure) to create instrumental variables (IVs) in a second GWAS of RLS (outcome).
11.	For the construction of our IVs, we selected the summary statistics from the GWAS meta-analysis of PD using the R package “MRInstruments” which is a part of TwoSampleMR R package.
12.	We calculated F-statistics and R2 in exposure.
13.	For the outcome, we used RLS GWAS full summary statistics.
14.	We performed clumping for the exposure and harmonization of datasets.
15.	Using the output from the 'mr' function, we generated a report containing tables and graphs summarising the results. 
16.	Steiger filtering was performed to exclude SNPs that explain more variance in the outcome than in the exposure.
17.	Report included an inverse-variance weighted (IVW) method, MR-Egger method, and Weighted median (WM).
18.	Heterogeneity was tested using Cochran’s Q test in the IVW and MR Egger methods.
19.	We performed MR PRESSO test to detect horizontal pleiotropy and detect possible outliers
20.	For each method, we constructed funnel plots to be able to manually inspect results and notice pleiotropic outliers.
The MR_PD(exposure)-RLS(outcome).R is a script containing all the code we utilized in our study for MR.

All scripts related to MR are intended for Rstudio. Script related to LDSC is intended for python. Table formatting was performed either in R or in Linux.

**LDSC_RLS-PD**
To assess whether there is genetic correlation between restless legs syndrome (RLS) and Parkinson’s disease (PD), we performed Linkage disequilibrium (LD) score regression (LDSC). 
The LDSR_RLS_PD.py is a script containing all the code we utilized in our study for LDSR.
