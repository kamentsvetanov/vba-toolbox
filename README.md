## Voxel-Based Analysis Toolbox: GLM, Mixed-Effects, and Commonality Analysis for Neuroimaging

A flexible MATLAB toolbox for voxel-wise statistical analysis, supporting GLM, mixed-effects models, and commonality analysis across multi-modal brain images.

### Background

This toolbox extends the original commonality analysis approach for neuroimaging, now providing a flexible framework for voxel-wise statistical analyses. Originally developed to support Wu et al., "Cerebral blood flow predicts multiple demand network activity and fluid intelligence across the lifespan", the toolbox now supports:

- General Linear Models (GLM) for voxel-wise regression analyses
- Mixed-effects models for hierarchical or repeated-measures designs
- Commonality analysis to decompose variance explained by multiple predictors

The toolbox works on multi-modal brain images, including ASL, T1w, and functional data. Like the original version, it leverages MATLAB’s fitlm and fitlme functions to compute voxel-wise coefficients, p-values, and residuals, enabling researchers to model both voxel-specific covariates and complex experimental designs.

### User case examples:

- Estimating variance explained by regional and systemic effects in RSFA maps [Tsvetanov et al., 2021](https://onlinelibrary.wiley.com/doi/full/10.1111/psyp.13714)
- Analyzing multi-modal imaging data across subjects and sessions
- Performing advanced voxel-level decomposition of predictor contributions

This new version is designed as a general-purpose, flexible voxel-based analysis toolbox, retaining the original capabilities while supporting a wider range of statistical modeling.

We extended this voxel-wise approach to commonality analysis in [Wu et al 2022](https://www.sciencedirect.com/science/article/pii/S0197458022002044).

![image](./figures/Figure_1.png)


### Dependencies
- [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) 
- [PALM](https://github.com/andersonwinkler/PALM)

Dependencies can be accessed at my [external](https://github.com/kamentsvetanov/CommonalityAnalysis/tree/main/code/external) repo from the subfolders 'spm12' and 'palm'.

The use of other external code in [.../code/external/](https://github.com/kamentsvetanov/CommonalityAnalysis/tree/main/code/external) of this package. 
