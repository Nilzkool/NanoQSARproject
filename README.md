# NanoQSARproject
Repo containing data analysis scripts used in the article **'Identifying nanodescriptors to predict the toxicity of nanomaterials: a case study on titanium dioxide'**

**Instructions on R scripts**
1.	R scripts ‘data_preprocessing_supplementary_files.R’ and ‘data_preprocessing_supplementary_files.R’ works on supplementary data files 1,2 and 3

2.	R scripts ‘data_preprocessing_combined.R’ and ‘data_preprocessing_combined.R’ works on all case studies_1

3.	The pre-processing scripts contain basic pre-processing steps like data-cleaning  and train-test splitting

4.	The data_analysis scripts contain MRMR ranking of input features, stepwise machine-learning modelling with linear regression and random forests, and Y-randomization tests on the best model in order

**Description of case studies**

**Supplementary data file 1** - This case study includes the in vitro data from Murugadoss et al.2020 (https://doi.org/10.1186/s12989-020-00341-7). In this study, cytotoxic effects induced by two TiO2 NMs of identical chemical composition and shape (near-spherical) but with different constituent (primary) particle sizes (17 and 117 nm) was evaluated, and compared their in vitro toxicity in different agglomeration states (small and large agglomerates). Case study 1 includes the measurements of the effects in multiple cell types and at different concentrations on the cell metabolic activity (measured by WST-1 assay), cell viability (lactate dehydrogenase assay), barrier integrity (transepithelial/transendothelial electrical resistance [TEER]), oxidative stress (changes in total glutathione level), pro-inflammatory mediators (interleukin 8 [IL-8], IL-6, tumor necrosis factor α [TNF-α] and IL-1β proteins measured by ELISA) and DNA damage (alkaline comet assay). Cells were exposed to the TiO2 NMs for 24 h in serum-free exposure conditions. Number of data rows (N) = 144.

**Supplementary data file 2** - This case study includes data extracted from 11 articles that were systematically selected from EMBASE database (references are provided in Table S1 in the manuscript). Case study 2 includes the measurements of cell viability (different colorimetric assays) and DNA damage (alkaline comet assay) induced by TiO2 NMs with different constituent particle sizes, crystal phases and surface charge, in multiple cell types and at different concentrations (µg/mL). In these studies, cells were exposed to the TiO2 NMs for 24 h. This case study also includes data from experiments performed with and without serum conditions. N= 49 (DNA damage) and 66 (cell viability).

**Supplementary data file 3** - This case study includes data extracted from a published study Thongkam et al.2017 (10.1093/mutage/gew056) performed in the ENPRA project (RISK ASSESSMENT OF ENGINEERED NANOPARTICLES, European Framework 7). The case study 3 includes measurements of cell viability and DNA damage (alkaline comet assay) induced by TiO2 NMs with different constituent particle sizes and crystal phases, in multiple cell types and at different concentrations (µg/mL). The duration of exposure was 4 h and experiments were performed in the presence of serum (10%) conditions. N = 201 (DNA damage) and 135 (cell viability).

**all case studies_1** - This file includes the measurements of DNA damage and cell viability combined from all the case studies.

