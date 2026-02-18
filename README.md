# Scale Dependence and Early Warning Signals in the Trophic Structure of Large Marine Ecosystems
This repository contains R scripts to reproduce empirical analyses and figures related to the manuscript "Scale Dependence and Early Warning Signals in the Trophic Structure Large Marine Ecosystems".

# R (v 4.5.1)
See the R folder for relevant analyses used to reproduce figures, analyses, scripts should be run in order.

Independent Scientific Trawl Data for Northwest Atlantic Large Marine Ecosystems provided by K. Frank.

# R Studio Session Info for all packages used throughout R scripts - sessionInfo()
R version 4.5.1 (2025-06-13)
Platform: aarch64-apple-darwin20
Running under: macOS Sequoia 15.6.1

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/Toronto
tzcode source: internal

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] ggrepel_0.9.6      spaa_0.2.5         janitor_2.2.1      readxl_1.4.5       mvtnorm_1.3-3      caret_7.0-1       
 [7] lattice_0.22-7     forecast_8.24.0    DirichletReg_0.7-2 Formula_1.2-5      bbmle_1.0.25.1     see_0.11.0        
[13] report_0.6.1       parameters_0.26.0  performance_0.14.0 modelbased_0.11.2  insight_1.3.0      effectsize_1.0.1  
[19] datawizard_1.1.0   correlation_0.8.7  bayestestR_0.16.0  easystats_0.7.4    beepr_2.0          cowplot_1.1.3     
[25] reshape2_1.4.4     zoo_1.8-14         vegan_2.7-1        permute_0.9-7      lubridate_1.9.4    forcats_1.0.0     
[31] stringr_1.5.1      dplyr_1.1.4        purrr_1.0.4        readr_2.1.5        tidyr_1.3.1        tibble_3.3.0      
[37] ggplot2_3.5.2      tidyverse_2.0.0   

loaded via a namespace (and not attached):
 [1] pROC_1.19.0.1        sandwich_3.1-1       rlang_1.1.6          magrittr_2.0.3       multcomp_1.4-28     
 [6] snakecase_0.11.1     tseries_0.10-58      compiler_4.5.1       mgcv_1.9-3           systemfonts_1.2.3   
[11] vctrs_0.6.5          quadprog_1.5-8       pkgconfig_2.0.3      prodlim_2025.04.28   tzdb_0.5.0          
[16] ragg_1.4.0           miscTools_0.6-28     recipes_1.3.1        cluster_2.1.8.1      R6_2.6.1            
[21] stringi_1.8.7        RColorBrewer_1.1-3   parallelly_1.45.0    rpart_4.1.24         cellranger_1.1.0    
[26] lmtest_0.9-40        numDeriv_2016.8-1.1  estimability_1.5.1   Rcpp_1.0.14          iterators_1.0.14    
[31] future.apply_1.20.0  audio_0.1-11         Matrix_1.7-3         splines_4.5.1        nnet_7.3-20         
[36] timechange_0.3.0     tidyselect_1.2.1     rstudioapi_0.17.1    timeDate_4041.110    maxLik_1.5-2.1      
[41] codetools_0.2-20     curl_6.3.0           listenv_0.9.1        plyr_1.8.9           quantmod_0.4.27     
[46] withr_3.0.2          urca_1.3-4           coda_0.19-4.1        future_1.58.0        survival_3.8-3      
[51] xts_0.14.1           gam_1.22-5           pillar_1.10.2        foreach_1.5.2        TTR_0.24.4          
[56] generics_0.1.4       hms_1.1.3            scales_1.4.0         globals_0.18.0       xtable_1.8-4        
[61] class_7.3-23         glue_1.8.0           emmeans_1.11.1       tools_4.5.1          data.table_1.17.6   
[66] ModelMetrics_1.2.2.2 gower_1.0.2          grid_4.5.1           bdsmatrix_1.3-7      colorspace_2.1-1    
[71] ipred_0.9-15         nlme_3.1-168         fracdiff_1.5-3       cli_3.6.5            textshaping_1.0.1   
[76] viridisLite_0.4.2    lava_1.8.1           gtable_0.3.6         digest_0.6.37        TH.data_1.1-3       
[81] farver_2.1.2         lifecycle_1.0.4      hardhat_1.4.1        MASS_7.3-65         

