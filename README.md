## Longitudinal Tracking of the Cystic Fibrosis Airway Metagenome in Infancy <br>
Marie-Madlen Pust<sup>1,2</sup>, Isa Rudolf<sup>1,2</sup>, Anna-Maria Dittrich<sup>1,2</sup> and Burkhard Tümmler<sup>1,2#</sup> <br> <br>
<sup>1</sup>Department of Paediatric Pneumology, Allergology and Neonatology, Hannover Medical School (MHH), Germany <br>
<sup>2</sup>Biomedical Research in Endstage and Obstructive Lung Disease Hannover (BREATH), German Center for Lung Research, Hannover Medical School, Germany <br>

<br><br>

### Background contamination

### Reference database used for read alignment

### R session
```
R version 4.1.2 (2021-11-01)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.3 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.so.3

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
[10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] usedist_0.4.0        umap_0.2.7.0         RVAideMemoire_0.9-81 stringr_1.4.0        purrr_0.3.4          plyr_1.8.6           tidyr_1.1.3          ggrepel_0.9.1        matrixStats_0.59.0   factoextra_1.0.7     pheatmap_1.0.12     
[12] ggdendro_0.1.22      rcompanion_2.4.1     ggpubr_0.4.0         dplyr_1.0.7          scales_1.1.1         ggplot2_3.3.5        vegan_2.5-7          lattice_0.20-44      permute_0.9-5        readr_1.4.0         

loaded via a namespace (and not attached):
 [1] TH.data_1.0-10     colorspace_2.0-2   ggsignif_0.6.2     ellipsis_0.3.2     class_7.3-20       modeltools_0.2-23  gld_2.6.2          rstudioapi_0.13    proxy_0.4-26       farver_2.1.0       RSpectra_0.16-0    fansi_0.5.0        mvtnorm_1.1-2     
[14] coin_1.4-1         codetools_0.2-18   splines_4.1.2      rootSolve_1.8.2.2  libcoin_1.0-8      knitr_1.33         jsonlite_1.7.2     broom_0.7.8        cluster_2.1.2      png_0.1-7          compiler_4.1.2     backports_1.2.1    assertthat_0.2.1  
[27] Matrix_1.3-4       fastmap_1.1.0      cli_3.0.1          htmltools_0.5.2    tools_4.1.2        gtable_0.3.0       glue_1.4.2         lmom_2.8           Rcpp_1.0.7         carData_3.0-4      vctrs_0.3.8        nlme_3.1-152       lmtest_0.9-38     
[40] xfun_0.24          lifecycle_1.0.0    rstatix_0.7.0      MASS_7.3-54        zoo_1.8-9          hms_1.1.0          parallel_4.1.2     sandwich_3.0-1     expm_0.999-6       RColorBrewer_1.1-2 yaml_2.2.1         Exact_2.1          gridExtra_2.3     
[53] reticulate_1.23    stringi_1.7.3      nortest_1.0-4      e1071_1.7-7        boot_1.3-28        rlang_0.4.11       pkgconfig_2.0.3    evaluate_0.14      labeling_0.4.2     cowplot_1.1.1      tidyselect_1.1.1   ggsci_2.9          magrittr_2.0.1    
[66] R6_2.5.0           DescTools_0.99.42  generics_0.1.0     multcompView_0.1-8 multcomp_1.4-17    DBI_1.1.1          pillar_1.7.0       withr_2.4.2        mgcv_1.8-36        survival_3.2-13    abind_1.4-5        tibble_3.1.6       crayon_1.4.1      
[79] car_3.0-12         utf8_1.2.1         rmarkdown_2.11     grid_4.1.2         data.table_1.14.0  digest_0.6.27      openssl_1.4.4      stats4_4.1.2       munsell_0.5.0      askpass_1.1      
```
