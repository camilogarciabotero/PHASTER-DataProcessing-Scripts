Data procesing on PHASTER outputs dataset
================

# Introduction

This repo hosts the complete PHASTER dataset from the 59 strains’
genomes analysed on PHASTER and several scripts developed for its
processing. Analysis of this processed data is published on (paper).

# Libraries

``` r
library(tidyverse)
library(magrittr)
library(readxl) 
library(rmarkdown)
library(knitr)
```

# Data importing and processing

The first step in this analysis is transforming the XSLSX file into a
friendly data frame for

``` r
df <- read_excel("Data/07-09-2020_PHASTER-raw.xlsx")
df %<>% mutate_at(vars(e_value), ~ as.numeric(as.character(.)))
df %<>% mutate_at(vars(species,completeness, genome), ~ as_factor(.))

kable(df[1:5,])
```

| species                               | genome     | completeness | score | region | hit\_number | cds\_position                | blast\_hit                      | e\_value | Source  |
| :------------------------------------ | :--------- | :----------- | ----: | -----: | ----------: | :--------------------------- | :------------------------------ | -------: | :------ |
| Bacillus amyloliquefaciens ATCC 13952 | CP009748.1 | Questionable |    70 |      1 |           1 | 1114538..1115125             | spore coat protein; KS08\_05605 |       NA | Genomic |
| Bacillus amyloliquefaciens ATCC 13952 | CP009748.1 | Questionable |    70 |      1 |           2 | complement(1115183..1115626) | spore coat protein; KS08\_05610 |       NA | Genomic |
| Bacillus amyloliquefaciens ATCC 13952 | CP009748.1 | Questionable |    70 |      1 |           3 | complement(1115775..1116257) | spore coat protein; KS08\_05615 |       NA | Genomic |
| Bacillus amyloliquefaciens ATCC 13952 | CP009748.1 | Questionable |    70 |      1 |           4 | complement(1116407..1116907) | spore coat protein; KS08\_05620 |       NA | Genomic |
| Bacillus amyloliquefaciens ATCC 13952 | CP009748.1 | Questionable |    70 |      1 |           5 | complement(1117000..1117314) | spore coat protein; KS08\_05625 |       NA | Genomic |
