Data procesing on PHASTER outputs dataset
================

# Introduction

This repo hosts the complete PHASTER dataset from the 59 strains’
genomes analysed on PHASTER and several scripts developed for its
processing. Analysis of this processed data is published on (paper).

## Libraries

``` r
library(tidyverse)
library(magrittr)
library(readxl) 
library(rmarkdown)
library(knitr)
```

## Data importing and processing

The first step in this analysis is transforming the XSLSX file into a
friendly data frame for R.

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

## Generating candiadates prophage organized by completeness (Incomplete, Questionable, Intact)

``` r
by_completeness <- df %>%
  filter(hit_number == 1) %>% 
  group_by(species,genome) %>%
  summarise(
    Incomplete = sum(str_count(completeness, "Incomplete")),
    Questionable = sum(str_count(completeness, "Questionable")),
    Intact = sum(str_count(completeness, "Intact"))
  )

#write_tsv(x = by_completeness, path = "Data/completeness_raw-01.tsv")

kable(by_completeness[1:5,])
```

| species                                   | genome         | Incomplete | Questionable | Intact |
| :---------------------------------------- | :------------- | ---------: | -----------: | -----: |
| Bacillus amyloliquefaciens ATCC 13952     | CP009748.1     |          3 |            4 |      2 |
| Bacillus mycoides ATCC 6462               | CP009692.1     |          3 |            5 |      0 |
| Bacillus amyloliquefaciens K2             | MOEA01000001.1 |          4 |            1 |      3 |
| Bacillus anthracis Rock3-42               | CM000732.1     |          2 |            0 |      0 |
| Bacillus subtilis subsp. stercoris D7XPN1 | JHCA01000001.1 |          0 |            0 |      1 |

## Generating the total number of prophage proteins per bacteria species

``` r
by_species_total <- df %>%
  group_by(species, genome) %>%
  summarise("Total prophage proteins" = sum(str_count(blast_hit, "PHAGE")))
```

``` 
  `summarise()` regrouping output by 'species' (override with `.groups` argument)
```

``` r
#write_tsv(x = by_species_total, path = "Data/total-proteins_raw-01.tsv")

kable(by_species_total[1:5,])
```

| species                                   | genome         | Total prophage proteins |
| :---------------------------------------- | :------------- | ----------------------: |
| Bacillus amyloliquefaciens ATCC 13952     | CP009748.1     |                     268 |
| Bacillus mycoides ATCC 6462               | CP009692.1     |                     121 |
| Bacillus amyloliquefaciens K2             | MOEA01000001.1 |                     189 |
| Bacillus anthracis Rock3-42               | CM000732.1     |                      49 |
| Bacillus subtilis subsp. stercoris D7XPN1 | JHCA01000001.1 |                      36 |
