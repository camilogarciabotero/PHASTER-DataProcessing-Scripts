Data procesing on PHASTER outputs dataset
================

# Introduction

This repo hosts the complete PHASTER dataset from the 59 strains’
genomes analysed on PHASTER and several scripts developed for its
processing. Analysis of this processed data is published on (paper).

# Libraries

# Data importing and processing

The first step in this analysis is transforming the XSLSX file into a
friendly data frame for

``` r
df <- read_excel("Data/07-09-2020_PHASTER-raw.xlsx")
df %<>% mutate_at(vars(e_value), ~ as.numeric(as.character(.)))
df %<>% mutate_at(vars(species,completeness, genome), ~ as_factor(.))

head(df)
```

``` 
  # A tibble: 6 x 10
    species genome completeness score region hit_number cds_position blast_hit
    <fct>   <fct>  <fct>        <dbl>  <dbl>      <dbl> <chr>        <chr>    
  1 Bacill… CP009… Questionable    70      1          1 1114538..11… spore co…
  2 Bacill… CP009… Questionable    70      1          2 complement(… spore co…
  3 Bacill… CP009… Questionable    70      1          3 complement(… spore co…
  4 Bacill… CP009… Questionable    70      1          4 complement(… spore co…
  5 Bacill… CP009… Questionable    70      1          5 complement(… spore co…
  6 Bacill… CP009… Questionable    70      1          6 complement(… spore co…
  # … with 2 more variables: e_value <dbl>, Source <chr>
```
