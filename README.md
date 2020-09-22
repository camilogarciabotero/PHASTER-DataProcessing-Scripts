Data procesing on PHASTER output dataset
================

This repo hosts the complete PHASTER dataset from the 59 strains’
genomes analysed on PHASTER and several scripts developed for its
processing. Analysis of this processed data is published on (paper). The
primary packages used here are in the Tidyverse set of libraries, the
rmarkdown and knitr packages were used to generate this report (Xie
2020; Allaire et al. 2020; Wickham et al. 2019).

## Libraries

``` r
library(tidyverse) 
library(magrittr)
library(readxl) 
library(rmarkdown)
library(knitr)
library(DT)
```

## Data importing and processing

The first step in this analysis is transforming the XSLSX file into a
friendly data frame for R.

``` r
df <- read_excel("Data/2020-09-14_PHASTER-raw.xlsx")
df %<>% mutate_at(vars(e_value), ~ as.numeric(as.character(.)))
df %<>% mutate_at(vars(species,completeness, genome), ~ as_factor(.))

glimpse(df)
```

``` 
  Rows: 12,046
  Columns: 10
  $ species      <fct> Bacillus amyloliquefaciens ATCC 13952, Bacillus amyloliq…
  $ genome       <fct> CP009748.1, CP009748.1, CP009748.1, CP009748.1, CP009748…
  $ completeness <fct> Questionable, Questionable, Questionable, Questionable, …
  $ score        <dbl> 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, …
  $ region       <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2,…
  $ hit_number   <dbl> 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 1…
  $ cds_position <chr> "1114538..1115125", "complement(1115183..1115626)", "com…
  $ blast_hit    <chr> "spore coat protein; KS08_05605", "spore coat protein; K…
  $ e_value      <dbl> NA, NA, NA, NA, NA, NA, 3.77e-05, NA, NA, NA, NA, 1.86e-…
  $ source       <chr> "Genomic", "Genomic", "Genomic", "Genomic", "Genomic", "…
```

## Fig 1C. Generating candiadate prophages organized by completeness (Incomplete, Questionable, Intact)

``` r
by_completeness <- df %>%
  filter(hit_number == 1) %>% 
  group_by(species,genome) %>%
  summarise(
    Incomplete = sum(str_count(completeness, "Incomplete")),
    Questionable = sum(str_count(completeness, "Questionable")),
    Intact = sum(str_count(completeness, "Intact"))
  )

write_tsv(x = by_completeness, path = "Data/completeness_raw-01.tsv")

kable(by_completeness[1:5,])
```

| species                                   | genome         | Incomplete | Questionable | Intact |
| :---------------------------------------- | :------------- | ---------: | -----------: | -----: |
| Bacillus amyloliquefaciens ATCC 13952     | CP009748.1     |          3 |            4 |      2 |
| Bacillus mycoides ATCC 6462               | CP009692.1     |          3 |            5 |      0 |
| Bacillus amyloliquefaciens K2             | MOEA01000001.1 |          4 |            1 |      3 |
| Bacillus anthracis Rock3-42               | CM000732.1     |          2 |            0 |      0 |
| Bacillus subtilis subsp. stercoris D7XPN1 | JHCA01000001.1 |          0 |            0 |      1 |

## Fig 1D. Generating the total number of prophage proteins per bacteria species

``` r
by_species_total <- df %>%
  group_by(species, genome) %>%
  summarise("Total prophage proteins" = sum(str_count(blast_hit, "PHAGE"))) %>% 
  summarize(
    across("Total prophage proteins", mean)
  )

write_tsv(x = by_species_total, path = "Data/total-proteins_raw-01.tsv")

kable(by_species_total[1:5,])
```

| species                                   | Total prophage proteins |
| :---------------------------------------- | ----------------------: |
| Bacillus amyloliquefaciens ATCC 13952     |                     268 |
| Bacillus mycoides ATCC 6462               |                     121 |
| Bacillus amyloliquefaciens K2             |                     189 |
| Bacillus anthracis Rock3-42               |                      49 |
| Bacillus subtilis subsp. stercoris D7XPN1 |                      36 |

## Fig S1. Generating the dataset of the intact prophages in each bacterial species and vizualizing on a bubble-plot.

``` r
intact_phages <- df %>%
  filter(completeness %in% c("Intact") & blast_hit %in% str_subset(blast_hit, "^PHAGE")) %>%
  mutate(Candidate = as.factor(str_extract(blast_hit, "[:graph:]+(?=:)")), region = as.factor(region)) %>%
  mutate(Candidate = str_remove(Candidate, "PHAGE_"),
         Candidate = str_replace(Candidate, "_", " "),
         Candidate = str_replace(Candidate, "Bacill", "Bacillus"),
         Candidate = str_replace(Candidate, "Brevib", "Brevibacillus"),
         Candidate = str_replace(Candidate, "Clostr", "Clostridium"),
         Candidate = str_replace(Candidate, "Thermu", "Thermus phi"),
         Candidate = str_replace(Candidate, "Paenib", "Paenibacillus"),
         Candidate = str_replace(Candidate, "Bacillus 1", "Bacillus virus 1"),
         Candidate = str_replace(Candidate, "Lister", "Listeria"),
         Candidate = str_replace(Candidate, "Geobac", "Geobacillus"),
         Candidate = str_replace(Candidate, "Staphy", "Staphylococcus"),
         ) %>%
  mutate(species = str_replace(species,"Bacillus", "B.")) %>% 
  separate(Candidate, c("Candidate", "Phage_NC_ID"), sep = "_NC_") %>%
  mutate(Candidate = str_remove(Candidate, "_")) %>% 
  select(species, score, region, Candidate, Phage_NC_ID) %>% 
  group_by(species, Candidate, region, score, Phage_NC_ID) %>% 
  summarise(
    CDS_hits = length(Candidate)
  ) %>% 
  group_by(species, region) %>% 
  filter(CDS_hits == max(CDS_hits)) %>% 
  mutate(Epithet = as_factor(word(species,2)))


write_tsv(x = intact_phages, path = "Data/intact-phages_raw-01.tsv")

kable(intact_phages[1:5,])
```

| species                         | Candidate             | region | score | Phage\_NC\_ID | CDS\_hits | Epithet           |
| :------------------------------ | :-------------------- | :----- | ----: | :------------ | --------: | :---------------- |
| B. amyloliquefaciens ATCC 13952 | Bacillus phi105       | 7      |   150 | 004167        |         6 | amyloliquefaciens |
| B. amyloliquefaciens ATCC 13952 | Brevibacillus Jimmer1 | 4      |   150 | 029104        |         9 | amyloliquefaciens |
| B. amyloliquefaciens DSM 7      | Bacillus phi105       | 2      |   150 | 004167        |         3 | amyloliquefaciens |
| B. amyloliquefaciens DSM 7      | Bacillus phi105       | 6      |   140 | 004167        |         6 | amyloliquefaciens |
| B. amyloliquefaciens DSM 7      | Bacillus phi105       | 7      |   130 | 004167        |         9 | amyloliquefaciens |

-----

``` r
intact_phages %>%
  as_tibble() %>%
  mutate(Species = species) %>%
  select(Species, Candidate, CDS_hits, Epithet) %>%
  pivot_wider(names_from = "Candidate", values_from = "CDS_hits", values_fill = 0, values_fn = sum) %>%
  pivot_longer(cols = -c(Species, Epithet), names_to = "Candidates", values_to = "CDS hits") %>%
  filter(`CDS hits` > 0) %>%
  ggplot(aes(Candidates, Species, size = `CDS hits`)) +
  geom_point() +
  facet_grid(rows = vars(Epithet), scales = "free", space = "free") +
  theme_bw() +
  scale_size_area(max_size = 8) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14, face = "italic"),
    axis.text.y = element_text(size = 14, face = "italic"),
    axis.title.x = element_text(size = 17, face = "bold"),
    axis.title.y = element_text(size = 17, face = "bold"),
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.position = "right"
  ) +
  labs(
    x = "Prophage candidates",
    y = "Bacterial species"
  ) +
  ggsave("Figs/bubble-plot-02.pdf", width = 12, height = 15)
```

<img src="README_files/figure-gfm/bubble-plot-1.png" style="display: block; margin: auto;" />

## License

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This
work is licensed under a
<a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative
Commons Attribution 4.0 International License</a>.

## References

<div id="refs" class="references">

<div id="ref-R-rmarkdown">

Allaire, JJ, Yihui Xie, Jonathan McPherson, Javier Luraschi, Kevin
Ushey, Aron Atkins, Hadley Wickham, Joe Cheng, Winston Chang, and
Richard Iannone. 2020. *Rmarkdown: Dynamic Documents for R*.
<https://CRAN.R-project.org/package=rmarkdown>.

</div>

<div id="ref-tidyverse2019">

Wickham, Hadley, Mara Averick, Jennifer Bryan, Winston Chang, Lucy
D’Agostino McGowan, Romain François, Garrett Grolemund, et al. 2019.
“Welcome to the tidyverse.” *Journal of Open Source Software* 4 (43):
1686. <https://doi.org/10.21105/joss.01686>.

</div>

<div id="ref-R-knitr">

Xie, Yihui. 2020. *Knitr: A General-Purpose Package for Dynamic Report
Generation in R*. <https://CRAN.R-project.org/package=knitr>.

</div>

</div>
