Data procesing on PHASTER output dataset
================

This repo hosts the complete PHASTER dataset from the 59 strains’
genomes analysed on PHASTER and several scripts developed for its
processing. Analysis of this processed data is published on (paper). The
primary packages used here are in the Tidyverse set of libraries, the
rmarkdown and knitr packages were used to generate this report (Xie
2021; Allaire et al. 2020; Wickham et al. 2019).

## Libraries

``` r
library(tidyverse)
library(magrittr)
library(readxl)
library(rmarkdown)
library(knitr)
library(DT)
library(patchwork)
```

## Data importing and processing

The first step in this analysis is transforming the XSLSX file into a
friendly data frame for R.

``` r
# Data set from all records from PHASTER detailed data
df_detailed <- read_excel("Data/2021-02-19_PHASTER-raw.xlsx")
df_detailed %<>% mutate_at(vars(e_value), ~ as.numeric(as.character(.)))
df_detailed %<>% mutate_at(vars(species,completeness, genome), ~ as_factor(.))

# Data set from the summary records from PHASTER summary data
df_summaries <- read_tsv("Data/2021-02-19_PHASTER-summaries.tsv")
```

## Fig 1C. Generating candiadate prophages organized by completeness (Incomplete, Questionable, Intact)

``` r
by_completeness <- df_detailed %>%
  filter(hit_number == 1) %>% 
  group_by(species,genome) %>%
  summarise(
    Incomplete = sum(str_count(completeness, "Incomplete")),
    Questionable = sum(str_count(completeness, "Questionable")),
    Intact = sum(str_count(completeness, "Intact"))
  )

write_tsv(x = by_completeness, path = "Data/completeness_raw-02.tsv")
```

      Warning: The `path` argument of `write_tsv()` is deprecated as of readr 1.4.0.
      Please use the `file` argument instead.

``` r
kable(by_completeness[1:5,])
```

| species                               | genome         | Incomplete | Questionable | Intact |
|:--------------------------------------|:---------------|-----------:|-------------:|-------:|
| Bacillus amyloliquefaciens ATCC 13952 | CP009748.1     |          3 |            4 |      2 |
| Bacillus mycoides ATCC 6462           | CP009692.1     |          2 |            0 |      0 |
| Bacillus mycoides ATCC 6470           | CP009692.9     |          8 |            2 |      0 |
| Bacillus amyloliquefaciens K2         | MOEA01000001.1 |          4 |            1 |      3 |
| Bacillus anthracis Rock3-42           | CM000732.1     |          2 |            0 |      0 |

## Fig 1D. Generating the total number of prophage proteins per bacteria species

``` r
by_species_total <- df_detailed %>%
  group_by(species, genome) %>%
  summarise("Total prophage proteins" = sum(str_count(blast_hit, "PHAGE"))) %>% 
  summarize(
    across("Total prophage proteins", mean)
  )

write_tsv(x = by_species_total, path = "Data/total-proteins_raw-02.tsv")

kable(by_species_total[1:5,])
```

| species                               | Total prophage proteins |
|:--------------------------------------|------------------------:|
| Bacillus amyloliquefaciens ATCC 13952 |                     268 |
| Bacillus mycoides ATCC 6462           |                      10 |
| Bacillus mycoides ATCC 6463           |                       0 |
| Bacillus mycoides ATCC 6464           |                       1 |
| Bacillus mycoides ATCC 6465           |                       0 |

## Fig S1. Generating the dataset of the intact prophages in each bacterial species and vizualizing on a bubble-plot.

``` r
phage_family <- tribble(
  ~most_common_phage, ~phage_family,
  'Brevibacillus Jimmer1' , 'Myoviridae',
  'Brevibacillus Jimmer2' , 'Myoviridae',
  'Brevibacillus Abouo', 'Myoviridae',
  'Brevibacillus Jenst' , 'Siphoviridae',
  'Brevibacillus Osiris', 'Myoviridae',
  'Paenibacillus HB10c2' , 'Siphoviridae',
  'Paenibacillus Harrison' , 'Siphoviridae', 
  'Paenibacillus Tripp' , 'Siphoviridae',
  'Deep s_D6E' , 'Myoviridae',
  'Listeria 2389' , 'Siphoviridae',
  'Listeria A006' , 'Siphoviridae',
  'Listeria vB_LmoS_188', 'Siphoviridae',
  'Listeria B054', 'Myoviridae',
  'Thermus phi OH2' , 'Myoviridae',
  'Geobacillus E3', 'Siphoviridae',
  'Geobacillus GBSV1' , 'Myoviridae',
  'Staphylococcus vB_SepS_SEP9' , 'Siphoviridae',
  'Staphylococcus SPbeta_like' , 'Siphoviridae',
  'Clostridium phiCT9441A' , 'Myoviridae',
  'Clostridium phiCD27' , 'Myoviridae',
  'Clostridium phiMMP02' , 'Myoviridae',
  'Clostridium phiCD505', 'Myoviridae',
  'Clostridium phiCT453A', 'Myoviridae',
  'Clostridium phiSM101' , 'Siphoviridae',
  'Clostridium c_st', 'Siphoviridae',
  'Clostridium phiCD111', 'Siphoviridae',
  'Clostridium phiCD211', 'Siphoviridae',
  'Bacillus SPP1' , 'Siphoviridae',
  'Bacillus phi105' , 'Siphoviridae',
  'Bacillus virus 1' , 'Myoviridae',
  'Bacillus Bam35c' , 'Tectiviridae',
  'Bacillus phBC6A51' , 'Siphoviridae',
  'Bacillus phBC6A52' , 'Siphoviridae',
  'Bacillus PM1' , 'Siphoviridae',
  'Bacillus WBeta' , 'Siphoviridae',
  'Bacillus phIS3501' , 'Siphoviridae',
  'Bacillus BtCS33' , 'Siphoviridae',
  'Bacillus SPbeta' , 'Siphoviridae',
  'Bacillus phi4J1' , 'Siphoviridae',
  'Bacillus BalMu1' , 'Myoviridae',
  'Bacillus AR9', 'Myoviridae',
  'Bacillus Bam35c', 'Tectiviridae',
  'Bacillus Bobb' , 'Herelleveridae',
  'Bacillus BtCS33', 'Siphoviridae',
  'Bacillus Eyuki' , 'Herelleveridae',
  'Bacillus G', 'Myoviridae',
  'Bacillus IEBH', 'Siphoviridae',
  'Bacillus JBP901', 'Herelleveridae',
  'Bacillus JL', 'Herelleveridae',
  'Bacillus Palmer', 'Podoviridae',
  'Bacillus PBS1', 'Myoviridae',
  'Bacillus PfEFR_5', 'Siphoviridae',
  'Bacillus phiCM3', 'Siphoviridae',
  'Bacillus phiNIT1', 'Herelleveridae',
  'Bacillus Pony', 'Herelleveridae',
  'Bacillus SP_10', 'Herelleveridae',
  'Bacillus vB_BhaS_171', 'Siphoviridae',
  'Bacillus Waukesha92' , 'Siphoviridae',
  'Bacillus Finn', 'Siphoviridae',
  'Bacillus Fah', 'Siphoviridae',
  'Enterobacteria phi92', 'Myoviridae',
  'Enterococcus phiEf11', 'Siphoviridae',
  'Enterobacteria vB_KleM_RaK2', 'Myoviridae',
  'Lactobacillus JCL1032', 'Siphoviridae',
  'Lactobacillus phiAT3', 'Siphoviridae',
  'Planktothrix PaV_LD', 'Siphoviridae',
  'Sphingomonas PAU', 'Myoviridae',
  'Cellulophaga phiSM', 'Myoviridae',
  'Staphylococcus SpaA1' , 'Siphoviridae',
  'Synechococcus S_SKS1' , 'Myoviridae',
  'Escherichia RCS47' , 'Myoviridae'
)
```

``` r
intact_phages_detailed <- df_detailed %>%
  filter(completeness %in% c("Intact") & blast_hit %in% str_subset(blast_hit, "^PHAGE")) %>%
  mutate(most_common_phage = as.factor(str_extract(blast_hit, "[:graph:]+(?=:)")), region = as.factor(region)) %>%
  mutate(most_common_phage = str_remove(most_common_phage, "PHAGE_"),
         most_common_phage = str_replace(most_common_phage, "_", " "),
         most_common_phage = str_replace(most_common_phage, "Bacill", "Bacillus"),
         most_common_phage = str_replace(most_common_phage, "Brevib", "Brevibacillus"),
         most_common_phage = str_replace(most_common_phage, "Clostr", "Clostridium"),
         most_common_phage = str_replace(most_common_phage, "Thermu", "Thermus phi"),
         most_common_phage = str_replace(most_common_phage, "Paenib", "Paenibacillus"),
         most_common_phage = str_replace(most_common_phage, "Bacillus 1", "Bacillus virus 1"),
         most_common_phage = str_replace(most_common_phage, "Lister", "Listeria"),
         most_common_phage = str_replace(most_common_phage, "Geobac", "Geobacillus"),
         most_common_phage = str_replace(most_common_phage, "Staphy", "Staphylococcus"),
         ) %>%
  mutate(species = str_replace(species,"Bacillus", "B.")) %>% 
  separate(most_common_phage, c("most_common_phage", "phage_nc_id"), sep = "_NC_") %>%
  mutate(most_common_phage = str_remove(most_common_phage, "_")) %>% 
  select(species, score, region, most_common_phage, phage_nc_id) %>% 
  group_by(species, most_common_phage, region, score, phage_nc_id) %>% 
  summarise(
    cds_hits = length(most_common_phage)
  ) %>% 
  group_by(species, region) %>% 
  filter(cds_hits == max(cds_hits)) %>% 
  mutate(epithet = as_factor(word(species,2)))%>%
  inner_join(phage_family)

# write_tsv(x = intact_phages_detailed, path = "Data/intact-phages-detailed_raw-01.tsv")

kable(intact_phages_detailed[1:5,])
```

| species                         | most\_common\_phage   | region | score | phage\_nc\_id | cds\_hits | epithet           | phage\_family |
|:--------------------------------|:----------------------|:-------|------:|:--------------|----------:|:------------------|:--------------|
| B. amyloliquefaciens ATCC 13952 | Bacillus phi105       | 7      |   150 | 004167        |         6 | amyloliquefaciens | Siphoviridae  |
| B. amyloliquefaciens ATCC 13952 | Brevibacillus Jimmer1 | 4      |   150 | 029104        |         9 | amyloliquefaciens | Myoviridae    |
| B. amyloliquefaciens DSM 7      | Bacillus phi105       | 2      |   150 | 004167        |         3 | amyloliquefaciens | Siphoviridae  |
| B. amyloliquefaciens DSM 7      | Bacillus phi105       | 6      |   140 | 004167        |         6 | amyloliquefaciens | Siphoviridae  |
| B. amyloliquefaciens DSM 7      | Bacillus phi105       | 7      |   130 | 004167        |         9 | amyloliquefaciens | Siphoviridae  |

``` r
bubbleplot_plot_detailed <- intact_phages_detailed %>%
  as_tibble() %>%
  mutate(species = species) %>%
  select(species, most_common_phage, cds_hits, epithet, phage_family) %>%
  pivot_wider(names_from = "most_common_phage", values_from = "cds_hits", values_fill = 0, values_fn = sum) %>%
  pivot_longer(cols = -c(species, epithet, phage_family), names_to = "most_common_phages", values_to = "cds_hits") %>%
  filter(cds_hits > 0) %>%
  ggplot(aes(most_common_phages, species, size = cds_hits, color = phage_family)) +
  geom_point() +
  # facet_grid(rows = vars(epithet), scales = "free", space = "free") +
  facet_grid(rows = vars(epithet), cols = vars(phage_family), scales = "free", space = "free") +
  theme_bw() +
  scale_color_brewer(palette = "Dark2") +
  scale_size_area(max_size = 8) +
  guides(color = guide_legend(override.aes = list(size=9))) +
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
    y = "Bacterial species",
    size = "CDS hits",
    color = "Phage family"
  ) +
  ggsave("Figs/bubbleplot-detailed-01.pdf", width = 12, height = 15)

bubbleplot_plot_detailed
```

<img src="README_files/figure-gfm/bubble_plot_detailed-1.png" style="display: block; margin: auto;" />

------------------------------------------------------------------------

Another apporach to find the prophage candidates from PHASTER data
evaluate the summary information. In this apporach the *most common
phage* is directly detected from PHASTER and its criteria is based on
adding unidentified prophage proteins to the PHAGE with most proteins
already identified using sequence similarity of the non identified
proteins.

``` r
intact_phages_summary_cleaned <- df_summaries %>%
  separate(most_common_phage, c("most_common_phage", "phage_proteins"), sep = "\\(") %>%
  relocate(species) %>%
  arrange(species) %>%
  mutate(
    phage_proteins = str_replace(phage_proteins, "\\)", ""),
    most_common_phage = str_remove(most_common_phage, "PHAGE_"),
    species = str_replace(species, "Bacillus", "B."),
    species = as.factor(species),
    phage_proteins = as.numeric(phage_proteins),
    most_common_phage = str_replace(most_common_phage, "_", " "),
    most_common_phage = str_replace(most_common_phage, "Bacill", "Bacillus"),
    most_common_phage = str_replace(most_common_phage, "Brevib", "Brevibacillus"),
    most_common_phage = str_replace(most_common_phage, "Clostr", "Clostridium"),
    most_common_phage = str_replace(most_common_phage, "Thermu", "Thermus phi"),
    most_common_phage = str_replace(most_common_phage, "Paenib", "Paenibacillus"),
    most_common_phage = str_replace(most_common_phage, "Bacillus 1", "Bacillus virus 1"),
    most_common_phage = str_replace(most_common_phage, "Lister", "Listeria"),
    most_common_phage = str_replace(most_common_phage, "Geobac", "Geobacillus"),
    most_common_phage = str_replace(most_common_phage, "Staphy", "Staphylococcus"),
    most_common_phage = str_replace(most_common_phage, "Entero phi92", "Enterobacteria phi92"),
    most_common_phage = str_replace(most_common_phage, "Entero phiEf11", "Enterococcus phiEf11"),
    most_common_phage = str_replace(most_common_phage, "Entero vB_KleM_RaK2", "Enterobacteria vB_KleM_RaK2"),
    most_common_phage = str_replace(most_common_phage, "Lactob", "Lactobacillus"),
    most_common_phage = str_replace(most_common_phage, "Sphing", "Sphingomonas"),
    most_common_phage = str_replace(most_common_phage, "Cellul", "Cellulophaga"),
    most_common_phage = str_replace(most_common_phage, "Plankt", "Planktothrix"),
    most_common_phage = str_replace(most_common_phage, "Synech", "Synechococcus"),
    most_common_phage = str_replace(most_common_phage, "Escher", "Escherichia")
  ) %>%
  separate(most_common_phage, c("most_common_phage", "phage_nc_id"), sep = "_NC_") %>%
  inner_join(phage_family)

bubbleplot_plot_summary <- intact_phages_summary_cleaned %>% 
  filter(completeness == "intact") %>%
  select(species, most_common_phage, phage_proteins, phage_family) %>%
  mutate(epithet = as_factor(word(species,2))) %>% 
  pivot_wider(names_from = "most_common_phage", values_from = "phage_proteins", values_fill = 0, values_fn = sum) %>% 
  pivot_longer(cols = -c(species, epithet, phage_family), names_to = "most_common_phage", values_to = "phage_proteins") %>%
  filter(phage_proteins > 0) %>%
  ggplot(aes(most_common_phage, species, size = phage_proteins, color = phage_family)) +
  geom_point() +
  facet_grid(rows = vars(epithet), scales = "free", space = "free") +
  theme_bw() +
  scale_size_area(max_size = 8) +
  scale_color_brewer(palette = "Dark2") +
  facet_grid(rows = vars(epithet), cols = vars(phage_family), scales = "free", space = "free") +
  guides(color = guide_legend(override.aes = list(size=9))) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 13, face = "italic"),
    axis.text.y = element_text(size = 13, face = "italic"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    strip.text.x = element_blank(),
    legend.position = "right"
  ) +
  labs(
    x = "Prophage candidates",
    y = "Bacterial species",
    size = "CDS hits",
    color = "Phage family"
) +
  ggsave("Figs/bubbleplot-summary-01.pdf", width = 14, height = 16)

bubbleplot_plot_summary
```

![](README_files/figure-gfm/bubble_plot_summary-1.png)<!-- -->

## License

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This
work is licensed under a
<a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative
Commons Attribution 4.0 International License</a>.

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-R-rmarkdown" class="csl-entry">

Allaire, JJ, Yihui Xie, Jonathan McPherson, Javier Luraschi, Kevin
Ushey, Aron Atkins, Hadley Wickham, Joe Cheng, Winston Chang, and
Richard Iannone. 2020. *Rmarkdown: Dynamic Documents for r*.
<https://github.com/rstudio/rmarkdown>.

</div>

<div id="ref-tidyverse2019" class="csl-entry">

Wickham, Hadley, Mara Averick, Jennifer Bryan, Winston Chang, Lucy
D’Agostino McGowan, Romain François, Garrett Grolemund, et al. 2019.
“Welcome to the <span class="nocase">tidyverse</span>.” *Journal of Open
Source Software* 4 (43): 1686. <https://doi.org/10.21105/joss.01686>.

</div>

<div id="ref-R-knitr" class="csl-entry">

Xie, Yihui. 2021. *Knitr: A General-Purpose Package for Dynamic Report
Generation in r*. <https://yihui.org/knitr/>.

</div>

</div>
