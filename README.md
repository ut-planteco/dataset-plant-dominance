# dataset-plant-dominance
Python and R scripts to analyse how abundance of herbaceous plant species can be explained by their above- and belowground traits and estimates of mycorrhizal dependence and flexibility.

# Title: Functional traits as predictors of global dominance and prevalence in herbaceous plants
## Short title: Predicting dominance and prevalence from traits

### Aim
Understanding why some plant species are abundant or widely distributed is a long-standing aim of plant ecology. This research investigated whether the position of herbaceous plant species within the plant economics spectrum (aboveground and belowground) and along the mycorrhizal collaboration gradient could explain their local abundance, geographic occupancy, and global abundance. 

### Methods
We used two published sources of global plant community data (sPlotOpen and the global biodiversity initiative facility; GBIF) to determine whether the local abundance, geographic occupancy and global abundance of herbaceous plant species can be explained by their above- and belowground traits and estimates of mycorrhizal dependence and flexibility.

### Results
Both above- and belowground traits were only weakly associated with local abundance, while geographic occupancy was associated with small plant size and traits indicating fast return on investment in aboveground tissues. Geographic occupancy was also related to belowground traits, being positively associated with specific root length and negatively associated with root diameter, and weakly positively associated with mycorrhizal flexibility. The traits associated with global abundance largely mirrored those associated with occupancy.

### Conclusions
Our analysis suggests that the local success of herbaceous plants is context-specific and there are no universal traits globally underlying high local abundance. By contrast, geographic occupancy aligns strongly with belowground traits, including specific root length and root diameter. The wide success of fine-rooted species may be related to nutrient enrichment during biogeographic history and in recent increasingly anthropogenic conditions.

> This repository contains GBIF, sPlotOpen, belowground, aboveground traits, The GBIF, sPlotOpen, belowground traits and aboveground traits. These datasets can be downloaded from following links: https://doi.org/10.15468/dl.4nqoev, https://doi.org/10.1111/geb.13346, https://doi.org/10.1111/geb.13179, https://doi.org/10.1038/s41586-021-03871-y, https://neo.gsfc.nasa.gov/view.php?datasetId=SEDAC_POP.

> As Human Footprint Index dataset is too big to include, this can be downloaded separately: https://hii-v2-downloads.wcshumanfootprint.org/data/HFP2009.zip

# Running the scripts

Unpack tamme2021.GBIF.spe_tot.csv.zip.001 file, it is packed into multiple volumes, most zip software should handle the unpacking.

First you need to execute python script step1_data_filtering.py

`python3 step1_data_filtering.py`

It will generate local abundance, geographic occpuancy and global abundance based on sPlotOpen dataset and geographic occupancy based on GBIF dataset. These generated text files can be used with the provided R script: `step2_analysis.R`.

To convert lm and glm output into table format, use `step3_convert_lmoutput_to_table.R` script.

## References:
> Carmona CP, Bueno CG, Toussaint A, Träger S, Diaz S, Moora M, Munson AD, Pärtel M, Zobel M, Tamme R. 2021. Fine-root traits in the global spectrum of plant form and function. Nature 597(7878): 683-687.

> Tamme R, Pärtel M, Kõljalg U, Laanisto L, Liira J, Mander Ü, Moora M, Niinemets Ü, Öpik M, Ostonen I, et al. 2021. Global macroecology of nitrogen-fixing plants. Global Ecology and Biogeography 30(2): 514-526.

> Guerrero-Ramirez NR, Mommer L, Freschet GT, Iversen CM, McCormack ML, Kattge J et al. 2021. Global root traits (GRooT) database. Global Ecology and Biogeography 30: 25-37.

> Sabatini FM, Lenoir J, Hattab T, Arnst EA, Chytrý M, Dengler J, De Ruffray P, Hennekens SM, Jandt U, Jansen F, et al. 2021. sPlotOpen - An environmentally balanced, open-access, global dataset of vegetation plots. Global Ecology and Biogeography 30(9): 1740-1764.

> Venter O, Sanderson EW, Magrach A, Allan JR, Beher J, Jones KR, Watson JE, et al. 2016. Global terrestrial Human Footprint maps for 1993 and 2009. Scientific data 3: 1-10.

> Center for International Earth Science Information Network - CIESIN - Columbia University, and Centro Internacional de Agricultura Tropical - CIAT. 2005. Gridded Population of the World, Version 3 (GPWv3): Population Density Grid. Palisades, NY: NASA Socioeconomic Data and Applications Center (SEDAC).
 
