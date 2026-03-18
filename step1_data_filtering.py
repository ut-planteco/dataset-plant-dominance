# Title: Functional traits as predictors of global dominance and prevalence in herbaceous plants
#
# Step 1: Filtering the input data for removeal non-woody and large plants, preparing published datasets for R
#
# This script uses sPlotOpen vegetation plots and GBIF occurrence data to generate input files for the R scripts.
# It gathers information from aboveground and belowground trait datasets. Exclude genus and species from the 
# aboveground and belowground list if they are woody. Woody classification is taken from (i) Carmona et al. 2021 
# woody/non-woody table, (ii) gRooT database (Guerrero-Ramأ­rez et al. 2021) when it is marked as woody, shrub or tree, 
# and (iii) hand picked list of genus that contain only woody species.
# Aboveground trait file will be used to determine the height of the plants, only select shorter than 2 meters.
#
# References:
# Carmona CP, Bueno CG, Toussaint A, Träger S, Diaz S, Moora M, Munson AD, Pärtel M, Zobel M, Tamme R. 2021. Fine-root traits in the global spectrum of plant form and function. Nature 597(7878): 683-687.
# Tamme R, Pärtel M, Kõljalg U, Laanisto L, Liira J, Mander Ü, Moora M, Niinemets Ü, Öpik M, Ostonen I, et al. 2021. Global macroecology of nitrogen-fixing plants. Global Ecology and Biogeography 30(2): 514-526.
# Guerrero-Ramirez NR, Mommer L, Freschet GT, Iversen CM, McCormack ML, Kattge J et al. 2021. Global root traits (GRooT) database. Global Ecology and Biogeography 30: 25-37
# Sabatini FM, Lenoir J, Hattab T, Arnst EA, Chytrý M, Dengler J, De Ruffray P, Hennekens SM, Jandt U, Jansen F, et al. 2021. sPlotOpen - An environmentally balanced, open-access, global dataset of vegetation plots. Global Ecology and Biogeography 30(9): 1740-1764.

from collections import defaultdict
import gzip

# list of plants to include, this initial list includes exclusion of the ones that are larger than 2m
short_plants = {"Heracleum mantegazzianum": True, "Asplenium polyodon": True, "Schoenoplectus acutus": True, "Impatiens pallida": True, "Spartina pectinata": True, "Zea mays": True, "Cladium mariscus": True, "Schoenoplectus tabernaemontani": True, "Phragmites australis": True, "Helianthus tuberosus": True, "Typha angustifolia": True}

header = True
with gzip.open('aboveground_traits.txt.gz', 'rt') as f: # aboveground traits from Carmona et al. 2021
	for line in f:
		if header:
			header = False
		else:
			col = line.strip().replace("\"", "").split(" ")
			if len(col) > 1:
				col[0] = col[0].replace("_", " ")
				if col[3] != "NA" and float(col[3]) <= 2:  # if plant height is provided and is equal or under 2m
					short_plants[col[0]] = True

tree_list = {}
with gzip.open('GRooTFullVersion.tsv.gz', 'rt', encoding = 'ISO-8859-1') as f: # gRooT database from Guerrero-Ramأ­rez et al. 2021
	for line in f:
		col = line.strip().replace("\"", "").split("\t")
		if len(col) > 22:
			if col[20] in ["woody", "shrub", "shrub/tree", "tree"] or col[22] == "woody": # find species that have woody, shrub or tree marked
				name = f"{col[13]} {col[14]}"
				tree_list[name] = True

with gzip.open('woodiness_1218sp.txt.gz', 'rt') as f: # non-woody/woody list of species from Carmona et al. 2021
	for line in f:
		col = line.strip().replace("\"", "").split(" ")
		if len(col) > 1:
			if col[1] == "woody": # find species that have woody marked - these will be igonred for trait analysis
				name = col[0].replace("_", " ")		
				tree_list[name] = True

with gzip.open('tree_genra.txt.gz', 'rt') as f:  # hand picked list of genus that have only woody species - these will be igonred for trait analysis
	for line in f:
		name = line.strip()
		tree_list[name] = True

sp = defaultdict(dict)
header = True
with gzip.open('tamme2021.GBIF.spe_tot.csv.gz', 'rt') as f: # GBIF dataset taken from Tamme et al. 2021, where species names are 
	for line in f:
		if header:
			header = False
			continue
		col = line.strip().replace("\"", "").split(",")
		if len(col) > 1:
			genus = col[1].split(" ")[0]
			if col[1] not in tree_list and genus not in tree_list and col[1] in short_plants:  # only select if not woody and shorter than 2m
				if col[1] in sp and col[2] in sp[col[1]]:
					sp[col[1]][col[2]] += float(col[3])
				else:
					sp[col[1]][col[2]] = float(col[3])

with open("data.gbif.geographic_occupancy.txt", "w") as file: # write for R script filtered GBIF geographic occupancy dataset
	for plant, plots in sp.items():
		count = len(plots)
		file.write(f"{plant}\t{count}\n")

splot_dggs = {}
with gzip.open('splot.to.dggs.txt.gz', 'rt') as f:  # sPlotOpen dataset from Sabatini et al. 2021, contains metadata
	for line in f:
		col = line.strip().replace("\"", "").split("\t")
		if len(col) > 1:
			splot_dggs[col[0]] = col[1]

splot_filter = {}
header = True
with gzip.open('sPlotOpen_header.txt.gz', 'rt') as f:  # sPlotOpen dataset from Sabatini et al. 2021, contains metadata
	for line in f:
		if header:
			header = False
			continue
		col = line.strip().replace("\"", "").split("\t")
		if len(col) > 20:
			# ecosystem type needs to be defined by the sPlotOpen dataset: 18 = forest; 19 = shrubland; 20 = grassland
			# additionally we only pick 1x1-10x10 plots - 1m^2 up to 100m^2
			if (col[18] == "TRUE" or col[19] == "TRUE" or col[20] == "TRUE") and col[10] != "NA" and float(col[10]) >= 1 and float(col[10]) <= 100:
				splot_filter[col[0]] = True

sp = defaultdict(dict)
sp_dggs = defaultdict(dict)
dggs_sp = defaultdict(dict)
header = True
with gzip.open('sPlotOpen_DT.txt.gz', 'rt') as f: # sPlotOpen dataset from Sabatini et al. 2021, contains abundance data
	for line in f:
		if header:
			header = False
			continue
		col = line.strip().replace("\"", "").split("\t")
		if len(col) > 1:
			genus = col[1].split(" ")[0]
			if col[1] not in tree_list and genus not in tree_list and col[0] in splot_filter and col[1] in short_plants and col[4] == "CoverPerc": # use only Cover Percentage values
				if col[3] == "0": # if somehow wrongly added information, default to very small value
					col[3] = 0.01
				if col[1] in sp and col[0] in sp[col[1]]:
					sp[col[1]][col[0]] += float(col[3]) # raw is colmun 3, relative is column 5
				else:
					sp[col[1]][col[0]] = float(col[3])
				# compare GBIF dggridR coordinates with sPlotOpen converted dggridR coordinates
				dggs = splot_dggs[col[0]]
				if col[1] in sp_dggs and dggs in sp_dggs[col[1]]:
					sp_dggs[col[1]][dggs] += float(col[3])
				else:
					sp_dggs[col[1]][dggs] = float(col[3])	
				# Revision: compare GBIF dggridR coordinates with sPlotOpen converted dggridR coordinates, this will add also local abundance and total abundance to the geographic occupancy list
				if dggs not in dggs_sp or col[1] not in dggs_sp[dggs]:
					dggs_sp[dggs][col[1]] = {}
				if col[0] in dggs_sp[dggs][col[1]]:
					dggs_sp[dggs][col[1]][col[0]] += float(col[3])
				else:
					dggs_sp[dggs][col[1]][col[0]] = float(col[3])				

with open("data.splotopen.dggs.geographic_occupancy.txt", "w") as file:	# write for R script filtered sPlotOpen binned hexad geographic occupancy
	for plant, plots in sp_dggs.items():
		count = len(plots)
		file.write(f"{plant}\t{count}\n")

with open("data.splotopen.geographic_occupancy.txt", "w") as file:  # write for R script filtered sPlotOpen geographic occupancy dataset
	for plant, plots in sp.items():
		count = len(plots)
		file.write(f"{plant}\t{count}\n")

with open("data.splotopen.local_abundance.txt", "w") as file:  # write for R script filtered sPlotOpen local abundance dataset	
	for plant, plots in sp.items():
		count = sum(plots.values()) / len(plots)
		file.write(f"{plant}\t{count}\n")

with open("data.splotopen.global_abundance.txt", "w") as file:  # write for R script filtered sPlotOpen global abundance dataset
	for plant, plots in sp.items():
		count = sum(plots.values())
		file.write(f"{plant}\t{count}\n")

res_occ = {}
res_loc = {}
res_glo = {}

for dggs in dggs_sp.keys():
	tot = {}
	# first find all the plots and then calculate the species occurrences
	for species in dggs_sp[dggs].keys():
		for plot in dggs_sp[dggs][species].keys():
			tot[plot] = True
	for species in dggs_sp[dggs].keys():
		cnt = 0
		total = 0
		for plot in dggs_sp[dggs][species].keys():
			cnt += 1
			total += dggs_sp[dggs][species][plot]
		#print(species, " ", cnt, " ", len(tot), " ", cnt / len(tot), " ", total, " ", total / len(tot))
		if species in res_occ:
			res_occ[species] += cnt / len(tot)
		else:
			res_occ[species] = cnt / len(tot)
		if species in res_loc:
			res_loc[species].append(total / len(tot))
		else:
			res_loc[species] = [total / len(tot)]
		if species in res_glo:
			res_glo[species].append((cnt / len(tot)) * (total / len(tot)))
		else:
			res_glo[species] = [(cnt / len(tot)) * (total / len(tot))]

with open("data.splotopen.hexad.geographic_occupancy.txt", "w") as file:# write for R script filtered sPlotOpen binned hexad geographic occupancy 
	for plant, count in res_occ.items():
		file.write(f"{plant}\t{count}\n")

with open("data.splotopen.hexad.local_abundance.txt", "w") as file:# write for R script filtered sPlotOpen binned hexad local abundance
	for plant, count in res_loc.items():
		val = sum(count) / len(count)
		file.write(f"{plant}\t{val}\n")

with open("data.splotopen.hexad.global_abundance.txt", "w") as file:# write for R script filtered sPlotOpen binned hexad global abundance 
	for plant, count in res_glo.items():
		val = sum(count)
		file.write(f"{plant}\t{val}\n")
