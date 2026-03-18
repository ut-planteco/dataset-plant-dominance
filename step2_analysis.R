# Title: Functional traits as predictors of global dominance and prevalence in herbaceous plants
#
# Step 2: Uses pre-filtered input data from Python script from published datasets (look references)
#
# This script uses sPlotOpen vegetation plots and GBIF occurrence data that is divided into
# local abundance, geographic occupancy and global abundance. For GBIF, only geographic occupancy
# is used. The abundance and occupancy datasets will be merged with pre-filtered belowground 
# (Guerrero-Ramأ­rez et al. 2021) and aboveground traits datasets. This script includes variety of 
# PCA and univariant analysis along with Human Impact Factor and density graphs.
#
# References:
# Carmona CP, Bueno CG, Toussaint A, Trأ¤ger S, Dأ­az S, Moora M, Munson AD, Pأ¤rtel M, Zobel M, Tamme R. 2021. Fine-root traits in the global spectrum of plant form and function. Nature 597(7878): 683-687.
# Tamme R, Pأ¤rtel M, Kأµljalg U, Laanisto L, Liira J, Mander أœ, Moora M, Niinemets أœ, أ–pik M, Ostonen I, et al. 2021. Global macroecology of nitrogen-fixing plants. Global Ecology and Biogeography 30(2): 514-526.
# Guerrero-Ramأ­rez NR, Mommer L, Freschet GT, Iversen CM, McCormack ML, Kattge J et al. 2021. Global root traits (GRooT) database. Global Ecology and Biogeography 30: 25-37
# Sabatini FM, Lenoir J, Hattab T, Arnst EA, Chytrأ½ M, Dengler J, De Ruffray P, Hennekens SM, Jandt U, Jansen F, et al. 2021. sPlotOpen - An environmentally balanced, open-access, global dataset of vegetation plots. Global Ecology and Biogeography 30(9): 1740-1764
#
# Tested on Windows 10 64BIT (22H2), AMD Ryzen 7 5700X, R 4.1.2

# set the folder to the working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

############################## READ IN LIBRARIES ###############################
library(ggplot2)
library(ggpmisc)
library(iNEXT)
library(PerformanceAnalytics)
library(car)
library(DarkDiv)
library(maptools)
library(maps)
library(npreg)
library(betapart)
library(vegan)
library(factoextra)
library(corrplot)
library(nlme)
library(U.PhyloMaker)
library(psych)
library(betapair)
library(mice)
library(RColorBrewer)
library(dggridR)

############################### DATA PREPARATION ###############################

# set up color palettes
cols = brewer.pal(n = 8, name = "Set1")
cols_alpha = brewer.pal(n = 8, name = "Pastel1")

# color palette with alpha
colorRampAlpha <- function(..., n, alpha) {
    colors <- colorRampPalette(...)(n)
    paste(colors, sprintf("%x", ceiling(255*alpha)), sep="")
}

# convert and hexad bin sPlotOpen coordinates to dggs grids
meta <- read.table("sPlotOpen_header.txt.gz", sep = "\t", header = T)
coords <- cbind.data.frame("plot" = meta$PlotObservationID, "lon" = meta$Longitude, "lat" = meta$Latitude)
dggs <- dgconstruct(res = 8)
coords$cell <- dgGEO_to_SEQNUM(dggs, coords$lon, coords$lat)$seqnum
write.table(coords, "sPlotOpen.dggsgrids.txt", sep = "\t", row.names = F, quote = FALSE)

# Load species list for sPlotOpen and GBIF
species.splotopen <- read.table("data.splotopen.local_abundance.txt", sep = "\t")
colnames(species.splotopen) <- c("species", "frequency")
species.gbif <- read.table("data.gbif.geographic_occupancy.txt", sep = "\t")
colnames(species.gbif) <- c("species", "frequency")

# Generate PCA space including only species present in sPlotOpen or GBIF datasets
# sPlotOpen + belowground traits
pr.data <- read.table("belowground_traits.txt.gz", sep="\t", header = T, row.names = 1)
pr.data <- pr.data[intersect(rownames(pr.data), species.splotopen$species), ]
# show how many species retained for the set
dim(pr.data)
# only retain complete cases and log transform the traits
pr.data.filtered1.splot <- log(na.omit(pr.data) + 1)
# calculate PCA use scaling
pr1.splot <- prcomp(pr.data.filtered1.splot, scale. = TRUE)
biplot(pr1.splot)
# rotate the PCA to be unison with the rest of PCAs
# Some R versions and operating systems can produce different orientation (e.g. one or both
# axes might be mirrored). To keep it consistant with the paper and other PCA outputs, orient 
# belowground RTD to be left and D to be bottom. For aboveground SLA should be top left, LA 
# top right and SSD bottom right.
swap <- pr1.splot$x[,1]
pr1.splot$x[,1] <- pr1.splot$x[,2]
pr1.splot$x[,2] <- swap
swap <- pr1.splot$rotation[,1]
pr1.splot$rotation[,1] <- pr1.splot$rotation[,2]
pr1.splot$rotation[,2] <- swap
biplot(pr1.splot)

# sPlotOpen + aboveground traits
pr.data <- read.table("aboveground_traits.txt.gz", sep=" ", header = T, row.names = 1)
rownames(pr.data) <- gsub("_", " ", rownames(pr.data))
pr.data <- pr.data[intersect(rownames(pr.data), species.splotopen$species), ]
# show how many species retained for the set
dim(pr.data)
# only retain complete cases and log transform the traits
pr.data.filtered2.splot <- log(na.omit(pr.data) + 1)
# calculate PCA use scaling
pr2.splot <- prcomp(pr.data.filtered2.splot, scale. = TRUE)
biplot(pr2.splot)
swap <- pr2.splot$x[,2]
pr2.splot$x[,2] <- -pr2.splot$x[,1]
pr2.splot$x[,1] <- swap
swap <- pr2.splot$rotation[,2]
pr2.splot$rotation[,2] <- -pr2.splot$rotation[,1]
pr2.splot$rotation[,1] <- swap
biplot(pr2.splot)

# GBIF + belowground traits
pr.data <- read.table("belowground_traits.txt.gz", sep="\t", header = T, row.names = 1)
pr.data <- pr.data[intersect(rownames(pr.data), species.gbif$species), ]
# show how many species retained for the set
dim(pr.data)
# only retain complete cases and log transform the traits
pr.data.filtered1.gbif <- log(na.omit(pr.data) + 1)
# calculate PCA use scaling
pr1.gbif <- prcomp(pr.data.filtered1.gbif, scale. = TRUE)
biplot(pr1.gbif)

# GBIF + aboveground traits
pr.data <- read.table("aboveground_traits.txt.gz", sep=" ", header = T, row.names = 1)
rownames(pr.data) <- gsub("_", " ", rownames(pr.data))
pr.data <- pr.data[intersect(rownames(pr.data), species.gbif$species), ]
# show how many species retained for the set
dim(pr.data)
# only retain complete cases and log transform the traits
pr.data.filtered2.gbif <- log(na.omit(pr.data) + 1)
# calculate PCA use scaling
pr2.gbif <- prcomp(pr.data.filtered2.gbif, scale. = TRUE)
biplot(pr2.gbif)
swap <- pr2.gbif$x[,2]
pr2.gbif$x[,2] <- -pr2.gbif$x[,1]
pr2.gbif$x[,1] <- swap
swap <- pr2.gbif$rotation[,2]
pr2.gbif$rotation[,2] <- -pr2.gbif$rotation[,1]
pr2.gbif$rotation[,1] <- swap
biplot(pr2.gbif)

# Combine both PCAs together
pr.data1 <- read.table("belowground_traits.txt.gz", sep="\t", header = T, row.names = 1)
pr.data1 <- pr.data1[intersect(rownames(pr.data1), species.splotopen$species), ]
pr.data2 <- read.table("aboveground_traits.txt.gz", sep=" ", header = T, row.names = 1)
rownames(pr.data2) <- gsub("_", " ", rownames(pr.data2))
pr.data3 <- na.omit(merge(pr.data1, pr.data2, by = 0))
rownames(pr.data3) <- pr.data3$Row.names
pr.data3 <- pr.data3[,-1]
pr.data.filtered3 = log(pr.data3 + 1)
pr3.splot = prcomp(pr.data.filtered3, scale. = TRUE)
biplot(pr3.splot)

################################### TABLE S2 ###################################
# Variance explained for PCAs
message("PCA info: belowground + sPlotOpen")
summary(pr1.splot)  # belowground + sPlotOpen
pr1.splot
pr1.splot$sdev^2
message("PCA info: belowground + GBIF")
summary(pr1.gbif)   # belowground + GBIF
pr1.gbif
pr1.gbif$sdev^2
message("PCA info: aboveground + sPlotOpen")
summary(pr2.splot)  # aboveground + sPlotOpen
pr2.splot
pr2.splot$sdev^2
message("PCA info: belowground + GBIF")
summary(pr2.gbif)   # aboveground + GBIF
pr2.gbif
pr2.gbif$sdev^2

############################# FIGURES 1, 2 and S5 ##############################
############################ TABLE 1, 2, S4 and S5 #############################
sink(file = "lm_output.txt") # write linear model statistics output to file
for(l in 1:3){
  ### Figure 1 configuration ###
  if(l == 1){
    files_frequency <- c(
      "data.splotopen.local_abundance.txt",
      "data.splotopen.local_abundance.txt",
      "data.splotopen.geographic_occupancy.txt",
      "data.splotopen.geographic_occupancy.txt",
      "data.splotopen.global_abundance.txt",
      "data.splotopen.global_abundance.txt"
    )
    files_traits <- c(
      "aboveground_traits.txt.gz",
      "belowground_traits.txt.gz",
      "aboveground_traits.txt.gz",
      "belowground_traits.txt.gz",
      "aboveground_traits.txt.gz",
      "belowground_traits.txt.gz"
    )
    titles <- c("a)", "b)", "c)", "d)", "e)", "f)")
  } else if(l == 2){
    ### Figure 1 configuration ###
    
    ### Figure 2 configuration ###
    files_frequency <- c(
      "data.gbif.geographic_occupancy.txt",
      "data.gbif.geographic_occupancy.txt"
    )
    files_traits <- c(
      "aboveground_traits.txt.gz",
      "belowground_traits.txt.gz"
    )
    titles <- c("a)", "b)")
  } else if(l == 3){
    ### Figure 2 configuration ###
    
    ### Figure S5 configuration ###
    files_frequency <- c(
      "data.splotopen.hexad.local_abundance.txt",
      "data.splotopen.hexad.local_abundance.txt",
      "data.splotopen.hexad.geographic_occupancy.txt",
      "data.splotopen.hexad.geographic_occupancy.txt",
      "data.splotopen.hexad.global_abundance.txt",
      "data.splotopen.hexad.global_abundance.txt"
    )
    files_traits <- c(
      "aboveground_traits.txt.gz",
      "belowground_traits.txt.gz",
      "aboveground_traits.txt.gz",
      "belowground_traits.txt.gz",
      "aboveground_traits.txt.gz",
      "belowground_traits.txt.gz"
    )
    titles <- c("a)", "b)", "c)", "d)", "e)", "f)")
    ### Figure S5 configuration ###
  }
    
  # Generates Figure 1, 2 and S5, showing PCA for local abundance, geographic occupancy and global abundance 
  # using above- and belowground traits for sPlotOpen and GBIF datasets
  if(l == 1){
    ### FIGURE 1 CONFIGURATION ###
    message("Generating graph Fig1")
    pdf("Fig1.pdf", width=12, height=18, pointsize=16)
    par(mfrow=c(3,2))
  } else if(l == 2){
    ### FIGURE 1 CONFIGURATION ###
    message("Generating graph Fig2")
    pdf("Fig2.pdf", width=12, height=6, pointsize=12)
    par(mfrow=c(1,2))
  } else if(l == 3){
    ### FIGURE 2 and S5 CONFIGURATION
    message("Generating graph FigS5")
    pdf("FigS5.pdf", width=12, height=18, pointsize=16)
    par(mfrow=c(3,2))
  }
  par(mar=c(6,6,2,2))
  for(k in 1:length(titles)){
    fn1 <- files_frequency[k]
    fn2 <- files_traits[k]
    
    cat("::", titles[k], "\n") # outputs which graph is calculated 
    message("::", titles[k], " generating graph based on ", fn1, " and ", fn2)
    freq <- read.table(fn1, sep = "\t", row.names = 1)
    colnames(freq) <- c("frequency")
    if(grepl("aboveground", fn2)){
      trait <- read.table(fn2, sep = " ", header = T, row.names = 1)
      rownames(trait) <- gsub("_", " ", rownames(trait))
    } else {
      trait <- read.table(fn2, sep = "\t", header = T, row.names = 1)
    }
    combined <- merge(freq, trait, by = 0)
    rownames(combined) <- combined$Row.names
    combined <- combined[,-1]
    # log transform the frequency and traits
    combined <- log(combined + 1)
    print(dim(combined)) # show the dimensions of frequency and trait intersecting dataset
    combined.imputed <- mice(combined[,-1], m=5, maxit=50, meth='pmm', seed=500, printFlag = FALSE)
    combined.complete <- complete(combined.imputed) # impute missing values, these will be drawn on the figure
  
    # use already existing trait space, which is calculated based on complete cases
    if(grepl("splotopen", fn1) & grepl("below", fn2)){
      pr <- pr1.splot
      data <- pr.data.filtered1.splot
    } else if(grepl("splotopen", fn1) & grepl("above", fn2)){
      pr <- pr2.splot
      data <- pr.data.filtered2.splot
    } else if(grepl("gbif", fn1) & grepl("below", fn2)){
      pr <- pr1.gbif
      data <- pr.data.filtered1.gbif
    } else if(grepl("gbif", fn1) & grepl("above", fn2)){
      pr <- pr2.gbif
      data <- pr.data.filtered2.gbif
    }
    p.values <- c()
    for(predictor in names(trait)){
      l <- lm(as.formula(paste("frequency ~", predictor)), data = combined)
      p.values <- c(p.values, summary(l)[[4]][8])
      print(summary(l))
    }
    message("p-values:", paste(round(p.values, 3), collapse = ",")) # show original p-values and calculate corrected p-values
    message("Corrected p-values (fdr):", paste(round(p.adjust(p.values, "fdr"), 3), collapse = ","))
    
    pts = predict(pr, combined.complete)
    pr.summary <- summary(pr)
    
    # draw PCA plot, do not fill yet
    plot(pr$x[ ,1], pr$x[ ,2],
        xlab = paste("PCA 1 (", round(pr.summary$importance[2,1] * 100, 1), "%)", sep = ""),
        ylab = paste("PCA 2 (", round(pr.summary$importance[2,2] * 100, 1), "%)", sep = ""),
        pch = 16, cex = 2.0, type = "n", cex.lab = 1.5, cex.axis = 1.2, axes = FALSE, main = "") 
    mtext(titles[k], 3, 0.5, adj = 0, cex = 1)
    
    # calculate ordisurf object
    obj = ordisurf(pts ~ combined$frequency, method = "REML", family="gaussian", select = TRUE, add = TRUE, col = "black", lwd.cl = 1, bubble = TRUE)
    
    data.range <- summary(as.vector(obj$grid$z))
    colfunc <- colorRampPalette(c(cols[2], "white", cols[1]))
    data.range.colors = colfunc(20)
    # color the area based on ordisurf object
    for(i in c(1:(length(obj$grid$x) - 1))){
        for(j in c(1:(length(obj$grid$y) - 1))){
            if(is.na(obj$grid$z[i, j])){

            } else {
                data.col <- data.range.colors[floor((obj$grid$z[i, j] - data.range[1]) / (data.range[6] - data.range[1]) * 19) + 1]
                rect(obj$grid$x[i], obj$grid$y[j], obj$grid$x[i + 1], obj$grid$y[j + 1], col = data.col, border = data.col) # coloured
            }
        }
    }
    tmp.pts <- pts
    tmp.prx <- rowSums(round(pr$x, 3))
    tmp.ptsx <- c()
    tmp.ptsy <- c()
    for(m in c(1:(dim(tmp.pts)[1]))){
        if(!(sum(round(tmp.pts[m,], 3)) %in% tmp.prx)){
            tmp.ptsx <- c(tmp.ptsx, tmp.pts[m, 1])
            tmp.ptsy <- c(tmp.ptsy, tmp.pts[m, 2])
        }
    }

    # draw imputed points (change PCH to 4 to distringuish them, otherwise use 16 as the main point)
    points(tmp.ptsx, tmp.ptsy, pch = 4, cex = 0.35, bg = "gray", col = "#777777")
    # draw complete points
    points(pr$x[ ,1], pr$x[ ,2], pch = 16, cex = 0.35, bg = "gray", col = "#777777")
    
    obj <- ordisurf(pts ~ combined$frequency, method = "REML", family="gaussian", select = TRUE, add = TRUE, col = "black", lwd.cl = 1, bubble = TRUE)
    rect(obj$grid$x[1] - 0.5, obj$grid$y[1] - 0.5, obj$grid$x[length(obj$grid$x)] + 0.5, obj$grid$y[length(obj$grid$y)] + 0.5, col = "#FFFFFF77", border = "#FFFFFF77") # coloured

    # draw AXIS
    axis(side = 1, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
    axis(side = 2, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
    #abline(h = 0, v = 0, lty = 3)
    box(lwd = 2)

    envfit <- envfit(pr, data, permutations = 999)
    plot(envfit, cex = 0.9, col="black") 
    
    data.pca <- merge(pts, combined, by = 0)
    #data.pca <- merge(pr$x, combined, by = 0) # if we want to test with imputed values
    l <- lm(frequency ~ PC1, data = data.pca)
    print(summary(l))
    l <- lm(frequency ~ PC2, data = data.pca)
    print(summary(l))
    l <- lm(frequency ~ PC1 * PC2, data = data.pca)
    print(summary(l))
  }
  dev.off()
}
sink(file = NULL) # disable writing the console output to file

########################## TABLE 1, 2 col% and colSD ###########################
files_frequency <- c(
  "data.splotopen.local_abundance.txt",
  "data.splotopen.geographic_occupancy.txt",
  "data.splotopen.global_abundance.txt",
  "data.gbif.geographic_occupancy.txt",
  "data.splotopen.hexad.local_abundance.txt",
  "data.splotopen.hexad.geographic_occupancy.txt",
  "data.splotopen.hexad.global_abundance.txt"
)
fr.data <- read.table("FR_species.csv", sep = ",", header = T)
sink(file = "lm_output-2.txt")
for(k in 1:length(files_frequency)){
  fn1 <- files_frequency[k]
  cat("::", fn1, "\n")
  freq <- read.table(fn1, sep = "\t", row.names = 1)
  colnames(freq) <- c("frequency")
  freq$frequency <- log(freq$frequency + 1)
  combined <- merge(freq, fr.data, by.x = 0, by.y = 2)
  combined
  p.values <- c()
  l <- lm(frequency ~ sd_col, combined)
  print(summary(l))
  p.values <- c(p.values, summary(l)[[4]][8])
  l <- lm(frequency ~ mean_col, combined)
  print(summary(l))
  p.values <- c(p.values, summary(l)[[4]][8])
  message("p-values:", paste(round(p.values, 3), collapse = ","))
  message("Corrected p-values (fdr):", paste(round(p.adjust(p.values, "fdr"), 3), collapse = ","))
}
sink(file = NULL)

################################## FIGURE S2 ###################################
# select all species from sPlotOpen and GBIF
all.species <- unique(c(species.splotopen$species, species.gbif$species))
# FIGURE S2
pr.data <- read.table("aboveground_traits.txt.gz", sep=" ", header = T, row.names = 1)
rownames(pr.data) <- gsub("_", " ", rownames(pr.data))
pr.data <- log(pr.data + 1)
pr.col <- read.table("FR_species.csv", sep=",", header = T, row.names = 2)
pr.col <- pr.col[,which(colnames(pr.col) %in% c("sd_col", "mean_col"))]
colnames(pr.col) <- c("colSD", "col%")
pr.combined <- merge(pr.data, pr.col, by = 0, all = T)
rownames(pr.combined) <- pr.combined$Row.names
pr.combined <- pr.combined[,-1]
pr.combined <- pr.combined[intersect(rownames(pr.combined), all.species), ]
pdf("FigS2.pdf", width=12, height=12, pointsize=16)
chart.Correlation(pr.combined, histogram=TRUE, pch=19)
dev.off()

################################## FIGURE S3 ###################################
pr.data <- read.table("belowground_traits.txt.gz", sep="\t", header = T, row.names = 1)
pr.data <- log(pr.data + 1)
pr.col <- read.table("FR_species.csv", sep=",", header = T, row.names = 2)
pr.col <- pr.col[,which(colnames(pr.col) %in% c("sd_col", "mean_col"))]
colnames(pr.col) <- c("colSD", "col%")
pr.combined <- merge(pr.data, pr.col, by = 0, all = T)
rownames(pr.combined) <- pr.combined$Row.names
pr.combined <- pr.combined[,-1]
pr.combined <- pr.combined[intersect(rownames(pr.combined), all.species), ]
pdf("FigS3.pdf", width=12, height=12, pointsize=16)
chart.Correlation(pr.combined, histogram=TRUE, pch=19)
dev.off()

################################## FIGURE S4 ###################################
pdf("FigS4.pdf", width=12, height=12, pointsize=16)
titles = c("a)", "b)", "c)", "d)")
par(mfrow = c(1, 4))
for(k in 1:4){
  if(k == 2){ # sPlotOpen + belowground
    pr <- pr1.splot
  } else if(k == 1){ # sPlotOpen + aboveground
    pr <- pr2.splot
  } else if(k == 4){ # GBIF + belowground
    pr <- pr1.gbif
  } else if(k == 3){ # GBIF + aboveground
    pr <- pr2.gbif
  }
  colFunc = colorRampPalette(c(cols[2], "#ffffff", cols[1]))
  corrplot(
    pr$rotation[,c(1:2)],
    method = "color", # Correlation plot method
    type = "full", 
    diag = TRUE,     
    tl.col = "black", 
    col = colFunc(40),
    bg = "white", 
    addshade = NULL,
    cl.pos = "n")
  text(x = par("usr")[1], y = par("usr")[4] - 1.5, labels = titles[k], pos = 4, xpd = NA, cex = 2)
}
dev.off()
  

sink(file = "lm_output-3.txt")
################################### TABLE S3 ###################################
files_frequency <- c(
  "data.splotopen.local_abundance.txt",
  "data.splotopen.geographic_occupancy.txt",
  "data.splotopen.global_abundance.txt",
  "data.gbif.geographic_occupancy.txt"
)
for(k in 1:length(files_frequency)){
  fn1 <- files_frequency[k]
  message("::", fn1)
  freq <- read.table(fn1, sep = "\t", row.names = 1)
  colnames(freq) <- c("frequency")
  freq$frequency <- log(freq$frequency + 1)
  if(grepl("splotopen", fn1)){
    pr1 <- pr1.splot
    data1 <- pr.data.filtered1.splot
    pr2 <- pr2.splot
    data2 <- pr.data.filtered2.splot
  } else if(grepl("gbif", fn1)){
    pr1 <- pr1.gbif
    data1 <- pr.data.filtered1.gbif
    pr2 <- pr2.gbif
    data2 <- pr.data.filtered2.gbif
  }
  pts1 <- predict(pr1, data1)
  pts2 <- predict(pr2, data2)
  
  data1.full <- cbind(data1, pts1)
  data2.full <- cbind(data2, pts2)
  data3.full <- merge(data1.full, data2.full, by = 0)
  data3.full <- merge(data3.full, freq, by.x = 1, by.y = 0)
  
  aPC1 = data3.full$PC1.y
  aPC2 = data3.full$PC2.y
  bPC1 = data3.full$PC1.x
  bPC2 = data3.full$PC2.x
  
  rvalues <- c()
  l <- lm(data3.full$frequency ~ aPC1 + aPC2 + bPC1 + bPC2)
  print(summary(l))
  rvalues <- c(rvalues, summary(l)[[8]])
  l <- lm(data3.full$frequency ~ aPC2 + bPC1 + bPC2)
  #print(summary(l))
  rvalues <- c(rvalues, summary(l)[[8]])
  l <- lm(data3.full$frequency ~ aPC1 + bPC1 + bPC2)
  #print(summary(l))
  rvalues <- c(rvalues, summary(l)[[8]])
  l <- lm(data3.full$frequency ~ aPC1 + aPC2 + bPC2)
  #print(summary(l))
  rvalues <- c(rvalues, summary(l)[[8]])
  l <- lm(data3.full$frequency ~ aPC1 + aPC2 + bPC1)
  #print(summary(l))
  rvalues <- c(rvalues, summary(l)[[8]])
  message("Full R^2:", round(rvalues[1], 4))
  names <- c("Full", "aPC1", "aPC2", "bPC1", "bPC2")
  for(i in 2:5){
    message(names[i], " R^2:", round(rvalues[1] - rvalues[i], 4))
  }
}
sink(file = NULL)

################################### TABLE S6 ###################################
files_frequency <- c(
  "data.splotopen.geographic_occupancy.txt",
  "data.splotopen.geographic_occupancy.txt",
  "data.gbif.geographic_occupancy.txt",
  "data.gbif.geographic_occupancy.txt"
)
files_traits <- c(
  "aboveground_traits.txt.gz",
  "belowground_traits.txt.gz",
  "aboveground_traits.txt.gz",
  "belowground_traits.txt.gz"
)
titles <- c("a)", "b)", "c)", "d)")

megatree <- read.tree("megatree.tre")

sink(file = "glm_output.txt")
for(k in 1:4){
  fn1 <- files_frequency[k]
  fn2 <- files_traits[k]
  
  cat("::", titles[k], "\n")
  message("::", titles[k], " calculating phylogeny corrected glm models ", fn1, " and ", fn2)
  freq <- read.table(fn1, sep = "\t", row.names = 1)
  colnames(freq) <- c("frequency")
  if(grepl("aboveground", fn2)){
    trait <- read.table(fn2, sep = " ", header = T, row.names = 1)
    rownames(trait) <- gsub("_", " ", rownames(trait))
  } else {
    trait <- read.table(fn2, sep = "\t", header = T, row.names = 1)
  }
  
  combined <- merge(freq, trait, by = 0)
  rownames(combined) <- gsub(" ", "_", combined$Row.names)
  combined <- combined[,-1]
  # log transform the frequency and traits
  combined <- log(combined + 1)
  
  spp.available <- intersect(rownames(combined), megatree$tip.label)
  
  combined <- combined[spp.available, ]
  
  p.values <- c()
  for(predictor in names(trait)){
    tmp <- combined[!is.na(combined[,predictor]), ]
    spp <- rownames(tmp)
    bm <- corPagel(1, megatree, form = ~spp);
    l <- gls(as.formula(paste("frequency ~", predictor)), data = tmp, correlation = bm)
    p.values <- c(p.values, summary(l)$tTable[8])
    print(summary(l))
  }
  message("p-values:", paste(round(p.values, 3), collapse = ","))
  message("Corrected p-values (fdr):", paste(round(p.adjust(p.values, "fdr"), 3), collapse = ","))
}
sink(file = NULL)

########################### TABLE S3 colSD and col% ############################
files_frequency <- c(
  "data.splotopen.geographic_occupancy.txt",
  "data.gbif.geographic_occupancy.txt"
)
fr.data <- read.table("FR_species.csv", sep = ",", header = T)
sink(file = "glm_output-2.txt")
for(k in 1:2){
  fn1 <- files_frequency[k]
  cat("::", fn1, "\n")
  freq <- read.table(fn1, sep = "\t", row.names = 1)
  colnames(freq) <- c("frequency")
  freq$frequency <- log(freq$frequency + 1)
  combined <- merge(freq, fr.data, by.x = 0, by.y = 2)
  rownames(combined) <- gsub(" ", "_", combined$Row.names)
  combined <- combined[,-1]
  
  spp.available <- intersect(rownames(combined), megatree$tip.label)
  
  combined <- combined[spp.available, ]
  
  for(predictor in c("sd_col", "mean_col")){
    tmp <- combined[!is.na(combined[,predictor]), ]
    spp <- rownames(tmp)
    bm <- corPagel(1, megatree, form = ~spp)
    l <- gls(as.formula(paste("frequency ~", predictor)), data = tmp, correlation = bm)
    print(summary(l))
  }
}
sink(file = NULL)

################################## FIGURE S1 ###################################
library(raster)
library(sf)
library(maps)
library(viridis)

pdf("FigS1.pdf", width=16, height=12)
par(mfcol = c(2, 2))

map("world", border = "grey", fill = TRUE, col = "grey", bg = "white", ylim = c(-60, 90), mar = c(0,0,0,0))
dat = read.table("population.splot.txt.gz", sep="\t", row.names = 1, header = TRUE)
points(dat$longitude, dat$latitude, pch = 16, col = "blue", cex = 0.7)
text(x = par("usr")[1], y = par("usr")[4] - 1.5, labels = "a)", pos = 4, xpd = NA, cex = 2)

map("world", border = "grey", fill = TRUE, col = "grey", bg = "white", ylim = c(-60, 90), mar = c(0,0,0,0))
dat = read.table("population.gbif.txt.gz", sep="\t", row.names = 1, header = TRUE)
points(dat$longitude, dat$latitude, pch = 16, col = "blue", cex = 0.7)
text(x = par("usr")[1], y = par("usr")[4] - 1.5, labels = "b)", pos = 4, xpd = NA, cex = 2)

par(mar = c(7.1, 7.1, 7.1, 6.1))
# Download HumanImpact data: https://wcshumanfootprint.org/v2/
data.humanimpact <- raster("HumanImpact/HFP2009.tif") # this needs to be downloaded
hi.values = na.omit(getValues(data.humanimpact))
plot(density(hi.values[hi.values >= 0.0], adjust = 3), xlab = "Human Footprint Index", main = "", xaxs = "i", yaxs = "i", xlim=c(0, 50), lwd = 2)
data.coords = read.table("population.gbif.txt.gz", sep="\t", row.names = 1, header = TRUE)
lines(density(na.omit(data.coords$hi[data.coords$hi > 0.0])), col = "red", lwd = 2)
data.coords = read.table("population.splot.txt.gz", sep="\t", row.names = 1, header = TRUE)
lines(density(na.omit(data.coords$hi[data.coords$hi > 0.0])), col = "blue", lwd = 2)
legend("topright", legend = c("globally", "GBIF", "sPlotOpen"), col = c("black", "red", "blue"), pch = 15)
text(x = par("usr")[1], y = par("usr")[4] - 1.5, labels = "c)", pos = 4, xpd = NA, cex = 2)

dat <- read.table("population_density_1degree.csv.gz", sep = ",")
dat <- dat[dat < 99999]
plot(density(unlist(log(dat)), adjust = 3), xlab = "Human population per square metre (log)", main = "", xaxs = "i", yaxs = "i", xlim=c(0, 10), lwd = 2)
dat = read.table("population.gbif.txt.gz", sep="\t", row.names = 1, header = TRUE)
lines(density(log(dat$population)), col = "red", lwd = 2)
dat = read.table("population.splot.txt.gz", sep="\t", row.names = 1, header = TRUE)
lines(density(log(dat$population)), col = "blue", lwd = 2)
legend("topright", legend = c("globally", "GBIF", "sPlotOpen"), col = c("black", "red", "blue"), pch = 15)
text(x = par("usr")[1], y = par("usr")[4] - 1.5, labels = "d)", pos = 4, xpd = NA, cex = 2)
dev.off()

################################### TABLE S1 ###################################
dat1 <- read.table("data.splotopen.geographic_occupancy.txt", sep = "\t", row.names = 1)
dat2 <- read.table("data.gbif.geographic_occupancy.txt", sep = "\t", row.names = 1)

splot.species <- rownames(dat1)
gbif.species <- rownames(dat2)
all.species <- unique(c(splot.species, gbif.species))
shared.species <- intersect(splot.species, gbif.species)

pr.data <- read.table("aboveground_traits.txt.gz", sep=" ", header = T, row.names = 1)
rownames(pr.data) <- gsub("_", " ", rownames(pr.data))
above.species <- rownames(pr.data)
nrow(pr.data)
pr.data <- log(pr.data + 1)
pr.col <- read.table("FR_species.csv", sep=",", header = T, row.names = 2)
pr.col <- pr.col[,which(colnames(pr.col) %in% c("sd_col", "mean_col"))]
colnames(pr.col) <- c("colSD", "col%")
pr.combined <- merge(pr.data, pr.col, by = 0, all = T)
rownames(pr.combined) <- pr.combined$Row.names
pr.combined <- pr.combined[,-1]
pr.combined <- pr.combined[intersect(rownames(pr.combined), all.species), ]
colSums(!is.na(pr.combined[splot.species,]))
colSums(!is.na(pr.combined[gbif.species,]))
colSums(!is.na(pr.combined[shared.species,]))
nrow(pr.data[intersect(rownames(pr.data), splot.species),])
nrow(pr.data[intersect(rownames(pr.data), gbif.species),])
nrow(pr.data[intersect(rownames(pr.data), shared.species),])
nrow(na.omit(pr.data[intersect(rownames(pr.data), splot.species),]))
nrow(na.omit(pr.data[intersect(rownames(pr.data), gbif.species),]))
nrow(na.omit(pr.data[intersect(rownames(pr.data), shared.species),]))

pr.data <- read.table("belowground_traits.txt.gz", sep="\t", header = T, row.names = 1)
below.species <- rownames(pr.data)
length(intersect(above.species, below.species))
nrow(pr.data)
pr.col <- read.table("FR_species.csv", sep=",", header = T, row.names = 2)
pr.col <- pr.col[,which(colnames(pr.col) %in% c("sd_col", "mean_col"))]
colnames(pr.col) <- c("colSD", "col%")
pr.combined <- merge(pr.data, pr.col, by = 0, all = T)
rownames(pr.combined) <- pr.combined$Row.names
pr.combined <- pr.combined[,-1]
pr.combined <- pr.combined[intersect(rownames(pr.combined), all.species), ]
colSums(!is.na(pr.combined[splot.species,]))
colSums(!is.na(pr.combined[gbif.species,]))
colSums(!is.na(pr.combined[shared.species,]))
nrow(pr.data[intersect(rownames(pr.data), splot.species),])
nrow(pr.data[intersect(rownames(pr.data), gbif.species),])
nrow(pr.data[intersect(rownames(pr.data), shared.species),])
nrow(na.omit(pr.data[intersect(rownames(pr.data), splot.species),]))
nrow(na.omit(pr.data[intersect(rownames(pr.data), gbif.species),]))
nrow(na.omit(pr.data[intersect(rownames(pr.data), shared.species),]))