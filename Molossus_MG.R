#Loading the packages
library(geomorph)
library(Morpho)
library(ggplot2)
library(phytools)
library(RColorBrewer)
install.packages("ape")
library(ape)
install.packages("dendextend")
library(dendextend)
install.packages("tidytree")
library(treeio)
install.packages("ggplotify")
library(ggtree)

#Change working directory
#Select the path to the file
setwd("")

#Read landmarks file
Especies_Mol_I<-readland.tps("Especies_Molossus_inferior.tps", readcurves = F)
Especies_Mol_I.sobrepo<-procSym(Especies_Mol_I) #Morpho package function
Especies_Mol_I.sobrepo2<-gpagen(Especies_Mol_I) #geomorph

plotAllSpecimens(Especies_Mol_I)
plot(Especies_Mol_I.sobrepo2)
dim(Especies_Mol_I)

sizeEI<-Especies_Mol_I.sobrepo$size #see scale

#Categorize the groups
gruposICVA<-as.factor(c(rep("Mcur", 2), rep("Mmil", 3), rep("Mpar", 5),rep("Mmel", 3), 
                        rep("M.flu", 6), rep("M.azt", 12), rep("M.coi", 3), rep("M.pre", 7), 
                        rep("M.mol", 174), rep("M.ruf", 84)))

PCA.gls <- gm.prcomp(Especies_Mol_I.sobrepo$rotated)
scoresIR<-Especies_Mol_I.sobrepo$PCscores[,1:12]

#Separating PC1 and PC2 values ​​for all specimens
PC1e2<-Especies_Mol_I.sobrepo$PCscores[,1:2]

#Establishing colors for species
cores_especies <- c(brewer.pal(n = 9, name = "Set1"), "black") 
names(cores_especies) <- c("Mcur", "Mmil", "Mpar", "Mmel",
                           "M.flu", "M.azt", "M.coi", "M.pre", "M.mol", "M.ruf")

# Define different shapes for the species
simbolos_especies <- c(
  Mcur = 0,  # Triangle
  Mmil = 1,  # Circle
  Mpar = 2,  # Square
  Mmel = 3,  # Diamond
  M.flu = 4,  # Inverted triangle
  M.azt = 5,  # Plus
  M.coi = 6,  # Cross
  M.pre = 7,  # Star
  M.mol = 8,  # Asterisk
  M.ruf = 9   # Right triangle
)

# Plot the PCA graph
# Generate the plot
ggplot(PC1e2, aes(x = PC1, y = PC2, color = gruposICVA, shape = gruposICVA)) +
  geom_point(aes(color = gruposICVA, shape = gruposICVA), size = 5, stroke = 1.5, na.rm = TRUE) +
  scale_shape_manual(values = simbolos_especies) + 
  scale_color_manual(values = cores_especies) +
  labs(x = "PCA1", y = "PCA2") +
  theme_minimal() +
  theme(
    legend.key.size = unit(1.5, "cm"),  # Increase the size of symbols in the legend
    legend.text = element_text(size = 14),  # Larger text in the legend
    axis.text = element_text(size = 12),  # Larger text for axes
    axis.title = element_text(size = 16)  # Larger title for axes
  )

# Create a data frame with the data
data <- data.frame(
  Groups = gruposICVA,
  CentroidSize = sizeEI
)

# Generate violin boxplot
violin_plot <- ggplot(data, aes(x = Groups, y = CentroidSize, fill = Groups)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  scale_fill_manual(values = cores_especies) +
  theme_minimal() +
  labs(
    x = "Species Groups",
    y = "Centroid Size"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

print(violin_plot)

# CVA analysis between groups with 1000 permutation tests and cross-validation            
cvall <- CVA(scoresIR, gruposICVA, rounds = 1000, cv = TRUE) 
typprobs <- typprobClass(cvall$CVscores, groups = gruposICVA)
print(typprobs)

# Plot CVA
# CVA with ggplot
cva_scores <- as.data.frame(cvall$CVscores)  # Convert matrix to data frame
cva_scores$species <- gruposICVA  # Add groups to the data frame

# Create the plot with CVA data
ggplot(cva_scores, aes(x = `CV 1`, y = `CV 2`)) +
  geom_point(aes(color = species, shape = species), size = 5, stroke = 1.5, na.rm = TRUE) +
  scale_shape_manual(values = simbolos_especies) +
  scale_color_manual(values = cores_especies) +
  coord_fixed(ratio = 1) +
  labs(x = "CV 1", y = "CV 2", color = "Species", shape = "Species") +
  theme_minimal() +
  theme(
    legend.key.size = unit(1.5, "cm"),  # Increase the size of symbols in the legend
    legend.text = element_text(size = 14),  # Larger text in the legend
    axis.text = element_text(size = 12),  # Larger text for axes
    axis.title = element_text(size = 16)  # Larger title for axes
  )

# ANOVA testing for dimorphism
gdf_inferior <- geomorph.data.frame(coords = scoresIR, groups = gruposICVA, size = sizeEI)

fit1_inferior <- procD.lm(coords ~ size, 
                          data = gdf_inferior, iter = 9999, 
                          RRPP = FALSE, print.progress = FALSE)  # Randomize raw values
summary(fit1_inferior)

# Testing centroid significance
fit2_inferior <- procD.lm(size ~ groups, 
                          data = gdf_inferior, iter = 9999, 
                          RRPP = FALSE, print.progress = FALSE)
summary(fit2_inferior)

# Extracting residuals
resI <- fit1_inferior$residuals
cvallCR <- CVA(resI, gruposICVA, rounds = 1000, cv = TRUE)  # With residuals

# CVA with residuals
cva_scoresCR <- as.data.frame(cvallCR$CVscores)  # Convert matrix to data frame
cva_scoresCR$species <- gruposICVA  # Add groups to the data frame

# Create the plot with CVA data with residuals
ggplot(cva_scoresCR, aes(x = `CV 1`, y = `CV 2`)) + 
  geom_point(aes(color = species, shape = species), size = 5, stroke = 1.5) + 
  scale_shape_manual(values = simbolos_especies) + 
  scale_color_manual(values = cores_especies) +
  coord_fixed(ratio = 1) + 
  labs(x = "CV 1", y = "CV 2", color = "Species", shape = "Species") + 
  theme_minimal() +
  theme(
    legend.key.size = unit(1.5, "cm"),  # Increase the size of symbols in the legend
    legend.text = element_text(size = 14),  # Larger text in the legend
    axis.text = element_text(size = 12),  # Larger text for axes
    axis.title = element_text(size = 16)  # Larger title for axes
  )

# Integration analysis of skull views
two.b.pls(scoresIR, scoresSR, iter = 999, seed = NULL, print.progress = TRUE)

# Regression scores
M <- Especies_Mol_I.sobrepo2$consensus
plethAllometry <- procD.lm(coords ~ groups, data = gdf_inferior)

# Create reference data frame for species
species_info <- data.frame(
  group = names(cores_especies),  # Species names
  color = unname(cores_especies),   # Associated colors
  shape = unname(simbolos_especies)  # Associated shapes
)

# Map colors and shapes based on gruposICVA
plot_colors <- species_info$color[match(gruposICVA, species_info$group)]
plot_shapes <- species_info$shape[match(gruposICVA, species_info$group)]

# Verify if the mapping is correct
head(data.frame(gruposICVA, plot_colors, plot_shapes))

allom.plot <- plot(
  plethAllometry,
  type = "regression",
  predictor = log(gdf_inferior$size),
  reg.type = "RegScore",
  pch = plot_shapes,  # Shapes
  col = plot_colors, 
  lwd = 2             # Line thickness
)

# Shape analysis for species groups. Smaller vs. larger species  
PC <- PCA.gls$x[, 1]
preds <- shape.predictor(Especies_Mol_I.sobrepo2$coords, x = cvall$CVscores[, 1], Intercept = FALSE, 
                         pred1 = min(cvall$CVscores[, 1]), pred2 = max(cvall$CVscores[, 1]))  # PC1 extremes

plotRefToTarget(M, preds$pred1, mag = 1)
plotRefToTarget(M, preds$pred2, mag = 1)

# Plot the landmarks for pred1
plotRefToTarget(M, preds$pred1, mag=1)
segments(preds$pred1[1, 1], preds$pred1[1, 2], preds$pred1[2, 1], preds$pred1[2, 2], lwd = 5, col = "blue")  # 1 to 2  # 1 to 2
segments(preds$pred1[2, 1], preds$pred1[2, 2], preds$pred1[7, 1], preds$pred1[7, 2], lwd = 5, col = "blue")  # 2 to 7
segments(preds$pred1[7, 1], preds$pred1[7, 2], preds$pred1[11, 1], preds$pred1[11, 2], lwd = 5, col = "blue")  # 7 to 11
segments(preds$pred1[11, 1], preds$pred1[11, 2], preds$pred1[17, 1], preds$pred1[17, 2], lwd = 5, col = "blue")  # 11 to 17
segments(preds$pred1[17, 1], preds$pred1[17, 2], preds$pred1[15, 1], preds$pred1[15, 2], lwd = 5, col = "blue")  # 17 to 15
segments(preds$pred1[15, 1], preds$pred1[15, 2], preds$pred1[16, 1], preds$pred1[16, 2], lwd = 5, col = "blue")  # 15 to 16
segments(preds$pred1[16, 1], preds$pred1[16, 2], preds$pred1[14, 1], preds$pred1[14, 2], lwd = 5, col = "blue")  # 16 to 14
segments(preds$pred1[14, 1], preds$pred1[14, 2], preds$pred1[18, 1], preds$pred1[18, 2], lwd = 5, col = "blue")  # 14 to 18
segments(preds$pred1[18, 1], preds$pred1[18, 2], preds$pred1[12, 1], preds$pred1[12, 2], lwd = 5, col = "blue")  # 18 to 12
segments(preds$pred1[12, 1], preds$pred1[12, 2], preds$pred1[8, 1], preds$pred1[8, 2], lwd = 5, col = "blue")  # 12 to 8
segments(preds$pred1[8, 1], preds$pred1[8, 2], preds$pred1[3, 1], preds$pred1[3, 2], lwd = 5, col = "blue")  # 8 to 3
segments(preds$pred1[3, 1], preds$pred1[3, 2], preds$pred1[1, 1], preds$pred1[1, 2], lwd = 5, col = "blue")  # 3 to 1
segments(preds$pred1[15, 1], preds$pred1[15, 2], preds$pred1[13, 1], preds$pred1[13, 2], lwd = 5, col = "blue")  # 15 to 13
segments(preds$pred1[13, 1], preds$pred1[13, 2], preds$pred1[14, 1], preds$pred1[14, 2], lwd = 5, col = "blue")  # 13 to 14
segments(preds$pred1[4, 1], preds$pred1[4, 2], preds$pred1[5, 1], preds$pred1[5, 2], lwd = 5, col = "blue")  # 4 to 5
segments(preds$pred1[5, 1], preds$pred1[5, 2], preds$pred1[6, 1], preds$pred1[6, 2], lwd = 5, col = "blue")  # 5 to 6
segments(preds$pred1[4, 1], preds$pred1[4, 2], preds$pred1[9, 1], preds$pred1[9, 2], lwd = 5, col = "blue")  # 4 to 9
segments(preds$pred1[6, 1], preds$pred1[6, 2], preds$pred1[10, 1], preds$pred1[10, 2], lwd = 5, col = "blue")  # 6 to 10

# Plot the landmarks for preds2
plotRefToTarget(M, preds$pred2, mag=1)
segments(preds$pred2[1, 1], preds$pred2[1, 2], preds$pred2[2, 1], preds$pred2[2, 2], lwd = 5, col = "red")  # 1 to 2  # 1 to 2
segments(preds$pred2[2, 1], preds$pred2[2, 2], preds$pred2[7, 1], preds$pred2[7, 2], lwd = 5, col = "red")  # 2 to 7
segments(preds$pred2[7, 1], preds$pred2[7, 2], preds$pred2[11, 1], preds$pred2[11, 2], lwd = 5, col = "red")  # 7 to 11
segments(preds$pred2[11, 1], preds$pred2[11, 2], preds$pred2[17, 1], preds$pred2[17, 2], lwd = 5, col = "red")  # 11 to 17
segments(preds$pred2[17, 1], preds$pred2[17, 2], preds$pred2[15, 1], preds$pred2[15, 2], lwd = 5, col = "red")  # 17 to 15
segments(preds$pred2[15, 1], preds$pred2[15, 2], preds$pred2[16, 1], preds$pred2[16, 2], lwd = 5, col = "red")  # 15 to 16
segments(preds$pred2[16, 1], preds$pred2[16, 2], preds$pred2[14, 1], preds$pred2[14, 2], lwd = 5, col = "red")  # 16 to 14
segments(preds$pred2[14, 1], preds$pred2[14, 2], preds$pred2[18, 1], preds$pred2[18, 2], lwd = 5, col = "red")  # 14 to 18
segments(preds$pred2[18, 1], preds$pred2[18, 2], preds$pred2[12, 1], preds$pred2[12, 2], lwd = 5, col = "red")  # 18 to 12
segments(preds$pred2[12, 1], preds$pred2[12, 2], preds$pred2[8, 1], preds$pred2[8, 2], lwd = 5, col = "red")  # 12 to 8
segments(preds$pred2[8, 1], preds$pred2[8, 2], preds$pred2[3, 1], preds$pred2[3, 2], lwd = 5, col = "red")  # 8 to 3
segments(preds$pred2[3, 1], preds$pred2[3, 2], preds$pred2[1, 1], preds$pred2[1, 2], lwd = 5, col = "red")  # 3 to 1
segments(preds$pred2[15, 1], preds$pred2[15, 2], preds$pred2[13, 1], preds$pred2[13, 2], lwd = 5, col = "red")  # 15 to 13
segments(preds$pred2[13, 1], preds$pred2[13, 2], preds$pred2[14, 1], preds$pred2[14, 2], lwd = 5, col = "red")  # 13 to 14
segments(preds$pred2[4, 1], preds$pred2[4, 2], preds$pred2[5, 1], preds$pred2[5, 2], lwd = 5, col = "red")  # 4 to 5
segments(preds$pred2[5, 1], preds$pred2[5, 2], preds$pred2[6, 1], preds$pred2[6, 2], lwd = 5, col = "red")  # 5 to 6
segments(preds$pred2[4, 1], preds$pred2[4, 2], preds$pred2[9, 1], preds$pred2[9, 2], lwd = 5, col = "red")  # 4 to 9
segments(preds$pred2[6, 1], preds$pred2[6, 2], preds$pred2[10, 1], preds$pred2[10, 2], lwd = 5, col = "red")  # 6 to 10

# Check Mahalanobis distances
# Plot Mahalanobis distances for Cvall
dendroS = hclust(cvall$Dist$GroupdistMaha)
dendroS$labels = levels(gruposICVA)
par(mar = c(4, 4.5, 1, 1))
dendroS = as.dendrogram(dendroS)
plot(dendroS, main = '', sub = '', xlab = "Species",
     ylab = 'Mahalanobis distance')

# Plot Mahalanobis distances for residuals in Cvall
dendroS = hclust(cvallCR$Dist$GroupdistMaha)
dendroS$labels = levels(gruposICVA)
par(mar = c(4, 4.5, 1, 1))
dendroS = as.dendrogram(dendroS)
plot(dendroS, main = '', sub = '', xlab = "Species",
     ylab = 'Mahalanobis distance')

# install.packages("dendextend") to compare previous dendrograms
# library(dendextend)
# Generate the first dendrogram with size effect
dendroS1 <- hclust(cvall$Dist$GroupdistMaha)
dendroS1$labels <- levels(gruposICVA)
dendroS1 <- as.dendrogram(dendroS1)

# Generate the second dendrogram without size effect
dendroS2 <- hclust(cvallCR$Dist$GroupdistMaha)
dendroS2$labels <- levels(gruposICVA)
dendroS2 <- as.dendrogram(dendroS2)

# Compare the two dendrograms
dend_list <- dendlist(dendroS1, dendroS2)

# Plot the two dendrograms side by side with connecting lines
tanglegram(dend_list, 
           highlight_distinct_edges = FALSE, 
           common_subtrees_color_lines = TRUE, 
           common_subtrees_color_branches = TRUE, 
           lab.cex = 1.2, edge.lwd = 2)

# Read the phylogenetic tree
# Install and load the treeio and ggplotify packages, if not already installed
# Load the tree
# Add the path to the file on your PC
tree <- read.beast("")

# Check the tree structure
str(tree)

# Check node ages
tree@phylo$node.label
tip <- c("MH185125_MH410726_Malvarezi", "EF080483_MH058051_Mfentoni", "MG191810_KM387368_Mbondae",
         "MH185179_MH058094_Msinaloae", "MH185186_MH058091_P.centralis", "KX355065_MH058057_Mverrilli", "JF454657_MH058046_E.auripendulus")

# Remove taxa not sampled in the morphological dataset
batstotal.tre <- drop.tip(tree, tip)

ggtree(batstotal.tre) +
  geom_tiplab(align = TRUE, linetype = 'dashed', linesize = .3)

# Boxplot for a target trait
df <- data.frame(x = gruposICVA, y = sizeEI)

# Calculate centroid means
sizecur <- mean(df[1:2, 2])
sizemil <- mean(df[3:5, 2])
sizepar <- mean(df[6:10, 2])
sizemel <- mean(df[11:13, 2])
sizeflu <- mean(df[14:19, 2])
sizeazt <- mean(df[20:31, 2])
sizecoi <- mean(df[32:34, 2])
sizepre <- mean(df[35:41, 2])
sizemol <- mean(df[42:215, 2])
sizeruf <- mean(df[216:299, 2])

# Calculate centroid means for each species and associate the value with each species
mediascent <- c(sizemil, sizepar, sizeazt,  sizemel, sizecoi, sizeflu, sizeruf, sizecur, sizepre, sizemol)
especiesmedia <- c("MH185146_MH058056_Mmilleri_Cuba", 
                   "129_Mparanaensis", 
                   "MH185133_MH058047_Maztecus_Mexico", 
                   "BDN1_Mmelini_BrasilPR",
                   "CA28_MH058078_Mcoibensis",  
                   "AMA114_Mfluminensis", 
                   "CESC64_Mrufus", 
                   "MH185139_MH058050_Mcurrentium_Paraguay",  
                   "MH185168_MH058079_Mpretiosus_Nicaragua", 
                   "CUMA73_Mmolossus_BrazilMA")

names(mediascent) <- especiesmedia
# Convert the 'treedata' tree to 'phylo'
bats.tree <- as.phylo(batstotal.tre)
# Check if the tree is of class phylo
class(bats.tree)

# Calculate ancestral states
ancestral_states <- fastAnc(bats.tree, mediascent, vars = FALSE, CI = FALSE)

# Simple visualization of ancestral states with phytools
plot(bats.tree)
nodelabels(ancestral_states, frame = "n", cex = 0.8)

# Plot contMap
obj <- contMap(bats.tree, mediascent)
ln.size <- log(setNames(mediascent, rownames(especiesmedia)))

names(ln.size) <- especiesmedia
mammal.contMap <- contMap(bats.tree, ln.size, plot = FALSE, res = 200)

# Change color scheme
mammal.contMap <- setMap(mammal.contMap,
                         c("white", "#FFFFB2", "#FECC5C", "#FD8D3C",
                           "#E31A1C"))
plot(mammal.contMap, fsize = c(0.7, 0.8),
     leg.txt = "log(centroid size)")
par(mar = c(5.1, 4.1, 4.1, 2.1))

# Testing phylogenetic signal
ps <- physignal(mediascent, bats.tree, iter = 9999)
ps

# Estimate effect size
ps_z <- physignal.z(mediascent, bats.tree)
ps_z

# Read the file with PC data
PC <- read.table("Especies-PC1-PC2.txt", header = TRUE)

# Check if all taxon names in PC are present in tree$tip.label
if (!all(PC$species %in% bats.tree$tip.label)) {
  stop("There are inconsistencies between the taxon names in the PC data and the phylogenetic tree.")
}

# Set column 1 as row names of the data frame
rownames(PC) <- PC$Especies

# Reorder PC data to match the order of tree taxon labels
PC2 <- PC[match(bats.tree$tip.label, PC$Especies), ]

# Phylomorphospace
phylomorphospace(bats.tree, PC[, 2:3], colors = setNames(c("blue", "red"), c(0, 1)))

# Comparing net rates of evolution among traits on phylogenies
# Creating group
gps <- c("129_Mparanaensis",
         "AMA114_Mfluminensis",
         "BDN1_Mmelini_BrasilPR",
         "CA28_MH058078_Mcoibensis",
         "CESC64_Mrufus",
         "CUMA73_Mmolossus_BrazilMA",
         "MH185133_MH058047_Maztecus_Mexico",
         "MH185139_MH058050_Mcurrentium_Paraguay",
         "MH185146_MH058056_Mmilleri_Cuba",
         "MH185168_MH058079_Mpretiosus_Nicaragua"
)

# Check if the tree is of class phylo
class(bats.tree)
bats.tree <- as.phylo(bats.tree) 
bats.tree$tip.label # Validate tree tip names
names(PC2) # Confirm vector names

PCs <- as.matrix(PC2[, 2:3])  # Select only PC1 and PC2 and transform into a matrix
str(PCs)  # Check if PCs is now a matrix
rownames(PCs) <- PC2$Especies  # Assign species names as row names
rownames(PCs)  # Confirm row names match species
all(bats.tree$tip.label %in% rownames(PCs))  # Check if all tree tips are in PCs
all(bats.tree$tip.label %in% gps) 
all(gps %in% PCs) # Compare 'mediascent' names with 'gps' to check if they match
str(gps)       # Check structure of gps
str(PCs)       # Check structure of PCs
all(gps %in% rownames(PCs))  # Compare gps with row names (species) of PCs
gp <- factor(gps)  # Ensure 'gps' is correctly converted to a factor

gps <- factor(gps)  # Convert 'gps' to a factor
# Reorganizing 'gps' according to rownames of 'PCs'
gps <- gps[match(rownames(PCs), gps)]
# Confirm species in gps are aligned with rownames(PCs)
all(gps == rownames(PCs))  # This should return TRUE if names are correctly aligned
rownames(PCs) <- gps
# Check if gps names are in rownames of PCs
setdiff(gps, rownames(PCs))  # This should return NULL or an empty vector
result <- compare.multi.evol.rates(PCs[rownames(PCs)], gps, bats.tree$tip.label, iter = 999)


PCs_filtrado <- PCs[rownames(PCs) %in% bats.tree$tip.label, , drop = FALSE]
gps_filtrado <- gps[rownames(PCs_filtrado)]
all(rownames(PCs_filtrado) %in% bats.tree$tip.label)  # Should return TRUE
result <- compare.multi.evol.rates(PCs[row.names(PCs)], gp = gps_filtrado, phy = bats.tree, iter = 999)
