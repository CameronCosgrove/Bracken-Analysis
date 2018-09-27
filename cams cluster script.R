# Loading libraries
library(recluster)
library(phytools)
library(maps)
library(stats)
library(cluster)

#Patch cluster####

#load data

patch_df <- read.csv("patch_cluster.csv", sep = ",", header = TRUE)
patch_df
patch_sub <- unique(patch_df$patch)  # Making a vector with sites in our dataset
species_sub <- unique(patch_df$species)  # Making a vector with species in our dataset

p_a_matrix_patch <- matrix(0, length(patch_sub), length(species_sub))
for (i in 1:nrow(p_a_matrix_patch)){
  temp_sites_3 <- patch_df[which( patch_df$patch == patch_sub[i]),]
  p_a_matrix_patch[i, which(species_sub%in%temp_sites_3$species)] <- 1
  print(i)
}
# Now let's name our rows and columns with the codes for the sites and the codes for the species.
rownames(p_a_matrix_patch) <- as.character(patch_df$patch[match(patch_sub, patch_df$patch)])
colnames(p_a_matrix_patch) <- as.character(patch_df$species[match(species_sub, patch_df$species)])
dim(p_a_matrix_patch)

# Check if the loop function worked alright and did its job
p_a_matrix_patch[1:10,1:10]

#Now I am going to cluster the data using UPGMA unweighted 
patch_upgma_sorensen <- recluster.cons(p_a_matrix_patch, dist = "sorensen", tr=100, p=0.5, method = "average")
cams_patch_sorensen_upgma_cons <- patch_upgma_sorensen$cons
plot(cams_patch_sorensen_upgma_cons, direction = "downwards", cex=0.5)
write.tree(cams_patch_sorensen_upgma_cons, "cams_patch_sorensen_upgma_cons.tre")


#quadrat cluster####
#import data

quad_data <- read.csv("cluster_quad.csv", sep=",", head=TRUE)
head(quad_data)  # View the first few columns

dim(quad_data)  # How many rows and columns are there?


quad_sub <- unique(quad_data$quadID)  # Making a vector with sites in our dataset
spp_sub <- unique(quad_data$sppID)  # Making a vector with species in our dataset


# First we'll create an empty matrix with our sites in the rows and our species in the columns. The loop function will place a `1` on a given cell when the species is present in an area and will fill out the remaining cells with a `0`.

p_a_matrix <- matrix(0, length(quad_sub), length(spp_sub))
for (i in 1:nrow(p_a_matrix)){
  temp_sites <- quad_data[which(quad_data$quadID == quad_sub[i]),]
  p_a_matrix[i, which(spp_sub%in%temp_sites$sppID)] <- 1
  print(i)
}

# Now let's name our rows and columns with the codes for the sites and the codes for the species.
rownames(p_a_matrix) <- as.character(quad_data$quadID[match(quad_sub, quad_data$quadID)])
colnames(p_a_matrix) <- as.character(quad_data$sppID[match(spp_sub, quad_data$sppID)])
dim(p_a_matrix)

# Check if the loop function worked alright and did its job
p_a_matrix[1:10,1:10]

#Remove unique values cause they just add noise
p_a_matrix_trim <- p_a_matrix[,which(!colSums(p_a_matrix) == 1)] 
dim(p_a_matrix_trim)

#single linkage####
cams_cluster <- recluster.cons(p_a_matrix, tr = 1, p = 0.5, dist = "simpson", method = "single")
cams_cluster_tmp <-cams_cluster$cons  # Selecting the consensus tree (we'll discuss it later)
plot(cams_cluster_tmp, direction = "downwards", cex = 0.5)  # You can change the "direction" argument to your liking.

write.tree(cams_cluster_tmp, "cams_cluster_trim_tmp.tre")

#complete linkage#### 
bol_completelink <- recluster.cons(p_a_matrix_trim, tr = 100, p = 0.5, method = "complete")
bol_completelink_cons <- bol_completelink$cons
plot(bol_completelink_cons, direction = "downwards", cex = 0.5)

#upgma####

# Full species presence/abscence matrix

quad_upgma <- recluster.cons(p_a_matrix, dist = "simpson", tr=100, p=0.5, method = "average")
cams_upgma_cons <- quad_upgma$cons
write.tree(cams_upgma_cons, "cams_upgma_cons.tre")
plot(cams_upgma_cons, direction = "downwards", cex=0.5)

# trimmed upgma

quad_upgma_trim <- recluster.cons(p_a_matrix_trim, tr=100, p=0.5, method = "average")
cams_upgma_cons_trim <- quad_upgma_trim$cons
write.tree(cams_upgma_cons_trim, "cams_upgma_cons_trim.tre")
plot(cams_upgma_cons_trim, direction = "downwards", cex=0.5)

recluster.dist()


#sac

spa_total <- specaccum(p_a_matrix)
plot(spa_total) #plots the species accumulation curve and the confidence intervals for sites.
plot(spa_total, ci.type="poly", col="#af8dc3", lwd=2, ci.lty=0, ci.col="#7fbf7b") #males a prettier plot

