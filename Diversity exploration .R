library(vegan)



#species area curve for bracken quadrats####



#Load data #
bracken_data <- read.csv("bracken_quad_data.csv", sep=",", head=TRUE)

bracken_quad_sub <- unique(bracken_data$quad_code)  # Making a vector with sites in our dataset
bracken_spp_sub <- unique(bracken_data$spp)  # Making a vector with species in our dataset

p_a_matrix_bracken <- matrix(0, length(bracken_quad_sub), length(bracken_spp_sub))
for (i in 1:nrow(p_a_matrix_bracken)){
  temp_sites3 <- bracken_data[which(bracken_data$quad_code == bracken_quad_sub[i]),]
  p_a_matrix_bracken[i, which(bracken_spp_sub%in%temp_sites3$spp)] <- 1
  print(i)
}

# Now let's name our rows and columns with the codes for the sites and the codes for the species.
rownames(p_a_matrix_bracken) <- as.character(bracken_data$quad_code[match(bracken_quad_sub, bracken_data$quad_code)])
colnames(p_a_matrix_bracken) <- as.character(bracken_data$spp[match(bracken_spp_sub, bracken_data$spp)])
dim(p_a_matrix_bracken)


spa_bracken <- specaccum(p_a_matrix_bracken)
plot(spa_bracken) #plots the species accumulation curve and the confidence intervals for sites.
plot(spa_bracken, ci.type="poly", col="#7fbf7b", lwd=2, ci.lty=0, ci.col="#af8dc3") #males a prettier plot



#species area curve for grassland quadrats####
#load data

grass_data <- read.csv("grass_quad_data.csv", sep=",", head=TRUE)


grass_quad_sub <- unique(grass_data$quad_ID)  # Making a vector with sites in our dataset
grass_spp_sub <- unique(grass_data$spp_ID)  # Making a vector with species in our dataset

#making the pa matrix

p_a_matrix_grass <- matrix(0, length(grass_quad_sub), length(grass_spp_sub))
for (i in 1:nrow(p_a_matrix_grass)){
  temp_sites4 <- grass_data[which(grass_data$quad_ID == grass_quad_sub[i]),]
  p_a_matrix_grass[i, which(grass_spp_sub%in%temp_sites4$spp_ID)] <- 1
  print(i)
}

# Now let's name our rows and columns with the codes for the sites and the codes for the species.
rownames(p_a_matrix_grass) <- as.character(grass_data$quad_ID[match(grass_quad_sub, grass_data$quad_ID)])
colnames(p_a_matrix_grass) <- as.character(grass_data$spp_ID[match(grass_spp_sub, grass_data$spp_ID)])
dim(p_a_matrix_grass)

#sac time! 
spa_grass <- specaccum(p_a_matrix_grass)
plot(spa_grass)  #plots the species accumulation curve and the confidence intervals for sites.
plot(spa_grass, ci.type="poly", col="#af8dc3", lwd=2, ci.lty=0, ci.col="#7fbf7b") #males a prettier plot




#speceis indivdual curve for bracken ####

#just to explore the differnece, this will not be used in may analysis

spi_bracken<-specaccum(p_a_matrix_bracken, method="rarefaction")
plot(spi_bracken)
plot(spi_bracken, ci.type="poly", col="red", lwd=2, ci.lty=0, ci.col="seagreen2") #males a prettier plot



#species area curve for all quadrats ####
#sac

spa_total1 <- specaccum(p_a_matrix)
plot(spa_total1) #plots the species accumulation curve and the confidence intervals for sites.
plot(spa_total1, ci.type="poly", col="skyblue3", lwd=2, ci.lty=0, ci.col="lightgreen", xlab = "Area", ylab = "Species found ")



#makes a prettier plot

#similarity comparison ####
p_a_matrix
sim_comp <- vegdist(p_a_matrix, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE)
sim_comp # can't do because I need a community matrix

#this I think assumes that only one individual of a speceis was found a each quadrat
#Actually I think I workd of p/a data so we good

species_pool_bracken <- specpool(p_a_matrix_bracken) #chao = 48.6 +/- 6.54
specpool_pool_grass <- specpool(p_a_matrix_grass) #chao = 78.5 +/- 36.9   could be an outcome of missidentifying grassed atrifically raising the number of unique spps, or a real quality of the grass spp pool.
species_pool_total <- specpool(p_a_matrix)

rarefy(p_a_matrix, se = TRUE, MARGIN = 1)



H<-diversity(p_a_matrix, "shannon", base = exp(1)) # can't do because I need a community matric with number of individuals, and I only have presence absecense. 

#jaccards simmilartiy ####
library(jaccard)

#load data
jaccard_df <- read.csv("jaccard.csv", sep=",", header = TRUE)

#Each vector is the same lenght and the order of species is the same in  both so they can be properly compared as two like binary sets.
bracken_jaccard <- as.vector(jaccard_df$bracken)  # Making a vector with binary speceis info in
grass_jaccard <- as.vector(jaccard_df$grass)  # Making a vector with binary species in
bracken_jaccard

jaccard(jaccard_df[[1]], jaccard_df[[2]])

jaccard(bracken_jaccard, grass_jaccard)
jaccard.test.bootstrap(jaccard_df$bracken, jaccard_df$grass, B=500)

vector

#Sorensens similarity #####
install.packages("OmicsMarkeR")

#plot both habitat intersection. Can tell yu about dispersal charateristics. Grassland has a larger speceis pool and plants with larger dispersal charateristics 
plot1 <- plot(spa_grass, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightgreen") 
  
plot2 <- plot(spa_bracken, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")

plot3 <- plot1 + plot2

spa_grass$habitat <- "grass"
spa_bracken$habitat <- "bracken"
spas <- rbind(spa_grass, spa_bracken)
plot(spas, ci.type="poly", col=habitat,lwd=2, ci.lty=0, ci.col=habitat)

#general beta tests
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("qvalue")
install.packages("jaccard", dependencies = TRUE, source = TRUE) # need to update R
library(jaccard)
