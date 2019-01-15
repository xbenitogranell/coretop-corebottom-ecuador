#Set WD
getwd("~/myfolder")

##loading libraries for functions used
library(vegan) #to perform nonmetric multidimensional analysis, simper and adonis2 
library(ade4) #to perfom PCA analysis
library(ggplot2) #to make some nice ordination plots
library(goeveg) #allow to select species for vegan ordination objects
library(tidyverse) #allow to manipulate tabulate data
library(spacemakeR) #allow to compute distance-based Moran eigenvectors. Installed using install.packages("spacemakeR", repos="http://R-Forge.R-project.org")
library(sp) #dependencies for spacemakeR
library(spdep) #dependencies for spacemakeR
library(adespatial) #dependencies for spacemakeR
library(scales) #allow to scale variables from 0 to 1 


######################################
#### Read and manipulate datasets ####
######################################

#Diatom datasets
diatoms <- read.csv("diatom_data.csv", row.names = 1)
# [,1:3] - lake, region and time points
# [,4:ncol(diatoms)] - diatom species relative abundances. Join between modern and fossil diatom datasets was performed prior analysis. 

site.time <- diatoms[,1:3]
diat <- diatoms[,4:ncol(diatoms)]

  #subset species data
  abund <- apply(diat, 2, max)
  n.occur <- apply(diat>0, 2, sum)
  diat <- diat[, abund>2 & n.occur>2] #relative abundances higher than 2% and present in more than 2 samples
  
  #Hellinger transformation
  diat <- decostand(diat, method = "hellinger")
  
  #Combine 
  spp <- cbind(diat, site.time)


#Environmental and climatic data
env_data <- read.csv("env_data.csv", row.names = 1)

#transform variables to meet assumptions of homogenity of variances
env_transformed <- transform(env_data, Elevation=sqrt(Elevation),Cond=log10(Cond+0.25), Ca=log10(Ca+0.25), Mg=log10(Mg+0.25), K=log10(K+0.25), TP=log10(TP+0.25), TN=log10(TN+0.25), NO3=log10(NO3+0.25), SO4=log10(SO4+0.25), 
                             MaxDepth=log10(MaxDepth+0.25), lake.area=log10(lake.area+0.25), MAT=log10(MAT+0.25), P.season=log10(P.season+0.25), MAP=log10(MAP+0.25), T.season=log10(T.season+0.25))


#panel correlation plots to assess data distribution
panel.hist <- function(x, ...) {     
  usr <- par("usr"); on.exit(par(usr))     
  par(usr = c(usr[1:2], 0, 1.5) )     
  h <- hist(x, plot = FALSE)     
  breaks <- h$breaks; nB <- length(breaks)     
  y <- h$counts; y <- y/max(y)     
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...) 
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {     
  usr <- par("usr"); on.exit(par(usr))     
  par(usr = c(0, 1, 0, 1))     
  r <- abs(cor(x, y, use = "complete"))   
  txt <- format(c(r, 0.123456789), digits = digits)[1]     
  txt <- paste0(prefix, txt)     
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)     
  text(0.5, 0.5, txt, cex = cex.cor * r) }


#Figure S1 of the manuscript
pairs(env_subset, diag.panel = panel.hist, upper.panel = panel.smooth, lower.panel = panel.cor, gap = 0, cex.labels = 1, cex=1.5, font.labels = 1) 

#Select variables with r Pearson <0.85
env_subset<- env_transformed[, which(names(env_data) %in% c("pH", "Cond", "Secchi", "Ca", "Mg", "K", "NO3", "SO4", "TN", "TP",
                                                            "MaxDepth", "Elevation", "lake.area",
                                                            "MAT", "P.season", "MAP", "T.season"))]

  
######################################
#### MORAN'S Eigenvector Maps  #######
######################################

#### Code  from Dray et al. 2012 #####

#Extract Latitude and Longitude
xy.coord <- data.frame(env_data$Lat, env_data$Long)

## create the Gabriel graph
xy.coord <- as.matrix(xy.coord)

nb1 <- graph2nb(gabrielneigh(xy.coord), sym = T)

lw1 <- nb2listw(nb1)

# create the MEMs
U <- scores.listw(lw1)

# test the significance of Moran's I for MEMs
mctests <- test.scores(U, lw1, 999)
colnames(U$vectors) <- paste("MEM", 1:ncol(U$vectors), sep="")
U.sel <- U
U.sel$vectors <- U.sel$vectors[, mctests[,2]<0.05] ## keep only MEMS with signif values

#Extract MEMS
mems <- U.sel$vectors


######################################
#### Principal Component Analysis ####
######################################
## Code adapted from Carles Alcaraz (carles.alcaraz@irta.cat)

#Perform PCA with NIPALS algorithm
PCA.data <- data.frame(scale(env_subset)) #scale variables
PCA.nipals <- nipals(PCA.data, nf = 2, rec = FALSE, niter = 1000, tol = 1e-09)

#Save summary matrix results  
PCA.summary <- NULL
PCA.summary <- matrix(data = NA, nrow = 6, ncol = 3, byrow = FALSE, dimnames = NULL)
colnames(PCA.summary) <- c("Parameter", "Value", "Explained variance")

PCA.summary[1, 1] <- c("Number of extracted factors")
PCA.summary[1, 2] <- PCA.nipals$nf

PCA.summary[2, 1] <- c("Number of variables")
PCA.summary[2, 2] <- nrow(PCA.nipals$co)

PCA.summary[3, 1] <- c("Eigenvalue 1")
PCA.summary[3, 2] <- PCA.nipals$eig[1]
PCA.summary[3, 3] <- (PCA.nipals$eig[1] / nrow(PCA.nipals$co)) * 100

PCA.summary[4, 1] <- c("Eigenvalue 2")
PCA.summary[4, 2] <- PCA.nipals$eig[2]
PCA.summary[4, 3] <- (PCA.nipals$eig[2] / nrow(PCA.nipals$co)) * 100

PCA.summary[5, 1] <- c("Number of iterations axis 1")
PCA.summary[5, 2] <- PCA.nipals$nb[1]

PCA.summary[6, 1] <- c("Number of iterations axis 2")
PCA.summary[6, 2] <- PCA.nipals$nb[2]

#Save the column coordinates (Component matrix) #ncol= n variables 
Component.matrix <- NULL
Component.matrix <- matrix(data = NA, nrow = ncol(PCA.data), ncol = 2, byrow = FALSE, dimnames = NULL)
Component.matrix[,1] <- rownames(PCA.nipals$co)
Component.matrix[,2] <- (((PCA.nipals$co[,1] ^ 2) * PCA.nipals$eig[1]) / PCA.nipals$eig[1]) + (((PCA.nipals$co[,2] ^ 2) * PCA.nipals$eig[2]) / PCA.nipals$eig[2])
Component.matrix <- cbind(Component.matrix[,1], PCA.nipals$co, Component.matrix[,2])
colnames(Component.matrix) <- c("Variable", "Component 1", "Component 2", "Power Importance")

#Save the column normed scores (Component scores coefficient matrix)
Component.coefficient.matrix <- NULL
Component.coefficient.matrix <- matrix(data = NA, nrow = ncol(PCA.data), ncol = 1, byrow = FALSE, dimnames = NULL)
Component.coefficient.matrix[,1] <- rownames(PCA.nipals$c1)
Component.coefficient.matrix <- cbind(Component.coefficient.matrix, PCA.nipals$c1)
colnames(Component.coefficient.matrix) <- c("Variable", "Component 1", "Component 2")

#Save the row coordinates (Factor Scores)
colnames(PCA.nipals$li) <- c("Component 1", "Component 2")
Factor.scores <- data.frame(cbind(PCA.data, PCA.nipals$li))


#Create data frame with site scores and regions
region <- site.time[1:21,] #lake regions of the modern diatom dataset
PCA.scores <- data.frame(component1=Factor.scores$Component.1, component2=Factor.scores$Component.2, lakes=site.time$region)

#extract factor scores (environmental variables)
comp1 <- as.numeric(Component.coefficient.matrix[,2])
comp2 <- as.numeric(Component.coefficient.matrix[,3])

 
#######################################
# Non-Metric Multidimensional Scaling #
#######################################

#Perform NMDS on diatom top-core and downcore datasets
diat.nmds <- metaMDS(diat)


#Perform Environmental Fitting procedure
diatTop <- subset(spp, time=="core top")
diatBottom <- subset(spp, time=="downcore")

diatTop <- diatTop[,!(names(diatTop) %in% c("site", "region", "time"))]
diatTop<- decostand(diatTop, method = "hellinger")

diatBottom <- diatBottom[,!(names(diatBottom) %in% c("site", "region", "time"))]
diatBottom<- decostand(diatBottom, method = "hellinger")

#NMDS with core-top diatom data
diat.nmds.top <- metaMDS(diatTop)

#Extract scores from diatom downcore NMDS ordination
diat.nmds.bottom <- metaMDS(diatBottom)
diat.nmds.scores.fossils <- as.data.frame(scores(diat.nmds.bottom, display = "sites"))
colnames(diat.nmds.scores.fossils) <- c("NMDS1.hist", "NMDS2.hist")

#combine environmental, MEMs and historical variables
explanatory <- cbind(env_subset, mems)
fit <- envfit(diat.nmds.top, explanatory, na.rm=TRUE, permutations = 999) 
fit


#Select the 20% most abundant species with 70% best environmental fit in NDMS for axes 1 & 2
selected_nmds <- ordiselect(diatTop, diat.nmds.top,  ablim = 0.2, fitlim = 0.7, method = "axes", env = fit)  


## Multivariate homogenity of group dispersions
mod.spp <- adonis2(diat ~ spp$time, method = "bray")
mod.spp


#######################################################
# Redundancy Analysis (RDA) and variance partitioning #
#######################################################

#merge datasets
df <- cbind(diatTop, mems, env_subset, diat.nmds.scores.fossils)

#remove rows with missing values
row.has.na <- apply(df, 1, function(x){any(is.na(x))})
sum(row.has.na)
df <- df[!row.has.na,]

#separate predictor matrices
spp <- df2[, names(df2) %in% names(diatTop)]
env <- df2[, names(df2) %in% names(env_subset)]
spatial <- df2[, names(df2) %in% names(data.frame(mems))]
hist <- df2[, names(df2) %in% names(diat.nmds.scores.fossils)]


#rda environmental
rda1 <- rda(spp, env)
rda1.R2a <- RsquareAdj(rda1)$adj.r.squared
rda1.fw.env <- forward.sel(spp, env, adjR2thresh = rda1.R2a)

rda1.fw.env$R2a.fullmodel <- rda1.R2a
write.csv(rda1.fw.env, "rda.env.fw.csv")

#rda spatial
rda1 <- rda(spp, spatial)
rda1.R2a <- RsquareAdj(rda1)$adj.r.squared
rda1.fw.spatial <- forward.sel(spp, spatial, adjR2thresh = rda1.R2a)

rda1.fw.spatial$R2a.fullmodel <- rda1.R2a
write.csv(rda1.fw.spatial, "rda.spatial.fw.csv")

#rda historical
rda1 <- rda(spp, hist)
rda1.R2a <- RsquareAdj(rda1)$adj.r.squared
rda1.fw.historical <- forward.sel(spp, hist, adjR2thresh = rda1.R2a)

rda1.fw.historical$R2a.fullmodel <- rda1.R2a
write.csv(rda1.fw.historical, "rda.hist.fw.csv")


#Extract subset of forward selected variables
#environmental
env.sel <- env[, which(names(env) %in% rda1.fw.env$variables)] 

#spatial
spatial.sel <- spatial[, which(names(spatial) %in% rda1.fw.spatial$variables)] 

#historical
hist.sel <- hist[, which(names(hist) %in% rda1.fw.historical$variables)] 


#Varpart
varpart <- varpart(spp, env.sel, spatial.sel, hist.sel)

plot(varpart, Xnames="")
text(locator(), c("Environmental","Spatial"), cex=1)


#Test pure effects
anova.cca(rda(spp, env.sel, cbind(spatial.sel, hist.sel)), perm.max = 999) ## test pure environmental (signif 0.05)
anova.cca(rda(spp, spatial.sel, cbind(env.sel, hist.sel)), perm.max = 999) ## test pure spatial (signif 0.005)  
anova.cca(rda(spp, hist.sel, cbind(env.sel, spatial.sel)), perm.max = 999) ## test pure spatial (signif 0.005)  


############################
### Plot grid PCA + NMDS ### 
############################
#Figure 2 of the manuscript#

par(mfrow=c(2,2))
par(mar=c(2,2,1,1), mgp=c(1.2,.5,0))

#Plot PCA site labels (=lakes)
  plot(PCA.scores$component1, PCA.scores$component2, type = "n", xlab = "PCA1", ylab = "PCA2")
  abline(h=0, col="grey")
  abline(v=0, col="grey")
  
  points(PCA.scores[(PCA.scores$lakes=="Andes"), 1:2], col="black", pch=21, bg="black")
  points(PCA.scores[(PCA.scores$lakes=="Inter Andean"), 1:2], col="black", pch=24, bg="black")

  lakelbls <- c("YAH", "YBO", "SPA", "CUN", "CUI", "LLA", "PIN", "COL", "KUY", "DCH", "HUA", "CHI", "CAR", "CUB", "EST", "YAN", "MAR", "JIG", "RIN", "FON", "PIC",
                    "YAH", "YBO", "SPA", "CUN", "CUI", "LLA", "PIN", "COL", "KUY", "DCH", "HUA", "CHI", "CAR", "CUB", "EST", "YAN", "MAR", "JIG", "RIN", "FON", "PIC")
  
  text(PCA.scores[,1:2], labels = lakelbls, pos = 2, cex = 0.9, offset = 0.3)
  
#Plot environmental variables    
  
  plot(comp1, comp2, pch=1, col="black", xlab = "PCA1", ylab = "PCA2")
  abline(h=0, col="grey")
  abline(v=0, col="grey")
  
  #Labels
  labels <- as.character(Component.coefficient.matrix[,1])
  text(comp1, comp2, labels = labels, pos = 1, cex = 0.9, offset = 0.3)
  
#Plot NMDS
  scrs <- scores(diat.nmds, display = "sites", choices = 1:2)
  plot(diat.nmds, type = "n", xlim=range(scrs[,1]), ylim=range(scrs[,2]))
  
  points(diat.nmds, display="sites", pch=21, col="#E69F00", bg="#E69F00", select=spp$time =="core top" & spp$region=="Andes")
  points(diat.nmds, display="sites", pch=24, col="#E69F00", bg="#E69F00", select=spp$time =="core top" & spp$region=="Inter Andean")
  
  points(diat.nmds, display="sites", pch=21, col="grey", bg="grey", select=spp$time=="downcore" & spp$region=="Andes")
  points(diat.nmds, display="sites", pch=24, col="grey", bg="grey", select=spp$time=="downcore" & spp$region=="Inter Andean")
  
  with(spp, ordiellipse(diat.nmds, spp$time, kind = "se", conf = 0.95, cex=1, col = c("#E69F00", "grey")))
 
  text(diat.nmds, labels = lakelbls, pos = 3, cex = 0.8, offset = 0.3)
  
  
  legend("bottomright",legend=c("core-top", "downcore"), pch=c(20), 
         col=c("#E69F00", "grey"), pt.bg =c("#E69F00", "grey" ), cex=0.7)
  
  legend("bottomleft", legend=c("High-elevation Andes", "Inter Andean plateau"), pch=c(21,24), cex=0.7)
  
  
#Plot species labels
  plot(diat.nmds, type = "n", xlim=range(scrs[,1]), ylim=range(scrs[,2]))
  
  #Plot envfit procedure
  plot(fit, cex=0.8, p.max=0.04) #plot only statistically significant variables
  
  #Plot selected species with best environmental fitting
  points(diat.nmds.top, display="species", select = selected_nmds, pch=3, col="red", cex=0.7)
  ordipointlabel(diat.nmds.top, display="species", select = selected_nmds, col="black", cex=0.6, add = TRUE)
  

#######################
### SIMPER analysis ###
#######################

sim <- simper(diat, site.time$time, permutations = 999)
summary(sim) #species identified are marked with * in Figure 3 using Inkscape


###################################################################################################
# Stacked species plots to visualize the amount of change between historical and modern time points
###################################################################################################

# Extract diatom species with best fit on the NMDS ordination to plot
data <- cbind(site.time, spp[, names(spp) %in% selected_nmds], depth=env_subset$MaxDepth, elevation=env_subset$Elevation)

# Make long format and calculate scaled relative abundances
spp_thin <- data %>% 
  gather(key = taxa, value = count, -time, -region, -site, -elevation, -depth) %>%
  mutate(time = as.character(time)) %>%
  mutate(region = as.character(region)) %>%
  arrange(depth) %>%
  group_by(region) %>%
  mutate(scaled_count = rescale(count)) %>%
  ungroup()


#Make the plot (Figure 3 of the manuscript)
plot <- ggplot(spp_thin, aes(time, scaled_count, fill = time)) +
  geom_col() +
  theme_bw() +
  facet_grid(region~taxa, space = "free") +
  scale_y_continuous(limits = c(0,1)) +
  labs(x = "", y = "Scaled relative abundance") +
  # customize the appearance
  theme(
    axis.text.x = element_blank(),
    # rotate the facet labels
    strip.text.x = element_text(angle = 40, hjust = 0, vjust = 0), 
    strip.text.y = element_text(angle = 0), #label text of lake water depth
    # turn off the label background
    strip.background = element_blank(),
    axis.ticks.x=element_blank()
  ) +
  scale_fill_manual(values=c("#999999", "#E69F00"), 
                    name="",
                    breaks=c("core top", "downcore"),
                    labels=c("core-top", "downcore"))

plot


#Make the plot (Figures S2 and S3 of the manuscript)
  #Make long format and arrange lakes by water depth (replace arrange() with elevation for make the Figure S3)
  spp_thin <- data %>% 
    gather(key = taxa, value = count, -time, -region, -site, -elevation, -depth) %>%
    mutate(time = as.character(time)) %>%
    mutate(region = as.character(region)) %>%
    arrange(depth) %>%
    mutate(depth_f = factor(depth)) 

  plot <- ggplot(spp_thin, aes(time, count, fill = time)) +
    geom_col() +
    theme_bw() +
    coord_flip() +
    facet_grid(depth_f~taxa) +
    #facet_grid(elevation_f~taxa) +
    labs(x = "", y = "Relative abundance (%)") +
    #scale_x_discrete(breaks = c(0,20)) +
    # customize the appearance
    theme(
      axis.text.y = element_blank(),
      axis.text.x = element_text(size  = 10,
                   angle = 45,
                   hjust = 1,
                   vjust = 1),
      # rotate the facet labels
      strip.text.x = element_text(angle = 40, hjust = 0, vjust = 0), 
      strip.text.y = element_text(angle = 0), #label text of lake water depth
      # turn off the label background
      strip.background = element_blank(),
      axis.ticks.x=element_blank()
    ) +
    scale_fill_manual(values=c("#999999", "#E69F00"), 
                      name="",
                      breaks=c("core top", "downcore"),
                      labels=c("core-top", "downcore"))
  
  plot

# voodoo that makes it so that facet labels can overlap
# https://stackoverflow.com/questions/49740215/ggplot-facet-grid-label-cut-off
plot_grob <- ggplotGrob(plot)
for(i in which(grepl("strip-t", plot_grob$layout$name))){
  plot_grob$grobs[[i]]$layout$clip <- "off"
}

# needed to draw the modified plot_grob
grid::grid.draw(plot_grob)

