setwd("~/Dropbox/Coral")

require(ape)
require(vegan)
require(picante)
require(ade4)

# Load list of coral-phylotype associaitons
Associations <- read.table("All_Associations.txt",header=TRUE,sep="\t")

# Convert associaitons into matrix form
Coral_Matrix <- as.data.frame.matrix(table(Associations[,2:3]))

# Load phylogenetic tree of all phylotypes
Tree <- read.tree("Tree_AllPhylotypes.tre")

# Note: There are 4 phylotypes that are in the dataset but not represented in the tree.
# Each of these phylotypes has only one observation. This code removes those 4 phylotypes. 
Coral_Matrix <- Coral_Matrix[,colnames(Coral_Matrix) %in% intersect(colnames(Coral_Matrix),Tree$tip.label)]

#RUN UNIFRAC CALCUATIONS
nreps <- 10                    ### Change number of bootstrap reps
rarefaction_depth <- 60         ### Change subsample depth 

Coral_Matrix_Keep_For_UniFrac <- Coral_Matrix[rowSums(Coral_Matrix) >= rarefaction_depth,]

Coral_Matrices <- list()
for (i in 1:nreps) {
  Coral_Matrices[[i]] <- rrarefy(Coral_Matrix_Keep_For_UniFrac,rarefaction_depth)
}

Coral_Distances <- list()
for (i in 1:nreps) {
  Pruned_Tree <- prune.sample(Coral_Matrices[[i]],Tree)
  Coral_Distances[[i]] <- unifrac(Coral_Matrices[[i]],Pruned_Tree)
  print(i) 
}

Sums = Reduce('+', Coral_Distances)
Average_Distance = Sums/length(Coral_Distances)
write.table(as.matrix(Average_Distance),"Average_UniFrac_Distance_R60_ALL_PHYLOTYPES.txt",sep="\t",col.names=NA)

#RUN UNIFRAC CALCUATIONS - FULL TREE
nreps <- 10                    ### Change number of bootstrap reps
rarefaction_depth <- 10         ### Change subsample depth 
Pruned_Tree <- prune.sample(Coral_Matrix,Tree)

Coral_Matrix_Keep_For_UniFrac <- Coral_Matrix[rowSums(Coral_Matrix) >= rarefaction_depth,]


Coral_Matrices <- list()
for (i in 1:nreps) {
  Coral_Matrices[[i]] <- rrarefy(Coral_Matrix_Keep_For_UniFrac,rarefaction_depth)
}

Coral_Distances <- list()
for (i in 1:nreps) {
  Coral_Distances[[i]] <- unifrac(Coral_Matrices[[i]],Pruned_Tree)
  print(i) 
}

Sums = Reduce('+', Coral_Distances)
Average_Distance = Sums/length(Coral_Distances)
write.table(as.matrix(Average_Distance),"Average_FT_UniFrac_Distance_R10_ALL_PHYLOTYPES.txt",sep="\t",col.names=NA)

#RUN BRAY-CURTIS CALCULATIONS
nreps <- 1000
rarefaction_depth <- 60

Coral_Matrix_Keep <- Coral_Matrix[rowSums(Coral_Matrix) >= rarefaction_depth,]

Coral_Matrices <- list()
for (i in 1:nreps) {
  Coral_Matrices[[i]] <- rrarefy(Coral_Matrix_Keep,rarefaction_depth)
}

Coral_Distances <- list()
for (i in 1:nreps) {
  Coral_Distances[[i]] <- vegdist(Coral_Matrices[[i]])
}

Sums = Reduce('+', Coral_Distances)
Average_Distance = Sums/length(Coral_Distances)
write.table(as.matrix(Average_Distance),"Average_BC_Distance_R60_ALL_PHYLOTYPES.txt",sep="\t",col.names=NA)

#Matrix of BRI Difs
BRI <- read.table("Coral_BRI.txt",header=TRUE,sep="\t")

rownames(BRI) <- BRI[,1]  
BRIs_vec <- as.vector(as.numeric(as.character(BRI[,2])))  
BRI_Dif <- outer(BRIs_vec,BRIs_vec,"-")
colnames(BRI_Dif) <- rownames(BRI)
rownames(BRI_Dif) <- rownames(BRI)
BRI_Dif <- abs(BRI_Dif)

write.table(BRI_Dif,"BRI_Dif.txt",sep="\t",col.names=NA)

#Get Matrix of Phylogenetic Distances

TreeDistances <- cophenetic.phylo(Tree)
write.table(TreeDistances,"TreeDistances.txt",sep="\t",col.names=NA)


#Write Function to Prepare Distance Matricies for Analysis
PrepMatrix <- function(x) {
  
  Matrix <- x
  Matrix[lower.tri(Matrix)] <- NA
  MeltedMatrix <- melt(Matrix)
  colnames(MeltedMatrix) <- c("Species1","Species2","Value")
  MeltedMatrix$Species2 <- gsub("\\."," ",MeltedMatrix$Species2)
  MeltedMatrix <- MeltedMatrix[MeltedMatrix$Species1 != MeltedMatrix$Species2,]
  MeltedMatrixComplete <- MeltedMatrix[complete.cases(MeltedMatrix),]
  return(MeltedMatrixComplete)
  
}

#Load Beta Div Calcs
require(reshape2)

BC5 <-  PrepMatrix(as.matrix(read.table("Average_BC_Distance_R5_ALL_PHYLOTYPES.txt",header=TRUE,row.names=1,sep="\t")))
BC10 <- PrepMatrix(as.matrix(read.table("Average_BC_Distance_R10_ALL_PHYLOTYPES.txt",header=TRUE,row.names=1,sep="\t")))
BC15 <- PrepMatrix(as.matrix(read.table("Average_BC_Distance_R15_ALL_PHYLOTYPES.txt",header=TRUE,row.names=1,sep="\t")))
BC20 <- PrepMatrix(as.matrix(read.table("Average_BC_Distance_R20_ALL_PHYLOTYPES.txt",header=TRUE,row.names=1,sep="\t")))
BC30 <- PrepMatrix(as.matrix(read.table("Average_BC_Distance_R30_ALL_PHYLOTYPES.txt",header=TRUE,row.names=1,sep="\t")))
BC40 <- PrepMatrix(as.matrix(read.table("Average_BC_Distance_R40_ALL_PHYLOTYPES.txt",header=TRUE,row.names=1,sep="\t")))
BC50 <- PrepMatrix(as.matrix(read.table("Average_BC_Distance_R50_ALL_PHYLOTYPES.txt",header=TRUE,row.names=1,sep="\t")))
BC60 <- PrepMatrix(as.matrix(read.table("Average_BC_Distance_R60_ALL_PHYLOTYPES.txt",header=TRUE,row.names=1,sep="\t")))

U5 <-  PrepMatrix(as.matrix(read.table("Average_UniFrac_Distance_R5_ALL_PHYLOTYPES.txt",header=TRUE,row.names=1,sep="\t")))
U10 <- PrepMatrix(as.matrix(read.table("Average_UniFrac_Distance_R10_ALL_PHYLOTYPES.txt",header=TRUE,row.names=1,sep="\t")))
U15 <- PrepMatrix(as.matrix(read.table("Average_UniFrac_Distance_R15_ALL_PHYLOTYPES.txt",header=TRUE,row.names=1,sep="\t")))
U20 <- PrepMatrix(as.matrix(read.table("Average_UniFrac_Distance_R20_ALL_PHYLOTYPES.txt",header=TRUE,row.names=1,sep="\t")))
U30 <- PrepMatrix(as.matrix(read.table("Average_UniFrac_Distance_R30_ALL_PHYLOTYPES.txt",header=TRUE,row.names=1,sep="\t")))
U40 <- PrepMatrix(as.matrix(read.table("Average_UniFrac_Distance_R40_ALL_PHYLOTYPES.txt",header=TRUE,row.names=1,sep="\t")))
U50 <- PrepMatrix(as.matrix(read.table("Average_UniFrac_Distance_R50_ALL_PHYLOTYPES.txt",header=TRUE,row.names=1,sep="\t")))
U60 <- PrepMatrix(as.matrix(read.table("Average_UniFrac_Distance_R60_ALL_PHYLOTYPES.txt",header=TRUE,row.names=1,sep="\t")))

FTU5 <-  PrepMatrix(as.matrix(read.table("Average_FT_UniFrac_Distance_R5_ALL_PHYLOTYPES.txt",header=TRUE,row.names=1,sep="\t")))
FTU10 <- PrepMatrix(as.matrix(read.table("Average_FT_UniFrac_Distance_R10_ALL_PHYLOTYPES.txt",header=TRUE,row.names=1,sep="\t")))
FTU15 <- PrepMatrix(as.matrix(read.table("Average_FT_UniFrac_Distance_R15_ALL_PHYLOTYPES.txt",header=TRUE,row.names=1,sep="\t")))
FTU20 <- PrepMatrix(as.matrix(read.table("Average_FT_UniFrac_Distance_R20_ALL_PHYLOTYPES.txt",header=TRUE,row.names=1,sep="\t")))
FTU30 <- PrepMatrix(as.matrix(read.table("Average_FT_UniFrac_Distance_R30_ALL_PHYLOTYPES.txt",header=TRUE,row.names=1,sep="\t")))
FTU40 <- PrepMatrix(as.matrix(read.table("Average_FT_UniFrac_Distance_R40_ALL_PHYLOTYPES.txt",header=TRUE,row.names=1,sep="\t")))
FTU50 <- PrepMatrix(as.matrix(read.table("Average_FT_UniFrac_Distance_R50_ALL_PHYLOTYPES.txt",header=TRUE,row.names=1,sep="\t")))
FTU60 <- PrepMatrix(as.matrix(read.table("Average_FT_UniFrac_Distance_R60_ALL_PHYLOTYPES.txt",header=TRUE,row.names=1,sep="\t")))

#IMPORT AND FORMAT BRI DIF AND PHYLO DIF
BRI_Dif <- PrepMatrix(as.matrix(read.table("BRI_Dif.txt",header=TRUE,row.names=1,sep="\t")))
  Comparisons_In_BRI_Dif <- paste(BRI_Dif$Species1,"_",BRI_Dif$Species2)
  BRI_Dif_With_Names <- cbind(Comparisons_In_BRI_Dif,BRI_Dif)

CoralTree <- read.nexus("CoralShortPhylo320.nex")  
CoralTreeDistances <- cophenetic.phylo(CoralTree)
CoralTreeDist <- PrepMatrix(CoralTreeDistances)
  CoralTreeDist$Species1 <- gsub("_"," ",CoralTreeDist$Species1)
  CoralTreeDist$Species2 <- gsub("_"," ",CoralTreeDist$Species2)
  
  Comparisons_In_CoralTreeDist <- paste(CoralTreeDist$Species1,"_",CoralTreeDist$Species2)
  Coral_Tree_Dist_With_Names <- cbind(Comparisons_In_CoralTreeDist,CoralTreeDist)
  
  
#WRITE CORRELATION FUNCTION
Beta_Div_Cor <- function(x) {

Beta <- x
Comparisions_In_Beta <- paste(Beta$Species1,"_",Beta$Species2)
Beta_With_Names <- cbind(Comparisions_In_Beta,Beta)
  Beta_With_Names <- Beta_With_Names[order(Beta_With_Names[,1]),]
  
BRI_Keep <- BRI_Dif_With_Names[BRI_Dif_With_Names[,1] %in% Beta_With_Names[,1],]
  BRI_Keep <- BRI_Keep[order(BRI_Keep[,1]),]

Coral_Keep <- Coral_Tree_Dist_With_Names[Coral_Tree_Dist_With_Names[,1] %in% Beta_With_Names[,1],]   
  Coral_Keep <- Coral_Keep[order(Coral_Keep[,1]),]
  
FullData <- cbind(Beta_With_Names,BRI_Keep,Coral_Keep)[,c(1,2,3,4,8,12)]
colnames(FullData) <- c("Comparison","Species1","Species2","BetaDiv","BRI_Dif","Phylo_Dif")

return(cor.test(FullData$BetaDiv,FullData$BRI_Dif)) #NOTE: Change to BRU_Dif for corr with BRI

}

#WRITE CORRELATION FUNCTION FOR INDO-PACIFIC ONLY
Corals_by_Ocean <- read.table("Corals_by_Ocean.txt",header=TRUE,sep="\t")
IP_Corals <- Corals_by_Ocean[Corals_by_Ocean$Ocean == "Indo-Pacific",]

Beta_Div_Cor_IP <- function(x) {
  
  Beta <- x
  Comparisions_In_Beta <- paste(Beta$Species1,"_",Beta$Species2)
  Beta_With_Names <- cbind(Comparisions_In_Beta,Beta)
  Beta_With_Names <- Beta_With_Names[order(Beta_With_Names[,1]),]
  Beta_With_Names <- Beta_With_Names[Beta_With_Names$Species1 %in% IP_Corals$Coral_Species & Beta_With_Names$Species2 %in% IP_Corals$Coral_Species,]
  
  BRI_Keep <- BRI_Dif_With_Names[BRI_Dif_With_Names[,1] %in% Beta_With_Names[,1],]
  BRI_Keep <- BRI_Keep[order(BRI_Keep[,1]),]
  
  Coral_Keep <- Coral_Tree_Dist_With_Names[Coral_Tree_Dist_With_Names[,1] %in% Beta_With_Names[,1],]   
  Coral_Keep <- Coral_Keep[order(Coral_Keep[,1]),]
  
  FullData <- cbind(Beta_With_Names,BRI_Keep,Coral_Keep)[,c(1,2,3,4,8,12)]
  colnames(FullData) <- c("Comparison","Species1","Species2","BetaDiv","BRI_Dif","Phylo_Dif")
  
  return(cor.test(FullData$BetaDiv,FullData$Phylo_Dif)) #NOTE: Change to BRI_Dif for corr with BRI
  
}

#
Beta_Div_Cor_IP(FTU5)

#WRITE CORRELATION FUNCTION FOR HT/VT CORALS ONLY

Metadata <- read.table("Nodes_10PLUS_New.txt",header=TRUE,sep="\t")
Coral_Metadata <- Metadata[Metadata$Node_Type == "Coral",]
HT_Corals <- droplevels(Coral_Metadata[Coral_Metadata$Symbiodinium_in_propagules == "No",1])
VT_Corals <- droplevels(Coral_Metadata[Coral_Metadata$Symbiodinium_in_propagules == "Yes",1])

HT_Corals <- gsub("_"," ",HT_Corals)
VT_Corals <- gsub("_"," ",VT_Corals)

Beta_Div_Cor_HTVT <- function(Beta,name) {

  Comparisions_In_Beta <- paste(Beta$Species1,"_",Beta$Species2)
  Beta_With_Names <- cbind(Comparisions_In_Beta,Beta)
  Beta_With_Names <- Beta_With_Names[order(Beta_With_Names[,1]),]
  Beta_With_Names_HT <- Beta_With_Names[Beta_With_Names$Species1 %in% HT_Corals & Beta_With_Names$Species2 %in% HT_Corals,]
  Beta_With_Names_VT <- Beta_With_Names[Beta_With_Names$Species1 %in% VT_Corals & Beta_With_Names$Species2 %in% VT_Corals,]
  
  BRI_Keep_HT <- BRI_Dif_With_Names[BRI_Dif_With_Names[,1] %in% Beta_With_Names_HT[,1],]
  BRI_Keep_HT <- BRI_Keep_HT[order(BRI_Keep_HT[,1]),]
  
  BRI_Keep_VT <- BRI_Dif_With_Names[BRI_Dif_With_Names[,1] %in% Beta_With_Names_VT[,1],]
  BRI_Keep_VT <- BRI_Keep_VT[order(BRI_Keep_VT[,1]),]
  
  Coral_Keep_HT <- Coral_Tree_Dist_With_Names[Coral_Tree_Dist_With_Names[,1] %in% Beta_With_Names_HT[,1],]   
  Coral_Keep_HT <- Coral_Keep_HT[order(Coral_Keep_HT[,1]),]
  
  Coral_Keep_VT <- Coral_Tree_Dist_With_Names[Coral_Tree_Dist_With_Names[,1] %in% Beta_With_Names_VT[,1],]   
  Coral_Keep_VT <- Coral_Keep_VT[order(Coral_Keep_VT[,1]),]
  
  FullData_HT <- cbind(Beta_With_Names_HT,BRI_Keep_HT,Coral_Keep_HT)[,c(1,2,3,4,8,12)]
  colnames(FullData_HT) <- c("Comparison","Species1","Species2","BetaDiv","BRI_Dif","Phylo_Dif")
  
  FullData_VT <- cbind(Beta_With_Names_VT,BRI_Keep_VT,Coral_Keep_VT)[,c(1,2,3,4,8,12)]
  colnames(FullData_VT) <- c("Comparison","Species1","Species2","BetaDiv","BRI_Dif","Phylo_Dif")
  
  HT_Phylo <- cor.test(FullData_HT$BetaDiv,FullData_HT$Phylo_Dif)
  VT_Phylo <- cor.test(FullData_VT$BetaDiv,FullData_VT$Phylo_Dif)
  
  HT_BRI <- cor.test(FullData_HT$BetaDiv,FullData_HT$BRI_Dif)
  VT_BRI <- cor.test(FullData_VT$BetaDiv,FullData_VT$BRI_Dif)
  
  Corrs <- c(HT_Phylo$estimate,VT_Phylo$estimate,HT_BRI$estimate,VT_BRI$estimate)
  Pvals <- c(HT_Phylo$p.val,VT_Phylo$p.val,HT_BRI$p.val,VT_BRI$p.val)
  Names <- c("Dataset","HT_Phylo","VT_Phylo","HT_BRI","VT_BRI","HT_Phylo_P","VT_Phylo_P","HT_BRI_P","VT_BRI_P")
  
  Output <- c(name,Corrs,Pvals)
  names(Output) <- Names
  
  return(Output)
  
}

HTVT_Corrs <- rbind(Beta_Div_Cor_HTVT(BC10,"BC10"),
Beta_Div_Cor_HTVT(BC15,"BC15"),
Beta_Div_Cor_HTVT(BC20,"BC20"),
Beta_Div_Cor_HTVT(BC30,"BC30"),
Beta_Div_Cor_HTVT(BC40,"BC40"),
Beta_Div_Cor_HTVT(BC50,"BC50"),
Beta_Div_Cor_HTVT(BC60,"BC60"),
Beta_Div_Cor_HTVT(U10,"U10"),
Beta_Div_Cor_HTVT(U15,"U15"),
Beta_Div_Cor_HTVT(U20,"U20"),
Beta_Div_Cor_HTVT(U30,"U30"),
Beta_Div_Cor_HTVT(U40,"U40"),
Beta_Div_Cor_HTVT(U50,"U50"),
Beta_Div_Cor_HTVT(U60,"U60"))

write.table(HTVT_Corrs,"BetaDiv_HTVT_Correlations.txt",sep="\t")

#

require(gplots)
HeatmapData <- as.matrix(read.table("Beta_Heatmap_Input.txt",header=TRUE,row.names=1,sep="\t"))

palette <- colorRampPalette(c("white","#a6bddb","#3690c0","#0570b0"))(n = 51)
heatmap.2(HeatmapData,scale="none",dendrogram="none",Rowv=FALSE,Colv=FALSE,trace="none",cellnote=HeatmapData,col=palette)


#

HeatmapData <- as.matrix(read.table("BetaDiv_HTVT_Correlations.txt",header=TRUE,row.names=1,sep="\t"))[,1:8]
HeatmapData <- round(HeatmapData,3)
HeatmapData[HeatmapData < 0] <- 0

palette <- colorRampPalette(c("white","#a6bddb","#3690c0","#0570b0"))(n = 51)
heatmap.2(HeatmapData,scale="none",dendrogram="none",Rowv=FALSE,Colv=FALSE,trace="none",cellnote=HeatmapData,col=palette,notecol="black")


###


Beta <- BC40
Name <- "Bray-Curtis RD40"

DensityPlots <- function(Beta,Name) {

VT <- Beta[Beta$Species1 %in% VT_Corals & Beta$Species2 %in% VT_Corals,]
HT <- Beta[Beta$Species1 %in% HT_Corals & Beta$Species2 %in% HT_Corals,]
    Between1 <- Beta[Beta$Species1 %in% VT_Corals & Beta$Species2 %in% HT_Corals,]
    Between2 <- Beta[Beta$Species1 %in% HT_Corals & Beta$Species2 %in% VT_Corals,]
Between <- rbind(Between1,Between2)

VT_Values <- cbind(VT$Value,"VT vs VT")
HT_Values <- cbind(HT$Value,"HT vs HT")
Bet_Values <- cbind(Between$Value,"HT vs VT")
AllValues <- as.data.frame(rbind(VT_Values,HT_Values,Bet_Values))
  colnames(AllValues) <- c("BetaDiv","Comparison")
  AllValues[,1] <- as.numeric(as.character(AllValues[,1]))

return(ggplot(AllValues,aes(x=BetaDiv,col=Comparison,fill=Comparison)) + geom_density(alpha=0.4) + theme_bw() + labs(x=Name,y="Density") + scale_fill_manual(values=c("#FD6467","#899DA4","#0B775E")) + scale_color_manual(values=c("#FD6467","#899DA4","#0B775E")) + theme(legend.position="none") + scale_x_continuous(limits=c(0,1)))

}

DensityPlotsLegend <- function(Beta,Name) {
  
  VT <- Beta[Beta$Species1 %in% VT_Corals & Beta$Species2 %in% VT_Corals,]
  HT <- Beta[Beta$Species1 %in% HT_Corals & Beta$Species2 %in% HT_Corals,]
  Between1 <- Beta[Beta$Species1 %in% VT_Corals & Beta$Species2 %in% HT_Corals,]
  Between2 <- Beta[Beta$Species1 %in% HT_Corals & Beta$Species2 %in% VT_Corals,]
  Between <- rbind(Between1,Between2)
  
  VT_Values <- cbind(VT$Value,"VT vs VT")
  HT_Values <- cbind(HT$Value,"HT vs HT")
  Bet_Values <- cbind(Between$Value,"HT vs VT")
  AllValues <- as.data.frame(rbind(VT_Values,HT_Values,Bet_Values))
  colnames(AllValues) <- c("BetaDiv","Comparison")
  AllValues[,1] <- as.numeric(as.character(AllValues[,1]))
  
  return(ggplot(AllValues,aes(x=BetaDiv,col=Comparison,fill=Comparison)) + geom_density(alpha=0.4) + theme_bw() + labs(x=Name,y="Density") + scale_fill_manual(values=c("#FD6467","#899DA4","#0B775E")) + scale_color_manual(values=c("#FD6467","#899DA4","#0B775E")) + theme(legend.position=c(0.3,0.8)) + scale_x_continuous(limits=c(0,1)))
  
}


grid.arrange(
  DensityPlotsLegend(BC10,"Bray-Curtis RD10"),
  DensityPlots(BC15,"Bray-Curtis RD15"),
  DensityPlots(BC20,"Bray-Curtis RD20"),
  DensityPlots(U10,"UniFrac RD10"),
  DensityPlots(U15,"UniFrac RD15"),
  DensityPlots(U20,"UniFrac RD20"),
  ncol=3)



grid.arrange(
DensityPlotsLegend(BC10,"Bray-Curtis RD10"),
DensityPlots(BC20,"Bray-Curtis RD20"),
DensityPlots(BC30,"Bray-Curtis RD30"),
DensityPlots(BC40,"Bray-Curtis RD40"),
DensityPlots(BC50,"Bray-Curtis RD50"),
DensityPlots(BC60,"Bray-Curtis RD60"),
DensityPlots(U10,"UniFrac RD10"),
DensityPlots(U20,"UniFrac RD20"),
DensityPlots(U30,"UniFrac RD30"),
DensityPlots(U40,"UniFrac RD40"),
DensityPlots(U50,"UniFrac RD50"),
DensityPlots(U60,"UniFrac RD60"),
ncol=6)
