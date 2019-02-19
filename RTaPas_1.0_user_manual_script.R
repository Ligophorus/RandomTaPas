 #############################################################################
# Script to implement — Random Tanglegram Partitions (Random TaPas):          #
# An Alexandrian approach to the cophylogenetic Gordian knot                  #
# J.A. Balbuena, O.A. Pérez-Escobar, C. Llopis-Belenguer, I. Blasco-Costa     #
# Submitted                                                                   #
# For questions/feedback: j.a.balbuena@uv.es                                  #
# LICENSE: MIT (https://opensource.org/licenses/MIT)                          #
# YEAR: 2018                                                                  #
# COPYRIGHT HOLDER: Symbiosis Ecol. & EVol. Lab @ U. Valencia                 #
 #############################################################################
# The Random TaPas algorithm is demonstrated with data of mitochondrial       #
# haplotypes of the cytochrome oxidase subunit 1 of the trematode, Coitocaecum#
# parvum, and those of its amphipod host, Paracalliope fluviatilis, in seven  #
# localities in South Island, New Zealand from Lagrue et al. (2016 Biol. J.   #
# Linn. Soc. 118:344–358)                                                     #
 #############################################################################
#                  Documentation and scripts available at                     #
#                https://github.com/Ligophorus/RandomTaPas/                   #
 #############################################################################
# Load libraries
library(paco)
library(phytools)
library(distory)
library(parallel)
# Set number of runs (N) for Random TaPas
N= 1e+4
# Functions (6):
# foo 1 of 6
trimHS.maxC <- function (N, HS, n, check.unique= FALSE) {
  # For N runs, chooses n unique one-to-one associations and trims
  # the h-s association matrix to include the n associations only.
  #
  # Args:
  #   N:  Number of runs.
  #   HS: Host-symbiont association matrix.
  #   n:  Number of unique associations
  #   check.unique: if TRUE discards duplicated trimmed matrices.
  # Returns:
  #   A list of trimmed matrices.
  trim.int <- function (x, HS, n) {
    HS.LUT <- which(HS == 1, arr.in= TRUE)
    HS.LUT <- cbind(HS.LUT, 1:nrow(HS.LUT))
    df <- as.data.frame(HS.LUT)
    hs.lut <- subset(df[sample(nrow(df)), ],
                     !duplicated(row) & !duplicated(col))
    if (nrow(hs.lut) < n) hs <- NULL else {
      hs.lut <- hs.lut[sample(nrow(hs.lut), n), ]
      hs <- diag(nrow(hs.lut))
      rownames(hs) <- rownames(HS[hs.lut[ ,1], ])
      colnames(hs) <- colnames(HS[ ,hs.lut[ ,2]])
      return(hs)
    }
  }
  trim.HS <- lapply(1:N, trim.int, HS= HS, n= n )
  if (check.unique == TRUE) trim.HS <- unique(trim.HS)
  if (length(trim.HS) < N)
                    warning("No. of trimmed H-S assoc. matrices < No. of runs")
  return(trim.HS)
}
# foo 2 of 6.
geo.D <- function (hs, treeH, treeS) {
  # For any trimmed matrix produced with trimHS.maxC, it prunes the host &
  # symbiont phylogenies to conform with the trimmed matrix and computes the
  # geodesic distance between the pruned trees
  # NOTE: This function can only be used with strictly bifurcating trees.
  #
  # Args.:
  #   hs: trimmed matrix
  #   treeH: host phylogeny
  #   treeS: symbiont phylogeny
  # Returns:
  #   A geodesic distance
  treeh <- ape::drop.tip(treeH, setdiff(treeH$tip.label, rownames(hs)))
  trees <- ape::drop.tip(treeS, setdiff(treeS$tip.label, colnames(hs)))
   # foo distory requires same labels in both trees. Dummy labels are produced.
         # 1st reorder hs as per tree labels:
  hs <- hs[treeh$tip.label, trees$tip.label]
         # 2nd swap trees labels with corresponding ones in treeh:
  hs.lut <- which(hs[treeh$tip.label, trees$tip.label]==1, arr.ind = TRUE)
  dummy.labels <- rownames(hs.lut)
  trees$tip.label <- dummy.labels
  combo.tree <- list(treeh, trees)
  gd <- distory::dist.multiPhylo(combo.tree)
  return(gd)
}
# foo 3 of 6
paco.ss <- function (hs, treeH, treeS, symmetric= FALSE,
                         proc.warns= FALSE, ei.correct= "sqrt.D") {
  # For any trimmed matrix produced with trimHS.maxC, it prunes the host &
  # symbiont phylogenies to conform with the trimmed matrix and computes with
  # PAco the squared sum of residuals of the Procrustes superimosition of the
  # host and symbiont configurations in Euclidean space.
  #
  # Args.:
  #   hs: trimmed matrix
  #   treeH: host phylogeny
  #   treeS: symbiont phylogeny
  #   symmetric: specifies the type of Procrustes superimposition
  #   proc.warns: switches on/off trivial warnings returned when treeH and
  #               treeS differ in size
  #   ei.correct: specifies how to correct potential negative eigenvalues
  # Returns:
  #   A sum of squared residuals
  eigen.choice <- c("none", "lingoes", "cailliez", "sqrt.D")
  if (ei.correct %in% eigen.choice == FALSE)
                    stop(writeLines("Invalid eigenvalue correction parameter.\r
               Correct choices are 'none', 'lingoes', 'cailliez' or 'sqrt.D'"))
  treeh <- ape::drop.tip(treeH, setdiff(treeH$tip.label, rownames(hs)))
  trees <- ape::drop.tip(treeS, setdiff(treeS$tip.label, colnames(hs)))
  # Reorder hs as per tree labels:
  hs <- hs[treeh$tip.label, trees$tip.label]
  DH <- cophenetic(treeh)
  DP <- cophenetic(trees)
  if (ei.correct == "sqrt.D"){DH<- sqrt(DH); DP<- sqrt(DP); ei.correct="none"}
  D <- paco::prepare_paco_data(DH, DP, hs)
  D <- paco::add_pcoord(D, correction= ei.correct)
  if (proc.warns == FALSE) D <- vegan::procrustes(D$H_PCo, D$P_PCo,
                                                    symmetric = symmetric) else
    D <- suppressWarnings(vegan::procrustes(D$H_PCo, D$P_PCo,
                                                    symmetric = symmetric))
  return(D$ss)
}
# foo 4 of 6
link.freq <- function (x, fx, HS, percentile= 0.05, sep= "-", below.p= TRUE) {
  # Determines the frequency of each host-symbiont association occurring in a
  # given percentile of cases that maximize phylogenetic congruence. 
  #
  # Args.:
  #   x: list of trimmed matrices produced by trimHS.maxC
  #   fx: vector of statistics produced with either geo.D or paco.ss
  #   percentile: percentile to evaluate
  #   sep: character that separates host and symbiont labels
  #   below.p: determines whether frequencies ar to be computed below or above
  #            the percentile set
  # Returns:
  #   Data frame with labels of hosts, symbionts and host-symbiont associations
  #   in three columns. Column 4 displays the frequency of occurrence of each
  #   host-symbiont association in p
  if (below.p == TRUE) percent <- which(fx <= quantile(fx, percentile)) else
    percent <- which(fx >= quantile(fx, percentile))
  trim.HS <- x[percent]
  paste.link.names <- function(X, sep) {
    X.bin <- which(X>0, arr.in=TRUE)
    Y <- diag(nrow(X.bin))
    Y <- diag(nrow(X.bin))
    rownames(Y) <- rownames(X)[X.bin[,1]]
    colnames(Y) <- colnames(X)[X.bin[,2]]
    pln <- paste(rownames(Y), colnames(Y), sep=sep)
    return(pln)
  }
  link.names <- t(sapply(trim.HS, paste.link.names, sep=sep))
  lf <- table(link.names)
  lf <- as.data.frame(lf*100 / length(trim.HS))
  HS.LUT <- which(HS ==1, arr.in=TRUE)
  linkf <-  as.data.frame(cbind(rownames(HS)[HS.LUT[,1]],
                                colnames(HS)[HS.LUT[,2]]))
  colnames(linkf) <- c('H', 'S')
  linkf$HS <- paste(linkf[,1], linkf[,2], sep=sep)
  linkf$Freq <- rep(0, nrow(linkf))
  linkf[match(lf[,1], linkf[,3]), 4] <- lf[,2]
  return(linkf)
}
# foo 5 of 6
One2one.f <- function (hs, reps= 1e+4) {
  # For a matrix of host-symbiont associations, it finds the maximum n for
  #  which one-to-one unique associations can be picked in trimHS.maxC over
  #  a number of runs.
  #
  # Args.:
  #   hs: matrix of host-symbiont associations
  #   reps: number of runs to evaluate
  # Returns:
  #   maximum n
  HS.LUT <- which(hs ==1, arr.in=TRUE)
  HS.LUT <- cbind(HS.LUT,1:nrow(HS.LUT))
  df <- as.data.frame(HS.LUT)
  V <- rep(NA,reps)
  for(i in 1:reps){
    hs.lut <- subset(df[sample(nrow(df)),],
                     !duplicated(row) & !duplicated(col))
    n <- sum(HS)
    while (n >0) {
      n <- n-1;
      if (nrow(hs.lut) == n) break
    }
    V[i]<- n
  }
  V <- min(V)
  return(V)
}
# foo 6 of 6
tangle.gram <- function(treeH, treeS, hs, colgrad, rcolgrad= TRUE, nbreaks=50,
                        fqtab, node.tag=TRUE, cexpt=1, ...) {
  # wrapper of cophylo.plot of package phytools is used for mapping as heatmap
  # the host-symbiont frequencies estimated by Random TaPas on a tanglegram. It
  # also plots the average frequency of occurrence of each terminal and
  # optionally, the fast maximum likelihood estimators of ancestral states of
  # each node.
  #
  # Args.:
  #   treeH: host phylogeny
  #   treeS: symbiont phylogeny
  #   hs: host-symbiont association matrix
  #   colgrad: vector defining the color gradient of the heatmap
  #   rcolgrad: color scale relative to max/min frequency of input
  #             or absolute (0-100)
  #   nbreaks: number of discrete values along color scale
  #   fqtab: dataframe produced with link.freq
  #   node.tag: specifies whether maximum likelihood estimators of ancestral
  #             states will be computed 
  #   cexpt: size of color points at terminals and nodes
  #   ...: any option available in cophylo.plot  
  # Returns:
  #   A tanglegram with quantitative information displayed as heatmap.				
  rescale.range <- function(x) {
    if (rcolgrad==TRUE) {
      x <- round(x)
      y <- range(x)
      col_lim <- (y[1]: y[2])-y[1]+1
      x <- x-y[1]+1
      new.range <- list(col_lim, x)
    } else { 
      new.range <- list(1:101, round(x)+1)
    }
    return(new.range)
  }
  NR <- rescale.range(fqtab[,4])
  rbPal <- colorRampPalette(colgrad) 
  linkcolor <- rbPal(nbreaks)[as.numeric(cut(NR[[1]],breaks = nbreaks))]
  HS.lut <- which(hs ==1, arr.ind=TRUE)
  linkhs <- cbind(rownames(hs)[HS.lut[,1]], colnames(hs)[HS.lut[,2]])
  obj <- phytools::cophylo(treeH,treeS, linkhs)
  phytools::plot.cophylo(obj, link.col=linkcolor[NR[[2]]], ...)
  
  Hfreq <- aggregate(fqtab[,4], by=list(freq = fqtab[,1]), FUN=mean)
  Sfreq <- aggregate(fqtab[,4], by=list(freq = fqtab[,2]), FUN=mean)
  
  Hfreq <- Hfreq[match(obj$trees[[1]]$tip.label, Hfreq$freq),]
  Sfreq <- Sfreq[match(obj$trees[[2]]$tip.label, Sfreq$freq),]
  
  if (node.tag==TRUE){
  fit.H<-fastAnc(obj$trees[[1]],Hfreq[,2])
  fit.S<-fastAnc(obj$trees[[2]],Sfreq[,2])
  NLH <- rescale.range (fit.H)
  NLS <- rescale.range (fit.S)
  nodelabels.cophylo(pch=16, col=linkcolor[NLH[[2]]], cex=cexpt)
  nodelabels.cophylo(pch=16, col=linkcolor[NLS[[2]]], cex=cexpt, which="right")
  }
  TLH <- rescale.range (Hfreq[,2])
  TLS <- rescale.range (Sfreq[,2])
  tiplabels.cophylo(pch=16, col=linkcolor[TLH[[2]]], cex=cexpt)
  tiplabels.cophylo(pch=16, col=linkcolor[TLS[[2]]], cex=cexpt, which="right")
}
############ end of function declaration ######################################
#
# Read data (It is assumed that input files are in the working directiory
# of the R session)
TreeH <- read.nexus("Pfluviatilis.tre")	      # Amphipod (host) consensus tree
TreeS <- read.nexus("Cparvum.tre")	      # Trematode (symbiont) consensus tree
mTreeH <- read.nexus("m_Pfluviatilis.tre") # 1,000 post. prob. amphipod trees
mTreeS <- read.nexus("m_Cparvum.tre")	   # 1,000 post. prob. trematode trees
HS <- as.matrix(read.table("PfxCp_assoc.txt", header = TRUE)) # Assoc. matrix
#
# Optional: Draw a tanglegram
HS.lut <- which(HS ==1, arr.ind=TRUE)
linkhs <- cbind(rownames(HS)[HS.lut[,1]], colnames(HS)[HS.lut[,2]])
obj <- cophylo(TreeH,TreeS, linkhs)
plot.cophylo(obj, link.lwd=1, link.lty=2, fsize=0.5, pts=FALSE,
             link.type="curved")
dev.off() #reset graphics device
# Strategy to set n
# 1. Find minimum value with which all runs can be completed
One2one.f(HS,N) #the output is 5 (sometimes 6, depending on the randomization)
# 2. Let's evaluate 5:10
X <- 5:10
Y <- rep(NA, length(X))
for(i in 1:length(X)) {
  THS <- trimHS.maxC(N, HS, n=X[i], check.unique=TRUE)
  Y[i] <- length(THS)
}
plot(X,Y, type="b", xlab="Number of unique H-S associations",
                    ylab="Number of runs accomplished")
# 3. We choose 8
n=8
THS <- trimHS.maxC(N, HS, n=n, check.unique=TRUE)
THS[sapply(THS, is.null)] <- NULL
#
####### Apply 2 global-fit methods: GD and PACo ##########
# NOTE PACo symmetric= TRUE in this example!!!!
GD <- sapply(THS, geo.D, treeH=TreeH, treeS= TreeS)
PACO <- sapply(THS, paco.ss, treeH=TreeH, treeS= TreeS, symmetric=TRUE)
#
######### extract h-s association frequencies
LFGD05 <- link.freq(THS, GD, HS, percentile=0.05)
LFPACO05 <- link.freq(THS, PACO, HS, percentile=0.05)
#
# Procedure to compute CIs
GD05 <- matrix(NA, length(mTreeH), nrow(LFGD05))
PACO05 <- matrix(NA, length(mTreeH), nrow(LFPACO05))
#
cores <- detectCores()
pb <- txtProgressBar(min = 0, max = length(mTreeH), style = 3)
cl <- makeCluster(cores-2)    # use nº CPUs - 2 for parallel computing
for(i in 1:length(mTreeH))		# CIs for GD
{
  GD.CI<-parallel::parSapply(cl, THS, geo.D, treeH=mTreeH[[i]],
                             treeS= mTreeS[[i]])
  LFGD05.CI <- link.freq(THS, GD.CI, HS, percentile=0.05)
  GD05[i,] <- LFGD05.CI[,4]
  setTxtProgressBar(pb, i)
}
close(pb)
#
pb <- txtProgressBar(min = 0, max = length(mTreeH), style = 3)
#
for(i in 1:length(mTreeH))		# CIs for PACo
{
  PA.CI<-parallel::parSapply(cl, THS, paco.ss, treeH=mTreeH[[i]],
                             treeS= mTreeS[[i]], symmetric=TRUE)
  LFPA05.CI <- link.freq(THS, PA.CI, HS, percentile=0.05)
  PACO05[i,] <- LFPA05.CI[,4]
  setTxtProgressBar(pb, i)
}
close(pb)
stopCluster(cl)
#
colnames(GD05) <- LFGD05[,3]
colnames(PACO05) <- LFPACO05[,3]
# GD - compute CIs and average
GD.LO <- apply(GD05, 2, quantile, 0.025)
GD.HI <- apply(GD05, 2, quantile, 0.975)
GD.AV <- apply(GD05, 2, mean)
# PACo - compute CIs and average
PACO.LO <- apply(PACO05, 2, quantile, 0.025)
PACO.HI <- apply(PACO05, 2, quantile, 0.975)
PACO.AV <- apply(PACO05, 2, mean)
#
# Barplot of frequencies per host-symbiont association
op <- par(mfrow=c(2,1),mgp = c(1, 0.2, 0), mar=c(2.8,2.5,0.2,1.5), tck = 0.02)
link.fq <-barplot(GD.AV, names.arg = LFGD05$HS,
       horiz=FALSE, cex.names = 0.4, las=2, cex.axis=0.5, ylab="Frequency (%)",
       ylim=c(0, max(GD.HI)), col="lightblue")
suppressWarnings(arrows(link.fq, GD.HI, link.fq, GD.LO, length= 0,
                        angle=90, code=3, col="darkblue"))
abline(h=n/sum(HS)*100, col="red")

link.fq <-barplot(PACO.AV, names.arg = LFPACO05$HS,
       horiz=FALSE, cex.names = 0.4, las=2, cex.axis=0.5, ylab="Frequency (%)",
       ylim=c(0, max(PACO.HI)),col="lightblue")
suppressWarnings(arrows(link.fq, PACO.HI, link.fq, PACO.LO, length= 0,
                        angle=90, code=3, col="darkblue"))
abline(h=n/sum(HS)*100, col="red")
par(op)
# Plot variance-to-mean-ratio values
V2M <- function(x) var(x)/mean(x)
vmrGD <- V2M(LFGD05[,4])
vmrPA <- V2M(LFPACO05[,4])
vmrMGD <- apply(GD05, 1, V2M)
vmrMPA <- apply(PACO05, 1, V2M)
boxplot(vmrMGD, vmrMPA, names = c("GD", "PACo"), ylab="Variance to mean ratio",
        col="lightblue")
text(1,vmrGD,"*",cex=2, col="darkblue")
text(2,vmrPA,"*",cex=2, col="darkblue")
dev.off()
# Tanglegram - trees are 1st rendered ultramteric to separate nodes
# (facilitates visualization)
trh <- compute.brtime(TreeH, TreeH$Nnode:1)
trs <- compute.brtime(TreeS, TreeS$Nnode:1)
# Set color scale - this one is supposed to be color-blind friendly
col.scale <- c("darkred", "lightblue","darkblue")
tangle.gram(trh, trs, HS, colgrad=col.scale, rcolgrad= TRUE,
            nbreaks=50, LFGD05, link.lwd=1, link.lty=1, fsize=0.5,
            pts=FALSE, link.type="curved", node.tag=TRUE,
            cexpt=1.2)
############################ END OF SCRIPT #####################################
