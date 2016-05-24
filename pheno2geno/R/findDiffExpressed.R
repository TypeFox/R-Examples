#
# find.diff.expressed.R
#
# Copyright (c) 2010-2012 GBIC: Danny Arends, Konrad Zych and Ritsert C. Jansen
# last modified May, 2012
# first written Mar, 2011
# Contains: find.diff.expressed, fake.population, plotRPpval
#

# find.diff.expressed
#
# DESCRIPTION:
#  Using Rank Product or student t-test analysis to select differentially expressed genes.
# PARAMETERS:
#   - population - Object of class population , must contain founders phenotypic data.
#   - verbose - Be verbose
#   - debugMode - 1: Print our checks, 2: print additional time information
#   - ... - parameters send to RP function
# OUTPUT:
#  An object of class population containing object of class RP in $founders$RP
#
find.diff.expressed <- function(population,use=c("ttest","rankprod"),verbose=FALSE,debugMode=0,...){
  #checks
  if(missing(population)) stop("provide population object\n")
  check.population(population)
  use <- match.arg(use)
  if(verbose && debugMode==1) cat("find.diff.expressed starting withour errors in checkpoints.\n")
  
  s<-proc.time()
  if(use=="rankprod"){
   if (requireNamespace("RankProd", quietly = TRUE)){
      rankProdRes <- RankProd::RP(population$founders$phenotypes,population$founders$groups,gene.names=rownames(population$founders$phenotypes),...)
      population$founders$RP <- rankProdRes
    }else{
      stop("Install RankProd package to use Rank Product analysis!\n")
    }
  }else{
    population$founders$RP$pval<- t(rbind(apply(population$founders$phenotypes,1,findUsingTTest.internal,population$founders$groups)))
  }
  e<-proc.time()
  if(verbose && debugMode==2)cat("Differentially expressed genes found in:",(e-s)[3],"seconds.\n")
  invisible(population)
}

############################################################################################################
#                  *** findUsingTTest.internal ***
#
# DESCRIPTION:
#  subfunction of find.diff.expressed using t-test to assess whether gene is differentially expressed
# 
# PARAMETERS:
#   phenoRow - single row of founders phenotype data
#   groupLabels - Specify which column of founders data belongs to group 0 and which to group 1.
#
# OUTPUT:
#  two p-values - for gene being up- and downregulated
#
############################################################################################################
findUsingTTest.internal <- function(phenoRow,groupLabels){
  a <- which(groupLabels==0)
  b <- which(groupLabels==1)
  if(mean(phenoRow[a],na.rm=T) < mean(phenoRow[b],na.rm=T)){
    what <- "less"
    return(c(0,t.test(phenoRow[a],phenoRow[b],alt=what)$p.value))
  }else{
    what <- "gre"
    return(c(t.test(phenoRow[a],phenoRow[b],alt=what)$p.value,0))
  }
}

############################################################################################################
#                  *** showRPpval ***
#
# DESCRIPTION:
#  showing pvals of RP for selected markers
# 
# PARAMETERS:
#   population - Object of class population , must contain founders phenotypic data.
#   markers - markers (specified by number) to be shown
# 
# OUTPUT:
#  none
#
############################################################################################################
showRPpval <- function(population,markers=1:10){
  #checks
  if(missing(population)) stop("provide population object\n")
  check.population(population)
  if(!is.numeric(markers)) stop("markers parameter must be numeric.\n")
  if(any(markers<1) || any(markers>nrow(population$founders$phenotypes))) stop("markers parameter must contain only values in between 1 and nr of markers (",nrow(population$founders$phenotypes),").\n")
  
  if(is.null(population$founders$RP$pval)) stop("Population object does not contain results of RP analysis run find.diff.expressed first.\n")
  
  toPrint <- matrix(0,length(markers),2)
  toPrint[,1] <- population$founders$RP$pval[markers,1]
  toPrint[,2] <- population$founders$RP$pval[markers,2]
  rownames(toPrint) <- rownames(population$founders$phenotypes)[markers]
  colnames(toPrint) <- c("up","down")
  print(toPrint)
}

############################################################################################################
#                  *** plotRPpval ***
#
# DESCRIPTION:
#  ploting pvals of RP for selected markers
# 
# PARAMETERS:
#   population - Object of class population , must contain founders phenotypic data.
#   markers - markers (specified by number) to be shown
#  treshold - treshold value, on which line is plotted (by default - 0.01)
# 
# OUTPUT:
#  none
#
############################################################################################################
plotRPpval <- function(population,thresholdRange=c(0.01,0.1,0.01)){
  #checks
  if(missing(population)) stop("provide population object\n")
  check.population(population)
  if(is.null(population$founders$RP$pval)) stop("Population object does not contain results of RP analysis run find.diff.expressed first.\n")
  thrRange <- seq(thresholdRange[1],thresholdRange[2],thresholdRange[3])
  n.upSelected <- NULL
  n.downSelected <- NULL
  for(threshold in thrRange){
    upNotNull <- population$founders$RP$pval[which(population$founders$RP$pval[,1] > 0),1]
    downNotNull <- population$founders$RP$pval[which(population$founders$RP$pval[,2] > 0),2]
    n.upSelected <- c(n.upSelected,length(which(upNotNull < threshold)))
    n.downSelected <- c(n.downSelected,length(which(downNotNull < threshold)))
  }
  xlim <- c(thresholdRange[1],thresholdRange[2])
  ylim <- c(min(c(n.upSelected, n.downSelected)),max(c(n.upSelected, n.downSelected)))
  plot(thrRange,n.upSelected ,main="RP analysis p-values",xlab="p-value",ylab="# markers selected", xlim=xlim, ylim=ylim, type="o")
  points(thrRange,n.downSelected,col="red",type="o")
  legend(x="topleft",legend=c("up regulated","down regulated"),col=c("black","red"),cex=0.8,pch=21,lwd=2,bg="white")
}
