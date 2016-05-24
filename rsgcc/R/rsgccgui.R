################################################
##gui.rsgcc
##Author: Chuang Ma
##Date: 2012-03-02


if( !require(gWidgetsRGtk2) ) install.packages("gWidgetsRGtk2", dependencies=TRUE)
if( !require(gWidgets) ) install.packages("gWidgets", dependencies=TRUE)
if( !require(fBasics)) install.packages("fBasics")
if( !require(gplots)) install.packages("gplots")
require(gWidgetsRGtk2)
require(gWidgets)
require(fBasics)
require(gplots)


options("guiToolkit"="RGtk2")


rsgcc.gui <- function(margins = c(1,1),
                      labRow = "", labCol = "",
                      lwid = c(0.5, 0.05, 0.01, 0.5, 0.01, 0.05, 0.5),
                      keynote = "FPKM"
                      ) {
  
availCor <- c("Gini correlation"="GCC", "Spearman correlation"="SCC", "Kendall correlation"="KCC", "Pearson correlation"="PCC", "Turkey's biweight"="BiWt")
availSim <- c("raw correlation(1-coef)" = "Raw", "absolute correlation(1-|coef|)" = "Abs", "squared correlation(1-|coef|^2)" = "Sqr")

availCluster <- c("complete linkage" = "complete", "average linkage" = "average", "median linkage" = "median", "centroid linkage" = "centroid", 
                   "mcquitty linkage" = "mcquitty", "single linkage" = "single", "ward linkage" = "ward")



data_file_dir <- NULL
data_matrix <- NULL
tsgene_matrix <- NULL
hm_dataframe <- NULL
CPUNum <- 1

DisplayFlag <- 0
TSGeneClusterFlag <- 0
CorMethodType <- "GCC"
SimilarityMethodType <- "Raw"
ClusterMethodType <- "complete"

OutputGPType <- FALSE
AbsCorThreshold <- 0.85
PermutationNum <- 2000
PValueThreshold <- 0.05

DisplayTimes <- 0

rowHC <- NULL
colHC <- NULL

ColMax <- 16711680
ColMean <- 122
ColMin <- 16776960

tsScore <- 0.95



  
##data to color index
coverDecToHex <- function(x) {
  
  ##tt <- x
  ##class(tt) <- "hexmode"
  tt <- .dec.to.hex(x)
  strlen <- str_length(tt)
  if( strlen > 6 ) {
    hh <- str_sub( tt, strlen-5, strlen)
  }else if( strlen < 6 ){
    hh <- paste(str_sub("0000000",1, 6-strlen), tt, sep="")
  }else {
    hh <- tt
  }
  paste("#", hh, sep="")
  
}


##start to run, update plot
updateRun <- function(h,...) {

  ##data_matrix
  x <- data_matrix
  if( is.null(x) ) {
    stop("Error: no GE data input.")
  }else {
    x <- data_matrix
  }
  cat("Starting to run...\n")
  print(SimilarityMethodType)
  if( TSGeneClusterFlag == 1 ) {
    tsgene_matrix <<- getsgene(x, tsThreshold = tsScore, Fraction= TRUE)$tsgene
    
    hm_dataframe <<- gcc.tsheatmap(tsgene_matrix, cpus = CPUNum, 
                            cormethod = CorMethodType, 
                            distancemethod  = SimilarityMethodType, 
                            clustermethod= ClusterMethodType, 
                            lwid = lwid, keynote = keynote,
                            margins = margins)
  }else {
    hm_dataframe <<- gcc.heatmap(x, cpus = CPUNum, 
                          cormethod = CorMethodType, 
                          distancemethod = SimilarityMethodType, 
                          clustermethod= ClusterMethodType, 
                          margins= margins, labRow = labRow, labCol = labCol)
  }
  cat("Finish!\n")
  
}


#file choose
fileChoose <- function(h, ...) {
 
   countLines <- function(filename) {
     cat("Line numer is:", length(readLines(filename)), "\n")
   }
   
   the_data <- gfile(text= "Select a file...", type="open", ..., action = "countLines", 
                      filter = list("All files" = list(patterns = c("*")),
                                    "text files" = list(patterns = c("*.txt")),
                                    "csv files" = list(patterns = c("*.csv"))),
                      handler = function(h,...) { do.call(h$action, list(h$file)) })
   
   data_file_dir <<- the_data
   cat("Loaded file dir: ", data_file_dir, "\n")
   
   
    lastdot <- sapply(gregexpr("\\.", the_data), tail, 1)
    if( lastdot < 0 )  { 
      data_matrix <<- as.matrix(read.table(the_data, sep="\t"))
    } else { 
        curfileName <- str_sub( the_data, lastdot+1, str_length(the_data)) 
        if(curfileName == "csv") {
          data_matrix <<- as.matrix(read.csv(the_data, sep="\t"))
        }else {
          data_matrix <<- as.matrix(read.table(the_data, sep="\t"))
        }
    }
}

#visualize data
visData <- function(h, ...) {
    
  if( svalue(h$obj) == TRUE ) {
    
    cat("Display GE Times: Yes.\n")
       
    if( DisplayTimes != 1 ) {
   
      if( is.null(data_matrix) ) {
        stop("Error: no GE data input")
      }
  
      cat("the dim info of data_matrix is:", dim(data_matrix), "\n")
  
      if( is.null( rownames(data_matrix)) ) { visx <- data_matrix }
      else {
        visx <- matrix(0, nrow = nrow(data_matrix), ncol = ncol(data_matrix)+1 )
        colnames(visx) <- c("rownames", colnames(data_matrix))
        visx[1:nrow(visx),2:ncol(visx)] <- data_matrix[1:nrow(data_matrix),1:ncol(data_matrix)]
        visx[,1] <- rownames(data_matrix)
      }
      
      DisplayTimes <<- 1
  
      win <- gwindow("Loaded gene expression data)")
      gp <- ggroup(container=win)
      group <- ggroup(horizontal=FALSE, border = TRUE, container=gp)
      odm <- gtable(visx, multiple = TRUE,  container= gp, expand = TRUE)
     
      cat("Finish to display GE data: Yes.\n")
    }

  }else {
    #do nothing
    DisplayTimes <<- 0
  }
}

updateHeatmapMax <- function(h,...) {
  
  ColMax <<- svalue(h$obj)

  cat("Color Index for Heatmap:", ColMax, ColMean, ColMin, "\n")
  
  #for TS genes, only two colors
  if( TSGeneClusterFlag == 1 ) {
      hmdata <- gcc.tsheatmap(tsgene_matrix, cpus = CPUNum, 
                           method = CorMethodType, 
                           distancemethod = SimilarityMethodType, 
                           clustermethod= ClusterMethodType,
                           rowhcdata = hm_dataframe$hcr,
                           colhcdata = hm_dataframe$hcc,
                           colrange = c(coverDecToHex(ColMin), coverDecToHex(ColMax)),
                           lwid = lwid)
    
  }else {
      hmdata <- gcc.heatmap(data_matrix, cpus = CPUNum, 
                            method = CorMethodType, 
                            distancemethod = SimilarityMethodType, 
                            clustermethod= ClusterMethodType,
                            rowhcdata = hm_dataframe$hcr,
                            colhcdata = hm_dataframe$hcc,
                            colrange = c(coverDecToHex(ColMin), coverDecToHex(ColMean), coverDecToHex(ColMax)),
                            margins= margins, labRow = "")
  }
}


##needed to be modified
updateHeatmapMean <- function(h,...) {
  
  ColMean <<- svalue(h$obj)
  cat("Color Index for Heatmap:", ColMax, ColMean, ColMin, "\n")
  
  #for TS genes, only two colors
  if( TSGeneClusterFlag == 0 ) {
      hmdata <- gcc.heatmap(data_matrix, cpus = CPUNum, 
                            method = CorMethodType, 
                            distancemethod = SimilarityMethodType, 
                            clustermethod= ClusterMethodType,
                            rowhcdata = hm_dataframe$hcr,
                            colhcdata = hm_dataframe$hcc,
                            colrange = c(coverDecToHex(ColMin), coverDecToHex(ColMean), coverDecToHex(ColMax)),
                            margins= margins, labRow = "")
  }else {
    #do nothing    
  }
  

}

##needed to be modified
updateHeatmapMin <- function(h,...) {
  
  ColMin <<- svalue(h$obj)
  
  cat("Start to save data.\n")
  cat("Color Index for Heatmap:", ColMax, ColMean, ColMin, "\n")
  
  #for TS genes, only two colors
  if( TSGeneClusterFlag == 1 ) {
      hmdata <- gcc.tsheatmap(tsgene_matrix, cpus = CPUNum, 
                           method = CorMethodType, 
                           distancemethod = SimilarityMethodType, 
                           clustermethod= ClusterMethodType,
                           rowhcdata = hm_dataframe$hcr,
                           colhcdata = hm_dataframe$hcc,
                           colrange = c(coverDecToHex(ColMin), coverDecToHex(ColMax)),
                           lwid = lwid)
    
  }else {    
      hmdata <- gcc.heatmap(data_matrix, cpus = CPUNum, 
                            method = CorMethodType, 
                            distancemethod = SimilarityMethodType, 
                            clustermethod= ClusterMethodType,
                            rowhcdata = hm_dataframe$hcr,
                            colhcdata = hm_dataframe$hcc,
                            colrange = c(coverDecToHex(ColMin), coverDecToHex(ColMean), coverDecToHex(ColMax)),
                            margins= margins, labRow = "")
    
  }
}

saveData <- function(h,...) {
  pdf(paste(data_file_dir, "_rsgcc_heatmap.pdf", sep=""), width = 5, height= 5)
  if( TSGeneClusterFlag == 1 ) {
      hm_dataframe <<- gcc.tsheatmap(tsgene_matrix, cpus = CPUNum, 
                           method = CorMethodType, 
                           distancemethod = SimilarityMethodType, 
                           clustermethod= ClusterMethodType,
                           rowhcdata = hm_dataframe$hcr,
                           colhcdata = hm_dataframe$hcc,
                           colrange = c(coverDecToHex(ColMin), coverDecToHex(ColMax)),
                           lwid = lwid)
    
  }else {
      hm_dataframe <<- gcc.heatmap(data_matrix, cpus = CPUNum, 
                            method = CorMethodType, 
                            distancemethod = SimilarityMethodType, 
                            clustermethod= ClusterMethodType,
                            rowhcdata = hm_dataframe$hcr,
                            colhcdata = hm_dataframe$hcc,
                            colrange = c(coverDecToHex(ColMin), coverDecToHex(ColMean), coverDecToHex(ColMax)),
                            margins= margins, labRow = "")
  }
  dev.off()
  
  GEMatrix_Data <- NULL
  cor_matrix <- NULL
  if( TSGeneClusterFlag == 1 ) {
    GEMatrix_Data <- tsgene_matrix
    cor_matrix <- hm_dataframe$hcr$pairmatrix
  }else {
    GEMatrix_Data <- data_matrix
    cor_matrix <- hm_dataframe$hcr$pairmatrix
  }
  cor_file_dir <- paste(data_file_dir, "_rsgcc_correlations.txt", sep="")
  genenum <- nrow(cor_matrix)
  cor_pairs <- matrix("NA", nrow = genenum*(genenum-1), ncol = 3)
  k <- 0
  for( i in 1:(genenum-1)) {
    for( j in (i+1):genenum){
      k <- k + 1
      cor_pairs[k,1] <- rownames(cor_matrix)[i]
      cor_pairs[k,2] <- colnames(cor_matrix)[j]
      cor_pairs[k,3] <- cor_matrix[i,j]
    }#end for j
  }#end for i
  
  write.table( cor_pairs[1:k,], cor_file_dir, sep="\t",  row.names = FALSE, col.names = FALSE, quote = FALSE)

  if( !require(ctc) ) { 
      biocLite <- NULL
      rm(biocLite)
      source("http://bioconductor.org/biocLite.R")
      biocLite("ctc")
  }
  require(ctc)
  
  cdt_file_dir <- paste(data_file_dir, "_rsgcc_cluster.cdt", sep="")
  r2atr( hm_dataframe$hcc$hc, distance = "GCC", file = paste(data_file_dir, "_rsgcc_cluster.atr", sep=""))
  r2gtr( hm_dataframe$hcr$hc, distance = "GCC", file = paste(data_file_dir, "_rsgcc_cluster.gtr", sep=""))
  r2cdt(hm_dataframe$hcr$hc, hm_dataframe$hcc$hc, GEMatrix_Data, file = paste(data_file_dir, "_rsgcc_cluster.cdt", sep=""))
  cat("Success in saving data.\n")
   
}

###########################################################################

## now layout
window <- gwindow("rsgcc(correlation and clustering analysis of gene expression data)")
BigGroup <- ggroup(container=window)
group <- ggroup(horizontal=FALSE, border = TRUE, container=BigGroup)


#add a button 
tmp <- gframe("Step 1: Load gene exp data", container = group)
Selectedfiles <- gbutton(c("Click here to load"), handler= fileChoose, container=group)

#read and display data
tmp <- gcheckbox("Display loaded data", checked = FALSE, handler = visData, container = group)


#find and cluster tissue/condition specific genes
tmp <- gcheckbox("Find ts-genes for clustering analysis", checked = FALSE, 
                 handler = function(h,...) { 
                   if( svalue(h$obj) == TRUE ) { 
                     TSGeneClusterFlag <<- 1
                     cat("rsgcc will find and cluster ts-genes.\n")
                    }else {
                      TSGeneClusterFlag <<- 0
                    } }, 
                 container = group)

#threshold for tsScore
tmp <- gframe("Threshold for tissue specificity score", container=group)
tsScoreAdjust <- gedit(text = "0.95", width = 30, coerce.with=as.numeric, horizontal=FALSE, 
                   handler=function(h,...){tsScore <<- svalue(h$obj); cat("Threshold for tissue specificity score:", tsScore, "\n")})
add(tmp, tsScoreAdjust)




tmp <- gframe("Step 2: Select a correlation method", container=group)
Correlations <- gradio(names(availCor), horizontal=FALSE, 
                       handler= function(h,...){ CorMethodType <<- availCor[svalue(Correlations)]; cat("Correlation method: ", CorMethodType, "\n")})
add(tmp, Correlations)


tmp <- gframe("Step 3: Specify a distance measure", container=group)
Similarities <- gradio(names(availSim), horizontal=FALSE, 
                       handler=function(h,...){SimilarityMethodType <<- availSim[svalue(Similarities)]; cat("Similarity method: ", SimilarityMethodType, "\n" )})
add(tmp, Similarities)

tmp <- gframe("Step 4: Choose a cluster method", container=group)
Clusters <- gcombobox(names(availCluster), horizontal=FALSE, 
                      handler=function(h,...) { ClusterMethodType <<- availCluster[svalue(Clusters)]; cat("Cluster method:", ClusterMethodType, "\n")} )
add(tmp, Clusters)

#number of cpus for computation
tmp <- gframe("Step 5: CPUs for correlation calcuation", container=group)
cpuAdjust <- gedit(text = "1", width = 30, coerce.with=as.numeric, horizontal=FALSE, 
                   handler=function(h,...){CPUNum <<- svalue(h$obj); cat("CPU number:", CPUNum, "\n")})
add(tmp, cpuAdjust)

StartRun <- gbutton(c("Start to run"), handler= updateRun, container=group)


#########################################################################
##set color for cluster
tmp <- gframe("Adjust colors for heat map", container=group)

tmp <- gframe("Color for max GE value", container=group)
ColorAdjust1 <- gslider(from=1, to=256^3-1, by=2000, value=16711680, handler=updateHeatmapMax)
add(tmp, ColorAdjust1, expand = TRUE)

tmp <- gframe("Color for median GE value", container=group)
ColorAdjust2 <- gslider(from=1, to=256^3-1, by=2000, value=1, handler=updateHeatmapMean)
add(tmp, ColorAdjust2, expand = TRUE)

tmp <- gframe("Color for min GE value", container=group)
ColorAdjust3 <- gslider(from=1, to=256^3-1, by=2000, value=16776960, handler=updateHeatmapMin)
add(tmp, ColorAdjust3, expand = TRUE)

##########################################################################


SaveData <- gbutton(c("Save correlations and cluster data"), handler= saveData, container=group)

add(BigGroup, ggraphics())


}

