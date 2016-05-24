# ldheatmap - Plots measures of pairwise linkage disequilibria for SNPs
# Copyright (C) 2004  J.Shin, S. Blay, N. Lewin-Koh, J.Graham, B.McNeney

# To cite LDheatmap in publications use:
# Shin J-H, Blay S, McNeney B and Graham J (2006). LDheatmap: An R
# Function for Graphical Display of Pairwise Linkage Disequilibria
# Between Single Nucleotide Polymorphisms. J Stat Soft, 16 Code Snippet 3

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

###########################################################################

"LDheatmap"<-
   function (gdat, genetic.distances=NULL, 
             distances="physical", LDmeasure="r", title="Pairwise LD",
             add.map=TRUE, add.key=TRUE, geneMapLocation=0.15, 
             geneMapLabelX=NULL, geneMapLabelY=NULL, 
             SNP.name=NULL, color=NULL,
             newpage=TRUE, name="ldheatmap", vp.name=NULL,
             pop=FALSE, flip=NULL, text=FALSE)
{

  makeImageRect <- function(nrow, ncol, cols, name, byrow=TRUE) {
    xx <- (1:ncol)/ncol   
    yy <- (1:nrow)/nrow
    if(byrow) {
      right <- rep(xx, nrow)
      top <- rep(yy, each=ncol)
    } else {
      right <- rep(xx, each=nrow)
      top <- rep(yy, ncol)
    }
    rectGrob(x=right, y=top, 
           width=1/ncol, height=1/nrow, 
           just=c("right", "top"), 
           gp=gpar(col=NA, fill=cols),
           name=name)
  }

  makeImageText <- function(nrow, ncol, cols, name) {
    cols <- as.character(cols)
    cols[is.na(cols)] <- ""
    cols <- paste("   ", cols)
    xx <- (1:ncol)/ncol   
    yy <- (1:nrow)/nrow
    right <- rep(xx, nrow)
    top <- rep(yy, each=ncol)
    textGrob(cols, x=right, y=top, 
           gp=gpar(cex=0.3),
           just=c("right", "top"), 
           name=name)
  }

  #_______________________Color Key__________________________________________##
  # Draw the Color Key
  if (is.null(color)) {
    if (inherits(gdat, "LDheatmap")) color <- gdat$color
    else  color <- grey.colors(20)
  }
  LDheatmap.Legend.add <- function(color, vp=heatmapVP){
    ImageRect<- makeImageRect(2,length(color), col=c(rep(NA,length(color)),color[length(color):1]),
                              "colorKey")
    keyVP <- viewport(x=1.1, y=-.10, height=.10, width=.5, just=c("right","bottom"), name="keyVP")
    #Adding the label 'Color key'
    if(LDmeasure=="r") {
       ttt<-expression(paste(R^2," Color Key"))
    } else {
      ttt<-"D' Color Key"
    }
    title<-textGrob(ttt, x=0.5, y=1.25, name="title", gp=gpar(cex=0.8))

    #Adding labels to the color key
    labels<-textGrob(paste(0.2*0:5), x=0.2*0:5,y=0.25, gp=gpar(cex=0.6), name="labels")

    #Drawing ticks at the bottom axis of the color key
    ticks<-segmentsGrob(x0=c(0:5)*0.2 , y0=rep(0.4,6), x1=c(0:5)*0.2 , y1=rep(0.5,6),name="ticks")

    #Drawing a box around the color key
    box <- linesGrob(x=c(0,0,1,1,0), y=c(0.5,1,1,0.5,0.5), name="box")

    key <- gTree(children=gList(ImageRect, title, labels, ticks, box), name = "Key", vp=keyVP)
    key
  }

  #_______________________Genetic Map________________________________________##
  # adds a genetic map to the heatmap plot along the diagonal
  if (is.null(flip)) {
    if (inherits(gdat, "LDheatmap") && !is.null(gdat$flipVP)) flip <- TRUE 
    else flip <- FALSE
  }
  LDheatmap.Map.add <- function(nsnps, add.map, genetic.distances, 
                       geneMapLocation=0.15,
                       geneMapLabelX=NULL, geneMapLabelY=NULL,
                       distances="physical", vp=NULL, 
                       SNP.name=NULL, ind=0, flip=FALSE){
    snp <- ((1:nsnps-1) + 0.5) / nsnps
#####################
    if(add.map){
    min.dist <- min(genetic.distances) 
    max.dist <- max(genetic.distances)
    total.dist <- max.dist - min.dist

    if(flip) geneMapLocation<- (-geneMapLocation)

    # Drawing the diagonal line 
    seq.x <- c(0.5*geneMapLocation + 1/(nsnps*2),
             1+0.5*geneMapLocation - 1/(nsnps*2))
    seq.y <- c(-0.5*geneMapLocation + 1/(nsnps*2),
             1-0.5*geneMapLocation - 1/(nsnps*2))
    diagonal<-linesGrob(seq.x, seq.y, gp=gpar(lty=1), name="diagonal", vp=vp)

    ## Adding line segments to the plot: (point1 <-> point2) 
    ## point1: relative position of a SNP on the scaled line
    ## point2: position of that SNP on the LD image  
    regionx <- seq.x[1] +
             ((genetic.distances-min.dist)/total.dist)*(seq.x[2]-seq.x[1])
    regiony <- seq.y[1] +
             ((genetic.distances-min.dist)/total.dist)*(seq.y[2]-seq.y[1]) 
    segments <- segmentsGrob(snp, snp, regionx, regiony, name="segments", vp=vp)

    ## Adding the text indicating Physical length of the region under study
    if (distances=="physical")
      mapLabel <- paste("Physical Length:", round((total.dist/1000),1),
                      "kb", sep="")
    else 
      mapLabel <- paste("Genetic Map Length:", round(total.dist,1),"cM",sep="")

    if (!flip) {
      if(is.null(geneMapLabelY)) geneMapLabelY <- 0.3
      if(is.null(geneMapLabelX)) geneMapLabelX <- 0.5
    }
    else {
      if(is.null(geneMapLabelY)) geneMapLabelY <- 0.8
      if(is.null(geneMapLabelX)) geneMapLabelX <- 0.4
    }
    title <- textGrob(mapLabel, geneMapLabelX, geneMapLabelY,
              gp=gpar(cex=0.9), just="left", name="title")

    geneMap <- gTree(children=gList(diagonal, segments, title), name="geneMap")

    ## Labelling some SNPs 
    if (!is.null(SNP.name) && (any(ind!=0))){
      symbols <- pointsGrob(snp[ind], snp[ind], pch="*",
                gp=gpar(cex=1.25, bg="blue", col="blue"), name="symbols", vp=vp)
      SNPnames <- textGrob(paste(" ", SNP.name), just="left", rot=-45,
              regionx[ind], regiony[ind], gp=gpar(cex=0.6, col="blue"), name="SNPnames", vp=vp)
      if (flip) {
        lenght_SNP_name <- max(nchar(SNP.name))
        long_SNP_name <- paste(rep(8,lenght_SNP_name), collapse="")
        name_gap <- convertWidth(grobWidth(textGrob(long_SNP_name)), "npc",valueOnly=TRUE)/sqrt(2)
        diagonal<-linesGrob(seq.x, seq.y, gp=gpar(lty=1), name="diagonal", vp=vp)
        #diagonal<-linesGrob(seq.x+name_gap, seq.y-name_gap, gp=gpar(lty=1), name="diagonal", vp=vp)
        segments <- segmentsGrob(snp, snp, regionx, regiony, name="segments", vp=vp)
        #segments <- segmentsGrob(snp+name_gap, snp-name_gap, regionx+name_gap, regiony-name_gap, name="segments", vp=vp)
        symbols <- NULL
        SNPnames <- textGrob(SNP.name, just="left", rot=-45,
              regionx[ind]-name_gap, regiony[ind]+name_gap, gp=gpar(cex=0.6, col="blue"), name="SNPnames", vp=vp)
              # snp[ind], snp[ind], gp=gpar(cex=0.6, col="blue"), name="SNPnames", vp=vp)
        title <- editGrob(title, y=unit(geneMapLabelY+name_gap, "npc"))
      }
      geneMap <- gTree(children=gList(diagonal, segments, title, symbols, SNPnames),name="geneMap")
    }} # if(add.map) end

    else if (!add.map && !is.null(SNP.name) && (any(ind!=0))){
      geneMap <- textGrob(paste(" ", SNP.name), just="left", rot=-45,
                          snp[ind], snp[ind], gp=gpar(cex=0.6, col="blue"),
                          name="SNPnames")
      if (flip) geneMap <- editGrob(geneMap, vp=vp)
    }
    else geneMap <- NULL

    geneMap
  }

  #____________________________________________________________________________#
  ## If genetic.distances is missing, calculate an equispaced default:
  if(is.null(genetic.distances)) {
     if (inherits(gdat,"data.frame"))
        genetic.distances=1000*(1:ncol(gdat))
     else if(inherits(gdat,"matrix"))
        genetic.distances=1000*(1:length(gdat[1,]))
     else    # gdat is of class LDheatmap
        genetic.distances = gdat$genetic.distances
  }
  #____________________________________________________________________________#
  ## Calculate or extract LDmatrix

    if(inherits(gdat,"snp.matrix")){
      require(chopsticks)
      ## Exclude SNPs with less than 2 alleles:
      # NOT YET IMPLEMENTED for snp.matrix
      #gvars <- unlist(sapply(gdat, function(x) genetics::nallele(x) == 2))
      #genetic.distances <- genetic.distances[gvars]
      #gdat <- gdat[gvars]

      ## Sort data in ascending order of SNPs map position:
      if(!is.vector(genetic.distances))
        {stop("Distance should be in the form of a vector")}
      o<-order(genetic.distances)
      genetic.distances<-genetic.distances[o]
      gdat<-gdat[,o]
      myLD <- ld.snp(gdat,depth=ncol(gdat))
      if(LDmeasure=="r")
        LDmatrix <- myLD[["rsq2"]]   
      else if (LDmeasure=="D")
        LDmatrix <- myLD[["dprime"]]  
      else 
        stop("Invalid LD measurement, choose r or D'.")      
      # LDmatrix is upper-left-triangular, rather than the usual upper-right.
      nsnp<-length(genetic.distances)
      tem<-matrix(NA,nrow=nsnp,ncol=nsnp)
      for(i in 1:(nsnp-1)) { tem[i,(i+1):nsnp]<-LDmatrix[i,1:(nsnp-i)] }
      LDmatrix<-tem # need something faster than the for loop
      row.names(LDmatrix)<-attr(myLD,"snp.names")
    }

    else if(inherits(gdat,"data.frame")){
      for(i in 1:ncol(gdat)) {
        if(!genetics::is.genotype(gdat[,i]))
          stop("column ",i," is not a genotype object\n")
      }

      ## Exclude SNPs with less than 2 alleles:
      gvars <- unlist(sapply(gdat, function(x) genetics::nallele(x) == 2))
      genetic.distances <- genetic.distances[gvars]
      gdat <- gdat[gvars]

      ## Sort data in ascending order of SNPs map position:
      if(!is.vector(genetic.distances))
        {stop("Distance should be in the form of a vector")}
      o<-order(genetic.distances)
      genetic.distances<-genetic.distances[o]
      gdat<-gdat[,o]
      myLD <- genetics::LD(gdat)
      if(LDmeasure=="r")
        LDmatrix <- myLD[[LDmeasure]]^2   
      else if (LDmeasure=="D'")
        LDmatrix <- abs(myLD[[LDmeasure]])  
      else 
        stop("Invalid LD measurement, choose r or D'.")      
    }
    else if(inherits(gdat,"LDheatmap")){
      LDmatrix <- gdat$LDmatrix
      distances <- gdat$distances
    }
    else if(inherits(gdat,"matrix")){
      if(nrow(gdat) != ncol(gdat))
        stop("The matrix of linkage disequilibrium measurements must be a square matrix")
      LDmatrix <- gdat
      LDmatrix[lower.tri(LDmatrix, diag=TRUE)] <- NA
    }
    else if(!missing(gdat))  
      stop(paste("No method for an object of class",class(gdat)))
    else
      stop("Need to supply LD matrix or genotypes")

  #____________________________________________________________________________#
  ## Draw the heat map
  heatmapVP <- viewport(width = unit(.8, "snpc"), height = unit(.8, "snpc"),
                       name=vp.name)
  flipVP <- viewport(width = unit(.8, "snpc"), height= unit(.8, "snpc"), y=0.6, angle=-45, name="flipVP")

  if(color[1]=="blueToRed") color = rainbow(20, start=4/6, end=0, s=.7)[20:1]
  if(newpage)grid.newpage()
  mybreak <- 0:length(color)/length(color)

  imgLDmatrix <- LDmatrix
  # if (flip) imgLDmatrix <-t(imgLDmatrix[,dim(LDmatrix)[2]:1])[,dim(LDmatrix)[2]:1]
  byrow<-ifelse(flip,FALSE,TRUE) #FALSE if flip=TRUE

  colcut <- as.character(cut(1-imgLDmatrix,mybreak,labels=as.character(color), include.lowest=TRUE))
  if(is.numeric(color)) colcut <- as.integer(colcut)
  ImageRect<-makeImageRect(dim(LDmatrix)[1],dim(LDmatrix)[2],colcut, name="heatmap",byrow)
  ImageText <- NULL
  if (text) ImageText<-makeImageText(dim(LDmatrix)[1],dim(LDmatrix)[2], round(imgLDmatrix, digits = 2), name="heatmaptext")
  title <- textGrob(title, 0.5, 1.05, gp=gpar(cex=1.0), name="title")
  if (flip) {
     ImageRect <- editGrob(ImageRect, vp=flipVP)
     if (text)
        ImageText <- editGrob(ImageText, vp=flipVP, rot=45, just="left")
  }
  heatMap <- gTree(children=gList(ImageRect, ImageText, title), name="heatMap")
  #____________________________________________________________________________#
  ## Draw a diagonal line indicating the physical or genetic map positions of the SNPs
  nsnps <- ncol(LDmatrix)
  step <- 1/(nsnps-1)
  ind <- match(SNP.name, row.names(LDmatrix), nomatch=0)
  geneMapVP <- NULL
  if (flip) geneMapVP <- flipVP
  geneMap <- LDheatmap.Map.add (nsnps, genetic.distances=genetic.distances,
                     geneMapLocation=geneMapLocation,add.map, 
                     geneMapLabelX=geneMapLabelX,
                     geneMapLabelY=geneMapLabelY,
                     distances=distances, vp=geneMapVP, 
                     SNP.name=SNP.name, ind=ind, flip=flip)
  #____________________________________________________________________________#
  ## Draw the Color Key
  if(add.key) Key <- LDheatmap.Legend.add(color, vp=heatmapVP)
  else Key <- NULL

  ## Assemble the heatmap, genetic map and color key into a grob and draw it
  LDheatmapGrob<-gTree(children=gList(heatMap, geneMap, Key),
                      vp=heatmapVP, name=name, cl="ldheatmap")
  grid.draw(LDheatmapGrob)
  if(pop){
    downViewport(heatmapVP$name)
    popViewport()} #pop the heat map viewport

  ldheatmap <- list(LDmatrix=LDmatrix, LDheatmapGrob=LDheatmapGrob, heatmapVP=heatmapVP, flipVP=geneMapVP,                
                    genetic.distances=genetic.distances, distances=distances, color=color)
  class(ldheatmap) <- "LDheatmap"
  invisible(ldheatmap)
} # function LDheatmap ends



preDrawDetails.ldheatmap <- function(x) {
  fontsize <- convertX(unit(1/20,"grobwidth", rectGrob()), "points")
  pushViewport(viewport(gp=gpar(fontsize=fontsize)))
}


postDrawDetails.ldheatmap <- function(x) {
  popViewport()
}


preDrawDetails.symbols <- function(x) {
  fontsize <- convertX(unit(1/20,"grobwidth", rectGrob()), "points")
  pushViewport(viewport(gp=gpar(fontsize=fontsize)))
}

postDrawDetails.symbols <- function(x) {
  popViewport()
}


