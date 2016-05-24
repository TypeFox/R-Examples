
###### CLASSES

# Basics, all that is needed for the bootstrap replicates

setClass("cabasicresults",
representation(
  Rprofile="matrix", Cprofile="matrix", Rweights="matrix", Cweights="matrix",
  Raxes="matrix", Caxes="matrix", r="numeric", mu="numeric") )

# Full, with variances etc for CRs
#  covariances are arrays with dim1=row/col, dim2/3=matrix (only upper triangle is non-zero)
#  axisvariances is number of axes for which vars/covs were calculated

setClass("cabootcrsresults", 
representation(
  br="cabasicresults", 
  DataMatrix="matrix", rows="numeric", columns="numeric", 
  rowlabels="character", collabels="character", 
  Rowprinccoord="matrix", Colprinccoord="matrix", Rowstdcoord="matrix", Colstdcoord="matrix",
  RowCTR="matrix", RowREP="matrix", ColCTR="matrix", ColREP="matrix", 
  RowVar="matrix", RowCov="array", ColVar="matrix", ColCov="array", 
  inertiasum="numeric", inertias="matrix",
  nboots="numeric", resampledistn="character", multinomialtype="character",
  sameaxisorder="numeric", 
  poissonzeronewmean="numeric", newzeroreset="numeric", 
  printdims="numeric", axisvariances="numeric" ) )


###### FUNCTIONS


#### Prints reasonably full results, including variances etc

printca <- function(x, datasetname="") {

#setGeneric("print", function(x,...) standardGeneric("print") )
#setMethod("print", signature(x="cabootcrsresults"), 
# function(x, datasetname="") {

## Printing macro
printwithaxes <- function(res, thenames) { 
names(res) <- thenames
print(res, digits=4)
}
## 

d <- min(x@printdims, x@br@r)
axnames <- character(length=d)
for (i in 1:d) { axnames[i] <- paste(" Axis",i) } 

cat("\n    RESULTS for Correspondence Analysis:", datasetname, "\n\n")
cat("Total inertia ", x@inertiasum, "\n\n")
cat("Inertias, percent inertias and cumulative percent inertias \n\n")
ins <- data.frame(x@inertias)
names(ins) <- c("Inertia","%  ","Cum. %")
print(ins, digits=6)
cat("\nRows in principal coordinates\n\n")
printwithaxes(data.frame(x@Rowprinccoord[ ,1:d], row.names=x@rowlabels), axnames)
cat("\nRow contributions (per mil)\n\n")
printwithaxes(data.frame(round(x@RowCTR[ ,1:d]*1000), row.names=x@rowlabels), axnames)
cat("\nRow representations (per mil)\n\n")
printwithaxes(data.frame(round(x@RowREP[ ,1:d]*1000), row.names=x@rowlabels), axnames)
cat("\nColumns in principal coordinates\n\n")
printwithaxes(data.frame(x@Colprinccoord[ ,1:d], row.names=x@collabels), axnames)
cat("\nColumn contributions (per mil)\n\n")
printwithaxes(data.frame(round(x@ColCTR[ ,1:d]*1000), row.names=x@collabels), axnames)
cat("\nColumn representations (per mil)\n\n")
printwithaxes(data.frame(round(x@ColREP[ ,1:d]*1000), row.names=x@collabels), axnames)
if (x@nboots>0) {
cat("\n\n  Results for Bootstrapping\n\n")
cat(x@nboots, "bootstrap replications with", x@resampledistn, "resampling\n")
if (x@resampledistn=="multinomial" & x@multinomialtype!="whole") 
 cat(paste("  ", 
  switch(x@multinomialtype,rowsfixed="with row sums constant",columnsfixed="with column sums constant"),
  "\n") )
cat("\nEstimated variances and covariances\n\n")
cat("Rows\n\n")
print(allvarscovs(x,"rows"),digits=4)
cat("\nColumns\n\n")
print(allvarscovs(x,"columns"),digits=4)
cat("\n\n")
} 
}


#### Prints brief 2-d results in similar style to ca package, with standard deviations

summaryca <- function(x, datasetname="") {

#setGeneric("summary", function(x,...) standardGeneric("summary") )
#setMethod("summary", signature(x="cabootcrsresults"), 
# function(x, datasetname="") {

colnames <- character(length=9)
colnames <- c("  Axis 1","StDev","Rep","Ctr","  Axis 2","StDev","Rep","Ctr","Quality")
colnamesnosd <- character(length=7)
colnamesnosd <- c("  Axis 1","Rep","Ctr","  Axis 2","Rep","Ctr","Quality")

cat("\n    SUMMARY RESULTS for Correspondence Analysis:", datasetname, "\n\n")
cat("Total inertia ", x@inertiasum, "\n\n")
cat("Inertias, percent inertias and cumulative percent inertias \n\n")
ins <- data.frame(x@inertias)
names(ins) <- c("Inertia","%  ","Cum. %")
print(ins, digits=4)
cat("\n")
if (x@nboots>0) { 
cat("Princ coords, std devs; rep and ctr (per mil); 2-d rep (per mil)\n\n")
} else {
cat("Princ coords; rep and ctr (per mil); 2-d rep (per mil)\n\n")
}
cat("Rows: \n")
rop <- data.frame( 
 round(x@Rowprinccoord[,1]*1000)/1000, round(sqrt(x@RowVar[,1])*1000)/1000, 
 round(x@RowREP[ ,1]*1000), round(x@RowCTR[ ,1]*1000),
 round(x@Rowprinccoord[,2]*1000)/1000, round(sqrt(x@RowVar[,2])*1000)/1000,
 round(x@RowREP[ ,2]*1000), round(x@RowCTR[ ,2]*1000),
 round(rowSums(x@RowREP[ ,1:2]*1000)), row.names=x@rowlabels )
if (x@nboots==0) {
  rop <- rop[,c(1,3,4,5,7,8,9)]
  names(rop) <- colnamesnosd
} else {
  names(rop) <- colnames
}
print(rop,digits=3)
cat("\n")

cat("Columns: \n")
cop <- data.frame(  
 round(x@Colprinccoord[,1]*1000)/1000, round(sqrt(x@ColVar[,1])*1000)/1000, 
 round(x@ColREP[ ,1]*1000), round(x@ColCTR[ ,1]*1000),
 round(x@Colprinccoord[,2]*1000)/1000, round(sqrt(x@ColVar[,2])*1000)/1000,
 round(x@ColREP[ ,2]*1000), round(x@ColCTR[ ,2]*1000),
 round(rowSums(x@ColREP[ ,1:2]*1000)), row.names=x@collabels )
if (x@nboots==0) {
  cop <- cop[,c(1,3,4,5,7,8,9)]
  names(cop) <- colnamesnosd
} else {
  names(cop) <- colnames
}
print(cop,digits=3)
cat("\n")

}


#### Plot results with confidence regions

# plot results (should be overloaded plot for cabootcrsresults class)
# two plots are produced, in each plot one set of points (rows or columns) is regarded as
# the primary set and is plotted in principal coordinates with confidence regions shown:
#   - one plot shows confidence regions for rows in principal coordinates
#   - one plot shows confidence regions for columns in principal coordinates
# the other set of points (columns or rows) is regarded as the secondary set and the plotting
# depends on the choice of biplot or french-style plot: 
#  biplot - secondary points shown as directions in standard coordinates
#  french - secondary points shown in principal coordinates
# 
# option of which axis pair to plot, default 1 and 2
# only bring in groups file now, allow for option of groups file or not,  
#  hence can change groups file for different runs, default is different colour
#  for each primary point and secondary points in black.
# 
# plots look better if saved as .eps or .pdf, rather than viewed in R or as .jpg or .png
# 
# use groupings and grouplabels (files or data frames) if groupings is non-null,
#  otherwise use plotsymbolscolours
# use othersmonochrome unless set to anything other than a colour, e.g. NULL

plotca <- function(x, datasetname="", showrowlabels=TRUE, showcolumnlabels=TRUE,
  groupings=NULL, grouplabels=NULL, plotsymbolscolours=c(19,"alldifferent",18,"alldifferent"), 
  othersmonochrome="black", crpercent=95, plottype="biplot", showrowcrs=TRUE, showcolumncrs=TRUE,
  firstaxis=1, lastaxis=2, plotallpairs=FALSE, picsize=c(-1,1)) {

#setGeneric("plot", function(x,...) standardGeneric("plot") )
#setMethod("plot", signature(x="cabootcrsresults"), 
# function(x, datasetname="", showrowlabels=TRUE, showcolumnlabels=TRUE,
#   groupings=NULL, grouplabels=NULL, plotsymbolscolours=c(19,"alldifferent",18,"alldifferent"), 
#   othersmonochrome="black", crpercent=95, plottype="biplot", showrowcrs=TRUE, showcolumncrs=TRUE,
#   firstaxis=1, lastaxis=2, plotallpairs=FALSE, picsize=c(-1,1)) {


## internal function to plot a  single picture

plotonepic <- function(a1,a2,plottype,things,nthings,nvars,Thingcoord,Varcoord,SBvar,SBcov,twoS,
                       inertiapc,resampledistn,multinomialtype,
                       thinggroup,thinggrlab,vargroup,vargrlab,thinglabels,varlabels,showcrs,picsizex,picsizey) {

eps <- 1e-15
critchisq2 <- qchisq(0.01*crpercent,2)
critchisq1 <- qchisq(0.01*crpercent,1)
theta <- seq(0,2*pi,0.001)
ellipsecoords <- rbind(sin(theta),cos(theta))

# horrible bodge to convert thinggrlab[[3]] to lists, mixing numbers and chars, also for vargrlab[[3]] in french plot
thinggrlab3 <- as.list(thinggrlab[[3]])
thinggrlab3int <- !is.na(as.integer(thinggrlab3))
for (i in  1:max(thinggroup[,2])) { if (thinggrlab3int[[i]]) { thinggrlab3[[i]] <- as.integer(thinggrlab3[[i]]) } }
if (plottype=="french") { 
vargrlab3 <- as.list(vargrlab[[3]])
vargrlab3int <- !is.na(as.integer(vargrlab3))
for (i in  1:max(vargroup[,2])) { if (vargrlab3int[[i]]) { vargrlab3[[i]] <- as.integer(vargrlab3[[i]]) } }
}

dev.new()
plot(Thingcoord[1,a1], Thingcoord[1,a2], xlim=picsizex, ylim=picsizey, 
    xlab=paste("Axis ", a1, "    ", inertiapc[a1], "%", sep=""), 
    ylab=paste("Axis ", a2, "    ", inertiapc[a2], "%", sep=""),  
    asp=1, pch=thinggrlab3[[thinggroup[1,2]]], col=thinggrlab[[4]][thinggroup[1,2]] )
for (i in 2:nthings) { points(Thingcoord[i,a1], Thingcoord[i,a2],
    asp=1, pch=thinggrlab3[[thinggroup[i,2]]], col=thinggrlab[[4]][thinggroup[i,2]] ) } 
abline(h=0,v=0)

if (!all(thinggrlab[[2]]=="")) { 
  labnum <- as.integer(thinggrlab3)
  labchar <- as.character(thinggrlab3)
  legend("topleft",thinggrlab[[2]],pch=labnum,col=thinggrlab[[4]],text.col=thinggrlab[[4]])
  for (i in 1:max(thinggroup[,2])) { if (is.na(labnum[[i]])) { labchar[[i]]<-thinggrlab3[[i]] } else { labchar[[i]]<-NA } }
  legend("topleft",thinggrlab[[2]],pch=labchar,col=thinggrlab[[4]],text.col=thinggrlab[[4]])
}

if (plottype=="biplot") { 
  if ((x@nboots>0)&(any(showcrs==TRUE))) {
  title(paste("Confidence regions for biplot of", things, "\n \n", datasetname ))
  title(paste("\n", resampledistn, "resampling,", 
   switch(multinomialtype, whole="", rowsfixed="row sums fixed,", columnsfixed="column sums fixed,"),
   x@nboots, "resamples \n"), font.main=1 )
  } else {
  title(paste("Biplot of", things, "\n", datasetname ))
  }
  for (i in 1:nvars) { lines(c(0,Varcoord[i,a1]), c(0,Varcoord[i,a2]), col=vargrlab[[4]][vargroup[[2]][i]]) } 
  grat <- cbind(Varcoord[,a1]/picsizex[1],Varcoord[,a1]/picsizex[2],Varcoord[,a2]/picsizey[1],Varcoord[,a2]/picsizey[2],0.95)/0.95
  cl <- 1.05/apply(grat,1,max)
  text(cl*Varcoord[ ,a1], cl*Varcoord[ ,a2], labels=varlabels, col=vargrlab[[4]][vargroup[[2]]], pos=4, cex=0.75 )
} else { # french
  if ((x@nboots>0)&(any(showcrs==TRUE))) {
  title(paste("Confidence regions for", things, "\n \n", datasetname ) )
  title(paste("\n", resampledistn, "resampling,", 
   switch(multinomialtype, whole="", rowsfixed="row sums fixed,", columnsfixed="column sums fixed,"),
   x@nboots, "resamples \n"), font.main=1 )
  } else {
  title(paste("Correspondence plot \n", datasetname ))
  }
  for (i in 1:nvars) { points(Varcoord[i,a1], Varcoord[i,a2], asp=1, 
         pch=vargrlab3[[vargroup[i,2]]], col=vargrlab[[4]][vargroup[i,2]] ) }
  text(Varcoord[ ,a1], Varcoord[ ,a2], labels=varlabels, col=vargrlab[[4]][vargroup[[2]]], pos=4, cex=0.75 )
  if (!all(vargrlab[[2]]=="")) {
    labnum <- as.integer(vargrlab3)
    labchar <- as.character(vargrlab3)
    legend("topright",vargrlab[[2]],pch=labnum,col=vargrlab[[4]],text.col=vargrlab[[4]])
    for (i in 1:max(vargroup[,2])) { if (is.na(labnum[[i]])) { labchar[[i]]<-vargrlab3[[i]] } else { labchar[[i]]<-NA } }
    legend("topright",vargrlab[[2]],pch=labchar,col=vargrlab[[4]],text.col=vargrlab[[4]])
  }  
}

for (i in 1:nthings) {
  if (thinggrlab[[5]][thinggroup[i,2]]) {
    text(Thingcoord[i,a1], Thingcoord[i,a2], labels=thinglabels[i], pos=4, cex=0.75, col=thinggrlab[[4]][thinggroup[i,2]] )
    if (showcrs[i]) {
      xbar <- Thingcoord[i,cbind(a1,a2)]
      V <- matrix(c(SBvar[i,a1],SBcov[i,min(a1,a2),max(a1,a2)],SBcov[i,min(a1,a2),max(a1,a2)],SBvar[i,a2]),2,2)
      E <- eigen(V, symmetric=TRUE)
      usec2 <- (1-twoS[i]) * (E$values[1]>eps)
      critchisq <- critchisq2 * usec2 + critchisq1 * (1-usec2)
      coords <- E$vectors %*% (critchisq*diag(E$values))^(1/2) %*% ellipsecoords 
      lines(xbar[1]+coords[1, ], xbar[2]+coords[2, ], pch=".", col=thinggrlab[[4]][thinggroup[i,2]] )
      # equivalent if using the ellipse package and chi-squared with 2 df: 
      # lines(ellipse(x=V, centre=xbar, npoints=1000), cex=1, pch=".", col=thinggrlab[[4]][thinggroup[i,2]] )
    }  
  }
}

if (any(Thingcoord[,a1]<picsizex[1])) {
 cat(paste("Warning: point outside plot limits, lowest x-value is ", min(Thingcoord[,a1]), "\n")) }
if (any(Thingcoord[,a1]>picsizex[2])) {
 cat(paste("Warning: point outside plot limits, largest x-value is ", max(Thingcoord[,a1]), "\n")) }
if (any(Thingcoord[,a2]<picsizey[1])) {
 cat(paste("Warning: point outside plot limits, lowest y-value is ", min(Thingcoord[,a2]), "\n")) }
if (any(Thingcoord[,a2]>picsizey[2])) {
 cat(paste("Warning: point outside plot limits, largest y-value is ", max(Thingcoord[,a2]), "\n")) }

} # plot one pic

##

if (!is.null(plotsymbolscolours)) { 
 if(!dim(array(plotsymbolscolours))==4) stop(paste("plotsymbolscolours must contain row symbol and colour, column symbol and colour\n\n"))
}
if (!any(plotsymbolscolours[c(2,4)]==c(colours(),"alldifferent","differentblues","differentreds"))) 
         stop(paste("colours must be alldifferent, differentblues, differentreds or R colour (type colours() for full list) \n\n"))
if ((crpercent<=0)|(crpercent>=100)) stop(paste("coverage percentage must be between 0 and 100 exclusive\n\n"))
if (!any(plottype==c("biplot","french"))) stop(paste("plotting must be biplot or french style\n\n"))
if (!any(class(showrowcrs)==c("integer","numeric","logical"))) stop(paste("showrowcrs must be logical or a vector of row numbers\n\n"))
if (!any(class(showcolumncrs)==c("integer","numeric","logical"))) stop(paste("showcolumncrs must be logical or a vector of row numbers\n\n"))
if ((firstaxis<1)|(firstaxis>x@axisvariances-1)) stop(paste("incorrect first axis =", firstaxis, "\n\n")) 
if (lastaxis>x@axisvariances) stop(paste("don't have variances for last axis =", lastaxis, "\n\n")) 
if (firstaxis>=lastaxis) stop(paste("last axis must be greater than first axis\n\n"))
if (!any(dim(array(picsize))==c(2,4))) stop(paste("picsize bounds are  lower,upper OR x lower,x upper,y lower,y upper \n\n"))
if (picsize[1]>=picsize[2]) stop(paste("incorrect axis scale picsize =", picsize[1], picsize[2], "\n\n"))
if (dim(array(picsize))==4) { 
 if (picsize[3]>=picsize[4]) stop(paste("incorrect y axis scale picsize =", picsize[3], picsize[4], "\n\n")) 
 if (abs((picsize[4]-picsize[3])-(picsize[2]-picsize[1]))>1e-10) stop(paste("x and y axes must be same length\n\n"))
}
options(warn=-1)
picsizey <- picsizex <- picsize[1:2]
if (dim(array(picsize))==4) picsizey <- picsize[3:4]

tworowS <- rowSums(x@DataMatrix>0)==2
twocolS <- colSums(x@DataMatrix>0)==2

# Groupings file or data frame contains 
#  row no.  row group
#  col no.  col group
# Group label file or data frame contains
#  row group number   group name   symbol   colour   plot ellipse? 
#  col group number   group name   symbol   colour   plot ellipse? 


if (is.null(groupings)) { # if no groupings of points supplied
  if (any(plotsymbolscolours[2]==c("alldifferent","differentreds","differentblues"))) { # row points different colours
    rowgroup <- as.data.frame(cbind(1:x@rows,1:x@rows))
    hv <- switch(plotsymbolscolours[2], "alldifferent"=c(0,0.85), "differentreds"=c(0,0.45), "differentblues"=c(0.5,0.85))
    rowgrlab <- as.data.frame(cbind(1:x@rows,"",plotsymbolscolours[1],rainbow(n=x@rows,start=hv[1],end=hv[2]),"T"),stringsAsFactors=FALSE)
  } else { # all same colour
    rowgroup <- as.data.frame(cbind(1:x@rows,rep(1,x@rows)))
    rowgrlab <- as.data.frame(cbind(1,"",plotsymbolscolours[1],plotsymbolscolours[2],"T"),stringsAsFactors=FALSE)
  } 
  class(rowgrlab[,1]) <- "integer"
  class(rowgrlab[,5]) <- "logical"
  if (any(plotsymbolscolours[4]==c("alldifferent","differentreds","differentblues"))) { # column points different colours
    colgroup <- as.data.frame(cbind(1:x@columns,1:x@columns))
    hv <- switch(plotsymbolscolours[4], "alldifferent"=c(0,0.85), "differentreds"=c(0,0.45), "differentblues"=c(0.5,0.85))
    colgrlab <- as.data.frame(cbind(1:x@columns,"",plotsymbolscolours[3],rainbow(n=x@columns,start=hv[1],end=hv[2]),"T"),stringsAsFactors=FALSE)
  } else { # all same colour
    colgroup <- as.data.frame(cbind(1:x@columns,rep(1,x@columns)))
    colgrlab <- as.data.frame(cbind(1,"",plotsymbolscolours[3],plotsymbolscolours[4],"T"),stringsAsFactors=FALSE)
  }
  class(colgrlab[,1]) <- "integer"
  class(colgrlab[,5]) <- "logical"
} else { # groupings in file  
  if (class(groupings)=="character") {
    rcgroup <- read.table(file=groupings,colClasses=c("integer","integer"))
  } else { # groupings in data frame
    rcgroup <- as.data.frame(groupings)
  } 
  rowgroup <- rcgroup[1:x@rows,]
  colgroup <- rcgroup[(x@rows+1):(x@rows+x@columns),]
  nrowgroups <- max(rowgroup[,2])
  ncolgroups <- max(colgroup[,2])
  if (class(grouplabels)=="character") { # group labels in file
    rcgrlab <- read.table(file=grouplabels,
                          colClasses=c("integer","character","character","character","logical"))
  } else { # group labels in data frame
    rcgrlab <- as.data.frame(grouplabels,stringsAsFactors=FALSE)
    class(rcgrlab[,1]) <- "integer"
    class(rcgrlab[,5]) <- "logical"
  }
  rowgrlab <- rcgrlab[1:nrowgroups,]
  colgrlab <- rcgrlab[(nrowgroups+1):(nrowgroups+ncolgroups),]
}

# quick option to plot only a few regions, overrides other options
rowcrs <- logical(length=x@rows)
columncrs <- logical(length=x@columns)
if (any(class(showrowcrs)==c("numeric","integer"))) { for (i in 1:length(showrowcrs)) { rowcrs[showrowcrs[i]]<-TRUE }
                           } else { rowcrs <- rowcrs | showrowcrs }
if (any(class(showcolumncrs)==c("numeric","integer"))) { for (i in 1:length(showcolumncrs)) { columncrs[showcolumncrs[i]]<-TRUE }
                              } else { columncrs <- columncrs | showcolumncrs }

# option for secondary set of points, with CRs not shown, to be monochrome
vrowgrlab <- rowgrlab
vcolgrlab <- colgrlab
if (any(othersmonochrome==colours())) {
  vrowgrlab[[4]] <- othersmonochrome
  vcolgrlab[[4]] <- othersmonochrome
}

# Plot row and col pictures for pairs of axes

if (showrowlabels==TRUE) { rowptlabels <- x@rowlabels } else { rowptlabels <- NULL }
if (showcolumnlabels==TRUE) { colptlabels <- x@collabels } else { colptlabels <- NULL }

for (a1 in firstaxis:(lastaxis-1)) {
  for (a2 in (a1+1):lastaxis) { 
    if ( (plotallpairs==TRUE) | ((a1==firstaxis)&(a2==lastaxis)) ) { 
if (plottype=="biplot") {
plotonepic(a1, a2, plottype, "rows", x@rows, x@columns, x@Rowprinccoord, x@Colstdcoord, x@RowVar, x@RowCov, tworowS, 
            x@inertias[,2], x@resampledistn, x@multinomialtype, 
            rowgroup, rowgrlab, colgroup, vcolgrlab, rowptlabels, colptlabels, rowcrs, picsizex, picsizey) 
plotonepic(a1, a2, plottype, "columns", x@columns, x@rows, x@Colprinccoord, x@Rowstdcoord, x@ColVar, x@ColCov, twocolS, 
            x@inertias[,2], x@resampledistn, x@multinomialtype, 
            colgroup, colgrlab, rowgroup, vrowgrlab, colptlabels, rowptlabels, columncrs, picsizex, picsizey) 
} else { # french
plotonepic(a1, a2, plottype, "rows", x@rows, x@columns, x@Rowprinccoord, x@Colprinccoord, x@RowVar, x@RowCov, tworowS, 
            x@inertias[,2], x@resampledistn, x@multinomialtype, 
            rowgroup, rowgrlab, colgroup, vcolgrlab, rowptlabels, colptlabels, rowcrs, picsizex, picsizey) 
plotonepic(a1, a2, plottype, "columns", x@columns, x@rows, x@Colprinccoord, x@Rowprinccoord, x@ColVar, x@ColCov, twocolS, 
            x@inertias[,2], x@resampledistn, x@multinomialtype, 
            colgroup, colgrlab, rowgroup, vrowgrlab, colptlabels, rowptlabels, columncrs, picsizex, picsizey)
} # plot type

} } } # pairs of pictures

options(warn=0)

} 


#### Extract a single 2 by 2 covariance matrix

covmat <- function(x, i, thing="row", axis1=1, axis2=2) {

## Printing macro
printwithaxes <- function(res, thenames) { 
names(res) <- thenames
print(res, digits=4)
}
## 

if (!(class(x)=="cabootcrsresults")) stop(paste("Must be of type cabootcrsresults\n\n"))
if (!any(thing==c("row","column"))) stop(paste("Must be row or column\n\n"))
if (axis1==axis2) stop(paste("What are you playing at?\n\n"))
if (!any(axis1==seq(1,x@axisvariances))) stop(paste("Covariance not available for these axes\n\n")) 
if (!any(axis2==seq(1,x@axisvariances))) stop(paste("Covariance not available for these axes\n\n"))
if ((thing=="row") & !any(i==seq(1,x@rows))) stop(paste("Invalid row number\n\n"))
if ((thing=="column") & !any(i==seq(1,x@columns))) stop(paste("Invalid column number\n\n"))

a1 <- min(axis1,axis2)
a2 <- max(axis1,axis2)
tname <- ""
if (thing=="row") {
V <- matrix(c(x@RowVar[i,axis1],x@RowCov[i,a1,a2],x@RowCov[i,a1,a2],x@RowVar[i,axis2]),2,2)
if (!is.null(x@rowlabels)) { tname <- paste("(",x@rowlabels[[i]],")") }
} else { # column
V <- matrix(c(x@ColVar[i,axis1],x@ColCov[i,a1,a2],x@ColCov[i,a1,a2],x@ColVar[i,axis2]),2,2)
if (!is.null(x@collabels)) { tname <- paste("(",x@collabels[[i]],")") }
}

cat(paste("Covariance matrix of", switch(thing,"row"="row","column"="column"), i, tname, "for axes", axis1,axis2,"\n\n"))
rcnames <- c(paste("Axis",axis1),paste("Axis",axis2))
printwithaxes(data.frame(V,row.names=rcnames),rcnames)

invisible(V)

}



#### Extract all vars and covs in readable form as a data frame

allvarscovs <- function(x, thing="rows") { 

## function to extract upper triangle

getcovs <- function(allC,n,ncovs) {
V <- matrix(0,n,ncovs)
for (i in 1:n) {
y <- allC[i,,]
V[i,] <- y[upper.tri(y)]
}
invisible(V)
} # getcovs

## 

if (!(class(x)=="cabootcrsresults")) stop(paste("Must be of type cabootcrsresults\n\n"))
if (!any(thing==c("rows","columns"))) stop(paste("Must be rows or columns\n\n"))

ncovs <- x@axisvariances*(x@axisvariances-1)/2
vcnames <- character(length=x@axisvariances+ncovs)
k <- 1
for (i in 1:x@axisvariances) { 
  vcnames[i] <- paste(" Var Axis",i) 
  if (i<x@axisvariances) {
    for (j in (i+1):x@axisvariances) {
      vcnames[x@axisvariances+k] <- paste(" Cov axes",i,j)
      k <- k+1
} } } 

if (thing=="rows") { 
Covs <- getcovs(x@RowCov,x@rows,ncovs)
allV <- data.frame(cbind(x@RowVar,Covs), row.names=x@rowlabels)
} else { # columns
Covs <- getcovs(x@ColCov,x@columns,ncovs)
allV <- data.frame(cbind(x@ColVar,Covs), row.names=x@collabels)
}

names(allV) <- vcnames

allV

}




