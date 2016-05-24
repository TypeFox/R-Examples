#-----------------------------------------------------------------------
#     Copyright (C) 2012-2014  Serge Iovleff, University Lille 1, Inria
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as
#    published by the Free Software Foundation; either version 2 of the
#    License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public
#    License along with this program; if not, write to the
#    Free Software Foundation, Inc.,
#    59 Temple Place,
#    Suite 330,
#    Boston, MA 02111-1307
#    USA
#
#    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
#
#-----------------------------------------------------------------------
# Plotting of a class  deriving from[\code{\linkS4class{IClusterModelBase}}]
##############################
# Adapted from Rmixmod package
##############################
# x is a class deriving from IClusterModel
# y is a list of variable
# ddensity : the density to plot on the histograms
.clusterPlot <- function(model, y, ddensity,...)
{
  # total number of variable in the data set
  nbVariable = ncol(model@component@data);
  # no y => display all variables
  if (missing(y)) { y=1:nbVariable; }
  else # perform some check
  {
    if (is.numeric(y)) # numbers of the columns to plot are given
    {
      if (max(y)>nbVariable)
        stop("In .clusterPlot, y indices mismatch the data dimension")
    }
    else # names of the variables to plot are given
    {
      if ( sum(y %in% colnames(model@component@data))!= length(y) )
      { stop(cat("In .clusterPlot, unknown variable: ", paste(y[which(!(y %in% colnames(model@component@data)))]),"\n"))}
    }
  }
  # get old par
  op <- par(no.readonly = TRUE) # the whole list of settable par's.
  # cluster parameters
  par(mar = rep(2.5,4), cex = .75, oma = c(0, 0, 3, 0))        # margin and font size
  nbCol = length(y)                       # size of the matrix screen
  split.screen(c(nbCol, nbCol))           # create layout matrix screens
  col = model@zi+2;                       # color for each group
  pch = rep(1, length.out = length(col)); # circles
  pch[model@component@missing[,1]] = 3;   # + for missing values
  # create histograms on the diagonal
  for ( i in 1:nbCol )
  {
    screen(i+((i-1)*nbCol))   # sreen(i,i)
    xValues<-seq( min(model@component@data[,y[i]]), max(model@component@data[,y[i]]), length.out = 200)
    density<-matrix(nrow=model@nbCluster, ncol=length(xValues))
    # loop over the clusters to generate densities
    for( k in 1:model@nbCluster )
    {  density[k,]<- ddensity(xValues, y[i], k, model);}
    # generate mixture density
    mixture<-apply(density,2,sum)
    if (is.numeric(y)) { xlab=colnames(model@component@data)[y[i]];}
    else               { xlab= y[i];}
    # TODO: check if xlab is empty
    if (is.null(xlab)) { xlab = paste("dimension ", i)}
    main=paste("Histogram of",xlab)
    h<-hist(model@component@data[,y[i]], xlab=xlab, main=main, ...)
    # add on the histogram the estimated densities
    ratio<-max(h$counts)/max(mixture)
    density<-density*ratio
    mixture<-mixture*ratio
    lines(xValues,mixture,col="azure4", lty=1, lwd=4)
    for( k in 1:model@nbCluster )
    { lines(xValues, density[k,], col=k+1, lty=2, lwd=2)}
  }
  # add biplots
  if (nbCol>1)
  {
    for ( i in 2:nbCol )
    {
      if (is.numeric(y)) { xlab=colnames(model@component@data)[y[i]];}
      else               { xlab= y[i];}
      if (is.null(xlab)) { xlab = paste("dimension ", i)}
      main=paste("Histogram of",xlab)
      for( j in 1:(i-1) )
      {
        screen(j+((i-1)*nbCol)) # screen(i,j)
        if (is.numeric(y)) {ylab=colnames(model@component@data)[y[j]];}
        else {ylab= y[j];}
        if (is.null(ylab)) { ylab = paste("dimension ", j)}
        main=paste("Histogram of",xlab)
        plot(model@component@data[,y[j]], model@component@data[,y[i]], col=col, pch=pch, xlab=xlab, ylab=ylab, ...)
      }
    }
  }
  #  mtext("Visualisation using latent logistic representation", outer = TRUE, cex = 1.5)
  close.screen(all.screens = TRUE)
  # restore plotting parameters
  par(op)
  invisible()
}
