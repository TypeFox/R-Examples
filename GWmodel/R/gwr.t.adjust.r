#The spgwr function gwr returns an object of class gwr which contains the parameter estimate and their standard errors. gwr.t.adjust returns the t values for the parameter
#estimates, unadjusted p values, and p values adjusted using (a) Bonferroni, (b) Benjamini-Hochberg, (c) Benjamini-Yekutieli and
#(d) Fotheringham-Byrne procedures.
#Author: MC
#Edited by BL

gwr.t.adjust <- function(gwm.Obj)
{
  hatmatrix <-gwm.Obj$GW.arguments$hatmatrix
  if(!hatmatrix)
     stop("No p-values to be adjusted")
  gwmx <- as.data.frame(gwm.Obj$SDF)
  #colnames(gmx)
  n <- dim(gwmx)[1]
  m <- dim(gwmx)[2]
  vnames<-all.vars(gwm.Obj$GW.arguments$formula)
  if(length(vnames)==length(names(gwm.Obj$lm$coefficients)))
     vnames[1]<-"Intercept"
  nv <- length(vnames)
  np <- nv
  ntests <- n * np
  enp <- gwm.Obj$GW.diagnostic$enp
  SDFnms <- names(gwmx)
  idx <- c()
  for(i in 1:m)
  {
    if(grepl("_TV", SDFnms[i]))
      idx <- c(idx, i)
  }
  tvals <- as.matrix(gwmx[, idx])
  pvals <- round(2 * (1 - pt(abs(tvals), ntests)), 3)

  bey_pvals <- round(p.adjust(pvals, "BY", n = ntests))
  beh_pvals <- round(p.adjust(pvals, "BH", n = ntests))
  bon_pvals <- round(p.adjust(pvals, "bonferroni", n = ntests))
  dim(bey_pvals) <- c(n,nv)
  dim(beh_pvals) <- c(n,nv)
  dim(bon_pvals) <- c(n,nv)
  #print(bey_pvals)
  colnames(tvals) <- paste(vnames, "_t", sep = "")
  colnames(pvals) <- paste(vnames, "_p", sep = "")
  colnames(bey_pvals) <- paste(colnames(pvals), "_by", sep = "")
  colnames(beh_pvals) <- paste(colnames(pvals), "_bh", sep = "")
  colnames(bon_pvals) <- paste(colnames(pvals), "_bo", sep = "")
  asf_pvals <- round(pvals * (1 + enp - (enp/ntests)), 3)
  asf_pvals[asf_pvals > 1] <- 1
  colnames(asf_pvals) <- paste(colnames(pvals), "_fb", sep = "")
  results <- list(t = tvals, p = pvals, by = bey_pvals, fb = asf_pvals,
      bo = bon_pvals, bh = beh_pvals)
  df.res<-data.frame(tvals, pvals, bey_pvals, asf_pvals,bon_pvals, beh_pvals)
  p4s <- proj4string(gwm.Obj$SDF)
  if(is(gwm.Obj$SDF, "SpatialPolygonsDataFrame"))
     polygons<-polygons(gwm.Obj$SDF)
  else
  {
     locat <- coordinates(gwm.Obj$SDF)
     rownames(locat)<-rownames(df.res)
  }
  griddedObj <- F
  if (is(gwm.Obj$SDF, "Spatial"))
  { 
      if (is(gwm.Obj$SDF, "SpatialPolygonsDataFrame"))
      {
         polygons<-polygons(gwm.Obj$SDF)
         #SpatialPolygons(regression.points)
         #rownames(gwres.df) <- sapply(slot(polygons, "polygons"),
                            #  function(i) slot(i, "ID"))
         SDF <-SpatialPolygonsDataFrame(Sr=polygons, data=df.res)
      }
      else
      {
         griddedObj <- gridded(gwm.Obj$SDF)
         SDF <- SpatialPointsDataFrame(coords=locat, data=df.res, proj4string=CRS(p4s), match.ID=F)
         gridded(SDF) <- griddedObj 
      }
  }
  else
      SDF <- SpatialPointsDataFrame(coords=locat, data=df.res, proj4string=CRS(p4s), match.ID=F)
   
#  if (is(gwm.Obj$SDF, "SpatialPolygonsDataFrame"))
#  {
#     SDF <-SpatialPolygonsDataFrame(Sr=polygons, data=df.res)
#  }
#  else
#     SDF <- SpatialPointsDataFrame(coords=locat, data=df.res, proj4string=CRS(p4s))
  
  res<-list(results=results, SDF=SDF)
  res
}