  ###################################################################
  #
  # This function is part of WACSgen V1.0
  # Copyright © 2013,2014,2015, D. Allard, BioSP,
  # and Ronan Trépos MIA-T, INRA
  #
  # This program is free software; you can redistribute it and/or
  # modify it under the terms of the GNU General Public License
  # as published by the Free Software Foundation; either version 2
  # of the License, or (at your option) any later version.
  #
  # This program is distributed in the hope that it will be useful,
  # but WITHOUT ANY WARRANTY; without even the implied warranty of
  # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  # GNU General Public License for more details. http://www.gnu.org
  #
  ###################################################################


  wacs.estimWT = function(yData,season,Vsel=NULL,Nclusters=NULL,plot.it=TRUE,DIR="./"){
  ###################################################################
  #
  #
  # wacs.estimWT     : estimates the Weather Type clustering for one given season
  #
  # ARGUMENTS
  #    yData     : contains the residuals of selected data;
  #                colnames(yData) = "NS.rain" "resid.V1" ... "resid.Vp"
  #    season    : season (used if plot.it)
  #    Vsel      : selection of variables (other than rain) used for clustering
  #    Nclusters : number of clusters that will be considered
  #    plot.it   : whether plots are produced or not
  #    DIR       : directory in which plots are printed
  #            
  # VALUE
  #   A list. First  element is the number of dray and wet types
  #           Second element is the vector of cluster labels
  #           Third  element is the matrix of cluster probabilities 
  ###################################################################
  
  Nd   = dim(yData)[1]
  Nv   = dim(yData)[2] - 1
 
  # some checks

  if (length(Nclusters)==0){
    Nclusters = 1:3
  }else{}
 
  FIND = TRUE
  if (length(Nclusters)==1){
        if (Nclusters==1){ 
          FIND = FALSE
        }else{}
  }else{}
  
# Building the array with variables on which clustering is performed  
  select=NULL
  if(length(Vsel)==0){ 
    select = rep(TRUE,Nv) # No selection => take all
  }else{}
  
  if(length(Vsel)>0){ # some selection
    if (max(Vsel) > Nv) stop ("[estimWT] variables selected out of range")
    for (v in 1:Nv){
      Vselect = FALSE
      if (any(Vsel==v)) Vselect = TRUE
      select = c(select,Vselect)
    }  
  }else{}
  V  = yData[,-1]
  V  = V[,select]
  Nv = dim(V)[2] - 1
  
  # Classification on dry days  
  V.dry=V[yData[,1] == -999,]
 
  EMdry       = Mclust(V.dry,G=Nclusters,modelNames="VVV",noise=FALSE)
  Nwt.dry     = EMdry$G
  clust.dry   = EMdry$classification
 
  if (plot.it){
    pdf(paste(DIR,"mclust_",season,"_dry.pdf",sep=""))  
    if (FIND){
      plot(EMdry,what="BIC")
      plot(EMdry,what="classification")
      toto = densityMclust(V.dry,G=Nwt.dry)
      plotDensityMclustd(toto,data=V.dry,points.col="green",points.pch=".")
      dev.off();rm(toto)
  }else{}
  }
  
  
  # Classification on wet days
  V =cbind(yData[,1],V)
  colnames(V)[1] = "NSrain"
  V.wet=V[V[,1] > -999,]
  
  EMwet     = Mclust(V.wet,G=Nclusters,modelNames="VVV",noise=FALSE)
  Nwt.wet   = EMwet$G
  clust.wet = Nwt.dry + EMwet$classification
  
  if (plot.it){
    pdf(paste(DIR,"mclust_",season,"_wet.pdf",sep=""))
    if (FIND){ 
      plot(EMwet,what="BIC")
      plot(EMwet,what="classification")
      toto = densityMclust(V.wet,G=Nwt.wet)
      plotDensityMclustd(toto,data=V.wet,points.col="blue",points.pch=".")
      dev.off();rm(toto)
  }else{}
  }
  
  # Merging both classifications
  Nwt        =  Nwt.dry + Nwt.wet 
  clust.all  =  NULL
  if (FIND){
    WT.z      = array(0,c(Nd,Nwt))
    clust.all = rep(1,dim(yData)[1])
    clust.all[ ( V[,1] == -999 ) ] = clust.dry
    clust.all[ ( V[,1] >  -999 ) ] = clust.wet  
    if (Nwt.dry >1){
      WT.z[ ( V[,1] == -999 ), 1:Nwt.dry ] = EMdry$z
    }
    if (Nwt.dry==1){
      WT.z[ ( V[,1] == -999 ), 1:Nwt.dry ] = rep(1,sum(V[,1] == -999))
    }
    if (Nwt.wet > 1){
      WT.z[ ( V[,1]  > -999 ), (Nwt.dry+1):(Nwt.dry+Nwt.wet) ] = EMwet$z
    }
    if (Nwt.wet==1){
      WT.z[ ( V[,1]  > -999 ), (Nwt.dry+1):(Nwt.dry+Nwt.wet) ] = rep(1,sum(V[,1]  > -999))
    }
  }
  else{
    WT.z  = array(0,c(Nd,2))
    WT.z[ ( V[,1] == -999 ), 1 ]   = 1
    WT.z[ ( V[,1]  > -999 ), 2 ]   = 1
    clust.all[ ( V[,1] == -999 ) ] = 1
    clust.all[ ( V[,1] >  -999 ) ] = 2
  }
  
  list.return        =  list(c(Nwt.dry,Nwt.wet),WT.z,clust.all)
  names(list.return) =  c("NumbWT","Prob","WT")
  return(list.return)
}
  
  
  
