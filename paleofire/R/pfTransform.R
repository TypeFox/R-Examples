#' Transform charcoal data for unique to multiple series
#' 
#' Charcoal data transformation, background estimation and homogenization for
#' unique to multiple series, accepts objects returned by
#' \code{\link{pfSiteSel}}.
#' 
#' 
#' @param ID An object returned by \code{\link{pfSiteSel}} or
#' \code{\link{pfTransform}}
#' @param add An object returned by \code{\link{pfAddData}}
#' @param Interpolate Logical, idicates wether data should be interpolated or
#' not, default=FALSE
#' @param Age Numeric, If Interpolate=TRUE, Age is used to specified the ages
#' where the interpolation took place, If Age=NULL (default) the interpolated
#' ages are automatically specified using the median resolution of the
#' record(s). If Age is specified as a vector (e.g. Age=(from=0,to=10000,
#' by=10)) the interpolation took place at specified ages.
#' @param method A character indicating the transformation method: "Z-Score",
#' Z-Score, "LOESS", Locally weighted regression, "SmoothSpline", Smoothing
#' spline, "Box-Cox", Box-Cox transformation, "MinMax", Minimax transformation,
#' "RunMed", Running median, "RunMean", Running mean, "RunQuantile", Running
#' quantile, "RunMin", Running min, "RunMax", Running max, "stl", Decompose a
#' time series into seasonal, trend and irregular components using loess, based
#' on \code{\link[stats]{stl}} function.
#' @param BasePeriod Numeric, a parameter specifying the base period for
#' calculating Z-score given in years BP (e.g. BasePeriod=c(0, 4000)), if empty
#' or unspecified the base period corresponds to record length.
#' @param span Numeric, the span parameter for the LOESS or Smoothing spline
#' methods
#' @param RunWidth Numeric, the width of the window for the"RunMed", "RunMean",
#' "RunQuantile", "RunMin", and "RunMax" methods in years.
#' @param RunQParam Numeric, the parameter specifying which quantile should be
#' calculated for the method "RunQuantile" (default=0.5 i.e. median).
#' @param stlYears Numeric, the bandwidth for stl decomposition, default=500
#' years.
#' @param alpha Numeric, alpha value to add before BoxCox calculation, see
#' \code{\link{pfBoxCox}}.
#' @param type Character, the type of Box-Cox transformation, see
#' \code{\link{pfBoxCox}} for details.
#' @param QuantType Character, by default QuantType="INFL" and influx are
#' automatically calculated, otherwise use QuantType="NONE" (not recommended).
#' @param MethodType Character, by default (MethodType=NULL) imply that when
#' for a specific site two charcoal unit exist the function pick the one define
#' by pref_unit. By passing different arguments to MethodType user can modify
#' the analysis to pick non preferred units by referring to more general
#' methods for instance MethodType = "POLS" will choose charcoal records from
#' pollen slides, or MethodType = "SIEV" sieved macro charcoal series. Type
#' (paleofiredata); levels(paleofiredata$METHOD) for available methods.
#' @param verbose Logical, verbose or not...
#' @return An object of the class "pfTransform".
#' @author O. Blarquez
#' @examples
#' 
#' ## Select the site Pas-de-Fond
#' ID=pfSiteSel(site_name=="Pas-de-Fond")
#' 
#' # Transform data sequentially using pfTransform function
#' tr=pfTransform(ID,method=c("MinMax","Box-Cox"))
#' 
#' ## Plot transformed data for the first site
#' plot(tr$Age[,1],tr$TransData[,1],type="l")
#' 
#' 
pfTransform=function(ID=NULL,
                     add=NULL,
                     Interpolate=FALSE,
                     Age=NULL,
                     method="Z-Score",
                     BasePeriod=c(-100,1e+09),
                     span=0.3,
                     RunWidth=500,
                     RunQParam=0.5,
                     stlYears=500,
                     type="BoxCox1964",
                     alpha=0.01,
                     QuantType="INFL",
                     MethodType=NULL,
                     verbose=TRUE
){
  
  
  ## TEST
#   ID=NULL;
#   add=NULL;
#   Interpolate=FALSE;
#   Age=NULL;
#   method="NULL";
#   BasePeriod=c(-100,1e+09);
#   span=0.3;
#   RunWidth=500;
#   RunQParam=0.5;
#   stlYears=500;
#   type="BoxCox1964";
#   alpha=0.01;
#   QuantType="INFL";
#   MethodType=NULL;
#   verbose=TRUE
  ## TEST
  
  
  ## Avoid no visible binding for global variable
  paleofiresites=NULL; rm(paleofiresites)
  
  # Value for warnings
  IDChar=ID
  
  # Check methods
  methods=c("stl", "Z-Score", "Box-Cox", "LOESS", "MinMax", "RunMed", "RunMean", "RunMin", "RunMax", "RunQuantile", "SmoothSpline", "Hurdle", 'NULL')
  warnmethod=method[(method %in% methods)==FALSE]
  if(length(warnmethod)!=0){stop(paste(warnmethod, "is not a valid method for pfTransform", sep=" "))}
  
  types=c("BoxCox1964", "JohnDraper")
  warntype=type[(type %in% types)==FALSE]
  if(length(warntype)!=0){stop(paste(warntype, "is not a valid type for pfBoxCox", sep=" "))}
  
  if(method=="RunMean" || method=="RunMin" || method=="RunMed" || method=="RunMax" || 
     method=="RunQuantile") {
    if("caTools" %in% rownames(installed.packages()) == FALSE) {install.packages("caTools")}
    if("gtools" %in% rownames(installed.packages()) == FALSE) {install.packages("gtools")}
  }
  
  if(identical(method, "Hurdle") ) {install.packages("pscl")}
  
  ## 0 Save parameters
  params=list(ID=ID,
              Interpolate=Interpolate,
              Age=Age,
              method=method,
              BasePeriod=BasePeriod,
              span=span,
              RunWidth=RunWidth,
              RunQParam=RunQParam,
              stlYears=stlYears,
              type=type,
              alpha=alpha)
  
  ## 1 Load charcoal paleofiredata
  if(verbose==TRUE){
    cat("Loading and preparing data...")
    cat("\n")}
  
  if(is.null(ID)==FALSE){
    if (is.list(ID) & length(ID)==2){
      
      data(paleofiredata,envir = environment())
      data(paleofiresites,envir = environment())
      
      #paleofiredata=na.omit(paleofiredata)
      # Sites are:
      ID=ID$id_site
      # Use only paleofiredata corresponding to ID
      paleofiredata=paleofiredata[paleofiredata[,1] %in% ID,]
      
      ## 1 Use Pref_Units 
      if(is.null(MethodType)){
        # Drop Non pref Units  
        for(i in ID){
          paleofiredata[paleofiredata[,1] %in% i &
                          !(paleofiredata[,5] %in% paleofiresites$pref_units[paleofiresites[,1]==i]),
                        7]=NA
        }
        paleofiredata=paleofiredata[!is.na(paleofiredata$TYPE),]
        ## Convert data to influx------
        if(QuantType=="INFL"){
          for(i in ID)  
            if( !(unique(paleofiredata[paleofiredata[,1]==i,7]) %in% "INFL") & 
                is.na(sum(paleofiredata[paleofiredata[,1]==i,2]))==FALSE){
              infl=influx(paleofiredata[paleofiredata[,1]==i,])
              paleofiredata[paleofiredata[,1]==i,4]=c(infl)
            }
        }
      } else {
        ## 2 User defined Method
        ## Drop duplicate units and keep only pref_unit in desired method
        paleofiredata=paleofiredata[paleofiredata[,6] %in% MethodType,]
        for(i in ID){
          if (length(unique(paleofiredata[paleofiredata[,1] %in% i,5]))>=2){
            paleofiredata[paleofiredata[,1] %in% i &
                            !(paleofiredata[,5] %in% paleofiresites$pref_units[paleofiresites[,1]==i]),7]=NA
          }
        }
        paleofiredata=paleofiredata[!is.na(paleofiredata$TYPE),]
        # Print which sites are dropped from the analysis
        cat(IDChar$site_name[!(IDChar$id_site %in% unique(paleofiredata[,1]))],"\n")
        cat(length(IDChar$site_name)-length(unique(paleofiredata[,1])), 
            " sites were excluded from the analysis \n")
        ID=unique(paleofiredata[,1])
        if(QuantType=="INFL"){
          for(i in ID)  
            if( !(unique(paleofiredata[paleofiredata[,1]==i,7]) %in% "INFL") & 
                is.na(sum(paleofiredata[paleofiredata[,1]==i,2]))==FALSE){
              infl=influx(paleofiredata[paleofiredata[,1]==i,])
              paleofiredata[paleofiredata[,1]==i,4]=c(infl)
            }
        }
        
      }
      ##-----
      ## Add users data
      if(is.null(add)==FALSE){
        # Add columns to match paleofiredata
        add$data=cbind(add$data,UNIT=NA,METHOD=NA,TYPE="INFL")
        # Then...
        paleofiredata=rbind(paleofiredata,add$data)
        ID=c(ID,unique(add$data[,1]))
      }
    }
  }
  if(is.null(ID)){
    paleofiredata=add$data
    ID=c(unique(add$data[,1]))
  }
  
  
  
  if (is.character(ID)){
    paleofiredata = read.csv(ID)
    ID=unique(paleofiredata[,1])
  }
  if (is.list(ID) & length(ID)>2)  {
    temp=ID$TransData
    depths=ID$IntDepths
    age=ID$Age
    sites=as.numeric(colnames(temp))
    ids=matrix(nrow=length(temp[,1]),ncol=length(temp[1,]))
    for (i in 1:length(temp[,1])){
      ids[i,]=sites
    }
    ids=c(ids)
    data=c(temp)
    age=c(age) 
    depths=c(depths)
    if(length(depths)==0) depths=rep(NA,length(age))
    paleofiredata=cbind(ids,depths,age,data)
    ID=unique( paleofiredata[,1])
  }
  if (is.matrix(ID)){
    paleofiredata=ID
    ID=unique(paleofiredata[,1])
  }
  
  # 2 Interpolate TRUE
  if (Interpolate==TRUE){
    # Interpolation procedure
    if (is.null(Age)) {
      res=matrix(ncol=1,nrow=length(ID))
      # Find the median time resolution for each paleofiredataset
      for (k in 1:length(ID)){
        resT=diff(paleofiredata[paleofiredata[,1]==ID[k],3])
        # Sometimes the last age is a copy of the previous one (why?)
        res[k]=c(median(resT[resT>0]))
      }
      # Find the median resolution for resampling
      step=round(median(res))
      minA=round(min(paleofiredata[,3]))
      maxA=round(max(paleofiredata[,3]))
      AgeN=seq(minA,maxA,step)
    }
    if (is.null(Age)==FALSE) {
      AgeN=Age
      #paleofiredata=paleofiredata[paleofiredata[,3]>min(AgeN),]
      #paleofiredata=paleofiredata[paleofiredata[,3]<max(AgeN),]
      #paleofiredata=paleofiredata[paleofiredata[,1] %in% ID,]
      ID=unique(paleofiredata[,1])
    }
    
    # Use linear interpolation to reconstruct a matrix of raw paleofiredata
    rawI=matrix(nrow=length(AgeN),ncol=length(ID))
    
    for (k in 1:length(ID)){
      if(length(paleofiredata[paleofiredata[,1]==ID[k],3])>=3){
        rawI[,k]=approx(paleofiredata[paleofiredata[,1]==ID[k],3],paleofiredata[paleofiredata[,1]==ID[k],4],AgeN, method = "linear")$y} else print(paste(IDChar$site_name[k], "has < 3 charcoal values and was excluded", sep=" "))
    }
    
    # Calculates Interpolated depths
    depthI=matrix(nrow=length(AgeN),ncol=length(ID))
    for (k in 1:length(ID)){
      if (is.na(sum(paleofiredata[paleofiredata[,1]==ID[k],2]))==F){
        if(length(paleofiredata[paleofiredata[,1]==ID[k],3])>=3){
          depthI[,k]=approx(paleofiredata[paleofiredata[,1]==ID[k],3],paleofiredata[paleofiredata[,1]==ID[k],2],AgeN,
                            method = "linear")$y
        }
      } else{depthI[,k]=NA}
    }
    
    ## Remove sites with less than 3 data values
    supp=c()
    for(i in 1:length(rawI[1,])){
      if(sum(!is.na(rawI[,i]))<3){supp[i]=1}else{supp[i]=0}
    }
    
    rawI=rawI[,supp==0]
    SuppSites=ID[supp==1]
    ID=ID[supp==0]  
    # Space for data
    transI=matrix(nrow=length(AgeN),ncol=length(ID))
    # Matrix of Ages (just a repeat)
    Ages=matrix(ncol=length(ID),nrow=length(AgeN))
    for (k in 1:length(ID)){
      Ages[,k]=c(AgeN)
    }
  }
  
  ## 3 No Interpolation:
  if (Interpolate==FALSE){
    # Which is the longest record?
    lengths=matrix(ncol=1,nrow=length(ID))
    for (k in 1:length(ID)){
      lengths[k]=c(length(paleofiredata[paleofiredata[,1] %in% ID[k],1]))
    }
    m=max(lengths)
    
    # Space for paleofiredata
    transI=matrix(nrow=m,ncol=length(ID))
    rawI=matrix(nrow=m,ncol=length(ID))
    depthI=matrix(nrow=m,ncol=length(ID))
    Ages=matrix(ncol=length(ID),nrow=m)
    
    # Matrix of Ages, rawData and depths
    for (k in 1:length(ID)){
      forNA=m-length(paleofiredata[paleofiredata[,1] %in% ID[k],3])
      AgeTemp=c(paleofiredata[paleofiredata[,1] %in% ID[k],3],rep(NA,forNA))
      Ages[,k]=c(AgeTemp)
      rawTemp=c(paleofiredata[paleofiredata[,1] %in% ID[k],4],rep(NA,forNA))
      rawI[,k]=c(rawTemp)
      depthTemp=c(paleofiredata[paleofiredata[,1] %in% ID[k],2],rep(NA,forNA))
      depthI[,k]=c(depthTemp)
    }
    ## End No Int
  }
  
  ## % Cat to see where we are
  if(verbose==TRUE){
    percent=seq(10,100,by=10)
    values=round(percent*length(ID)/100)
    cat("Transforming...")
    cat("\n")
    cat("Percentage done: ")
  }
  
  # Play with transformations!
  for (j in 1:length(method)){
    methodj=method[j]
    if (j>=2){rawI=transI}
    
    # Transformations
    for (k in 1:length(ID)){
      
      tmp=cbind(Ages[,k],rawI[,k])
      tmp=na.omit(tmp)
      ## At least 3 data values! 
      if(sum(tmp[,1])>3 & ID[k]!=882){
        # Not Tamagaucia site (882)!
        if(methodj=="stl") {
          agesI=seq(tmp[1,1],tmp[length(tmp[,1]),1],1)
          # Stl requires evenly spaced data
          forTS=approx(tmp[,1],tmp[,2],agesI)$y
          x=ts(forTS,start=1,frequency=stlYears)
          dim(x)=NULL
          stlResult=stl(x,"per")$time.series[,2]
          transI[,k]=approx(agesI,stlResult,Ages[,k])$y
        }
        if(methodj=='NULL'){
          transI[,k]=approx(tmp[,1],tmp[,2],Ages[,k])$y
        }
        if (methodj=="Z-Score") {
          mu=mean(tmp[tmp[,1]>=BasePeriod[1] & tmp[,1]<=BasePeriod[2],2])
          sigma=sd(tmp[tmp[,1]>=BasePeriod[1] & tmp[,1]<=BasePeriod[2],2])
          # No data in BasePeriod return scale:
          # Single value in BasePeriod return scale:
          if (is.na(mu) | is.na(sigma) | sigma==0){transI[,k]=approx(tmp[,1],scale(tmp[,2]),Ages[,k])$y}
          # Z-Score otherwise
          else {transI[,k]=approx(tmp[,1],(tmp[,2]-mu)/sigma,Ages[,k])$y}            
        }
        if (methodj=="Box-Cox") {
          transI[,k]=approx(tmp[,1],pfBoxCox(tmp[,2],alpha=alpha,type=type),Ages[,k])$y
        }
        if (methodj=="LOESS") {
          transI[,k]=approx(tmp[,1],predict(loess(tmp[,2]~ tmp[,1], span=span)),Ages[,k])$y
        }
        if (methodj=="MinMax") {
          transI[,k]=approx(tmp[,1],pfMinMax(tmp[,2]),Ages[,k])$y
        }
        if (methodj=="RunMed") {
          w=round(RunWidth/((max(tmp[,1])-min(tmp[,1]))/length(tmp[,1])))
          if (gtools::odd(w)) w=w else w=w+1
          transI[,k]=approx(tmp[,1],runmed(tmp[,2],w),Ages[,k])$y
        }
        if (methodj=="RunMean") {
          w=round(RunWidth/((max(tmp[,1])-min(tmp[,1]))/length(tmp[,1])))
          if (gtools::odd(w)) w=w else w=w+1
          transI[,k]=approx(tmp[,1],caTools::runmean(tmp[,2],w),Ages[,k])$y
        }
        if (methodj=="RunMin") {
          w=round(RunWidth/((max(tmp[,1])-min(tmp[,1]))/length(tmp[,1])))
          if (gtools::odd(w)) w=w else w=w+1
          transI[,k]=approx(tmp[,1],caTools::runmin(tmp[,2],w),Ages[,k])$y
        }
        if (methodj=="RunMax") {
          w=round(RunWidth/((max(tmp[,1])-min(tmp[,1]))/length(tmp[,1])))
          if (gtools::odd(w)) w=w else w=w+1
          transI[,k]=approx(tmp[,1],caTools::runmax(tmp[,2],w),Ages[,k])$y
        }
        if (methodj=="RunQuantile") {
          w=round(RunWidth/((max(tmp[,1])-min(tmp[,1]))/length(tmp[,1])))
          if (gtools::odd(w)) w=w else w=w+1
          transI[,k]=approx(tmp[,1],caTools::runquantile(tmp[,2],w,RunQParam),Ages[,k])$y
        }
        if (methodj=="SmoothSpline") {
          transI[,k]=approx(tmp[,1],smooth.spline(tmp[,1],tmp[,2],spar=span)$y,Ages[,k])$y
        }
        #         if (methodj=="GAM"){
        #           transI[,k]=approx(tmp[,1],gam(tmp[,2]~s(tmp[,1]))$fitted.values,Ages[,k])$y
        #         }
        if (methodj=="Hurdle"){
          # Transform data to count using pfMinMax
          tmp[,2]=round(pfMinMax(tmp[,2])*100)
          transI[,k]=approx(tmp[,1],pscl::hurdle(tmp[,2]~tmp[,1])$fitted.values,Ages[,k])$y
        }
      }
      if(k %in% values & verbose==TRUE & j==length(method)){
        cat(percent[values==k])
        cat(" ")
      }
    }
    ## j loop end
  }
  if(verbose==TRUE) cat("\n")
  
  ### End Return Results
  colnames(transI)=ID
  output=structure(list(Age=structure(Ages,col.names=as.character(ID),class="matrix" ),
                        IntDepths=structure(depthI,col.names=as.character(ID),class="matrix" ),
                        IntData=structure(rawI,col.names=as.character(ID),class="matrix" ),
                        TransData=structure(transI,col.names=as.character(ID),class="matrix"),
                        Method=method,
                        params=params
  ))
  class(output)="pfTransform"  
  return(output)
  
  ###
}


influx=function(x){
  ## Calculate Sed Acc
  d1=c()
  t1=c()
  for(k in 2:(length(x[,1])-1)){
    d1[k]=x[k+1,2]-x[k-1,2]
    t1[k]=x[k+1,3]-x[k-1,3] 
  }
  sedacc=(d1*100)/t1
  sedacc[1]=sedacc[2]
  sedacc=c(sedacc,sedacc[length(sedacc)])
  ## Calculate Influx
  infl=(x[,4]*sedacc)
  return(infl)
}

