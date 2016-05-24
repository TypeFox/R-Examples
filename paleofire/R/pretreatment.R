#' Calculate particules accumulation rates for sediment records
#' 
#' This is the R version of the CharAnalysis CharPretreatment.m function
#' originally develloped by P. Higuera and available at
#' https://sites.google.com/site/charanalysis
#' 
#' 
#' @param serie A proxy record to be transformed in accumulation rates, could
#' be particule counts, surfaces, volumes, etc.
#' @param params A matrix with the following colums: CmTop, CmBot, AgeTop,
#' AgeBot, Volume, in the same order.
#' @param Int Logical specifying whether the function interpolates particle
#' zero counts, default TRUE
#' @param first,last Date of the first, last sample for accumulation rate
#' calculation, if NULL first, last are automatically specified as the the
#' minimum and maximum ages of the record respectively
#' @param yrInterp Temporal resolution of the interpolated accumulation rates,
#' if NULL, yrInterp is automatically specified as the median resolution of the
#' record
#' @return Return an output structure with the following:
#' \item{cmI}{interpolated depths} \item{ybpI}{interpolated ages}
#' \item{accI}{accumulation rates}
#' @author O. Blarquez translated from P. Higuera CharPretreatment.m function
#' @examples
#' 
#' \dontrun{
#' # In this example we will use the charcoal record of the Lac du Loup from Blarquez et al. (2010).
#' # Blarquez, O., C. Carcaillet, B. Mourier, L. Bremond, and O. Radakovitch. 2010. Trees in the 
#' # subalpine belt since 11 700 cal. BP: origin, expansion and alteration of the modern forest. 
#' # The Holocene 20:139-146.
#' 
#' # Load raw charcoal data in mm^2
#' A=read.csv("http://blarquez.com/public/code/loupchar.csv")
#' C_=A[,6] # charcoal areas
#' P_=A[,1:5] # CmTop, CmBot, AgeTop, AgeBot, Volume 
#' 
#' # Calculates charcoal accumulation rate (CHAR, mm2.cm-2.yr-1)
#' CHAR=pretreatment(params=P_,serie=C_,Int=TRUE)
#' plot(CHAR)
#' }
#' 
pretreatment=function(params,serie,Int=TRUE,first=NULL,last=NULL,yrInterp=NULL){
  
  ## This is the R version of the CharAnalysis CharPretreatment.m function
  ## originally develloped by P. Higuera and available at https://sites.google.com/site/charanalysis
  ## Requires a matrix as input with the following columns 
  ## CmTop CmBot AgeTop AgeBot Volume 
  ## And a serie from which to calculate accumulation rates 
  
  A=cbind(params,serie)
  
  if(is.null(first)) first=min(A[,3])
  if(is.null(last)) last=max(A[,4])
  
  ## Redefine the matrix and values
  cm=A[,1]
  cmB=A[,2]
  count=A[,6]
  vol=A[,5]
  con=count/vol
  ybp=A[,3]
  
  ## Interpolate or not zero counts
  if (Int==TRUE){
    ## SCREEEN RECORD FOR MISSING VALUES:
    missingValuesIndex = which(count<=0)
    nMissingValues = length(missingValuesIndex)
    if (nMissingValues > 0){  # if some levels were not sampled...
      startIn = missingValuesIndex[which(c(99, diff(missingValuesIndex)) > 1)] # Index start of gaps
      # created by missing values.
      endIn = missingValuesIndex[which(diff(missingValuesIndex) > 1)] # Index end of gaps   
      # created by missing values. 
      endIn =c(endIn, missingValuesIndex[length(missingValuesIndex)]) # Add on last sample as last
      # end point for gaps.
      gapIn = cbind(startIn, endIn) # Index values for start (j = 1) and end 
      # (j = 2) of each gap.
      
      nGaps = length(startIn)  # Number of gaps in the record  
      cmGaps = sum(na.omit(cm[gapIn[,2]]-cm[gapIn[,1]])) # Sum 
      # of all cm of gap(s).
      yrGaps = sum(na.omit(ybp[gapIn[,2]]-ybp[gapIn[,1]],1)) # Sum 
      # of all cm of gap(s).
      # disp(['NOTE: ' num2str(nMissingValues) ' missing value(s) '...
      #       num2str(length(gapIn)) ' gap(s) totaling '...
      #       num2str(cmGaps) ' cm and ' num2str(round(yrGaps)) ' years.'])
      # disp('      Values created via interpolation.')
      
      for (i in 1:nGaps ){ # Fill in gaps via interpolation.
        # Last level is missing
        if(is.na(cm[gapIn[i,2]+1])) lag2=0 else lag2=1
        x = cbind(cm[gapIn[i,1]-1] , cm[gapIn[i,2]+1] )# [cm] End 
        # First level is missing
        if(length(x)==1) lag1=0 else lag1=1
        x = cbind(cm[gapIn[i,1]-lag1] , cm[gapIn[i,2]+lag2] )# [cm] End 
        # depths for interpolation.
        cmInterp = diff(as.numeric(A[gapIn[i,2],1:2]))
        xi = seq(A[gapIn[i,1]-lag1,1] ,A[gapIn[i,2]+lag2,1], cmInterp)
        # [cm] Desired interpolated depths. 
        y = c(A[gapIn[i,1]-lag1,5] , A[gapIn[i,2]+lag2,5]) # [cm^3] 
        # End volume for interpolation.
        y2 = c(A[gapIn[i,1]-lag1,6] , A[gapIn[i,2]+lag2,6]) # [cm^3] 
        # [pieces] 
        # End charcoal cound for interpolation.    
        yi = approx(x,y,xi)$y
        y2i = approx (x,y2,xi)$y
        if (length(yi[2:length(yi)-1]) != length(gapIn[i,1]:gapIn[i,2])){
          yi =yi[-2]      # Trim inerpolated values if there are more 
          y2i = y2i[-2]   # interpolated values than needed, given 
          # variable sampling intervals around gap. 
        }
        A[gapIn[i,1]:gapIn[i,2],5] = yi[2:length(yi)-1]    # Fill in 
        # interpolated sample volumes.
        A[gapIn[i,1]:gapIn[i,2],6] = y2i[2:length(y2i)-1]   # Fill in 
        # interpolated sample counts.
      } 
    } else  gapIn = NA # If no missing values, create empty 
    # variable to pass to CharPeakAnalysisResults.
  }
  
  ## RETRIEVE VARIABLES FROM INPUT FILES:
  cm = A[,1]    # [cm] sample depths.
  count =  A[,6] # [#] charcoal counts
  vol =  A[,5] # [cm^3] sample volumes.
  con = count / vol # [# cm-3] charcoal con.
  ybp =  A[,3]  # [cal ybp] age at top of sample
  
  ## Calculate sediment accumulation rate
  sedacc=c()
  for(i in 1:length(ybp)-1){
    sedacc[i]=c((cm[i+1]-cm[i])/(ybp[i+1]-ybp[i]))
  }
  sedacc=c(sedacc,0)
  
  ## Calculate yrInterp
  if(is.null(yrInterp)){
    yrInterp=round(median(diff(ybp)))
  }
  
  ## INTERPOLATE CHAROCAL DATA TO yrInterp INTERVALS:
  cmTop = A[,1]  # [cm] Depth at top of sample.
  cmBot = A[,2]  # [cm] Depth at bottom of sample.
  ageTop = A[,3] # [yr BP] Age at top of sample.
  ageBot = A[,4] # [yr BP] Age at bottom of sample.
  
  ybpI =seq(first, last, yrInterp) # [yr BP] Years to 
  # resample record to.
  
  propMatrix = matrix(nrow=length(ybpI),ncol=length(ybp))
  
  for (i in 1:length(ybpI)){ # For each resampled sample.
    rsAgeTop = ybpI[i]
    rsAgeBot = rsAgeTop+yrInterp
    
    for (j in 1:length(ybp)){# For each raw sample.
      # If raw sample straddles rsAgeBot
      if (ageTop[j] >= rsAgeTop && ageTop[j] < rsAgeBot && ageBot[j] > rsAgeBot){
        propMatrix[i,j] = c(rsAgeBot - ageTop[j])
      }
      # If raw sample straddles rsAgeTop
      if (ageTop[j] < rsAgeTop && ageBot[j] <= rsAgeBot && ageBot[j] > rsAgeTop){
        propMatrix[i,j] = c(ageBot[j] - rsAgeTop)
      }
      # If raw sample is entirely within resampled sample.
      if (ageTop[j] >= rsAgeTop && ageBot[j] <= rsAgeBot){
        propMatrix[i,j] =c( ageBot[j] - ageTop[j])
      }
      # If raw sample is entirely outside of resamples sample (i.e.
      # resampling finer than actual sample).
      if (ageTop[j] < rsAgeTop && ageBot[j] > rsAgeBot){
        propMatrix[i,j] = c(rsAgeBot - rsAgeTop)
      }
    }
  } # End making porportion matrix
  propMatrix = propMatrix / yrInterp
  
  # Determine values for each resampled interval.
  cmI = c()
  countI =  c()
  volI =  c()
  conI =  c()
  sedAccI =  c()
  
  for (i in 1:length(ybpI)){  # For each resampled sample
    inc = which(propMatrix[i,]>0) # Index for raw samples contributing to 
    # resampled sample.
    
    # Charcoal.cmI(i) = sum(Charcoal.cm(in) .* propMatrix(i,in)')
    countI[i] = sum(count[inc] * t(propMatrix[i,inc]))
    volI[i] = sum(vol[inc] * t(propMatrix[i,inc]))
    conI[i] = sum(con[inc] * t(propMatrix[i,inc]))
    sedAccI[i] = sum(sedacc[inc] * t(propMatrix[i,inc]))
  }
  cmI = approx(ybp,cm,ybpI)$y # [cm] 
  # interpolated depths.
  
  ## DERIVE CHARCOAL ACCUMULATOIN RATE:
  acc=c()
  acc = (con * sedacc) # [#/cm2/yr] take sedAcc and 
  # multiply by Charcoal.con to get  Charcoal.acc
  accI = (conI * sedAccI) # [#/cm2/yr] take sedAccI and 
  # multiply by Charcoal.conI to get Charcoal.accI 
  
  ## Return values 
  output=structure(list(cmI=cmI,ybpI=ybpI,accI=accI,ageTop=ageTop,ageBot=ageBot,yrInterp=yrInterp,acc=acc))
  class(output)="CHAR"
  return(output)
  ## Et Hop
}





#' Plot CHAR
#' 
#' Plot an object of the class "CHAR" returned by the pretreatment function.
#' Original accumulation rates are presented using grey bars, accumulation
#' rates interpolated at equal time steps are presented by a black curve.
#' 
#' @method plot CHAR
#' @export
#' @param x An object of the class "CHAR".
#' @param \dots \dots{}
#' @author O. Blarquez
#' @examples
#' 
#' ## In this example we will use the charcoal record of the Lac du Loup (Blarquez et al. 2010)
#' ## Load raw charcoal data in mm^2
#' A=read.csv("http://blarquez.com/public/code/loupchar.csv")
#' C_=A[,6] # charcoal areas
#' P_=A[,1:5] # CmTop, CmBot, AgeTop, AgeBot, Volume 
#' 
#' 
#' ## Calculates charcoal accumulation rate (CHAR, mm2.cm-2.yr-1)
#' CHAR=pretreatment(params=P_,serie=C_,Int=TRUE)
#' plot(CHAR)
#' 
#' 
plot.CHAR=function(x,...){
  
  ## PLOT
  # Values for plotting
  age=c(matrix(c(x$ageTop, x$ageBot), 2, byrow = TRUE)) 
  accInit<-rep(x$acc,each=2)
  ageInt=c(matrix(c(x$ybpI,x$ybpI+x$yrInterp), 2, byrow = TRUE)) 
  accInt<-rep(x$accI,each=2)
  
  # plot
  plot(age,accInit,type="l",col="grey",ylim=c(min(na.omit(accInit)),max(na.omit(accInit))))
  polygon(c(age,age[length(age)]), c(accInit,-100), col='grey') 
  lines(age,accInit,type="l",col="grey")
  lines(ageInt,accInt,type="l")
}
