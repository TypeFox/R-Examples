
################################################################################
##  compareCrimes.R
##  Various functions for creating evidence variables
################################################################################

#-- Required R Packages
#require(geosphere) # for compareSpatial() if lat/long coordinates used

##  compareCrimes
##==============================================================================
#'  Creates evidence variables by calculating \sQuote{distance} between crime pairs
#'
#'  Calculates spatial and temporal distance, difference in categorical, and
#'   absolute value of numerical crime variables
##  Inputs:
#'  @param Pairs (n x 2) matrix of crimeIDs
#'  @param crimedata data.frame of crime incident data. There must be a column 
#'   named \code{crimedata} that refers to the crimeIDs given in \code{Pairs}. 
#'   Other column names must correspond to what is given in \code{varlist} list. 
#'  @param varlist a list with elements named: crimeID, spatial, temporal,
#'    categorical, and numerical. Each element should be a vector of the column 
#'    names of \code{crimedata} corresponding to that feature:
#'   \itemize{  
#'    \item crimeID: crime ID for the crimedata that is matched to \code{Pairs}
#'    \item spatial: X,Y coordinates (in long,lat or Cartesian) of crimes
#'    \item temporal: DT.FROM, DT.TO of crimes. If times are uncensored, then 
#'     only DT.FROM needs to be provided.
#'    \item categorical: (optional) categorical crime variables
#'    \item numerical: (optional) numerical crime variables
#'   }
#'  @param binary (logical) match/no match or all combinations for categorical
#'    data
#'  @param longlat (logical) are spatial coordinates in (long,lat)?
#'  @param show.pb (logical) show the progress bar
#'  @param \ldots other arguments passed to hidden functions
##  Outputs:
#'  @return data.frame of various proximity measures between the two crimes
#'  \itemize{
#'   \item If \code{spatial} data is provided: the euclidean distance 
#'   (if \code{longlat = FALSE}) or Haversine great circle distance 
#'   (\code{\link[geosphere]{distHaversine}} if \code{longlat = TRUE}) is 
#'   returned (in kilometers).
#'   
#'   \item If \code{temporal} data is provided: the expected absolute time 
#'    difference is returned: 
#'    \itemize{
#'      \item temporal - overall difference (in days)  [0,max]
#'      \item tod - time of day difference (in hours)  [0,12]
#'      \item dow - fractional day of week difference (in days) [0,3.5]
#'    }
#'    
#'   \item If \code{categorical} data is provided: if \code{binary = TRUE} then a 
#'    1 if the categories of each crime match and a 0 if they do not match. If 
#'    \code{binary = FALSE}, then a factor of merged values (in form of f1:f2)
#'    
#'   \item If \code{numerical} data is provided: the absolute difference is 
#'    returned.
#'  }
#'  @examples
#'  data(crimes)
#'  pairs = t(combn(crimes$crimeID[1:4],m=2))   # make some crime pairs
#'  
#'  varlist = list(
#'    spatial = c("X", "Y"),
#'    temporal = c("DT.FROM","DT.TO"),
#'    categorical = c("MO1",  "MO2", "MO3"))    # crime variables list
#'      
#'  compareCrimes(pairs,crimes,varlist,binary=TRUE)    
#'  
#'  @references
#'  Porter, M. D. (2014). A Statistical Approach to Crime Linkage.
#'    \emph{arXiv preprint arXiv:1410.2285.}.
#'  \url{http://arxiv.org/abs/1410.2285}
#'  @export
##==============================================================================
compareCrimes <- function(Pairs,crimedata,varlist,binary=TRUE,longlat=FALSE,
                          show.pb=FALSE,...){  

  if(class(crimedata) != 'data.frame') crimedata = as.data.frame(crimedata)
  if(!all(unlist(varlist) %in% colnames(crimedata))) 
    stop("There are names in varlist that don't match column names in crimedata")  
  
  valid.types = c("crimeID","spatial","temporal","categorical","numerical")
  if( !all(names(varlist) %in% valid.types) ) 
    warning(paste(setdiff(names(varlist),valid.types), '(from varlist)',
                  'is not a valid type of crime data and will not be used to', 
                  'make evidence variables.'))
  
  crimeID = as.character(crimedata$crimeID)
  i1 = match(Pairs[,1],crimeID)
  i2 = match(Pairs[,2],crimeID)

  d.spat = d.temp = d.cat = d.num = NA

  if(!is.null(varlist$spatial)){
    spatial = crimedata[varlist$spatial]  
    if(ncol(spatial) != 2) stop("spatial must be a data.frame of spatial coordinates 
                                with two columns")
    d.spat = compareSpatial(spatial[i1,1:2],spatial[i2,1:2],longlat=longlat)
  }
  
  if(!is.null(varlist$temporal)){
    temporal = crimedata[varlist$temporal]  
    if(ncol(temporal) == 1) temporal = cbind(temporal,temporal) # no censoring
    d.temp = compareTemporal(temporal[i1,1:2],temporal[i2,1:2],show.pb=show.pb,...)  
  }

  catNames = varlist$categorical
  if(!is.null(catNames)){
    d.cat = data.frame(matrix(nrow=length(i1),ncol=length(catNames),
                                    dimnames=list(NULL,catNames)))
    for(cat in catNames){
      d.cat[,cat] = compareCategorical(crimedata[i1,cat],crimedata[i2,cat],
                                             binary=binary)
    }  
  }

  numNames = varlist$numerical
  if(!is.null(numNames)){
    d.num = data.frame(matrix(nrow=length(i1),ncol=length(numNames),
                                  dimnames=list(NULL,numNames)))
    for(num in numNames){
      d.num[,num] = compareNumeric(crimedata[i1,num],crimedata[i2,num])
    } 
  }  

  Dist = data.frame(spatial=d.spat,d.temp,d.cat,d.num)
  Dist[sapply(Dist,function(x) all(is.na(x)))] <- list(NULL)  
  E = cbind(Pairs,Dist,stringsAsFactors=FALSE)  # combine crime pair info with evidence variables
  rownames(E) <- NULL  # remove rownames  
return(E)
}


##  compareTemporal
##==============================================================================
#'  Make temporal evidence variable from (possibly uncertain) temporal info
#'
#'  Calculates the temporal distance between crimes
##  Inputs:
#'  @param DT1 (n x 2) data.frame of (DT.FROM,DT.TO) for the crimes
#'  @param DT2 (n x 2) data.frame of (DT.FROM,DT.TO) for the crimes
#'  @param show.pb (logical) show the progress bar
#'  @param \ldots other arguments passed to \code{\link{expAbsDiff.circ}}
##  Outputs:
#'  @return data.frame of expected absolute differences:
#'    \itemize{
#'      \item temporal - overall difference (in days)  [0,max]
#'      \item tod - time of day difference (in hours)  [0,12]
#'      \item dow - fractional day of week difference (in days) [0,3.5]
#'    }
#'  @keywords internal
##==============================================================================
compareTemporal <- function(DT1,DT2,show.pb=FALSE,...){
  #- interval length (in hours)  
  L1 = as.numeric(abs(difftime(DT1[,2],DT1[,1],units='hours')))  
  L2 = as.numeric(abs(difftime(DT2[,2],DT2[,1],units='hours')))  
  #- time (seconds since origen)
  day1 = as.numeric(julian(DT1[,1])) 
  day2 = as.numeric(julian(DT2[,1]))  
  #- time of day (at origen, not actual)
  tod1 = (as.numeric(DT1[,1])/3600) %% 24
  tod2 = (as.numeric(DT2[,1])/3600) %% 24
  #- time of day (at origen, not actual)
  dow1 = (as.numeric(DT1[,1])/(3600*24)) %% 7            
  dow2 = (as.numeric(DT2[,1])/(3600*24)) %% 7            

  #- calculate temporal absolute differences
  n = length(L1)
  temporal = tod = dow = numeric(n)
  if(show.pb) pb = txtProgressBar(style=3,max=n)  
  for(i in 1:n){
    temporal[i] = expAbsDiff(c(day1[i],day1[i]+L1[i]/24),c(day2[i],day2[i]+L2[i]/24))
    tod[i] = expAbsDiff.circ(c(tod1[i],tod1[i]+L1[i]),c(tod2[i],tod2[i]+L2[i]),mod=24,...)
    dow[i] = expAbsDiff.circ(c(dow1[i],dow1[i]+L1[i]/24),c(dow2[i],dow2[i]+L2[i]/24),mod=7,...)
    if(show.pb) setTxtProgressBar(pb,i)
  } 
  if(show.pb) close(pb)
return(data.frame(temporal,tod,dow) )             # expected absolute difference
}



##  expAbsDiff
##==============================================================================
#'  Expected absolute difference of two uniform RVs 
#'
#'  Calculates the expected absolute difference of two uniform rv's
##  Inputs:
#'  @param X c(min,max)
#'  @param Y c(min,max)
##  Outputs:
#'  @return the expected absolute difference
#'  @keywords internal
##==============================================================================
expAbsDiff <- function(X,Y){
  if(X[2]<X[1]) stop("X[2] < X[1]")
  if(Y[2]<Y[1]) stop("Y[2] < Y[1]")
  if(X[1]<=Y[1]){  # set Sx to have minimum lower bound
    Sx = X
    Sy = Y
  } else{
    Sx = Y
    Sy = X    
  }
  
  # Scenario 1 (no overlap)
  if(Sx[2] <= Sy[1]){ 
    return(mean(Sy) - mean(Sx))
  }

  bks = sort(c(Sx,Sy))
  sz = diff(bks)
  mids = bks[-1] - sz/2   
  
  # Scenario 2 (partial overlap)
  if(Sx[2] <= Sy[2]){ 
    px = sz*c(1,1,0) / diff(Sx)
    py = sz*c(0,1,1) / diff(Sy) 
    return( 
      (mids[2]-mids[1])*px[1]*py[2] +
      (mids[3]-mids[1])*px[1]*py[3]+
      (sz[2]/3)*px[2]*py[2] +
      (mids[3]-mids[2])*px[2]*py[3] )
  }
  
  # Scenario 3 (Y is subset X)
  if(Sx[2] > Sy[2]){
    px = sz*c(1,1,1) / diff(Sx)
    #py = sz*c(0,1,0) / diff(Sy)   
    return(
      (mids[2]-mids[1])*px[1] + 
      (sz[2]/3)*px[2] +
      (mids[3]-mids[2])*px[3]  
      )
  }
}  


##  expAbsDiff.circ
##==============================================================================
#'  Expected absolute difference of two circular uniform RVs
#'
#'  Estimates the expected circular temporal distance between crimes using discrete 
#'   FFT or numerical integration
##  Inputs:
#'  @param X c(min,min+length). X[1] must be >= 0 and X[2] >= X[1]. It is possible
#'   that X[2] can be > mod. I.e., do not do \code{X \%\% mod}
#'  @param Y c(min,min+length). Same conditions from X applies.
#'  @param mod the period of time. E.g., mod=24 for time of day (in hours), 
#'   mod=7 for day of week (in days)
#'  @param n number of bins for discretization of continuous time domain.   
#'  @param method use convolution or monte carlo integration (\code{montecarlo})
##  Outputs:
#'  @return the expected absolute difference
#'  @keywords internal
##  Note: There is 1440 min/day, so n = 2000 should give close to minute resolution 
##==============================================================================
expAbsDiff.circ <- function(X,Y,mod=24,n=2000,method="convolution"){
  if(X[1]<0 | Y[1]<0) stop("X and Y must be >0")
  if(diff(X)>=mod | diff(Y)>=mod) return(mod/4)  # uniform over mod
  if(diff(X)==0) return(getD(X[1],Y,mod=mod))           
  if(diff(Y)==0) return(getD(Y[1],X,mod=mod))    
  if(method=="convolution"){
    while((n %% 2) != 0) n = nextn(n+1)          # n must be even and nextn
    theta = seq(0,mod,length=n+1)
    delta = diff(theta[1:2])
    x = diff(punif(theta,X[1],X[2])) + diff(punif(theta+mod,X[1],X[2]))
    y = diff(punif(theta,Y[1],Y[2])) + diff(punif(theta+mod,Y[1],Y[2]))  
    conv = convolve(x,y,conj=TRUE,type="circular")
    tt = ifelse(theta<=mod/2,theta,mod-theta)[-(n+1)]
    d = sum(tt*conv)
  }
  if(method=="numerical"){
    if(diff(Y)<diff(X)){
      tt = seq(Y[1],Y[2],length=n)
      d = mean(getD(tt,X,mod=mod))
    } else{
      tt = seq(X[1],X[2],length=n)
      d = mean(getD(tt,Y,mod=mod))    
    }
  }
return(d)  
}


##  getD
##==============================================================================
#'  Expected absolute distance of a circular uniform RV to a point
#'
#'  Calculates the expected circular distance between a uniform rv and a point
##  Inputs:
#'  @param y vector of times in [0,mod)
#'  @param X c(min,min+length). X[1] must be >= 0 and X[2] >= X[1]. It is possible
#'   that X[2] can be > mod. I.e., do not do \code{X \%\% mod}
#'  @param mod the period of time. E.g., mod=24 for time of day (in hours), 
#'   mod=7 for day of week (in days)
##  Outputs:
#'  @return the expected absolute difference
#'  @keywords internal
##==============================================================================
getD <- function(y,X,mod=24){
  if(X[1] > mod | X[1] < 0) stop("Minimum X[1] not within limits [0,mod)")
  if(X[2] < X[1]) stop("X[2] must be >= X[1]")
  y = (y - X[1]) %% mod  
  B = X[2] - X[1]        # length of interval 
  if(B == 0) return(mod/2-abs(mod/2-abs(y)))  # For |X| = 0
  if(B >= mod) return(rep(mod/4,length(y)))   # For long intervals
  D = numeric(length(y))  
  if(diff(X)>=mod/2){    
    K = mod - B/2 - (mod/2)^2/B
    u = y-mod/2
    i1 = (y <= B-mod/2)
    D[i1] = y[i1]*(1-mod/B) + K
    i2 = (y > B-mod/2 & y <= mod/2)
    D[i2] = (y[i2]-B/2)^2/B  +B/4   # mod/2 - (u^2+(B-u)^2)/(2*B)
    i3 = (y > mod/2 & y <= B)
    D[i3] =  (B-y[i3])*(1-mod/B) + K
    i4 = (y > B)
    D[i4] = mod/2 - B/4 - ((u[i4]-B/2))^2/B
  }
  else{            
    u = y-B/2
    i1 = (y<B)
    D[i1] = u[i1]^2/B + B/4
    i2 = (y>=B & y<=mod/2)
    D[i2] = u[i2]
    i3 = (y>mod/2 & y<=B+mod/2)
    D[i3] = mod/2 - ((y[i3]-mod/2)^2 + (B-y[i3]+mod/2)^2) /(2*B)  
    i4 = (y > B+mod/2)  
    D[i4] = mod - u[i4]
  }
  return(D)  
}


##  compareSpatial
##==============================================================================
#'  Make spatial evidence variables
#'
#'  Calculates spatial distance between crimes (in km)
##  Inputs:
#'  @param C1 (n x 2) matrix of coordinates for the crimes
#'  @param C2 (n x 2) matrix of coordinates for the crimes
#'  @param longlat (logical) if true, the the coordinates are in (Long,Lat), else
##      assume a suitable project where euclidean distance can be applied
##  Outputs:
#'  @return numeric vector of distances between the crimes (in km)
#'  @keywords internal
##==============================================================================
#library(geosphere)
compareSpatial <- function(C1,C2,longlat=FALSE){
  if(longlat)  d = geosphere::distHaversine(C1,C2)
  else         d = sqrt(rowSums((C1 - C2)^2))
  return(d/1000)            # distance in km
}

##  compareNumeric
##==============================================================================
#' Make evidence variables from numeric crime data
#'
#' Calculates absolute difference between crimes variables
##  Inputs:
#'  @param C1 length n numerical values of crime attributes
#'  @param C2 length n numerical values of crime attributes
##  Outputs:
#'  @return numeric vector of absolute differences
#'  @keywords internal
##==============================================================================
compareNumeric <- function(C1,C2) abs(C1 - C2)


##  compareCategorical
##==============================================================================
#'  Make evidence variables from categorical crime data
#'
#'  Compares categorical crime data to check if they match.
##  Inputs:
#'  @param C1 length n categorical values of crime attributes
#'  @param C2 length n categorical values of crime attributes
#'  @param levs the levels of all possible values
#'  @param  binary (logical) match/no match or all combinations
##  Outputs:
#'  @return if binary=TRUE: 1 for match, 0 for non-matches;
#'  if binary=FALSE: factor vector of merged values (in form of f1:f2)
#'  @keywords internal
##==============================================================================
compareCategorical <- function(C1,C2,levs,binary=FALSE){
  if(binary){
    match.na = is.na(C1) & is.na(C2)   # counts matching NA's as a match
    match.value = as.character(C1)==as.character(C2)
    B = ifelse(match.value | match.na,1,0)
    B[is.na(B)] = 0    # set rest of NAs to non-match
    return(factor(B,levels=c(0,1)))
  }
  C1 = as.factor(C1); C2 = as.factor(C2)
  if(missing(levs)) levs = sort(unique(levels(C1),levels(C2)))
  A = data.frame(C1=factor(C1,levels=levs), C2=factor(C2,levels=levs))
  flip = which(is.na(A[,1]) | unclass(A[,1]) > unclass(A[,2]))
  A[flip,] = A[flip,2:1]    # Sort rows
  B = paste(A[,1],A[,2],sep=':')        # Merge values
  B = factor(B,levels=catLevels(levs))  # Make into factor
  return(B)
}

##  catLevels
##==============================================================================
#'  Make levels for merging category predictors
##  Inputs:
#'  @param levs levels of a catagorical variable (factor)
##  Outputs:
#'  @return levels for a new categorical variable of form f1:f2
#'  @keywords internal
##==============================================================================
catLevels <- function(levs){
  levs = unique(c(levs,NA))  # Add NA if not already included
  nlevs = length(levs)
  a = NULL
  for(i in 1:nlevs){
    b = cbind(levs[i],levs[i:nlevs])
    a = rbind(a,b)
  }
  levs2 = paste(a[,1],a[,2],sep=':')
  return(levs2)
}

