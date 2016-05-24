# A class for histogram-valued data ---------------------------------------
#' Class distributionH.
#' @name distributionH-class
#' @rdname distributionH-class
#' @exportClass distributionH
setClass(
  Class = "distributionH",
  representation = representation(x="numeric",
                                  p="numeric",
                                  m="numeric",
                                  s="numeric"),
  #distributionH=
  validity = function(object){
    
    if (length(object@x)<=0){
      x=NA}
    else {
      if(length(object@p)){
        nv=length(object@x)
        p=c(0,c(1:nv)/nv)}
      if (length(object@x)!=length(object@p)){
        stop("the x and p vectors must be of the same length")}
      if (length(object@x)==1){
        object@x=c(object@x,object@x)
        object@p=c(0,1)
      }
      nv=length(object@x)
      
      if ( !is.na(object@x) && min(object@x[2:nv]-object@x[1:(nv-1)])<0){
        return("the x must be in not descending order")
      }
      if (!is.na(object@p) && (min(object@p[2:nv]-object@p[1:(nv-1)])<0 || object@p[1]<0 || object@p[nv]>1 || object@p[1]!=0 || object@p[nv]!=1))
      {return("the p must be in not descending order from 0 to 1")}
      return(TRUE)
    }
  }
)
#' Constructor method of distributionH class
#' 
#' Class \code{distributionH} defines a histogram object
#' 
#' @name distributionH
#' @rdname distributionH-class
#' @aliases initialize,distributionH-method
#' @description Class \code{"distributionH"} desfines an histogram object
#' The class describes a histogram by means of its cumulative distribution
#' function. The methods are develoved accordingly to the L2 Wasserstein
#' distance between distributions.
#'  
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("distributionH", x, p, m, s)}.
#' @param .Object the type ("distributionH")
#' @param x a numeric vector. it is the domain of the distribution (i.e. the
#' extremes of bins).
#' @param p a numeric vector (of the same lenght of x). It is the cumulative distribution function CDF.
#' @param m (optional) a numeric value. Is the mean of the histogram.
#' @param s (optional) a numeric positive value. It is the standard deviation of a histogram.
#' @author Antonio Irpino
#' @references Irpino, A., Verde, R. (2015) \emph{Basic
#' statistics for distributional symbolic variables: a new metric-based
#' approach} Advances in Data Analysis and Classification, DOI
#' 10.1007/s11634-014-0176-4
#' @keywords classes
#' @examples  #---- initialize a distributionH object mydist
#' # from a simple histogram 
#' # ----------------------------
#' # | Bins    |  Prob  | cdf   |
#' # ----------------------------
#' # | [1,2)   |  0.4   | 0.4   |
#' # | [2,3]   |  0.6   | 1.0   |
#' # ----------------------------
#' # | Tot.    |  1.0   | -     |
#' # ----------------------------
#' mydist=new("distributionH",c(1,2,3),c(0, 0.4, 1))
#' str(mydist)
#' # OUTPUT
#' # Formal class 'distributionH' [package "HistDAWass"] with 4 slots
#' #   ..@@ x: num [1:3] 1 2 3 the quantiles
#' #   ..@@ p: num [1:3] 0 0.4 1 the cdf
#' #   ..@@ m: num 2.1 the mean
#' #   ..@@ s: num 0.569 the standard deviation
#' @seealso \code{\link{meanH}} computes the mean. \code{\link{stdH}} computes the standard deviation.
# showClass("distributionH")
setMethod("initialize", "distributionH",
          definition=function(.Object,x=numeric(0),p=numeric(0),m=numeric(0),s=numeric(0)){
            .Object@x=x
            .Object@p=p
            if (length(x)>0){
              validObject(.Object)
              if (length(m)==0) .Object@m=meanH(.Object) else .Object@m=m
              if (length(s)==0) .Object@s=stdH(.Object) else .Object@s=s
            }
            
            return(.Object)
          }
)

# A class of a matrix of histogram-valued data ----------------------------
#' Class MatH.
#' @name MatH-class
#' @rdname MatH-class
#' @exportClass MatH
#' @docType class
#' @description   Class \code{MatH} defines a matrix of \code{distributionH} objects
#' @author Antonio Irpino
#' @references Irpino, A., Verde, R. (2015) \emph{Basic
#' statistics for distributional symbolic variables: a new metric-based
#' approach} Advances in Data Analysis and Classification, DOI
#' 10.1007/s11634-014-0176-4
#' @keywords classes
#' @examples
#' 
#' ##---- create a list of six distributionH objects
#' ListOfDist<-vector("list",6)
#' ListOfDist[[1]]<-distributionH(c(1,2,3),c(0, 0.4, 1))
#' ListOfDist[[2]]<-distributionH(c(7,8,10,15),c(0, 0.2, 0.7, 1))
#' ListOfDist[[3]]<-distributionH(c(9,11,20),c(0, 0.5, 1))
#' ListOfDist[[4]]<-distributionH(c(2,5,8),c(0, 0.3, 1))
#' ListOfDist[[5]]<-distributionH(c(8,10,15),c(0,  0.75, 1))
#' ListOfDist[[6]]<-distributionH(c(20,22,24),c(0, 0.12, 1))
#' 
#' ## create a MatH object filling it by columns
#' MyMAT=new("MatH",nrows=3,ncols=2,ListOfDist=ListOfDist,
#'   names.rows=c("I1","I2","I3"), names.cols=c("Var1","Var2"),by.row=FALSE)
#' 
#' showClass("MatH")
#'
#' @param .Object the object type "MatH"
#' @param ListOfDist a vector or a list of \code{distributionH} objects
#' @param names.rows a vector or list of strings with thenames of the rows
#' @param names.cols a vector or list of strings with thenames of the columns (variables)

setClass(
  Class = "MatH",
  representation = representation(M="matrix"),
)
#' Constructor method for MatH class
#' @name MatH
#' @rdname MatH-class
#' @aliases initialize,MatH-method
setMethod("initialize", "MatH",
          definition=function(
            .Object,nrows=1,ncols=1,ListOfDist=NULL, names.rows=NULL, names.cols=NULL,by.row=FALSE)
            {
            tt=new("list");
            .Object@M=matrix(tt,nrows,ncols)
            
            if (length(ListOfDist)>0){
              nOBJ=length(ListOfDist);
              if (by.row){
                count=0;
                for (i in 1:nrows){
                  for (j in 1:ncols){                  
                    count=count+1
                    if (count==nOBJ) count=1
                    .Object@M[i,j][[1]]=ListOfDist[[count]]
                  }
                }
              }
              else
              {
                count=0;
                for (j in 1:ncols){
                  for (i in 1:nrows){                  
                    count=count+1
                    if (count>nOBJ) count=1
                    .Object@M[i,j][[1]]=ListOfDist[[count]]
                  }
                }
              }
            }
            else
            {
              for (i in 1:nrows){
                for (j in 1:ncols){
                  .Object@M[i,j][[1]]=new("distributionH")
                }
              }
            }
          
            if (length(names.rows)>0){
              
              count=0; 
              rnames=vector("list",nrows)
              
              for (i in 1:nrows){
                count=count+1
                
                if (count>length(names.rows)){
                  rnames[[count]]=paste("I",count,sep="")
                }
                else{
                  rnames[[count]]=names.rows[[count]]               }
              }
              rownames(.Object@M)=rnames
            }
            else{
              rownames(.Object@M)=paste("I",1:nrows,sep="")}
            
            if (length(names.cols)>0){
              count=0
              cnames=vector("list",ncols)
              for (i in 1:ncols){
                count=count+1
                if (count>length(names.cols)){
                  cnames[[count]]=paste("X",count,sep="")
                  }                
                else{
                  cnames[[count]]=names.cols[[count]] 
                  }
              colnames(.Object@M)=cnames
              
              }
            }else  
              colnames(.Object@M)=paste("X",1:ncols,sep="")
            
            return(.Object)
          }
)
## Classes for Histogram Time Series ------
## A single distribution with a time stamp ------
#' Class TdistributionH
#'
#' Class \code{TdistributionH} defines a histogram with a time (point or period)
#' 
#' @name TdistributionH-class
#' @rdname TdistributionH-class
#' @exportClass TdistributionH
setClass(
  Class = "TdistributionH",
  representation = representation(tstamp="numeric",
                                  period="list"),
  contains = "distributionH",
  validity = function(object){
    if (length(object@tstamp)>1){
      stop("time stamp slot may contain a single value")
    }
    if (length(object@period)>2){
      stop("period slot may contain no more that two ordered values")
    }
    if ((length(object@period)==2)&&(object@period[[1]]>object@period[[2]])){
      stop("period slot error: the end of period is lower than the starting one")
    }
  }
)
## Initialize ------
#' Constructor method of TdistributionH Class
#' 
#' @name TdistributionH
#' @rdname TdistributionH-class
#' @aliases initialize,TdistributionH-method
#' @param .Object the type of object ("TdistributionH") a \code{"distributionH"} object with a time reference
#' @param tstamp a numeric value related to  a timestamp
#' @param period a list of two values, the starting time and the ending time (alternative to tstamp if the
#' distribution is observed along a period and not on a timestamp)
#' @param x a vector of increasing values, the domain of the distribution (the same of \code{distributionH} object)
#' @param p a vector of increasing values from 0 to 1, 
#' the CDF of the distribution (the same of \code{distributionH} object)
#' @param m a number, the mean of the distribution (the same of \code{distributionH} object)
#' @param s a positive number, the standard deviation of the distribution (the same of \code{distributionH} object)
setMethod("initialize", "TdistributionH",
          definition=function(.Object,
                              tstamp=numeric(0),
                              period=list(start=-Inf,end=-Inf),
                              x=numeric(0),
                              p=numeric(0),
                              m=numeric(0),
                              s=numeric(0)){
            .Object@x=x
            .Object@p=p
            .Object@m=m
            .Object@s=s
            .Object@tstamp=tstamp
            .Object@period=period
            if (length(x)>0){
              validObject(.Object)
              if (length(m)==0) .Object@m=meanH(.Object) else .Object@m=m
              if (length(s)==0) .Object@s=stdH(.Object) else .Object@s=s
            }
            return(.Object)
          }
)
#Coerce a TdistributionH into a distributionH
setAs(from = "TdistributionH",to = "distributionH",
      function(from,to){
        to=new("distributionH",
               x=from@x,
               p=from@p,
               m=from@m,
               s=from@s)
        return(to)
      }
)

## A mulvariate MatH with a time stamp  ----
#' Class TMatH
#'
#' Class \code{TMatH} defines a matrix of histograms, a \code{TMatH} object, with a time (a timepoint or a time window).
#' 
#' @name TMatH-class
#' @rdname TMatH-class
#' @exportClass TMatH
setClass(
  Class = "TMatH",
  representation = representation(tstamp="numeric",
                                  period="list"),
  contains = "MatH",
  validity = function(object){
    if (length(tstamp)>1){
      stop("time stamp slot may contain a single value")
    }
    if (length(period)>2){
      stop("period slot may contain no more that two ordered values")
    }
    if ((length(period)==2)&&(period[1]>period[2])){
      stop("period slot error: the end of period is lower than the starting one")
    }
  }
)
## Initialize TMatH ------
#' Constructor method of TdistributionH Class
#' 
#' @name TMatH
#' @rdname TMatH-class
#' @aliases initialize,TMatH-method
#' @param .Object the type of object ("TMatH")
#' @param tstamp a vector of time stamps, numeric.
#' @param period a list of pairs with a vectorof starting time and a vector ofending time.
#' This parameter is used alternatively to \code{tstamp} if the distributions are related to time periods
#' instead of timestamps
#' @param mat a \code{MatH} object
setMethod("initialize", "TMatH",
          definition=function(.Object,
                              tstamp=numeric(0),
                              period=list(start=-Inf,end=-Inf),
                              mat=new("MatH")
          ) {
            .Object@M=mat
            .Object@tstamp=tstamp
            .Object@period=period
            
            return(.Object)
          }
)

## A Histogram Time Series HTS-----
#' Class HTS
#'
#' Class \code{HTS} defines a histogram time series, i.e. a set of histograms observed along time
#' 
#' @name HTS-class
#' @rdname HTS-class
#' @exportClass HTS
setClass(
  Class = "HTS",
  representation = representation(data="vector"),
  validity = function(object){
    if (length(object@data)>0){
      count=0
      while (count<length(object@data)){
        count=count+1
        if (count==1){
          type=is(object@data[[1]])[1]
          if (!(type=="TdistributionH")||(type=="TMatH")){
            stop("data must be  TdistributionH or TMatH objects")
          }
        }
        else{
          if (is(object@data[[count]])[1]!=type){
            stop("all data must be of the same type")
          }
          if (type=="TMatH"){
            if ((nrow(object@data[[i]]@M)!=nrow(object@data[[i-1]]@M))||
                  (ncol(object@data[[i]]@M)!=ncol(object@data[[i-1]]@M))){
              stop("all TMatH objects must have the same dimensions")
            }
          }
        }
      }
    }
  }
)
## Initialize HTS ------
#' Constructor method of HTS Class (Histogram Time Series)
#' 
#' @name HTS
#' @rdname HTS-class
#' @aliases initialize,HTS-method
#' @param .Object the object type ("HTS") a histogram time series
#' @param epocs the number of histograms (one for each timepoint or period)
#' @param ListOfTimedElements a vector of \code{TdistributionH} objects
setMethod("initialize", "HTS",
          definition=function(.Object,epocs=1, ListOfTimedElements=c(new("TdistributionH"))){
            counts=0
            for (i in 1:epocs){
              counts=counts+1
              if (counts>length(ListOfTimedElements)){counts=1}  
              .Object@data=c(.Object@data,ListOfTimedElements[[counts]])
            }
            validObject(.Object)
            return(.Object)
          }
)