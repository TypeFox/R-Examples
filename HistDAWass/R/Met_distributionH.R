#Constructor
#' Wrapper function distributionH
#' 
#' A histogram object can be created also with the function \code{distributionH(...)}, the costructor function for creating an object containing the description of
#' a histogram.
#' 
#' @name distributionH
#' @rdname distributionH-class
#' @export
#' @return A \code{distributionH} object
#' @examples
#' # or using
#' mydist=distributionH(x=c(1,2,3),p=c(0,0.4, 1))
#' @import methods
distributionH=function(x=numeric(0),p=numeric(0)){
  object=new("distributionH",x=x,p=p)
  return(object)
}

#Get

#' Method \code{get.m}: the mean of a distribution 
#' @name get.m
#' @rdname get.m-methods
#' @exportMethod get.m
setGeneric("get.m", function(object) standardGeneric("get.m"))
#' @rdname get.m-methods
#' @aliases get.m,distributionH-method
#' @description This functon return the mean of a \code{distributionH} object.
#' @param object a \code{distributionH} object
#' @return A numeric value
#' @examples
#' D=distributionH(x=c(1,2,3,4),p=c(0,0.2,0.6, 1))
#' get.m(D) #returns the mean of D
setMethod("get.m","distributionH",
          function(object){
            if (!is.null(object)){
              return(object@m)
            }  else
              return(NA)
          }
)

#' Method \code{get.s}: the standard deviation of a distribution
#' @name get.s
#' @rdname get.s-methods
#' @exportMethod get.s
setGeneric("get.s", function(object) standardGeneric("get.s"))

#' @rdname get.s-methods
#' @aliases get.s,distributionH-method
#' @description This functon return the standard deviation of a \code{distributionH} object.
#' @param object a \code{distributionH} object.
#' @return A numeric positive value, the standard deviation.
#' @examples
#' D=distributionH(x=c(1,2,3,4),p=c(0,0.2,0.6, 1))
#' get.s(D) # returns the standard deviation of D
setMethod("get.s","distributionH",
          function(object){
            if (!is.null(object)){
              return(object@s)
            }  else
              return(NA)
          }
)
#' Method \code{get.distr}: show the distribution
#' @name get.distr
#' @rdname get.distr-methods
#' @exportMethod get.distr
setGeneric("get.distr", function(object) standardGeneric("get.distr"))

#' @rdname get.distr-methods
#' @aliases get.distr,distributionH-method
#' @description This functon return the cumulative distribution function of a \code{distributionH} object.
#' @param object a \code{distributionH} object.
#' @return A data frame: the first column contains the domain the second the CDF values.
#' @examples
#' D=distributionH(x=c(1,2,3,4),p=c(0,0.2,0.6, 1))
#' get.distr(D) #a data.frame describing the CDF of D
setMethod("get.distr","distributionH",
          function(object){
            if (!is.null(object)){
              MAT=cbind(object@x,object@p)
              colnames(MAT)=c('x','p')
              return(MAT=as.data.frame(MAT))
            }  else
              return(NA)
          }
)
#' Method \code{get.histo}: show the distribution with bins
#' @name get.histo
#' @rdname get.histo-methods
#' @exportMethod get.histo
setGeneric("get.histo", function(object) standardGeneric("get.histo"))

# Get the distribution of a histogram
#' @rdname get.histo-methods
#' @aliases get.histo,distributionH-method
#' @description This functon return a data.frame describing the histogram of a \code{distributionH} object.
#' @param object a \code{distributionH} object.
#' @return A matrix: the two columns contains the bounds of the histogram the third contains the probablity (or the relative frequency) of the bin.
#' @examples
#' D=distributionH(x=c(1,2,3,4),p=c(0,0.2,0.6, 1))
#' get.histo(D) #returns the histogram representation of D by a data.frame
setMethod("get.histo","distributionH",
          function(object){
            if (!is.null(object)){
              M=get.distr(object)
              MAT=cbind(M$x[1:(length(M$x)-1)], 
                        M$x[2:length(M$x)],
                        M$p[2:length(M$p)]-M$p[1:(length(M$p)-1)])
              colnames(MAT)=c('min.x','max.x', 'p')
              return(MAT=as.data.frame(MAT))
            }  else
              return(NA)
          }
)

# Basic statistics of distributions --------
#' Method \code{meanH}: computes the mean of a distribution
#' @name meanH
#' @rdname meanH-methods
#' @exportMethod meanH
setGeneric("meanH", function(object) standardGeneric("meanH"))
#' Method \code{stdH}: computes the standard deviation of a distribution
#' @name stdH
#' @rdname stdH-methods
#' @exportMethod stdH
setGeneric( "stdH", function(object) standardGeneric("stdH"))
#' Method \code{skewH}: computes the skewness of a distribution
#' @name skewH
#' @rdname skewH-methods
#' @exportMethod skewH
setGeneric("skewH", function(object) standardGeneric("skewH"))
#' Method \code{kurtH}: computes the kurthosis of a distribution
#' @name kurtH
#' @rdname kurtH-methods
#' @exportMethod kurtH
setGeneric("kurtH", function(object) standardGeneric("kurtH"))
#' Method \code{crwtransform}: returns the centers and the radii of bins of a distribution
#' @name crwtransform
#' @rdname crwtransform-methods
#' @exportMethod crwtransform
setGeneric("crwtransform", function(object) standardGeneric("crwtransform"))

#' @rdname meanH-methods
#' @aliases meanH,distributionH-method 
#' @description Mean of a histogram (First moment of the distribution)
#' @param object a \code{distributionH} object
#' @return the mean of the distribution
#' @author Antonio Irpino
#' @keywords distribution
#' @examples
#' 
#' ##---- A mydist distribution ----
#' mydist<-distributionH(x=c(1,2,3,10), p=c(0,0.1,0.5,1))
#' ##---- Compute the mean of mydist ----
#' meanH(mydist) #---> 4.4
#'
setMethod("meanH","distributionH",
          function(object){
            if (!is.null(object@x)||!is.null(object@p)){
              resu=crwtransform(object)
              c=resu[[1]]
              w=resu[[3]]
              m=t(w)%*%c
              m=as.numeric(m)
              return(m)
            }  else
              stop('Something wrong, null domain or cdf')
          }
)
#' @rdname stdH-methods
#' @aliases stdH,distributionH-method 
#' @description Standard deviation of a histogram (i.e., the square root of the centered
#' second moment)
#'  
#' @param object a \code{distributionH} object
#' @return A value for the standard deviation
#' @author Antonio Irpino
#' @keywords distribution
#' @examples
#' 
#' ##---- A mydist distribution ----
#' mydist<-distributionH(x=c(1,2,3,10), p=c(0,0.1,0.5,1))
#' ##---- Compute the standard deviation of mydist ----
#' stdH(mydist) #---> 2.563851
#'
setMethod("stdH","distributionH",
          function(object){
            if (!is.null(object@x)||!is.null(object@p)){
              resu=crwtransform(object)
              c=resu[[1]]
              r=resu[[2]]
              w=resu[[3]]
              std=sqrt(abs(sum(w*c^2+1/3*w*r^2)-(sum(w*c))^2))
              std=as.numeric(std)
              return(std)
            }  else
              stop('Something wrong, null domain or cdf')
          }
)
#' @rdname skewH-methods
#' @aliases skewH,distributionH-method 
#' @description Skewness of a histogram (using the third standardized moment)
#'  
#' @param object a \code{distributionH} object
#' @return A value for the skewness index
#' @author Antonio Irpino
#' @keywords distribution
#' @examples
#' 
#' ##---- A mydist distribution ----
#' mydist<-distributionH(x=c(1,2,3,10), p=c(0,0.1,0.5,1))
#' ##---- Compute the skewness of mydist ----
#' skewH(mydist) #---> -1.186017
#'
setMethod("skewH","distributionH",
          function(object){
            if (!is.null(object@x)||!is.null(object@p)){
              resu=crwtransform(object)
              cs=(resu[[1]]-object@m)/(object@s)
              rs=resu[[2]]/(object@s)
              w=resu[[3]]
              sk=t(w)%*%(cs*(rs^2+cs^2));
              sk=as.numeric(sk)
              return(sk)
            }  else
              stop('Something wrong, null domain or cdf')
          }
)
#' @rdname kurtH-methods
#' @aliases kurtH,distributionH-method 
#' @description Kurtosis of a histogram (using the fourth standardized moment)
#' @param object a \code{distributionH} object
#' @return A value for the kurtosis index, 3 is the kurtosis of a Gaussian
#' distribution
#' @author Antonio Irpino
#' @keywords distribution
#' @examples
#' 
#' ##---- A mydist distribution ----
#' mydist<-distributionH(x=c(1,2,3,10), p=c(0,0.1,0.5,1))
#' ##---- Compute the kurtosis of mydist ----
#' kurtH(mydist) #---> 1.473242
#' 

setMethod("kurtH","distributionH",
          function(object){
            if (!is.null(object@x)||!is.null(object@p)){
              resu=crwtransform(object)
              cs=(resu[[1]]-object@m)/(object@s)
              rs=resu[[2]]/(object@s)
              w=resu[[3]]
              ku=0.2*t(w)%*%(5*cs^4+10*(cs^2)*rs^2+rs^4)
              ku=as.numeric(ku)
              return(ku)
            }  else
              stop('Something wrong, null domain or cdf')
          }
)
#' @rdname crwtransform-methods
#' @aliases crwtransform,distributionH-method 
#' 
#' @description Centers and ranges calculation for bins of a histogram. It is useful for a
#' very fast computation of statistics and methods based on the L2 Wassertein
#' distance between histograms. 
#' @param object a \code{distributionH} object
#' @return A list containing \item{$Centers }{The midpoints of the bins of the
#' histogram} \item{$Radii }{The half-lenghts of the bins of the histogram}
#' \item{$Weights }{The relative frequencies or the probailities associated with
#' each bin (the sum is equal to 1)}
#' @author Antonio Irpino
#' @references Irpino, A., Verde, R., Lechevallier, Y. (2006) \emph{Dynamic
#' clustering of histograms using Wasserstein metric}, In: Proceedings of
#' COMPSTAT 2006, Physica-Verlag, 869-876
#' @keywords distribution
#' @examples
#' 
#' ##---- A mydist distribution ----
#' mydist<-distributionH(x=c(1,2,3,10), p=c(0,0.1,0.5,1))
#' ##---- Compute the cfd value for q=5 (not observed) ----
#' crwtransform(mydist)

setMethod("crwtransform","distributionH",
          function(object){
            if (!is.null(object@x)||!is.null(object@p)){
              resu=list()
              nv=length(object@x)
              c=(object@x[2:nv]+object@x[1:(nv-1)])/2
              r=(object@x[2:nv]-object@x[1:(nv-1)])/2
              w=object@p[2:nv]-object@p[1:(nv-1)]
              resu=list(Centers=c,Radii=r,Weights=w)
              return(resu)
            }  else
              stop('Something wrong, null domain or cdf')
          }
)
# Overloading of the sum of two distribution according to the L2 w --------
#' Method +
#' @name +
#' @aliases +,distributionH,distributionH-method
#' @description the sum of two distribution according to the L2 Wasssertein
#' @param e1 a \code{distributionH} object or a number
#' @param e2 a \code{distributionH} object or a number
#' @return a \code{distributionH} object
#' @export
#' @docType methods
#' @rdname plus-methods
#' 
setMethod("+",
          signature(e1 = "distributionH",e2="distributionH"),
          function (e1, e2) 
          {
            tmp=register(e1,e2)
            x=callGeneric(tmp[[1]]@x,tmp[[2]]@x)
            OBJ_NEW=new("distributionH",x,tmp[[1]]@p,(tmp[[1]]@m+tmp[[2]]@m))
          }
)
#' Method +
#' @name +
#' @aliases +,numeric,distributionH-method
#' @description the sum of a number and a distribution according to the L2 Wasssertein
#' @export
#' @docType methods
#' @rdname plus-methods

setMethod("+",
          signature(e1 = "numeric",e2="distributionH"),
          function (e1, e2) 
          {
            x=callGeneric(rep(e1,length(e2@x)),e2@x)
            OBJ_NEW=new("distributionH",x,e2@p,(e1+e2@m),e2@s)
          }
)
#' Method +
#' @name +
#' @aliases +,distributionH,numeric-method
#' @description the sum of adistribution and a number according to the L2 Wasssertein
#' @export
#' @docType methods
#' @rdname plus-methods
setMethod("+",
          signature(e1 = "distributionH", e2="numeric"),
          function (e1, e2) 
          {
            x=callGeneric(rep(e2,length(e1@x)),e1@x)
            OBJ_NEW=new("distributionH",x,e1@p,(e2+e1@m),e1@s)
          }
)
# Overloading of the difference of two distributions according to the L2 wasserstein --------
#' Method -
#' @name minus
#' @aliases -,distributionH,distributionH-method 
#' @description the difference of two distribution according to the L2 Wasssertein
#' @export
setMethod("-",
          signature(e1 = "distributionH",e2="distributionH"),
          function (e1, e2) 
          {
            tmp=register(e1,e2)
            x=callGeneric(tmp[[1]]@x,tmp[[2]]@x)
            OBJ_NEW=new("distributionH",x,tmp[[1]]@p)
          }
)
#' Method -
#' @name minus
#' @aliases -,numeric,distributionH-method 
#' @param e1 a \code{distributionH} object or a number
#' @param e2 a \code{distributionH} object or a number
#' @description the difference of a number and a distribution according to the L2 Wasssertein
setMethod("-",
          signature(e1 = "numeric",e2="distributionH"),
          function (e1, e2) 
          {
            x=callGeneric(rep(e1,length(e2@x)),e2@x)
            OBJ_NEW=new("distributionH",x,e2@p)
          }
)
#' Method -
#' @name minus
#' @aliases -,distributionH,numeric-method 
#' @description the difference of a distribution and a number according to the L2 Wasssertein
#' @note it may not works properly if the difference is not a distribution
setMethod("-",
          signature(e1 = "distributionH", e2="numeric"),
          function (e1, e2) 
          {
            x=callGeneric(e1@x,rep(e2,length(e1@x)))
            OBJ_NEW=new("distributionH",x,e1@p,(e1@m-e2),e1@s)
          }
)

# Overloading of the product of a number by a distribution according to the L2 w --------
#' Method *
#' @name *-methods
#' @aliases *,distributionH,distributionH-method
#' @description the product of a number and a distribution according to the L2 Wasssertein
#' @param e1 a \code{distributionH} object or a number
#' @param e2 a \code{distributionH} object or a number
#' @export
setMethod("*",
          signature(e1 = "distributionH",e2="distributionH"),
          function (e1, e2) 
          {
            stop('please use dotpW function product between distributions')
          }
)
#' Method *
#' @name *-methods
#' @aliases *,numeric,distributionH-method
#' @description the product of a number and a distribution according to the L2 Wasssertein
setMethod("*",
          signature(e1 = "numeric",e2="distributionH"),
          function (e1, e2) 
          {
            x=callGeneric(rep(e1,length(e2@x)),e2@x)
            OBJ_NEW=new("distributionH",x,e2@p,e1*e2@m,e1*e2@s)
          }
)
#' Method *
#' @name *-methods
#' @aliases *,distributionH,numeric-method
#' @description the product of a number and a distribution according to the L2 Wasssertein
setMethod("*",
          signature(e1 = "distributionH", e2="numeric"),
          function (e1, e2) 
          {
            x=callGeneric(rep(e2,length(e1@x)),e1@x)
            OBJ_NEW=new("distributionH",x,e1@p,e2*e1@m,e2*e1@s)
          }
)

# Utilities for single or couples of distributionH --------------------------------
#' Method \code{checkEmptyBins}
#' @name checkEmptyBins
#' @rdname checkEmptyBins-methods
#' @exportMethod checkEmptyBins
setGeneric("checkEmptyBins", function(object) standardGeneric("checkEmptyBins"))
#' Method \code{compQ}
#' @name compQ
#' @rdname compQ-methods
#' @exportMethod compQ
setGeneric("compQ", function(object,p) standardGeneric("compQ"))
#' Method \code{compP}
#' @name compP
#' @rdname compP-methods
#' @exportMethod compP
setGeneric("compP", function(object,q) standardGeneric("compP"))
#' Method \code{register}
#' @name register
#' @rdname register-methods
#' @exportMethod register
setGeneric("register", function(object1,object2) standardGeneric("register"))

#' @rdname register-methods
#' @aliases register,distributionH-method 
#' @description Given two \code{distributionH} objects, it returns two equivalent distributions such that 
#' they share the same cdf values. This function is useful for computing basic statistics.
#' 
#' @param object1 A \code{distributionH} object
#' @param object2 A \code{distributionH} object
#' @return The two \code{distributionH} objects in input sharing the same cdf (the \code{p}
#' slot)
#' @author Antonio Irpino
#' @references Irpino, A., Lechevallier, Y. and Verde, R. (2006): \emph{Dynamic
#' clustering of histograms using Wasserstein metric} In: Rizzi, A., Vichi, M.
#' (eds.) COMPSTAT 2006. Physica-Verlag, Berlin, 869-876.\cr Irpino, A.,Verde,
#' R. (2006): \emph{A new Wasserstein based distance for the hierarchical
#' clustering of histogram symbolic data} In: Batanjeli, V., Bock, H.H.,
#' Ferligoj, A., Ziberna, A. (eds.) Data Science and Classification, IFCS 2006.
#' Springer, Berlin, 185-192.
#' @keywords distribution
#' @examples
#' 
#' ##---- initialize two distributionH objects mydist1 and mydist2
#'  mydist1=distributionH(c(1,2,3),c(0, 0.4, 1))
#'  mydist2=distributionH(c(7,8,10,15),c(0, 0.2, 0.7, 1))
#'  ## register the two distributions
#'  regDist=register(mydist1,mydist2)
#'  
#' ## OUTPUT:
#' ## regDist$[[1]]
#' ## An object of class "distributionH"
#' ## Slot "x": [1] 1.0 1.5 2.0 2.5 3.0
#' ## Slot "p": [1] 0.0 0.2 0.4 0.7 1.0
#' ## ...
#' ## regDist$[[2]] 
#' ## An object of class "distributionH"
#' ## Slot "x": [1] 7.0 8.0 8.8 10.0 15.0
#' ## Slot "p": [1] 0.0 0.2 0.4  0.7  1.0
#' ## ...
#' 
setMethod(f="register",signature=c(object1="distributionH",object2="distributionH"),
          function(object1,object2){
            if (!identical(object1@p,object2@p)){
              commoncdf=sort(unique(c(object1@p,object2@p)))
              nr=length(commoncdf)
              result=matrix(0,nr,3)
              result[,3]=commoncdf
              result[1,1:2]=c(object1@x[1],object2@x[1])
              result[nr,1:2]=c(object1@x[length(object1@x)],object2@x[length(object2@x)])
              if (nr>2) {
                for (i in c(2:(nr-1))){
                  result[i,1] =compQ(object1,result[i,3])
                  result[i,2] =compQ(object2,result[i,3])
                }
              } 
              o1=new("distributionH",as.vector(result[,1]),as.vector(commoncdf))
              o2=new("distributionH",as.vector(result[,2]),as.vector(commoncdf))
              return(c(o1,o2))
            }
            else{
              return(c(object1,object2))
            }
            
          }
)
#' @rdname checkEmptyBins-methods
#' @aliases checkEmptyBins,distributionH-method 
#' @description The method checking for empty bins in a distribution, i.e. if two cdf consecutive
#' values are equal. In that case a probability value of \code{1e-7} is
#' assigned to the empty bin and the cdf is recomputed. This methods is useful
#' for numerical reasons.
#' 
#' 
#' @param object a \code{distributionH} object
#' @return A \code{distributionH} object without empty bins
#' @author Antonio Irpino
#' @keywords distribution
#' @examples
#' 
#' ##---- A mydist distribution with an empty bin i.e. two consecutive values of p are equal----
#' mydist<-distributionH(x=c(1,2,3,10), p=c(0,0.5,0.5,1))
#' ##---- Checks for empty byns and returns the newdist object without empty bins ----
#' newdist<-checkEmptyBins(mydist)
#'
setMethod(f="checkEmptyBins",signature="distributionH",
            function(object){
            w=object@p[2:length(object@p)]-object@p[1:(length(object@p)-1)]
            TOL=1e-08
            if (length(which(w<=TOL))){
              w[which(w<TOL)]=10*TOL
              object@p=c(0,cumsum(w))/sum(w)
              
            }
            return(object)
          }
)
#' @rdname compQ-methods
#' @aliases compQ,distributionH-method 
#' @description Compute the quantile value of a histogram for a given probability. 
#' 
#' 
#' @param object an object of \env{distributionH} class
#' @param p a number between 0 and 1
#' @return \deqn{y= F^{-1}(p)=Q(p)} A number that is the quantile of the passed
#' histogram \env{object} at level \env{p}.
#' @author Antonio Irpino
#' @examples
#' 
#' ##---- A mydist distribution ----
#' mydist<-distributionH(x=c(1,2,3,10), p=c(0,0.1,0.5,1))
#' ##---- Compute the quantile of mydist for different values of p ----
#' y<-compQ(mydist,0.5) #the median
#' y<-compQ(mydist,0) #the minimum
#' y<-compQ(mydist,1) #the maximum
#' y<-compQ(mydist,0.25) #the first quartile
#' y<-compQ(mydist,0.9) #the ninth decile
#' 
setMethod(f="compQ",signature=c(object="distributionH",p="numeric"),
          function(object,p){
            # %Computes the p-th quantile p=[0,1] of the distribution o1
            # %INPUT  - p  a value in [0,1]
            # %           - o1 a distribution
            # %OUTPUT-  res the computed quantile
            # %example:q=compQ(o1,0.5) returns the median of the
            # % distribution o1
            
            #Check for errors
            if (p<0 || p>1) stop("p must be a value between 0 and 1")
            
            if (p<=0) return(q=object@x[1])
            if (p>=1) return (q=object@x[length(object@x)])
            
            ini=max(object@x[object@p<=p])
            pos1=which.max(object@x[object@p<=p])
            pos2=pos1+1
            fin=object@x[pos2];
            if (ini==fin){
              q=ini
            }
            else{
              q=ini+(object@x[pos2]-object@x[pos1])*(p-object@p[pos1])/(object@p[pos2]-object@p[pos1])
            }
            return(q)
          }
)
#' @rdname compP-methods
#' @aliases compP,distributionH-method 
#' @description Compute the cdf probability at a given value for a histogram
#' 
#' @param object is an object of \env{distributionH} class
#' @param q is a numeric value
#' @return Returns a value between 0 and 1.
#' @keywords distribution
#' @examples
#' 
#' ##---- A mydist distribution ----
#' mydist<-distributionH(x=c(1,2,3,10), p=c(0,0.1,0.5,1))
#' ##---- Compute the cfd value for q=5 (not observed) ----
#' p<-compP(mydist,5)
#'
setMethod(f="compP",signature=c(object="distributionH",q="numeric"),
          function(object,q){
            domain=object@x
            cumul=object@p
            if (q<=domain[1]) return(p=0)
            if (q>=domain[length(domain)]) return(p=1)
            
            ini=max(cumul[domain<=q])
            pos1=which.max(cumul[domain<=q])
            pos2=pos1+1
            fin=cumul[pos2]
            
            if (ini==fin){
              p=fin}
            else{
              p=ini+(cumul[pos2]-cumul[pos1])*(q-domain[pos1])/(domain[pos2]-domain[pos1])}
            return(as.numeric(p))
          }
)
# L2 Wasserstein distance between two distributions and related results ----
#' Method \code{WassSqDistH}
#' @name WassSqDistH
#' @rdname WassSqDistH-methods
#' @exportMethod WassSqDistH
setGeneric("WassSqDistH", function(object1,object2,...) standardGeneric("WassSqDistH"))# Wasserstein distance between two distributions
#' Method \code{rQQ}
#' @name rQQ
#' @rdname rQQ-methods
#' @exportMethod rQQ
setGeneric("rQQ", function(e1,e2) standardGeneric("rQQ"))# Quantile-Quantile correlation between twi distributions

#' @rdname WassSqDistH-methods
#' @aliases WassSqDistH,distributionH-method 
#' @description Computes the squared L2 Wasserstein distance between two \code{distributionH} objects.
#' @param object1 is an object of \env{distributionH} class
#' @param object2 is an object of \env{distributionH} class
#' @param ... optional parameters
#' @param details (optional, default=FALSE) is a logical value, if TRUE returns the decomposition of the distance
#' @return 
#' If \code{details=FALSE}, the function returns the squared L2 Wasserstein distance.\cr
#' If \code{details=TRUE}, the function returns list containing the squared distance, its 
#' decomposition in three parts (position, size and shape) and the correlation coefficient between the quantile functions.
#' @references
#' Irpino, A. and Romano, E. (2007): \emph{Optimal histogram representation of large data sets: 
#' Fisher vs piecewise linear approximations}. RNTI E-9, 99-110.\cr
#' Irpino, A., Verde, R. (2015) \emph{Basic
#' statistics for distributional symbolic variables: a new metric-based
#' approach} Advances in Data Analysis and Classification, DOI 10.1007/s11634-014-0176-4
#' @keywords distribution
#' @examples
#' ##---- create two distributionH objects ----
#'  mydist1=distributionH(x=c(1,2,3),p=c(0, 0.4, 1))
#'  mydist2=distributionH(x=c(7,8,10,15),p=c(0, 0.2, 0.7, 1))
#'# -- compute the squared L2 Waaserstein distance
#' WassSqDistH(mydist1,mydist2)
#'# -- compute the squared L2 Waaserstein distance with details
#' WassSqDistH(mydist1,mydist2,details=TRUE)
setMethod(f="WassSqDistH",signature=c(object1="distributionH",object2="distributionH"),
          #Computes the L2 Wasserstein squared distance between two distributions
          #INPUT: object1 and object2 - two distributionH objects
          #OUTPUT: A list containing the distance and its decomposition in three parts (position, size and shape)
          function(object1=object1,object2=object2,details=FALSE){
            tmp=register(object1,object2)
            nv=length(tmp[[1]]@p)
            w=tmp[[1]]@p[2:nv]-tmp[[1]]@p[1:(nv-1)]
            c1=0.5*(tmp[[1]]@x[2:nv]+tmp[[1]]@x[1:(nv-1)])
            c2=0.5*(tmp[[2]]@x[2:nv]+tmp[[2]]@x[1:(nv-1)])
            r1=0.5*(tmp[[1]]@x[2:nv]-tmp[[1]]@x[1:(nv-1)])
            r2=0.5*(tmp[[2]]@x[2:nv]-tmp[[2]]@x[1:(nv-1)])
            D=t(w)%*%((c1-c2)^2+1/3*(r1-r2)^2)
            if (details) {
              DC=(object1@m-object2@m)^2
              DS=(object1@s-object2@s)^2
              if(abs(D-DC-DS)<1e-10) DR=0 else DR=abs(D-DC-DS)
              
              rho=1-abs(D-DC-DS)/(2*object1@s*object2@s)
              if (rho<0) rho=0
              resu=c(D,DC,DS,DR,rho)
              names(resu)=c("SQ_W_dist", "POSITION", "SIZE", "SHAPE", "rQQ")
              return(resu)}
            else return(as.numeric(D))
          }
)

#' Method \code{dotpW}
#' @name dotpW
#' @rdname dotpW-methods
#' @exportMethod dotpW
setGeneric("dotpW", function(e1,e2) standardGeneric("dotpW"))#dotproduct from L2 Wasserstein
#' @rdname dotpW-methods
#' @aliases dotpW,distributionH-method 
#' @description The dot product of two distributions inducing the L2 Wasserstein metric
#' 
#' @param e1 a \code{distributionH} object or a number
#' @param e2 a \code{distributionH} object or a number
#' @return A numeric value
#' @author Antonio Irpino
#' @references  Irpino, A., Verde, R. (2015) \emph{Basic
#' statistics for distributional symbolic variables: a new metric-based
#' approach} Advances in Data Analysis and Classification, DOI 10.1007/s11634-014-0176-4
#' @keywords distribution
#' @examples
#' 
#' ## let's define two distributionH objects
#' mydist1<-distributionH(x=c(1,2,3,10), p=c(0,0.1,0.5,1))
#' mydist2<-distributionH(x=c(5,7,15), p=c(0,0.7,1))
#' 
#' ## the dot product between the distributions
#' dotpW(mydist1,mydist2) #---> 39.51429
#' 
#' ## the dot product between a distribution and a numeric
#' dotpW(mydist1,3)  #---> 13.2
#' dotpW(3,mydist1)  #---> 13.2
#' 
#' 

setMethod("dotpW",
          signature(e1 = "distributionH",e2="distributionH"),
          definition = function (e1, e2) 
          {
            tmp=register(e1,e2)
            w=tmp[[1]]@p[2:length(tmp[[1]]@p)]-tmp[[1]]@p[1:(length(tmp[[1]]@p)-1)]
            c1=(tmp[[1]]@x[2:length(tmp[[1]]@p)]+tmp[[1]]@x[1:(length(tmp[[1]]@p)-1)])/2
            c2=(tmp[[2]]@x[2:length(tmp[[1]]@p)]+tmp[[2]]@x[1:(length(tmp[[1]]@p)-1)])/2
            r1=(tmp[[1]]@x[2:length(tmp[[1]]@p)]-tmp[[1]]@x[1:(length(tmp[[1]]@p)-1)])/2
            r2=(tmp[[2]]@x[2:length(tmp[[1]]@p)]-tmp[[2]]@x[1:(length(tmp[[1]]@p)-1)])/2              
            dprod=sum(w*(c1*c2+1/3*r1*r2))
            return(dprod)
          }
)
#' @rdname dotpW-methods
#' @aliases dotpW,distributionH-method 
#' @description The dot product of a number (considered as an impulse distribution function) and a distribution
setMethod("dotpW",
          signature(e1 = "numeric",e2="distributionH"),
          function (e1, e2) 
          {
            dprod=e1*e2@m
            return(dprod)
          }
)
#' @rdname dotpW-methods
#' @aliases dotpW,distributionH-method 
#' @description The dot product of a distribution and a number (considered as an impulse distribution function).
setMethod("dotpW",
          signature(e1 = "distributionH", e2="numeric"),
          function (e1, e2) 
          {
            dprod=e2*e1@m
            return(dprod)
          }
)
#' @rdname rQQ-methods
#' @aliases rQQ,distributionH-method 
#' @description Quantile-Quantile correlation between two distributions
#' @param e1  A \code{distributionH} object
#' @param e2  A \code{distributionH} object
#' @return Pearson correlation index between quantiles
#' @author Antonio Irpino
#' @references  Irpino, A., Verde, R. (2015) \emph{Basic
#' statistics for distributional symbolic variables: a new metric-based
#' approach} Advances in Data Analysis and Classification, DOI 10.1007/s11634-014-0176-4
#' @examples
#' 
#' ##---- initialize two distributionH object mydist1 and mydist2
#'  mydist1<-distributionH(x=c(1,2,3),p=c(0, 0.4, 1))
#'  mydist2<-distributionH(x=c(7,8,10,15),p=c(0, 0.2, 0.7, 1))
#'  ## computes the rQQ
#'  rQQ(mydist1,mydist2)
#'  ## OUTPUT 0.916894
#' 
#' @export rQQ
setMethod("rQQ",
          signature(e1 = "distributionH",e2="distributionH"),
          definition = function (e1, e2) 
          {
            rQQ=(dotpW(e1,e2)-e1@m*e2@m)/(e1@s*e2@s)
            return(rQQ)
          }
)
## ---- Show overridding for distributionH and MatH ----
#' Method show for distributionH
#' @name show
#' @rdname show-distributionH-methods
#' @docType methods
#' @aliases show,distributionH-method 
#' @description An overriding show function for a \code{distributionH} object. The function returns a representation 
#' of the histogram, if the number of bins is high the central part of the histogram is truncated. 
#' @param object a \code{distributionH} object
#' @examples
#' ##---- initialize a distributionH
#'  mydist<-distributionH(x=c(7,8,10,15),p=c(0, 0.2, 0.7, 1))
#'  # show the histogram
#'  mydist
setMethod("show",
          signature(object="distributionH"),
          definition=function(object){
            if(length(object@p)>2){
            mymat=matrix(0,length(object@p)-1,2)
            if (length(object@p)>11){
              cat("Output shows the first five and the last five bins due to eccesive length \n")
              mymat=matrix(0,12,2)
              mymat[1,1]="X"
              mymat[1,2]="p"
              count=0
              for (i in 2:6){
                count=count+1
                mymat[count+1,1]=paste("[",format(object@x[(i-1)],digits=5),"-",format(object@x[i],digits=5),")",sep="") 
                mymat[count+1,2]=paste(format(object@p[i]-object@p[i-1],digits=4))
              }
              count=count+1
              mymat[count+1,1]=paste("...") 
              mymat[count+1,2]=paste("...")
              for (i in (length(object@p)-4):length(object@p)){
                count=count+1
                mymat[count+1,1]=paste("[",format(object@x[(i-1)],digits=5)," ; ",format(object@x[i],digits=5),")",sep="")  
                mymat[count+1,2]=paste(format(object@p[i]-object@p[i-1],digits=4))
              }
              rownames(mymat)=c("", paste("Bin",1:5,sep="_"),"...",
                                paste("Bin",(length(object@p)-5):(length(object@p)-1),sep="_"))
              write.table(format(mymat, justify="right",digits=5),
                          row.names=T, col.names=F, quote=F)  
            }
            else{
              mymat=matrix(0,length(object@p),2)
              mymat[1,1]="X"
              mymat[1,2]="p"
              for (i in 2:(length(object@p)-1)){
                mymat[i,1]=paste("[ ",format(object@x[(i-1)],digits=5)," ; ",format(object@x[i],digits=5)," )",sep="")  
                mymat[i,2]=paste(format(object@p[i]-object@p[i-1],digits=4))
              }
              mymat[length(object@p),1]=paste("[ ",format(object@x[(length(object@p)-1)],digits=5),
                                              " ; ",format(object@x[length(object@p)],digits=5)," ]",sep="")  
              mymat[length(object@p),2]=paste(format(object@p[length(object@p)]-object@p[length(object@p)-1],digits=4))
              
              rownames(mymat)=c(" ", paste("Bin",1:(length(object@p)-1),sep="_"))
              write.table(format(mymat, justify="right"),
                          row.names=T, col.names=F, quote=F)
            }
            cat(paste("\n mean = ",object@m, "  std  = ",object@s,"\n "))
            
            } else (cat("Empty distributionH\n"))
          }
)




## --- Plot  for distributionH  ----
#' plot for a distributionH object
#' @name plot-distributionH
#' @docType methods
#' @aliases plot,distributionH-method
#' @description A plot function for a \code{distributionH} object. The function returns a representation 
#' of the histogram.
#' @param x  a \code{distributionH} object
#' @param ... other optional parameters
#' @param type (optional) a string describing the type of plot, default="HISTO".\cr Other allowed types are 
#' \cr"CDF"=Cumulative distribution function, \cr"QF"= quantile function, \cr"DENS"=a density approximation, 
#' \cr"HBOXPLOT"=horizontal boxplot, \cr"VBOXPLOT"= vertical boxplot,
#' @param col (optional) a string the color of the plot, default="green".
#' @param border (optional) a string the color of the border of the plot, default="black".
#' @examples
#' ##---- initialize a distributionH
#'  mydist<-distributionH(x=c(7,8,10,15),p=c(0, 0.2, 0.7, 1))
#'  # show the histogram
#'  plot(mydist) #plots mydist
#'  plot(mydist, type="HISTO", col="red", border="blue") #plots mydist
#'  plot(mydist, type="DENS", col="red", border="blue") #plots a density approximation for mydist
#'  plot(mydist, type="HBOXPLOT") #plots a horizontal boxplot for mydist
#'  plot(mydist, type="VBOXPLOT") #plots a vertical boxplot for mydist
#'  plot(mydist, type="CDF") #plots the cumulative distribution function of mydist
#'  plot(mydist, type="QF") #plots the quantile function of mydist
#' @importFrom utils write.table
#' @export
setMethod("plot",
          signature(x = "distributionH" ),
          function (x,  type="HISTO",col="green",border="black") 
          { 
            plot.gg(x,  type=type,col=col,border=border)
          }
)