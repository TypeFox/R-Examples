#' Wrapper function of \code{MatH} class
#' 
#' This function create a matrix of histogram data, i.e. a \code{MatH}
#' object
#' 
#' @name MatH
#' @rdname MatH-class
#' @export
#' @param x (optional, default= an empty \code{distributionH} object) a list of
#' \code{distributionH} objects
#' @param nrows (optional, default=1)an integer, the number of rows.
#' @param ncols (optional, default=1) an integer, the number of columns (aka
#' variables).
#' @param rownames (optional, default=NULL) a list of strings containing the
#' names of the rows.
#' @param varnames (optional, default=NULL) a list of strings containing the
#' names of the columns (aka variables).
#' @param by.row (optional, default=FALSE) a logical value, TRUE the matrix is
#' row wise filled, FALSE the matrix is filled column wise.
#' @return A \code{matH} object
#' @examples
#' 
#' #bulding an empty 10 by 4 matrix of histograms
#' MAT=MatH(nrows=10,ncols=4)

MatH=function(x=list(new('distributionH')), nrows=1, ncols=1,rownames=NULL,varnames=NULL, by.row=FALSE ){
  MAT=new('MatH', 
          nrows=nrows,
          ncols=ncols,
          ListOfDist=x, 
          names.rows=rownames, 
          names.cols=varnames, by.row=by.row)
  return(MAT)
}

# overriding of "[" operator for MatH object ----
 #' extract from a MatH Method [
 #' @name [
 #' @rdname extract-methods
 #' @aliases [,MatH,ANY,ANY,ANY-method
 #' [,MatH-method
  #' @description This method overrides the "[" operator for a  \code{matH} object.
 #' @param x a \code{matH} object
 #' @param i  a set of integer values identifying the rows
 #' @param j  a set of integer values identifying the columns
 #' @param ... not useful
 #' @param drop a logical value inherited from the basic method "[" but not used (default=TRUE)
 #' @return A \code{matH} object
 #' @examples
 #' D=BLOOD #the BLOOD dataset
 #' SUB_D=BLOOD[c(1,2,5),c(1,2)]
 #' @importFrom stats variable.names
 #' @export
setMethod("[",
          signature(x = "MatH"),
          function (x, i, j, ..., drop=TRUE) 
          {
            if (missing(i) &&  missing(j)) {
              i=c(1:nrow(x@M))
              j=c(1:ncol(x@M))
            }
            else{
              if (missing(i)) i=c(1:nrow(x@M))
              if (missing(j)) j=c(1:ncol(x@M))
            }
            
            MAT=new("MatH",length(i),length(j))
            for (r in 1:length(i)){
              for (c in 1:length(j)){
                MAT@M[r,c][[1]]=x@M[i[r],j[c]][[1]]
              }
            }
            
            rownames(MAT@M)=row.names(x@M)[i]
            colnames(MAT@M)=variable.names(x@M)[j]
            return(MAT)
          }
)
# methods for getting information from a MatH
#' Method get.MatH.nrows
#' @name get.MatH.nrows
#' @rdname get.MatH.nrows-methods
#' @exportMethod get.MatH.nrows
setGeneric("get.MatH.nrows",function(object) standardGeneric("get.MatH.nrows"))
#' @rdname get.MatH.nrows-methods
#' @aliases get.MatH.nrows,MatH-method
#' @description It returns the number of rows of a \code{MatH} object
#' @param object  a \code{MatH} object
#' @return An integer, the number of rows.
setMethod(f="get.MatH.nrows",signature=c(object="MatH"),
          function(object){
            return(nrow(object@M))
          }
)
#' Method get.MatH.ncols
#' @name get.MatH.ncols
#' @rdname get.MatH.ncols-methods
#' @exportMethod get.MatH.ncols

setGeneric("get.MatH.ncols",function(object) standardGeneric("get.MatH.ncols"))
#' @rdname get.MatH.ncols-methods
#' @aliases get.MatH.ncols,MatH-method
#' @description It returns the number of columns of a \code{MatH} object
#' @param object  a \code{MatH} object
#' @return An integer, the number of columns.
setMethod(f="get.MatH.ncols",signature=c(object="MatH"),
          function(object){
            return(ncol(object@M))
          }
)
#' Method get.MatH.rownames
#' @name get.MatH.rownames
#' @rdname get.MatH.rownames-methods
#' @exportMethod get.MatH.rownames
setGeneric("get.MatH.rownames",function(object) standardGeneric("get.MatH.rownames"))
#' @rdname get.MatH.rownames-methods
#' @aliases get.MatH.rownames,MatH-method
#' @description It returns the labels of the rows of a \code{MatH} object
#' @param object  a \code{MatH} object
#' @return A vector of char, the label of the rows.
setMethod(f="get.MatH.rownames",signature=c(object="MatH"),
          function(object){
            return(rownames(object@M))
          }
)
#' Method get.MatH.varnames
#' @name get.MatH.varnames
#' @rdname get.MatH.varnames-methods
#' @exportMethod get.MatH.varnames
setGeneric("get.MatH.varnames",function(object) standardGeneric("get.MatH.varnames"))
#' @rdname get.MatH.varnames-methods
#' @aliases get.MatH.varnames,MatH-method
#' @description It returns the labels of the columns, or the names of the variables, of a \code{MatH} object
#' @param object  a \code{MatH} object
#' @return A vector of char, the labels of the columns, or the names of the variables.
setMethod(f="get.MatH.varnames",signature=c(object="MatH"),
          function(object){
            return(colnames(object@M))
          }
)
#' Method get.MatH.main.info
#' @name get.MatH.main.info
#' @rdname get.MatH.main.info-methods
#' @exportMethod get.MatH.varnames
setGeneric("get.MatH.main.info",function(object) standardGeneric("get.MatH.main.info"))
#' @rdname get.MatH.main.info-methods
#' @aliases get.MatH.main.info,MatH-method
#' @description It returns the number of rows, of columns the labels of rows and columns of a \code{MatH} object.
#' @param object  a \code{MatH} object
#' @return A list of char, the labels of the columns, or the names of the variables.
#' @slot nrows - the number of rows
#' @slot ncols - the number of columns
#' @slot rownames - a vector of char, the names of rows
#' @slot varnames - a vector of char, the names of columns
#' 
setMethod(f="get.MatH.main.info",signature=c(object="MatH"),
          function(object){
              return(list(nrows=get.MatH.nrows(object), ncols=get.MatH.ncols(object),
                        rownames=get.MatH.rownames(object),varnames=get.MatH.varnames(object)))
          }
)
#' Method get.MatH.stats
#' @name get.MatH.stats
#' @rdname get.MatH.stats-methods
#' @exportMethod get.MatH.stats
setGeneric("get.MatH.stats",function(object,...) standardGeneric("get.MatH.stats"))

#' @rdname get.MatH.stats-methods
#' @aliases get.MatH.stats,MatH-method
#' @description It returns statistics for each distribution contained in a \code{MatH} object.
#' @param object  a \code{MatH} object
#' @param ... a set of other parameters
#' @param stat (optional) a string containing the required statistic. Default='mean'\cr
#' - \code{stat='mean'} - for computing the mean of each histogram\cr
#' - \code{stat='median'} - for computing the median of each histogram\cr
#' - \code{stat='min'} - for computing the minimum of each histogram\cr
#' - \code{stat='max'} - for computing the maximum of each histogram\cr
#' - \code{stat='std'} - for computing the standard deviatio of each histogram\cr
#' - \code{stat='skewness'} - for computing the skewness of each histogram\cr
#' - \code{stat='kurtosis'} - for computing the kurtosis of each histogram\cr
#' - \code{stat='quantile'} - for computing the quantile ot level \code{prob} of each histogram\cr
#' @param prob (optional)a number between 0 and 1 for computing the value once choosen the \code{'quantile'} option for \code{stat}.
#' @return A list 
#' @slot stat - the chosen statistic
#' @slot prob - level of probability if stat='quantile'
#' @slot MAT - a matrix of values
#' @examples
#' get.MatH.stats(BLOOD) # the means of the distributions in BLOOD dataset
#' get.MatH.stats(BLOOD,stat='median') # the medians of the distributions in BLOOD dataset
#' get.MatH.stats(BLOOD,stat='quantile', prob=0.5) #the same as median
#' get.MatH.stats(BLOOD,stat='min') # minima of the distributions in BLOOD dataset
#' get.MatH.stats(BLOOD,stat='quantile', prob=0) #the same as min
#' get.MatH.stats(BLOOD,stat='max') # maxima of the distributions in BLOOD dataset
#' get.MatH.stats(BLOOD,stat='quantile', prob=1) #the same as max
#' get.MatH.stats(BLOOD,stat='std') # standard deviations of the distributions in BLOOD dataset
#' get.MatH.stats(BLOOD,stat='skewness') #skewness indices of the distributions in BLOOD dataset
#' get.MatH.stats(BLOOD,stat='kurtosis') #kurtosis indices of the distributions in BLOOD dataset
#' get.MatH.stats(BLOOD,stat='quantile',prob=0.05) 
#' #the fifth percentiles of distributions in BLOOD dataset

setMethod(f="get.MatH.stats",signature=c(object="MatH"),
          function(object, stat='mean', prob=0.5){
            r=get.MatH.nrows(object)
            c=get.MatH.ncols(object)
            MAT=matrix(NA,get.MatH.nrows(object),get.MatH.ncols(object))
            rownames(MAT)=get.MatH.rownames(object)
            colnames(MAT)=get.MatH.varnames(object)
            for (i in 1:r){
              for (j in 1:c){
                if (length(object@M[i,j][[1]]@x)>0){
                  if (stat=='mean'){
                    MAT[i,j]=object@M[i,j][[1]]@m
                  }
                  if (stat=='std'){
                    MAT[i,j]=object@M[i,j][[1]]@s
                  }
                  if (stat=='skewness'){
                    MAT[i,j]=skewH(object@M[i,j][[1]])
                  }
                  if (stat=='kurtosis'){
                    MAT[i,j]=kurtH(object@M[i,j][[1]])
                  }
                  if (stat=='median'){
                    MAT[i,j]=compQ(object = object@M[i,j][[1]],p=0.5)
                  }
                  if (stat=='quantile'){
                    MAT[i,j]=compQ(object = object@M[i,j][[1]],p=prob)
                  }
                  if (stat=='min'){
                    MAT[i,j]=compQ(object = object@M[i,j][[1]],p=0)
                  }
                  if (stat=='max'){
                    MAT[i,j]=compQ(object = object@M[i,j][[1]],p=1)
                  }
                  
                }
              }
            }
            if (stat=='quantile'){
              return(list(stat=stat, prob=prob, mat=MAT))
            } else{
            return(list(stat=stat, mat=MAT))
            }
          }
)




# methods for collating by row or by column two MatHs ----
#' Method WH.bind.row
#' @name WH.bind.row
#' @rdname WH.bind.row-methods
#' @exportMethod WH.bind.row
setGeneric("WH.bind.row",function(object1,object2) standardGeneric("WH.bind.row"))#
#' Method WH.bind.col
#' @name WH.bind.col
#' @rdname WH.bind.col-methods
#' @exportMethod WH.bind.col
setGeneric("WH.bind.col",function(object1,object2) standardGeneric("WH.bind.col"))#
#' Method WH.bind
#' @name WH.bind
#' @rdname WH.bind-methods
#' @exportMethod WH.bind
setGeneric("WH.bind",function(object1,object2,byrow) standardGeneric("WH.bind"))#
#' @rdname WH.bind.row-methods
#' @aliases WH.bind.row,MatH-method
#' @description It attaches two \code{MatH} objects with the same columns by row.
#' @param object1  a \code{MatH} object
#' @param object2  a \code{MatH} object
#' @return a \code{MatH} object,  
#' @examples
#' M1<-BLOOD[1:3,]
#' M2<-BLOOD[5:8,]
#' MAT<-WH.bind.row(M1,M2)

setMethod(f="WH.bind.row",signature=c(object1="MatH",object2="MatH"),
          function(object1,object2){
            ncol1=ncol(object1@M)
            ncol2=ncol(object2@M)
            nrow1=nrow(object1@M)
            nrow2=nrow(object2@M)
            if (ncol1!=ncol2){stop("The two matrix must have the same number of columns")}
            NewMat=new("MatH", nrows=nrow1+nrow2,ncols=ncol1)
            NewMat@M=rbind(object1@M,object2@M)
            return(NewMat)
          }
          )
#' @rdname WH.bind.col-methods
#' @aliases WH.bind.col,MatH-method
#' @description It attaches two \code{MatH} objects with the same rows by colums.
#' @param object1  a \code{MatH} object
#' @param object2  a \code{MatH} object
#' @return a \code{MatH} object,  
#' @examples
#' M1<-BLOOD[1:10,1]
#' M2<-BLOOD[1:10,3]
#' MAT<-WH.bind.col(M1,M2)
setMethod(f="WH.bind.col",signature=c(object1="MatH",object2="MatH"),
          function(object1,object2){
            ncol1=ncol(object1@M)
            ncol2=ncol(object2@M)
            nrow1=nrow(object1@M)
            nrow2=nrow(object2@M)
            if (nrow1!=nrow2){stop("The two matrix must have the same number of rows")}
            NewMat=new("MatH", nrows=nrow1,ncols=ncol1+ncol2)
            NewMat@M=cbind(object1@M,object2@M)
            return(NewMat)
          }
)
#' @rdname WH.bind-methods
#' @aliases WH.bind,MatH-method
#' @description It attaches two \code{MatH} objects with the same columns by row, or the same rows by colum.
#' @param object1  a \code{MatH} object
#' @param object2  a \code{MatH} object
#' @param byrow  a logical value (default=TRUE) attaches the objects by row
#' @return a \code{MatH} object,  
#' @examples
#' # binding by row 
#' M1<-BLOOD[1:10,1]
#' M2<-BLOOD[1:10,3]
#' MAT<-WH.bind(M1,M2, byrow=TRUE)
#' # binding by col
#' M1<-BLOOD[1:10,1]
#' M2<-BLOOD[1:10,3]
#' MAT<-WH.bind(M1,M2, byrow=FALSE)
#' @seealso \code{\link{WH.bind.row}} for binding by row, \code{\link{WH.bind.col}} for binding by column 
setMethod(f="WH.bind",signature=c(object1="MatH",object2="MatH"),
          function(object1,object2, byrow=TRUE){
            ncol1=ncol(object1@M)
            ncol2=ncol(object2@M)
            nrow1=nrow(object1@M)
            nrow2=nrow(object2@M)
            if (byrow==TRUE){
              NewMat=WH.bind.row(object1,object2)
            }
            else{
              NewMat=WH.bind.col(object1,object2)
            }
            return(NewMat)
          }
)

# methods for MatH based on the L2 Wasserstein distance between distributions ----
#' Method WH.mat.sum
#' @name WH.mat.sum
#' @rdname WH.mat.sum-methods
#' @exportMethod WH.mat.sum
setGeneric("WH.mat.sum",function(object1,object2) standardGeneric("WH.mat.sum"))#ok matrix sum
#' Method WH.mat.prod
#' @name WH.mat.prod
#' @rdname WH.mat.prod-methods
#' @exportMethod WH.mat.prod
setGeneric("WH.mat.prod",function(object1,object2,...) standardGeneric("WH.mat.prod"))#ok matrix product
#' @rdname WH.mat.sum-methods
#' @aliases WH.mat.sum,MatH-method
#' @description It sums two \code{MatH} objects, i.e. two matrices of distributions, 
#' by summing the quantile functions of histograms. This sum is consistent with 
#' a set of distributions equipped with a L2 wasserstein metric.
#' @param object1  a \code{MatH} object
#' @param object2  a \code{MatH} object
#' @return a \code{MatH} object,  
#' @examples
#' # binding by row 
#' M1<-BLOOD[1:5,]
#' M2<-BLOOD[6:10,]
#' MAT<-WH.mat.sum(M1,M2)
setMethod(f="WH.mat.sum",signature=c(object1="MatH",object2="MatH"),
          #sums two MatH, i.e.  two matrices of distributionsH
          #INPUT: 
          #OUTPUT: 
          
          function(object1,object2){
            nrows1=nrow(object1@M)
            ncols1=ncol(object1@M)
            nrows2=nrow(object1@M)
            ncols2=ncol(object1@M)
            if (!identical(dim(object1@M),dim(object2@M))){
              stop("the two matrices must be of the same dimension")}
            else
            {
              MATS=object1
              TMP=new("MatH",1,2)
              for (r in 1:nrows1){
                for (c in 1:ncols1){
                  TMP@M[1,1][[1]]=object1@M[r,c][[1]]
                  TMP@M[1,2][[1]]=object2@M[r,c][[1]]
                  TMP=registerMH(TMP)
                  MATS@M[r,c][[1]]=new("distributionH",
                                       (TMP@M[1,1][[1]]@x+TMP@M[1,2][[1]]@x),
                                       TMP@M[1,1][[1]]@p,
                                       (TMP@M[1,1][[1]]@m+TMP@M[1,2][[1]]@m))
                }
              }
            }
            return(MATS)
          }
          
)
#' @rdname WH.mat.prod-methods
#' @aliases WH.mat.prod,MatH-method
#' @description It is the matrix product of two \code{MatH} objects, i.e. two matrices of distributions, 
#' by using the dot product of two histograms that is consistent with 
#'  a set of distributions equipped with a L2 wasserstein metric.
#' @param object1  a \code{MatH} object
#' @param object2  a \code{MatH} object
#' @param ... other optional parameters
#' @param traspose1 a logical value, default=FALSE. If TRUE trasposes object1
#' @param traspose2 a logical value, default=FALSE. If TRUE trasposes object2
#' @return a matrix of numbers  
#' @examples
#' 
#' M1<-BLOOD[1:5,]
#' M2<-BLOOD[6:10,]
#' MAT<-WH.mat.prod(M1,M2,traspose1=TRUE, traspose2=FALSE)
setMethod(f="WH.mat.prod",signature=c(object1="MatH",object2="MatH"),
          #sums two MatH, i.e.  two matrics of distributionsH
          #INPUT: 
          #OUTPUT: 
          
          function(object1,object2,traspose1=FALSE,traspose2=FALSE){
            if (traspose1==TRUE)
            {
              #trasposing the first matrix
              object1@M=t(object1@M)
            }
            if (traspose2==TRUE)
            {
              #trasposing the second matrix
              object2@M=t(object2@M)
            }
            nrows1=nrow(object1@M)
            ncols1=ncol(object1@M)
            nrows2=nrow(object2@M)
            ncols2=ncol(object2@M)
            if (ncols1!=nrows2){
              cat("Fisrt matrix dimensions ", nrow(object1@M), "x", ncol(object1@M), "\n",
                  "Second matrix dimensions ", nrow(object2@M), "x", ncol(object2@M), "\n")
              stop("Dimensions of matrices are not compatible")}
            
            
            
            MAT=matrix(0,nrows1,ncols2)
#             cat("Fisrt matrix dimensions ", nrow(object1@M), "x", ncol(object1@M), "\n",
#                 "Second matrix dimensions ", nrow(object2@M), "x", ncol(object2@M), "\n")
            for (r in 1:nrows1){
              for (c in 1:ncols2){
                for (els in 1:ncols1){
                  MAT[r,c]=MAT[r,c]+dotpW(object1@M[r,els][[1]],object2@M[els,c][[1]])
                  
                }
              }
            }
            return(MAT)
          }
)

#L2 Wasserstein basic operations and basic statistics for matrices of distributionH ----
#' Method WH.vec.sum
#' @name WH.vec.sum
#' @rdname WH.vec.sum-methods
#' @exportMethod WH.vec.sum
setGeneric("WH.vec.sum",function(object,...) standardGeneric("WH.vec.sum"))#OK weighted sum of a vector of distributionH
#' Method WH.vec.mean
#' @name WH.vec.mean
#' @rdname WH.vec.mean-methods
#' @exportMethod WH.vec.mean
setGeneric("WH.vec.mean",function(object,...) standardGeneric("WH.vec.mean"))#OK weighted mean of a vector of distributionH
#' Method WH.SSQ
#' @name WH.SSQ
#' @rdname WH.SSQ-methods
#' @exportMethod WH.SSQ
setGeneric("WH.SSQ",function(object,...) standardGeneric("WH.SSQ"))#weighted de-codeviance matrix
#' Method WH.var.covar
#' @name WH.var.covar
#' @rdname WH.var.covar-methods
#' @exportMethod WH.var.covar
setGeneric("WH.var.covar",function(object,...) standardGeneric("WH.var.covar"))#weighted variance variance matrix
#' Method WH.correlation
#' @name WH.correlation
#' @rdname WH.correlation-methods
#' @exportMethod WH.correlation
setGeneric("WH.correlation",function(object,...) standardGeneric("WH.correlation"))#weighted corelation matrix
#' Method WH.SSQ2
#' @name WH.SSQ2
#' @rdname WH.SSQ2-methods
#' @exportMethod WH.SSQ2
setGeneric("WH.SSQ2",function(object1,object2,...) standardGeneric("WH.SSQ2"))#weighted de-codeviance matrix
#' Method WH.var.covar2
#' @name WH.var.covar2
#' @rdname WH.var.covar2-methods
#' @exportMethod WH.var.covar2
setGeneric("WH.var.covar2",function(object1,object2,...) standardGeneric("WH.var.covar2"))#weighted variance variance matrix
#' Method WH.correlation2
#' @name WH.correlation2
#' @rdname WH.correlation2-methods
#' @exportMethod WH.correlation2
setGeneric("WH.correlation2",function(object1,object2,...) standardGeneric("WH.correlation2"))#weighted corelation matrix

#' @rdname WH.vec.sum-methods
#' @aliases WH.vec.sum,MatH-method
#' @description Compute a histogram that is the weighted sum of the set of histograms contained
#' in a \code{MatH} object, i.e. a matrix of histograms, consistent with 
#' a set of distributions equipped with a L2 wasserstein metric.
#' @param object  a \code{MatH} object
#' @param ... optional arguments 
#' @param w it is possible to add a vector of weights (positive numbers) having the same size of the \code{MatH object}, 
#' default = equal weights for all cells 
#' @return a \code{distributionH} object, i.e. a histogram  
#' @examples
#' hsum<-WH.vec.sum(BLOOD)
#' # generate a set of random weights
#' RN<-runif(get.MatH.nrows(BLOOD)*get.MatH.ncols(BLOOD))
#' hsum<-WH.vec.sum(BLOOD,w=RN)

setMethod(f="WH.vec.sum",signature=c(object="MatH"),
          function(object,w=numeric(0)){
            nrows=nrow(object@M)
            ncols=ncol(object@M)
            nelem=nrows*ncols
            if (missing(w)) {
              w=rep(1,nelem)
            } 
            else {
              if (length(object@M)!=length(w)) 
                stop('Wheights must have the same dimensions of the input matrix of distributions')
              if (min(w)<0) 
                stop('Weights must be positive!!')
            }
            w=matrix(w,nrows,ncols)
            
            SUM=new("distributionH",c(0,0),c(0,1))
            for (c in 1:ncols){
              for (r in 1:nrows){
                SUM=SUM+w[r,c]*object@M[r,c][[1]]
              }
            }
            return(SUM)
          }
)
#' @rdname WH.vec.mean-methods
#' @aliases WH.vec.mean,MatH-method
#' @description Compute a histogram that is the weighted mean of the set of histograms contained
#' in a \code{MatH} object, i.e. a matrix of histograms, consistent with 
#' a set of distributions equipped with a L2 wasserstein metric.
#' @param object  a \code{MatH} object
#' @param ... optional arguments 
#' @param w it is possible to add a vector of weights (positive numbers) having the same size of
#'  the \code{MatH object}, default = equal weights for all 
#' @return a \code{distributionH} object, i.e. a histogram  
#' @examples
#' hmean<-WH.vec.mean(BLOOD)
#' # generate a set of random weights
#' RN<-runif(get.MatH.nrows(BLOOD)*get.MatH.ncols(BLOOD))
#' hmean<-WH.vec.mean(BLOOD,w=RN)

setMethod(f="WH.vec.mean",signature=c(object="MatH"),
          function(object,w=numeric(0)){
            #if (length(object@M)==1) return(object)
            nrows=nrow(object@M)
            ncols=ncol(object@M)
            nelem=nrows*ncols
            if (missing(w)) {
              w=rep(1/nelem,nelem)
            } 
            else {
              if (length(object@M)!=length(w)) 
                stop('Wheights must have the same dimensions of the input matrix of distributions')
              if (min(w)<0) 
                stop('Weights must be positive!!')
            }
            w=matrix(w,nrows,ncols)
            w=w/sum(w)
            MEAN=new("distributionH",c(0,0),c(0,1))
            for (c in 1:ncols){
              for (r in 1:nrows){
                MEAN=MEAN+w[r,c]*object@M[r,c][[1]]
              }
            }
            
            return(MEAN)
          }
)
#' @rdname WH.SSQ-methods
#' @aliases WH.SSQ,MatH-method
#' @description Compute the sum-of-squares-deviations (from the mean) matrix of a \code{MatH} object, i.e. 
#' a matrix of numbers, consistent with 
#' a set of distributions equipped with a L2 wasserstein metric.
#' @param object  a \code{MatH} object
#' @param ... some optional parameters 
#' @param w it is possible to add a vector of weights (positive numbers) 
#' having the same size of the rows of the \code{MatH object}, 
#' default = equal weight for each row
#' @return a squared \code{matrix} with the weighted sum of squares  
#' @examples
#' WH.SSQ(BLOOD)
#' # generate a set of random weights
#' RN<-runif(get.MatH.nrows(BLOOD))
#' WH.SSQ(BLOOD,w=RN)
setMethod(f="WH.SSQ",signature=c(object="MatH"),
          function(object,w=numeric(0)){
            nrows=nrow(object@M)
            ncols=ncol(object@M)
            nelem=nrows*ncols
            if (missing(w)) {
              w=rep(1,nrows)
            } 
            else {
              if (nrows!=length(w)) 
                stop('Wheights must have the same length of rows of the input matrix of distributions')
              if (min(w)<0) 
                stop('Weights must be positive!!')
            }
            w=matrix(w,nrows,1)
            #w=w/sum(w)
            DEV_MAT=matrix(0,ncols,ncols)
            colnames(DEV_MAT)=colnames(object@M)
            rownames(DEV_MAT)=colnames(object@M)
            #compute the means
            MEANS=new("MatH",1,ncols)
            for (v1 in 1:ncols){
              MEANS@M[1,v1][[1]]=WH.vec.mean(object[,v1],w)
            }            
            for (v1 in 1:ncols){
              for (v2 in v1:ncols){
                for (indiv in 1:nrows){
                   if (v1==v2){
                     DEV_MAT[v1,v2]=DEV_MAT[v1,v2]+w[indiv,1]*((object@M[indiv,v1][[1]]@s)^2+(object@M[indiv,v1][[1]]@m)^2)
                   }else{
                    DEV_MAT[v1,v2]=DEV_MAT[v1,v2]+w[indiv,1]*dotpW(object@M[indiv,v1][[1]],object@M[indiv,v2][[1]])
                   }
                }
                if (v2>v1){
                  DEV_MAT[v1,v2]=DEV_MAT[v1,v2]-sum(w)*dotpW(MEANS@M[1,v1][[1]],MEANS@M[1,v2][[1]])
                  DEV_MAT[v2,v1]=DEV_MAT[v1,v2]
                }else{
                  DEV_MAT[v1,v1]=DEV_MAT[v1,v1]-sum(w)*(MEANS@M[1,v1][[1]]@s^2+MEANS@M[1,v1][[1]]@m^2)
                }
              }
            }
            if(ncols==1){
              return(as.vector(DEV_MAT))
            }
            else return(DEV_MAT)
          }
)
#' @rdname WH.var.covar-methods
#' @aliases WH.var.covar,MatH-method
#' @description Compute the variance-covariance matrix of a \code{MatH} object, i.e. 
#' a matrix of values consistent with 
#' a set of distributions equipped with a L2 wasserstein metric.
#' @param object  a \code{MatH} object
#' @param ... some optional parameters 
#' @param w it is possible to add a vector of weights (positive numbers) 
#' having the same size of the rows of the \code{MatH object}, 
#' default = equal weight for each row
#' @return a squared \code{matrix} with the (weighted) variance-covariance values
#' @references Irpino, A., Verde, R. (2015) \emph{Basic
#' statistics for distributional symbolic variables: a new metric-based
#' approach} Advances in Data Analysis and Classification, DOI
#' 10.1007/s11634-014-0176-4  
#' @examples
#' WH.var.covar(BLOOD)
#' # generate a set of random weights
#' RN<-runif(get.MatH.nrows(BLOOD))
#' WH.var.covar(BLOOD,w=RN)
setMethod(f="WH.var.covar",signature=c(object="MatH"),
          function(object,w=numeric(0)){
            nrows=nrow(object@M)
            ncols=ncol(object@M)
            nelem=nrows*ncols
            if (missing(w)) {
              w=rep(1,nrows)
            } 
            else {
              if (nrows!=length(w)) 
                stop('Weights must have the same length of rows of the input matrix of distributions')
              if (min(w)<0) 
                stop('Weights must be positive!!')
            }
            w=matrix(w,nrows,1)
            w=w/sum(w)
            COV_MAT=WH.SSQ(object,w)
            return(COV_MAT)
            
          }
)
#' @rdname WH.correlation-methods
#' @aliases WH.correlation,MatH-method
#' @description Compute the correlation matrix of a \code{MatH} object, i.e. 
#' a matrix of values consistent with 
#' a set of distributions equipped with a L2 wasserstein metric.
#' @param object  a \code{MatH} object
#' @param ... some optional parameters 
#' @param w it is possible to add a vector of weights (positive numbers) 
#' having the same size of the rows of the \code{MatH object}, 
#' default = equal weight for each row
#' @return a squared \code{matrix} with the (weighted) correlations indices
#' @references Irpino, A., Verde, R. (2015) \emph{Basic
#' statistics for distributional symbolic variables: a new metric-based
#' approach} Advances in Data Analysis and Classification, DOI
#' 10.1007/s11634-014-0176-4  
#' @examples
#' WH.correlation(BLOOD)
#' # generate a set of random weights
#' RN<-runif(get.MatH.nrows(BLOOD))
#' WH.correlation(BLOOD,w=RN)
setMethod(f="WH.correlation",signature=c(object="MatH"),
          function(object,w=numeric(0)){
            nrows=nrow(object@M)
            ncols=ncol(object@M)
            nelem=nrows*ncols
            if (missing(w)) {
              w=rep(1,nrows)
            } 
            else {
              if (nrows!=length(w)) 
                stop('Wheights must have the same length of rows of the input matrix of distributions')
              if (min(w)<0) 
                stop('Weights must be positive!!')
            }
            w=matrix(w,nrows,1)
            w=w/sum(w)
            COV_MAT=WH.var.covar(object,w)
            CORR_MAT=as.matrix(COV_MAT)
            
            for (v1 in 1:ncols){
              for (v2 in v1:ncols){
                CORR_MAT[v1,v2]= COV_MAT[v1,v2]/sqrt((COV_MAT[v1,v1]*COV_MAT[v2,v2]))
                CORR_MAT[v2,v1]=CORR_MAT[v1,v2]
              }
            }
            return(CORR_MAT)
            
          }
)
#' @rdname WH.SSQ2-methods
#' @aliases WH.SSQ2,MatH-method
#' @description Compute the sum-of-squares-deviations (from the mean) matrix using two  \code{MatH} objects having the same number of rows,
#'  It returns a rectangular a matrix of numbers, consistent with 
#' a set of distributions equipped with a L2 wasserstein metric.
#' @param object1  a \code{MatH} object
#' @param object2  a \code{MatH} object
#' @param ... some optional parameters 
#' @param w it is possible to add a vector of weights (positive numbers) 
#' having the same size of the rows of the \code{MatH object}, 
#' default = equal weight for each row
#' @return a rectangular \code{matrix} with the weighted sum of squares  
#' @examples
#' M1<-BLOOD[,1]
#' M2<-BLOOD[,2:3]
#' WH.SSQ2(M1,M2)
#' # generate a set of random weights
#' RN<-runif(get.MatH.nrows(BLOOD))
#' WH.SSQ2(M1,M2,w=RN)
setMethod(f="WH.SSQ2",signature=c(object1="MatH",object2="MatH"),
          function(object1, object2, w=numeric(0)){
            nrows1=nrow(object1@M)
            ncols1=ncol(object1@M)
            
            nrows2=nrow(object2@M)
            ncols2=ncol(object2@M)
            
            if (nrows1!=nrows2){stop('The two matrices have a different number of rows')}
            
            if (missing(w)) {
              w=rep(1,nrows1)
            } 
            else {
              if (nrows1!=length(w)) 
                stop('Wheights must have the same length of rows of the input matrix of distributions')
              if (min(w)<0) 
                stop('Weights must be positive!!')
            }
            w=matrix(w,nrows1,1)
            #w=w/sum(w)
            DEV_MAT=matrix(0,ncols1,ncols2)
            rownames(DEV_MAT)=colnames(object1@M)
            colnames(DEV_MAT)=colnames(object2@M)
            
            #compute the means
            MEANS1=new("MatH",1,ncols1)
            for (v1 in 1:ncols1){
              MEANS1@M[1,v1][[1]]=WH.vec.mean(object1[,v1],w)
            }
            MEANS2=new("MatH",1,ncols2)
            for (v2 in 1:ncols2){
              MEANS2@M[1,v2][[1]]=WH.vec.mean(object2[,v2],w)
            }
            
            
            for (v1 in 1:ncols1){
              for (v2 in 1:ncols2){
                for (indiv in 1:nrows1){
                  DEV_MAT[v1,v2]=DEV_MAT[v1,v2]+w[indiv,1]*dotpW(object1@M[indiv,v1][[1]],object2@M[indiv,v2][[1]])
                }
                
                DEV_MAT[v1,v2]=DEV_MAT[v1,v2]-sum(w)*dotpW(MEANS1@M[1,v1][[1]],MEANS2@M[1,v2][[1]])
              }
            }
            if(ncols1==1&&ncols2==1){
              return(as.vector(DEV_MAT))
            }
            else return(DEV_MAT)
          }
)
#' @rdname WH.var.covar2-methods
#' @aliases WH.var.covar2,MatH-method
#' @description Compute the covariance matrix using two  \code{MatH} objects having the same number of rows,
#'  It returns a rectangular a matrix of numbers, consistent with 
#' a set of distributions equipped with a L2 wasserstein metric.
#' @param object1  a \code{MatH} object
#' @param object2  a \code{MatH} object
#' @param ... some optional parameters 
#' @param w it is possible to add a vector of weights (positive numbers) 
#' having the same size of the rows of the \code{MatH object}, 
#' default = equal weight for each row
#' @return a rectangular \code{matrix} with the weighted sum of squares  
#' @examples
#' M1<-BLOOD[,1]
#' M2<-BLOOD[,2:3]
#' WH.var.covar2(M1,M2)
#' # generate a set of random weights
#' RN<-runif(get.MatH.nrows(BLOOD))
#' WH.var.covar2(M1,M2,w=RN)
setMethod(f="WH.var.covar2",signature=c(object1="MatH",object2="MatH"),
          function(object1, object2, w=numeric(0)){
            nrows1=nrow(object1@M)
            ncols1=ncol(object1@M)
            
            nrows2=nrow(object2@M)
            ncols2=ncol(object2@M)
            
            if (nrows1!=nrows2){stop('The two matrices have a different number of rows')}
            
            if (missing(w)) {
              w=rep(1,nrows1)
            } 
            else {
              if (nrows1!=length(w)) 
                stop('Wheights must have the same length of rows of the input matrix of distributions')
              if (min(w)<0) 
                stop('Weights must be positive!!')
            }
            w=matrix(w,nrows1,1)
            w=w/sum(w)
            VAR_MAT=matrix(0,ncols1,ncols2)
            rownames(VAR_MAT)=colnames(object1@M)
            colnames(VAR_MAT)=colnames(object2@M)
            
            #compute the means
            MEANS1=new("MatH",1,ncols1)
            for (v1 in 1:ncols1){
              MEANS1@M[1,v1][[1]]=WH.vec.mean(object1[,v1],w)
            }
            MEANS2=new("MatH",1,ncols2)
            for (v2 in 1:ncols2){
              MEANS2@M[1,v2][[1]]=WH.vec.mean(object2[,v2],w)
            }
            
            
            for (v1 in 1:ncols1){
              for (v2 in 1:ncols2){
                for (indiv in 1:nrows1){
                  VAR_MAT[v1,v2]=VAR_MAT[v1,v2]+w[indiv,1]*dotpW(object1@M[indiv,v1][[1]],object2@M[indiv,v2][[1]])
                }
                
                VAR_MAT[v1,v2]=VAR_MAT[v1,v2]-sum(w)*dotpW(MEANS1@M[1,v1][[1]],MEANS2@M[1,v2][[1]])
              }
            }
            if(ncols1==1&&ncols2==1){
              return(as.vector(VAR_MAT))
            }
            else return(VAR_MAT)
          }
)
#' @rdname WH.correlation2-methods
#' @aliases WH.correlation2,MatH-method
#' @description Compute the correlation matrix using two  \code{MatH} objects having the same number of rows,
#'  It returns a rectangular a matrix of numbers, consistent with 
#' a set of distributions equipped with a L2 wasserstein metric.
#' @param object1  a \code{MatH} object
#' @param object2  a \code{MatH} object
#' @param ... some optional parameters 
#' @param w it is possible to add a vector of weights (positive numbers) 
#' having the same size of the rows of the \code{MatH object}, 
#' default = equal weight for each row
#' @return a rectangular \code{matrix} with the weighted sum of squares  
#' @examples
#' M1<-BLOOD[,1]
#' M2<-BLOOD[,2:3]
#' WH.correlation2(M1,M2)
#' # generate a set of random weights
#' RN<-runif(get.MatH.nrows(BLOOD))
#' WH.correlation2(M1,M2,w=RN)
setMethod(f="WH.correlation2",signature=c(object1="MatH",object2="MatH"),
          function(object1, object2, w=numeric(0)){
            nrows1=nrow(object1@M)
            ncols1=ncol(object1@M)
            
            nrows2=nrow(object2@M)
            ncols2=ncol(object2@M)
            
            if (nrows1!=nrows2){stop('The two matrices have a different number of rows')}
            
            if (missing(w)) {
              w=rep(1,nrows1)
            } 
            else {
              if (nrows1!=length(w)) 
                stop('Wheights must have the same length of rows of the input matrix of distributions')
              if (min(w)<0) 
                stop('Weights must be positive!!')
            }
            w=matrix(w,nrows1,1)
            w=w/sum(w)
            COV_MAT=WH.var.covar2(object1,object2,w)
            CORR_MAT=as.matrix(COV_MAT)
            
            for (v1 in 1:ncols1){
              for (v2 in 1:ncols2){
                CORR_MAT[v1,v2]= COV_MAT[v1,v2]/sqrt(WH.var.covar(object1[,v1],w)*WH.var.covar(object2[,v2],w))
              }
            }
            if (length(CORR_MAT)==1) return(as.vector(CORR_MAT)) 
            else return(CORR_MAT)
            
          }
)

# Utility methods for registration of distributions ----
#' Method is.registeredMH
#' @name is.registeredMH
#' @rdname is.registeredMH-methods
#' @exportMethod is.registeredMH
setGeneric("is.registeredMH",function(object) standardGeneric("is.registeredMH"))#OK
#' @rdname is.registeredMH-methods
#' @aliases is.registeredMH,MatH-method
#' @description Checks if a \code{MatH} contains histograms described by the same number of
#' bins and the same cdf.
#' 
#' @param object A \code{MatH} object
#' @return a \code{logical} value \code{TRUE} if the distributions share the
#' same cdf, \code{FALSE} otherwise.
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
#' ##---- initialize three distributionH objects mydist1 and mydist2
#'  mydist1=new("distributionH",c(1,2,3),c(0, 0.4, 1))
#'  mydist2=new("distributionH",c(7,8,10,15),c(0, 0.2, 0.7, 1))
#'  mydist3=new("distributionH",c(9,11,20),c(0, 0.8, 1))
#'  ## create a MatH object
#'  MyMAT=new("MatH",nrows=1,ncols=3,ListOfDist=c(mydist1,mydist2,mydist3), 1,3)
#'  is.registeredMH(MyMAT)
#'  ## [1] FALSE #the distributions do not share the same cdf 
#'  ## Hint: check with str(MyMAT)
#'  
#'  ## register the two distributions
#'  MATregistered=registerMH(MyMAT)
#'  is.registeredMH(MATregistered)
#'  ## TRUE #the distributions share the same cdf
#'  ## Hint: check with str(MATregistered)
#' 
setMethod(f="is.registeredMH",signature=c(object="MatH"),
          #check if all the distributions share the same cdf
          #INPUT: object11  - a vector or a matrix two distributions
          #OUTPUT: resu - a matrix of distributionH objects with
          #recomputed quantiles on a common cdf
          function(object){
            nrows=nrow(object@M)
            ncols=ncol(object@M)
            ndis=nrows*ncols
            #Check if the distribution are registered
            OK=1;
            count=1
            r=1
            tmpcdf=object@M[1,1][[1]]@p
            while (OK==1){
              count=count+1
              if (count<=ndis){
                if (!identical(tmpcdf,object@M[count][[1]]@p)) {
                  OK=0
                  return(FALSE)
                }
              }
              else {
                OK=0
                return(TRUE)
              }
            }
          }
)
#' Method registerMH
#' @name registerMH
#' @rdname registerMH-methods
#' @exportMethod registerMH
setGeneric("registerMH",function(object) standardGeneric("registerMH"))#OK
#' @rdname registerMH-methods
#' @aliases registerMH,MatH-method
#' @description \code{registerMH} method registers a set of distributions of a \code{MatH} object
#' All the
#' distribution are recomputed to obtain distributions sharing the same
#' \code{p} slot. This methods is useful for using fast computation of all
#' methods based on L2 Wasserstein metric. The distributions will have the same
#' number of element in the \code{x} slot without modifing their density
#' function.
#' 
#' 
#' @param object  A \code{MatH} object (a matrix of distributions)
#' @return A \code{MatH} object, a matrix of distributions sharing the same
#' \code{p} slot (i.e. the same cdf).
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
#'#initialize three distributionH objects mydist1 and mydist2 
#'  mydist1=new("distributionH",c(1,2,3),c(0, 0.4, 1))
#'  mydist2=new("distributionH",c(7,8,10,15),c(0, 0.2, 0.7, 1))
#'  mydist3=new("distributionH",c(9,11,20),c(0, 0.8, 1))
#'# create a MatH object
#'
#'  MyMAT=new("MatH",nrows=1,ncols=3,ListOfDist=c(mydist1,mydist2,mydist3), 1,3)
#'# register the two distributions
#'   MATregistered=registerMH(MyMAT)
#  
#'# OUTPUT the structure of MATregstered
#' str(MATregistered)
#'#   Formal class 'MatH' [package "HistDAWass"] with 1 slots
#'#   .. @@ M:List of 3
#'#   .. ..$ :Formal class 'distributionH' [package "HistDAWass"] with 4 slots
#'#   .. .. .. ..@@ x: num [1:6] 1 1.5 2 2.5 2.67 ...
#'#   .. .. .. ..@@ p: num [1:6] 0 0.2 0.4 0.7 0.8 1
#'#   ...
#'#   .. ..$ :Formal class 'distributionH' [package "HistDAWass"] with 4 slots
#'#   .. .. .. ..@@ x: num [1:6] 7 8 8.8 10 11.7 ...
#'#   .. .. .. ..@@ p: num [1:6] 0 0.2 0.4 0.7 0.8 1
#'#   ...
#'#   .. ..$ :Formal class 'distributionH' [package "HistDAWass"] with 4 slots
#'#   .. .. .. ..@@ x: num [1:6] 9 9.5 10 10.8 11 ...
#'#   .. .. .. ..@@ p: num [1:6] 0 0.2 0.4 0.7 0.8 1
#'#   ...
#'#   .. ..- attr(*, "dim")= int [1:2] 1 3
#'#   .. ..- attr(*, "dimnames")=List of 2
#'#   .. .. ..$ : chr "I1"
#'#   .. .. ..$ : chr [1:3] "X1" "X2" "X3"
# 
setMethod(f="registerMH",signature=c(object="MatH"),
          #register a row or a column vector of qfs of distributionH:
          #if the cdf are different a a matrix resu is returned with the quantiles of the two
          #distribution computed at the same levels of a common vector of cdfs.
          #INPUT: object11  - a vector or a matrix two distributions
          #OUTPUT: resu - a matrix of distributionH objects with
          #recomputed quantiles on a common cdf
          function(object){
            nrows=nrow(object@M)
            ncols=ncol(object@M)
            ndis=nrows*ncols
            #Check if the distributions are registered
            if (is.registeredMH(object)){return(object)}
            commoncdf=numeric(0)
            for (i in 1:nrows){
              for (j in 1:ncols){
                commoncdf=rbind(commoncdf,t(t(object@M[i,j][[1]]@p)))
              }
            }
            commoncdf=sort(unique(commoncdf))
            #check for tiny bins and for very long vectors of wheights
            diffs=commoncdf[2:length(commoncdf)]-commoncdf[1:(length(commoncdf)-1)]
            diffs[which(diffs<1e-8)]=0
            commoncdf=sort(unique(cumsum(x = c(0,diffs))))
            commoncdf=commoncdf/commoncdf[length(commoncdf)]
#             todelete=which(diffs<1e-8)
#             if (length(todelete)>0){
#               commoncdf=as.vector(commoncdf[-todelete,1])
#               if (coomoncdf[length(commoncdf)]<1){
#                 coomoncdf=c(commoncdf,1)
#               }
#             }
            #end of check
            nr=length(commoncdf)
            result=matrix(0,nr,(ndis+1))
            result[,(ndis+1)]=commoncdf
            NEWMAT=new("MatH",nrows,ncols)  
            for (r in 1:nrows){
              for (c in 1:ncols){
                x=numeric(0)
                for (rr in 1:nr){
                  x=c(x,compQ(object@M[r,c][[1]],commoncdf[rr]))  
                }
                NEWMAT@M[r,c][[1]]=new("distributionH",x,commoncdf)
              }
            }
            return(NEWMAT)
          }
)




#' Method Center.cell.MatH Centers all the cells of a matrix of distributions
#' @name Center.cell.MatH
#' @rdname Center.cell.MatH-methods
#' @exportMethod Center.cell.MatH
setGeneric("Center.cell.MatH",function(object) standardGeneric("Center.cell.MatH"))#OK
#' @rdname Center.cell.MatH-methods
#' @aliases Center.cell.MatH,MatH-method
#' @description The function transform a MatH object (i.e. a matrix of distributions), 
#' such that each distribution is shifted and has a mean equal to zero 
#' @param object a MatH object, a matrix of distributions.
#' @return A \code{MatH} object, having each distribution with a zero mean.
#' @examples
#' CEN_BLOOD=Center.cell.MatH(BLOOD)
#' get.MatH.stats(BLOOD, stat="mean")
setMethod(f="Center.cell.MatH",signature=c(object="MatH"),
          function(object){
            nr=get.MatH.nrows(object)
            nc=get.MatH.ncols(object)
            NM=object
            for (i in 1:nr){
              for (j in 1:nc){
                NM@M[i,j][[1]]@x=NM@M[i,j][[1]]@x-NM@M[i,j][[1]]@m
                NM@M[i,j][[1]]@m=0
              }
            }
            return(NM)
          }
)
## Show overridding ----
#' Method show for MatH
#' @name show-MatH
#' @rdname show-MatH-methods
#' @docType methods
# @aliases show,distributionH-method
# @name show
# @rdname show-MatH
#' @aliases show,MatH-method 
#' @description An overriding show method for a \code{MatH} object. The method returns a representation 
#' of the matrix using the mean and the standard deviation for each histogram. 
#' @param object  a \code{MatH} object
#' @examples
#' show(BLOOD)
#' print(BLOOD)
#' BLOOD
setMethod("show",
          signature(object="MatH"),
          definition = function(object){
            cat("a matrix of distributions \n", paste(ncol(object@M)," variables ",
                                                     nrow(object@M), " rows \n" ), "each distibution in the cell is represented by the mean and the standard deviation \n ")
            mymat=matrix(0,nrow(object@M)+1,ncol(object@M))
            for (i in 1:ncol(object@M)){mymat[1,i]=colnames(object@M)[i]}
            for (i in 1:nrow(object@M)){
              for (j in 1:ncol(object@M)){
                if(length(object@M[i,j][[1]]@x)==0){
                  mymat[i+1,j]=paste("Empty distribution")
                }
                else{
                  if ((abs(object@M[i,j][[1]]@m)>1e5 || abs(object@M[i,j][[1]]@m)<1e-5)&&
                        (object@M[i,j][[1]]@s>1e5 || object@M[i,j][[1]]@s<1e-5))  {
                    mymat[i+1,j]=paste("[m=",format(object@M[i,j][[1]]@m,digits=5,scientific=TRUE),
                                       " ,s=",format(object@M[i,j][[1]]@s,digits=5,scientific=TRUE),"]")
                  }
                  if ((abs(object@M[i,j][[1]]@m)<=1e5 && abs(object@M[i,j][[1]]@m)>=1e-5)&&
                        (object@M[i,j][[1]]@s<=1e5 || object@M[i,j][[1]]@s>=1e-5))  {
                    mymat[i+1,j]=paste("[m=",format(object@M[i,j][[1]]@m,digits=5),
                                       " ,s=",format(object@M[i,j][[1]]@s,digits=5),"]")
                  }
                  if ((abs(object@M[i,j][[1]]@m)>1e5 || abs(object@M[i,j][[1]]@m)<1e-5)&&
                        (object@M[i,j][[1]]@s<=1e5 && object@M[i,j][[1]]@s>=1e-5)) {
                    mymat[i+1,j]=paste("[m=",format(object@M[i,j][[1]]@m,digits=5,scientific=TRUE),
                                       " ,s=",format(object@M[i,j][[1]]@s,digits=5),"]")
                  }
                  if ((abs(object@M[i,j][[1]]@m)<=1e5 && abs(object@M[i,j][[1]]@m)>=1e-5)&&
                        (object@M[i,j][[1]]@s>1e5 || object@M[i,j][[1]]@s<1e-5))  {
                    mymat[i+1,j]=paste("[m=",format(object@M[i,j][[1]]@m,digits=5),
                                       " ,s=",format(object@M[i,j][[1]]@s,digits=5,scientific=TRUE),"]")
                  }
                }
              }
            }
            
            rownames(mymat)=c(
              paste(rep(" ",nchar(rownames(object@M)[1])),collapse=""),
              row.names(object@M))
            write.table(format(mymat,justify="centre"),row.names=T, col.names=F,quote=F)
          }
)
 if (!isGeneric("plot")) { 
   setGeneric("plot", 
              function(x, y, ...) standardGeneric("plot")) 
 } 
## --- Plot overloading
#' Method plot for a matrix  of histograms
#' @name plot-MatH
#' @docType methods
#' @rdname plot-MatH
#' @aliases plot,MatH-method
#' @description An overloading plot function for a \code{MatH} object. The method returns a graphical representation 
#' of the matrix of histograms. 
#' @param x a \code{distributionH} object
#' @param y not used in this implementation
#' @param type (optional) a string describing the type of plot, default="HISTO".\cr
#'  Other allowed types are \cr
#'  "DENS"=a density approximation, \cr
#'  "BOXPLOT"=l boxplot
#' @param border (optional) a string the color of the border of the plot, default="black".
#' @examples
#' plot(BLOOD) #plots BLOOD dataset
#' \dontrun{
#' plot(BLOOD, type="HISTO",  border="blue") #plots a matrix of histograms
#' plot(BLOOD, type="DENS",  border="blue") #plots a matrix of densities
#' plot(BLOOD, type="BOXPLOT") #plots a  boxplots
#' }
#' @importFrom utils write.table
#' @export

setMethod("plot",
          signature(x = "MatH"),
          function (x, y="missing", type="HISTO",border="black") 
          {plot.M(x, type=type, border=border)

          }
)
