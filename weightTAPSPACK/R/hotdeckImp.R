#' Impute data using hotdeck imputation method
#'
#' Impute missing data in a dataframe using hotdeck imputation method
#'
#' @param data A data frame with missing values to be imputed
#' @param outcome Character strings indicating outcome variable(s) of the user's interest 
#' @param covars Character strings indicating covariates of the user's interest
#'
#' @return A data frame with imputed values
#' @docType methods
#' @author David G. Carlson \email{carlson.david@@wustl.edu}  Michelle Torres: \email{smtorres@@wustl} Taeyong Park \email{t.park@@wustl.edu}
#' @seealso \code{\link{weightTAPSPACK}} \code{\link{weightTAPS}} \code{\link{variablesTAPS}} \code{\link{weightTAPSoutput}} \code{\link{simpleWeight}} \code{\link{subsetTAPS}} \code{\link{attritTAPS}} \code{\link{hotdeckImp}} \code{\link{wavesTAPS}}
#' @rdname hotdeckImp
#' @aliases hotdeckImp,ANY-method
#' @export
setGeneric(name="hotdeckImp",
           def=function(data, outcome, covars)
           {standardGeneric("hotdeckImp")}
)

setMethod(f="hotdeckImp",
          definition=function(data, outcome, covars){ 
            holddata<-data
            data<-cbind(data[ ,outcome],data[,covars])
            numberOfNA<-apply(data, 2, FUN=function(x) length(which(is.na(x))))
            data<-data[ , numberOfNA<nrow(data)] # Delete those variables that have no data
            names<-list()
            for (i in 1:ncol(data)){
              names[[i]]<-levels(data[,i])
            }
            data<-data.frame(data, stringsAsFactors=TRUE)
            factorVar<-logical(0)
            dataMat<-matrix(NA, nrow(data), ncol(data)) # Convert the input data into a matrix to use the impute.NN_HD function
            for (i in 1:ncol(data)){
              factorVar[i]<-is.factor(data[,i]) # Before we unclass the input data, we should be aware of what variables are factors
              dataMat[,i]<-unclass(data)[[i]] # Unclass the input data to use the impute.NN_HD function
            }
            lengthNA<-numeric(0)
            for (j in 1:nrow(dataMat)){
              lengthNA[j]<-length(which(apply(dataMat, 1, is.na)[,j])) # Hot deck imputation does not work for the observations that are missing all covariates of interest; we delete these observations
            }
            dataMat<-dataMat[!lengthNA==ncol(dataMat), ] # Delete the observations that are missing all covariates of interest
            imputedSubset<-impute.NN_HD(dataMat) # Run the impute.NN_HD function that performs nearest neighbor hot deck imputation 
            imputedSubset<-data.frame(imputedSubset) # Convert the output into a data frame
            if (length(which(factorVar))>2) {
              imputedSubset[,factorVar]<-data.frame(apply(imputedSubset[,factorVar], 2, as.factor)) # Convert the variables that were originally a factor into a factor
            } else if (length(which(factorVar))==1) {
              imputedSubset[,factorVar]<-data.frame(as.factor(imputedSubset[,factorVar]))
            } else {
              imputedSubset<-imputedSubset
            }
            for (i in 1:length(names)){
              if(!length(names[[i]])==0) levels(imputedSubset[,i])<-names[[i]]
            }
            colnames(imputedSubset)<-colnames(data)
            holddata[,names(imputedSubset)[-(1:length(outcome))]]<-imputedSubset[,names(imputedSubset)[-(1:length(outcome))]]
            return(holddata)
          }
)
