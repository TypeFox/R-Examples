#' Impute data using multiple imputation method
#'
#' Impute missing data in a dataframe using multiple imputation method
#'
#' @param data A data frame with missing values to be imputed
#' @param outcome Character strings indicating outcome variable(s) of the user's interest 
#' @param covars Character strings indicating covariates of the user's interest
#' @param m A nunmeric object indicating number of data sets to create - default is 5
#'
#' @return A data frame with imputed values
#' @docType methods
#' @author David G. Carlson \email{carlson.david@@wustl.edu} Michelle Torres: \email{smtorres@@wustl} Taeyong Park \email{t.park@@wustl.edu}
#' @seealso \code{\link{weightTAPSPACK}} \code{\link{weightTAPS}} \code{\link{variablesTAPS}} \code{\link{weightTAPSoutput}} \code{\link{simpleWeight}} \code{\link{subsetTAPS}} \code{\link{attritTAPS}} \code{\link{hotdeckImp}} \code{\link{wavesTAPS}}
#' @rdname multipleImp
#' @aliases multipleImp,ANY-method
#' @export
setGeneric(name="multipleImp",
           def=function(data, outcome, covars, m=5)
           {standardGeneric("multipleImp")}
)

setMethod(f="multipleImp",
          definition=function(data, outcome, covars, m=5){ # The input data come from subsetTAPS.R; m specifies the number of imputed data sets the user want to create. 
            holddata<-data
            holddataRep<-rep(holddata, m) 
            holddataList<-rep(list(NA), m) 
            for(j in 1:m){
              for(i in 1:ncol(holddata)){
                k<-i%%ncol(holddata)
                if(k==0) {
                  k <- ncol(holddata)
                }
                else k <- k
                holddataList[[j]][k]<-list(holddataRep[[i]]) 
              }
              names(holddataList[[j]])<-colnames(holddata)
            }
            data<-cbind(data[ ,outcome],data[,covars])
            numberOfNA<-apply(data, 2, FUN=function(x) length(which(is.na(x))))
            dataForImputing<-data[ , numberOfNA>0 & numberOfNA<nrow(data)] 
            data<-data[ , numberOfNA<nrow(data)] # Delete those variables that have no data
            data<-data.frame(data, stringsAsFactors=TRUE)
            test<-logical()
            for (i in 1:ncol(data)){
              test<-c(test,is.character(data[,i]))
            }
            test<-!test
            data<-data[,test]
            test<-!test
            imputed<-mice(data, m) # Use the mice function to create "m" imputed data sets
            dataRep<-rep(data, m) # To accommodate "m" imputed data sets, replicate the input data "m" times
            imputedSubset<-rep(list(NA), m) # imputedSubset will be the final output. Create an empty list that will contain the final results.
            for(j in 1:m){
              for(i in 1:ncol(data)){
                k<-i%%ncol(data)
                if(k==0) {
                  k <- ncol(data)
                }
                else k <- k
                imputedSubset[[j]][k]<-list(dataRep[[i]]) # This produces a list that contains the replicated input data set
              }
              names(imputedSubset[[j]])<-colnames(data)
            }
            # Now, we want to plug the imputed data yielded by the mice function into the missing values in imputedSubset 
            q1<-which(names(imputedSubset[[1]]) %in% names(dataForImputing))
            q2<-list()
            for (i in 1:length(q1)){
              q2[[i]]<-which(is.na(imputedSubset[[1]][[q1[i]]]))
            }
            for(j in 1:m){
              for (i in 1:length(q1)){
                if(!is.null(imputedSubset[[j]][[q1[i]]][q2[[i]]])&!is.null(imputed[[5]][[q1[i]]][,j])) imputedSubset[[j]][[q1[i]]][q2[[i]]]<-imputed[[5]][[q1[i]]][,j]
              }
            }
            for(i in 1:length(imputedSubset)){
              imputedSubset[[i]]<-rbind(imputedSubset[[i]],data[,test])
            }
            qNames<-names(imputedSubset[[1]])[-(1:length(outcome))]
            for(j in 1:m){
              for(i in 1:length(qNames)){
                holddataList[[j]][qNames[i]]<-imputedSubset[[j]][qNames[i]]
              }
            }
            for(i in 1:length(holddataList)){
              holddataList[[i]]<-as.data.frame(holddataList[[i]])
            }
            return(holddataList)
          }
)
