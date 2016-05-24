#' Subset TAPS data
#'
#' Subset TAPS data by outcome and covariates specified.
#'
#' @param outcome A character vector of the names of outcome variables of interest
#' @param covars A character vector of the names of the covariate variables of interest
#' @param weight A logical argument specifying whether to use TAPS base weights or not - default is FALSE
#' @param refusedasNA A logical argument specifying whether to consider the response refused as an NA - default and suggested value is TRUE
#'
#' @return A dataframe of subsetted TAPS data
#' @docType methods
#' @author David G. Carlson \email{carlson.david@@wustl.edu}  Michelle Torres: \email{smtorres@@wustl} Taeyong Park \email{t.park@@wustl.edu}
#' @seealso \code{\link{weightTAPSPACK}} \code{\link{weightTAPS}} \code{\link{variablesTAPS}} \code{\link{weightTAPSoutput}} \code{\link{simpleWeight}} \code{\link{attritTAPS}} \code{\link{multipleImp}} \code{\link{hotdeckImp}} \code{\link{wavesTAPS}}
#' @rdname subsetTAPS
#' @aliases subsetTAPS,ANY-method
#' @export
setGeneric(name="subsetTAPS",
           def=function(outcome,covars=NULL,weight=FALSE,refusedasNA=TRUE)
           {standardGeneric("subsetTAPS")}
)

setMethod(f="subsetTAPS",
          definition=function(outcome,covars=NULL,weight=FALSE,refusedasNA=TRUE){

            if(weight) weighter<-as.numeric(dd$basewt) else weighter<-1
            if(is.null(covars)){
              q<-numeric()
              for(i in 1:length(outcome)){
                
                if(!outcome[i]%in%colnames(TAPScum)){
                  holder<-readline(paste("Outcome variable",outcome[i],"is not a valid variable name\n and is being dropped, press enter to continue"))
                  q<-c(q,i)
                }
              }
              if(length(q)!=0) outcome<-outcome[-q]
              if(is.na(outcome[1])) stop("No outcome variables, function terminated")
              data<-data.frame(wustlid=as.numeric(TAPScum$wustlid),TAPScum[outcome],baseweights=weighter)
            }else{
              q<-numeric()
              for(i in 1:length(outcome)){
                
                if(!outcome[i]%in%colnames(TAPScum)){
                  holder<-readline(paste("Outcome variable",outcome[i],"is not a valid variable name\n and is being dropped, press enter to continue"))
                  q<-c(q,i)
                }
              }
              if(length(q)!=0) outcome<-outcome[-q]
              if(is.na(outcome[1])) stop("No outcome variables, function terminated")
              q<-numeric()
              for(i in 1:length(covars)){
                
                if(!covars[i]%in%colnames(TAPScum)){
                  holder<-readline(paste("Covariate variable",covars[i],"is not a valid variable name\n and is being dropped, press enter to continue"))
                  q<-c(q,i)
                }
              }
              if(length(q)!=0) covars<-covars[-q]
             if(!is.na(covars[1])){
               data<-data.frame(wustlid=as.numeric(TAPScum$wustlid),TAPScum[outcome],TAPScum[covars],baseweights=weighter)
             }else{
               data<-data.frame(wustlid=as.numeric(TAPScum$wustlid),TAPScum[outcome],baseweights=weighter)
             }
              
            }
            data[is.null(data)]<-NA
            if(refusedasNA){
              data[data=="Refused"]<-NA
              if(!is.null(covars)){
                for(i in 1:length(covars)){
                  if(is.factor(data[,covars[i]])) data[,covars[i]]<-factor(data[,covars[i]])
                  
                }
              }
              
            }
            for(i in 1:length(outcome)){
              data<-as.data.frame(data[!is.na(data[outcome[i]]),])
            }
            impdata <- merge(data, dd, by="wustlid")   
            impdata[is.null(impdata)]<-NA
            
            return(impdata)
          }
)
