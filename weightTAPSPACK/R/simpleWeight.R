#' Find weights for survey data
#'
#' Find weights for survey from a dataframe of observations and population margins. This function performs a raking process based on the population margins specified by the user.
#'
#' @param data A dataframe of observations and values for population categories
#' @param pop.margins A list of known population parameters
#' @param names The names of the population parameters. The dataframe MUST include variables with the same names as the population parameters specified.
#' @param trunc_at A numeric object specifying where to truncate the weights (what should the max weight be?) - default is 5
#'
#' @return A vector of weights
#' @docType methods
#' @author David G. Carlson \email{carlson.david@@wustl.edu}  Michelle Torres: \email{smtorres@@wustl} Taeyong Park \email{t.park@@wustl.edu}
#' @seealso \code{\link{weightTAPSPACK}} \code{\link{weightTAPS}} \code{\link{variablesTAPS}} \code{\link{weightTAPSoutput}} \code{\link{subsetTAPS}} \code{\link{attritTAPS}} \code{\link{multipleImp}} \code{\link{hotdeckImp}} \code{\link{wavesTAPS}}
#' @rdname simpleWeight
#' @aliases simpleWeight,ANY-method
#' @export
setGeneric(name="simpleWeight",
           def=function(data,pop.margins,names=c(~agegend,~ethm,~educat,~regmetro,~incomcat,~ppnet),trunc_at=5)
           {standardGeneric("simpleWeight")}
)

setMethod(f="simpleWeight",
          definition=function(data,pop.margins,names=c(~agegend,~ethm,~educat,~regmetro,~incomcat,~ppnet),trunc_at=5){
            if(all(data$baseweights==1)){
              tapsdesign <- svydesign(ids = ~1, strata = NULL, data=data, weights = NULL)
              
              
            }else{
              test1<-try(svydesign(ids = ~1, strata = ~strata, data=data, weights = ~baseweights))
              if(class(test1)[[1]]=="try-error"){
                print("Strata were missing in the subsetted data and had to be set to NULL.\n Consider not using baseweights.")
                c<-readline("\n Press enter to continue:")
                tapsdesign <- svydesign(ids = ~1, strata = NULL, data=data, weights = ~baseweights)
              }else{
                tapsdesign <- svydesign(ids = ~1, strata = ~strata, data=data, weights = ~baseweights)
                test<-try(rake(tapsdesign, sample.margins=names, population.margins = pop.margins, control=list(maxit=100,epsilon=1,verbose=FALSE)),silent=TRUE)
                if(class(test)[[1]]=="try-error"){
                  print("Strata were missing in the subsetted data and had to be set to NULL.\n Consider not using baseweights.")
                  c<-readline("\n Press enter to continue:")
                  tapsdesign <- svydesign(ids = ~1, strata = NULL, data=data, weights = ~baseweights)
                }
              }
            }
            taps.design.raked <- rake(tapsdesign, sample.margins=names, population.margins = pop.margins, control=list(maxit=100,epsilon=1,verbose=FALSE))
            n.obs<-length(weights(taps.design.raked))
            if(max(weights(taps.design.raked))>=(trunc_at/n.obs)){
              print(paste('Warning - the weights will be truncated. The highest pre-truncated weight is: ', max(weights(taps.design.raked))*n.obs))
              c<-readline('\n Press enter to continue:')
              test<-try(taps.design.raked<-trimWeights(taps.design.raked, upper=(trunc_at/n.obs), trim=TRUE, strict=TRUE), silent=TRUE)
              if(class(test)[[1]]=="try-error"){
                taps.design.raked<-trimWeights(taps.design.raked, upper=(trunc_at/n.obs), trim=TRUE)
              }else{
                taps.design.raked<-trimWeights(taps.design.raked, upper=(trunc_at/n.obs), trim=TRUE, strict=TRUE)
              }
            }
            final.weights<-weights(taps.design.raked)*n.obs
            return(final.weights)  
          }
)
