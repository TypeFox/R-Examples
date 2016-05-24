#'
#' Trims eiwild-object
#' 
#' @description
#' Trims eiwild-object with burn-in and thinning
#' 
#' @param obj object of type eiwild
#' @param burnin number of draws to be cut away from the beginning of the chain. default=0
#' @param thinning number specifying the thinning interval. default=1
#' @param sample specifies sample size after burn-in and thinning (default is \code{NULL})
#' 
#' @return eiwild-object
#' 
#' @examples
#' \dontrun{
#' # loading some fake election data
#' data(topleveldat)
#' form <- cbind(CSU_2, SPD_2, LINK_2, GRUN_2) ~ cbind(CSU_1, SPD_1, Link_1)
#' set.seed(1234)
#' res <- indAggEi(form=form, aggr=aggr, indi=indi, IDCols=c("ID","ID"),
#'                 sample=1000, thinning=2, burnin=100,verbose=100)
#' 
#' res
#' res2<- trimEiwild(res, burnin=100, thinning=3)
#' res2
#' }
#' 
#' @export
#' 

trimEiwild <- function(obj, burnin=0, thinning=1, sample=NULL){
   
   if(class(obj)[1]!="eiwild")
      stop("\"obj\" has to be of class eiwild!", call.=FALSE)
   
   draws <- nrow(obj$draws$alphaDraws)
   ret <- obj
   
   if(burnin < 0)
      stop("burnin has to be a positive number!", call.=FALSE)
   if(thinning < 1)
      stop("thinning has to be greater than 0!", call.=FALSE)
   
   alphaDraws <- obj$draws$alphaDraws[(burnin+1):draws,]
   alphaDraws <- alphaDraws[seq(1, nrow(alphaDraws), by=thinning),]
   if(!is.null(sample) && sample<nrow(alphaDraws))
      alphaDraws <- alphaDraws[1:sample,]
   alphaDraws <- mcmc(alphaDraws)
   ret$draws$alphaDraws <- alphaDraws
   
   if("betaDraws" %in% names(obj$draws)){
      betaDraws <- obj$draws$betaDraws[(burnin+1):draws,]
      betaDraws <- betaDraws[seq(1, nrow(betaDraws), by=thinning),]
      if(!is.null(sample) && sample<nrow(betaDraws))
         betaDraws <- betaDraws[1:sample,]
      betaDraws <- mcmc(betaDraws) 
      ret$draws$betaDraws <- betaDraws
   }
   
   cellCounts <- obj$draws$cellCounts[(burnin+1):draws,]
   cellCounts <- cellCounts[seq(1, nrow(cellCounts), by=thinning),]
   if(!is.null(sample) && sample<nrow(cellCounts))
      cellCounts <- cellCounts[1:sample,]
   cellCounts <- mcmc(cellCounts)
   ret$draws$cellCounts <- cellCounts      
   
   return(ret)
}

