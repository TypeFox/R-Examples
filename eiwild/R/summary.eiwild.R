#'
#' summary method for \code{eiwild}-object
#' 
#' @param object \code{eiwild}-object
#' @param cred percentage for credibility interval of alphas and cellCounts
#' @param ... no function. included for S3 generic/method consistency
#' 
#' @return tables and matrices
#' \itemize{
#'    \item \code{relative:} global beta values calculated with \code{cellCounts}
#'    \item \code{absolut:} \code{cellCounts} mean
#'    \item \code{alphaMeans:} Means of \code{alphaDraws}
#'    \item \code{relativeCol:} proportions with \code{colSum=1}
#'    \item \code{countsCred:} Credibility Interval of length \code{cred} for \code{cellCounts}
#'    \item \code{alphaCred:} Credibility Interval of length \code{cred} for \code{alphaDraws}
#'    \item \code{realtiveCred:} Credibility Interval of length \code{cred} for global beta values calculated with \code{cellCounts}
#' }
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
#' sumRes <- summary(res)
#' sumRes
#' str(sumRes)
#' }
#' 
#' 
#' @export
#' @exportClass summary.eiwild
#' @S3method summary eiwild

summary.eiwild <- function(object, cred=0.95, ...){
   
   r <- ncol(object$rowdf)
   c <- ncol(object$coldf)
   le <- nrow(object$draws$alphaDraws)
   
   # absolute and relative means for alphas and cellCounts
   alMeans <- matrix(colMeans(object$draws$alphaDraws),r,c)
   cellMeans <- matrix(colMeans(object$draws$cellCounts),r,c)
   rel <- cellMeans/rowSums(cellMeans)
   relCol <- t(t(cellMeans)/colSums(cellMeans))
   rownames(rel) <- colnames(object$rowdf)
   colnames(rel) <- colnames(object$coldf)
   dimnames(alMeans) <- dimnames(cellMeans) <- dimnames(relCol) <- dimnames(rel)
   
   # credibility interval
   if(cred>1|cred<0)
      stop("\"cred\" has to be between 0 and 1!" ,call.=FALSE)
   inter <- c((1-cred)/2, cred + (1-cred)/2)
   countsCred <- apply(object$draws$cellCounts, 2, "quantile", inter)
   alphaCred <- apply(object$draws$alphaDraws, 2, "quantile", inter)
   
   relAll <- t(sapply(1:le, function(k){ # make cell-wise probability per sample-iteration
      tmp <- matrix(object$draws$cellCounts[k,],r,c)
      c( tmp / rowSums(tmp))
      }))
   relativeCred <- apply(relAll, 2, "quantile", inter)
   colnames(relativeCred) <- gsub("alpha","betaGlob",colnames(alphaCred))
   
   ret <- list(relative = rel,
               absolut = round(cellMeans),
               alphaMeans = alMeans,
               relativeCol = relCol,
               countsCred = round(countsCred),
               alphaCred = alphaCred,
               relativeCred = relativeCred)
   class(ret) <- c("summary.eiwild", class(ret))
   return(ret)
}

#' @S3method print summary.eiwild

print.summary.eiwild <- function(x, ...){
   cat("relative cellMeans:\n")
   print(round(x$relative,4))
   cat("\n")
   cat("absolute cellMeans:\n")
   print(round(x$absolut,4))  
}


#' get profit and loss of partys
#' 
#' @param x table of ecological inference (see \code{\link{summary.eiwild}})
#' @param rnd rounding of values (default is \code{1})
#' @param zero replace negative values with \code{0} (default=\code{TRUE})
#' @param which if table isn't square it has to be a vector giving 2 arguments (see Details)
#' 
#' @details
#' if table isn't square the row or column not to be calculated has to be given in \code{wich}.
#' 
#' First element has to be \code{"r"} for row or \code{"c"} for column.
#' 
#' 2nd element has to give name of row or column.
#' 
#' @return
#' table with balance BUT (!!!) order of rows or cols maybe changed
#' 
#' @seealso
#' \code{\link{summary.eiwild}}
#' 
#' @export
#'   

## TO DO: make it work with more than one extra column/row
getBalance <- function(x, rnd=1, zero=TRUE, which=NULL){
   
   r <- nrow(x)
   c <- ncol(x)   
   if( !any(r==c | !is.null(which)) )
      stop("Either \"which\" has to be specified or \"x\" has to be square!" ,call.=FALSE)
   
      # Für den Fall dass which spezifiziert ist
   if(!is.null(which)){
      if( !which[1] %in% c("r","c")){
         stop("\"which[1]\" has to be \"r\" or \"l\" ", call.=FALSE)
      } else{
         if(which[1] =="r"){
            which2 <- which(rownames(x)==which[2])
            x2 <- x[-which2,]
         } else{
            which2 <- which(colnames(x)==which[2])
            x2 <- x[,-which2]
         }
      }
      if(ncol(x2)!=nrow(x2))
         stop("\"which\" wasn't specified enough. Table isn't square!", call.=FALSE)
      
      x3 <- x2 - t(x2)
      
      if(which[1] =="r"){
         ret <- rbind(x3,x[which2,])
         rownames(ret)[nrow(ret)] <- rownames(x)[which2]
      } else{
         ret <- cbind(x3,x[,which2])
         colnames(ret)[ncol(ret)] <- colnames(x)[which2]
      } 
   } else if(r==c){
      ret <- x - t(x)
   }
   
   if(zero==TRUE)
      ret[which(ret<0)] <- 0 
   
   ret <- round(ret/rnd)*rnd
   
   return(ret)
}



