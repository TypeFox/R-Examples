#' @export
#' @import lme4
summary.arfimaMLM <-
function(object, ...){
  d<- object$d
  result <- summary(object$result, ...)

  out <- list(d=d,result=result)
  
  if(!is.null(object$ecm)){
    out$ecm <- summary(object$ecm, ...)
  }

  if(!is.null(object$arma)){
      arma_print <- data.frame(Variable=NA)
      for(i in 1:length(object$arma)){
          se_tmp <- sqrt(diag(object$arma[[i]]$var.coef))
          names(se_tmp) <- paste0(names(se_tmp),".SE")
          arma_print <- merge(arma_print, all=TRUE,
              data.frame(Variable=ls(object$arma)[[i]]
                         , t(as.matrix(object$arma[[i]]$coef))
                         , t(as.matrix(se_tmp))
                         )
          )
      }
      out$arma <- arma_print[!is.na(arma_print$Variable)
                              , (order(names(arma_print)[-1])+1)]
      rownames(out$arma) <- arma_print$Variable[!is.na(arma_print$Variable)]
  }
  
  class(out) <- "summary.arfimaMLM"
  out
}
