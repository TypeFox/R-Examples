#' Marking product-limit plots at the censored times.
#' 
#' This function is invoked and controlled by \code{plot.prodlim}.
#' 
#' This function should not be called directly. The arguments can be specified
#' as \code{atRisk.arg} in the call to \code{plot.prodim}.
#' 
#' @param x The values of the curves at \code{times}.
#' @param times The times where there curves are plotted.
#' @param nlost The number of subjects lost to follow-up (censored) at
#' \code{times}.
#' @param pch The symbol used to mark the curves.
#' @param col The color of the symbols.
#' @param ...  Arguments passed to \code{points}.
#' @return Nil
#' @author Thomas Alexander Gerds <tag@@biostat.ku.dk>
#' @seealso \code{\link{plot.prodlim}}, \code{\link{confInt}},
#' \code{\link{atRisk}}
#' @keywords survival
#' @export
markTime <- function(x,times,nlost,pch,col,...){
  mtimeList=lapply(1:length(x),function(i){
    who=nlost[[i]]>0 & !is.na(nlost[[i]])
    mark.x=times[who]
    mark.y=x[[i]][who]
    if (length(col)<length(x)) mcol=col[1] else mcol=col[i]
    if (length(pch)<length(x)) mpch=pch[1] else mpch=pch[i]
    points(x=mark.x,y=mark.y,col=mcol,pch=mpch,...)
    ##       cbind(mark.x,mark.y)
    invisible(NULL)
  })
}


