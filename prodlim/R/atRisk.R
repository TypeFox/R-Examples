#' Drawing numbers of subjects at-risk of experiencing an event below
#' Kaplan-Meier and Aalen-Johansen plots.
#' 
#' This function is invoked and controlled by \code{plot.prodlim}.
#' 
#' This function should not be called directly. The arguments can be specified
#' as \code{atRisk.arg} in the call to \code{plot.prodim}.
#' 
#' @param x an object of class `prodlim' as returned by the
#' \code{prodlim} function.
#' @param newdata see \code{plot.prodim}
#' @param times Where to compute the atrisk numbers.
#' @param line Distance of the atrisk numbers from the inner plot.
#' @param col The color of the text.
#' @param labelcol The color for the labels. Defaults to col.
#' @param interspace Distance between rows of atrisk numbers.
#' @param cex Passed on to \code{mtext} for both atrisk numbers and
#' labels.
#' @param labels Labels for the at-risk rows.
#' @param title Title for the at-risk labels
#' @param titlecol The color for the title. Defaults to 1 (black).
#' @param pos The value is passed on to the \code{mtext} argument
#' \code{at} for the labels (not the atriks numbers).
#' @param adj Passed on to \code{mtext} for the labels (not the atriks
#' numbers).
#' @param dist If \code{line} is missing, the distance of the upper
#' most atrisk row from the inner plotting region: par()$mgp[2].
#' @param adjust.labels If \code{TRUE} the labels are left adjusted.
#' @param ... Further arguments that are passed to the function
#' \code{mtext}.
#' @return Nil
#' @author Thomas Alexander Gerds <tag@@biostat.ku.dk>
#' @seealso \code{\link{plot.prodlim}}, \code{\link{confInt}},
#' \code{\link{markTime}}
#' @keywords survival
#' @export
atRisk <- function(x,
                   newdata,
                   times,
                   line,
                   col,
                   labelcol=NULL,
                   interspace,
                   cex,
                   labels,
                   title="",
                   titlecol=NULL,
                   pos,
                   adj,
                   dist,
                   adjust.labels=TRUE,
                   ...){
    if (missing(times)) times <- seq(0,x$maxtime,x$maxtime/10)
    if (x$model=="competing.risks"){
        px <- lifeTab(object=x,times=times,cause=1,newdata=newdata,stats=NULL)[[1]]
    }
    else if (x$model=="survival"){
        px <- lifeTab(object=x,times=times,newdata=newdata,stats=NULL)
    }
    if (is.matrix(px) || is.data.frame(px))
        sumx <- lapply(data.frame(px)[,grep("n.risk",colnames(px)),drop=FALSE],function(x)x)
    else
        sumx <- lapply(px,function(v){
                           u <- v[,grep("n.risk",colnames(v)),drop=FALSE]
                           if (NCOL(u)>1){
                               ulist <- lapply(1:NCOL(u),function(i)u[,i])
                               names(ulist) <- colnames(u)
                               ulist
                           }
                           else
                               u
                       })
    if (is.list(sumx[[1]]))
        sumx <- unlist(sumx,recursive=FALSE)
    if (all(sapply(sumx,NCOL))==1)
        nlines <- length(sumx)
    if (missing(line)){
        line <- par()$mgp[2] + dist +
            (0:(2*nlines-1)) *interspace -(nlines-1)
    }
    if (missing(cex)) cex <- 1
    ## if (missing(pos)) pos <- min(times)
    if (missing(pos)) pos <- par()$usr[1]
    if (missing(adj)) adj <- 1
    if (missing(labels))
        if (length(names(sumx)==nlines))
            labels <- paste("",names(sumx),"",sep="")
        else
            labels <- rep("",nlines)
    ## c("No.   \nsubjects",rep("",nlines-1))
    # title for no. at-risk below plot
    # --------------------------------------------------------------------
    if (is.null(titlecol)){
        tcol <- 1
    } else {
          if (is.na(titlecol[1]))
              tcol <- 1
          else
              tcol <- titlecol[1]
      }
    ##
    if (!is.null(title))
        mtext(title,
                        side=1,
                        at=pos,
                        col=tcol,
                        line=line[1]-1,
                        adj=adj,
                        cex=cex,
                        outer=FALSE,
                        xpd=NA,
                        ...)
    # labeling the no. at-risk below plot
    # --------------------------------------------------------------------
    ## if (is.null(adjust.labels) || adjust.labels==TRUE){
    ## labels <- format(labels,justify="left")}
    if (length(col)==nlines/2) ## 1 cluster level
        col <- rep(col,rep(2,length(col)))
    lapply(1:nlines,function(y){
               mtext(text=as.character(sumx[[y]]),
                               side=1,
                               at=times,
                               line=rep(line[y],length(times)),
                               col=rep(col[y],length(times)),
                               cex=cex,
                               outer=FALSE,
                               xpd=NA,
                               ...)
               if (is.null(labelcol)){
                   lcol <- col[y]
               } else {
                     if (is.na(labelcol[y]))
                         lcol <- labelcol[1]
                     else
                         lcol <- labelcol[y]
                 }
               ## print(labels[y])
               mtext(text=labels[y],
                               side=1,
                               at=pos,
                               col=labelcol[y],
                               ## col=1,
                               line=line[y],
                               adj=adj,
                               cex=cex,
                               outer=FALSE,
                               xpd=NA,
                               ...)
           })
}
