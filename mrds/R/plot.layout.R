#' Layout for plot methods in mrds
#'
#' This function does the paging, using \code{devAskNewPage()}. This means we
#' can just call plots and R will make the prompt for us
#' Warning, this function has side effects! It modifies \code{devAskNewPage}!
#'
#' Code is stolen and modified from plot.R in mgcv by Simon Wood
#'
#' @param which which plots are to be created
#' @param pages number of pages to span the plots accross
#'
#' @author David L. Miller, based on code by Simon N. Wood
#' @importFrom grDevices dev.interactive
plot.layout <- function(which,pages){

  # how many plots are there
  n.plots <- length(which)

  # handle edge cases
  if (n.plots==0) stop("No terms to plot - nothing for plot() to do.")
  if (pages>n.plots) pages<-n.plots
  if (pages<0) pages<-0

  # figure out how to display things
  if(pages!=0){
    ppp <- n.plots%/%pages
    if (n.plots%%pages!=0){
      ppp<-ppp+1
      while (ppp*(pages-1)>=n.plots) pages<-pages-1
    }
    # now figure out number of rows and columns
    c <- r <- trunc(sqrt(ppp))
    if (c<1) r <- c <- 1
    if (c*r < ppp) c <- c + 1
    if (c*r < ppp) r <- r + 1
    oldpar<-par(mfrow=c(r,c))
  }else{
    ppp<-1
    oldpar<-par()
  }

  # should we ask for new pages?
  if((pages==0&&prod(par("mfcol"))<n.plots&&dev.interactive())||
      pages>1&&dev.interactive()){
    ask <- TRUE
  }else{
    ask <- FALSE
  }

  if(length(which)==1 | pages==1){
    ask <- FALSE
  }

  if(ask){
    oask <- devAskNewPage(TRUE)
  }else{
    oask <- devAskNewPage()
  }

  return(oask)
}
