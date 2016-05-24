plot.FRBmultireg <- function(x, expl, resp, confmethod = c("BCA","basic"), onepage = TRUE,...) {

FRBres <- x
confmethod <- match.arg(confmethod)
currentAsk <- devAskNewPage(ask = NULL)
if (!is.null(x$bootest)) {
  allexpl <- rownames(FRBres$coefficients)
  allresp <- colnames(FRBres$coefficients)
  p <- length(allexpl)
  q <- length(allresp)
  
  if (missing(expl)) {expl <- allexpl; explMissing <- TRUE} else explMissing <- FALSE
  if (missing(resp)) {resp <- allresp; respMissing <- TRUE} else respMissing <- FALSE
  nexpl <- length(expl)
  nresp <- length(resp)
  
  if (is.numeric(expl)) {
      if (any(!(expl %in% c(1:p)))) stop(paste("indices in 'expl' should be in 1:",p,sep=""))
      rinds <- expl
      expl <- allexpl[rinds]
  }
  if (is.numeric(resp)) {
      if (any(!(resp %in% c(1:q)))) stop(paste("indices in 'resp' should be in 1:",q,sep=""))
      cinds <- resp
      resp <- allresp[cinds]
  }
  if (any(!(expl %in% allexpl))|any(!(resp %in% allresp))) stop("One or more specified variable names does not match")
  rinds <- match(expl, allexpl)
  cinds <- match(resp, allresp)
  
    R <- ncol(FRBres$bootest$centered)
  
  ## either all plots are put on one page, but then there is a severe limit for nexpl!
  if (onepage) {
  #  # if expl or resp was not specified, take a subset in case total is too much
  #  giveWarning <- FALSE
  #  if (explMissing & nexpl>5) { 
  #    giveWarning <- TRUE
  #    expl <- expl[1:5] 
  #    nexpl <- 5 
  #  }
  #  if (respMissing & nresp>4) { 
  #    giveWarning <- TRUE
  #    resp <- resp[1:4] 
  #    nresp <- 4 
  #  }
  #  if (giveWarning) warning("Number of plots too large, subset was selected: consider specifying variables in
  #  arguments 'expl' and 'resp'; or set 'onepage=FALSE'")
  #
  #  # otherwise, let's see if it fits on one page anyway
  #  if (nexpl>5 || nresp>5) warning("Large number of plots to fit on one page, may fail: consider specifying (fewer) variables in
  #  arguments 'expl' and 'resp'; or set 'onepage=FALSE'")
  #  
   
    # try whether all plots can fit on the page, otherwise lower the number of coefficients one by one...
    par(mfrow=c(nexpl, nresp))
    a <- try(hist(1:R), silent=TRUE) ; devAskNewPage(ask = FALSE)
    if (class(a)=="try-error") {
      DoesNotFit <- TRUE
      giveWarning <- TRUE
    }
    else { 
      DoesNotFit <- FALSE
      giveWarning <- FALSE
    } 
    while (DoesNotFit & nresp>0) {
      if (nexpl > nresp) nexpl <- nexpl - 1
      else nresp <- nresp - 1
      par(mfrow=c(nexpl, max(nresp,1)))
      a <- try(hist(1:R), silent=TRUE); 
      if (class(a)!="try-error") DoesNotFit <- FALSE
    }
    if (nresp==0) stop("Something is wrong: plot margins too small?")
    if (giveWarning) warning("Number of plots too large to fit on the page, subset was selected: consider specifying (fewer) 
      variables in 'expl' and 'resp'; or enlarge graphics device; or set 'onepage=FALSE'")
    
    par(mfrow=c(nexpl, nresp))
    for(i in 1:nexpl) {
        for (j in 1:nresp) {
            vecposition <- (cinds[j] - 1)*p + rinds[i]
            if (confmethod=="basic") {
                lowerlim <- FRBres$bootest$CI.basic[vecposition,1]
                upperlim <- FRBres$bootest$CI.basic[vecposition,2]
            }
            else {
                lowerlim <- FRBres$bootest$CI.bca[vecposition,1]
                upperlim <- FRBres$bootest$CI.bca[vecposition,2]
            }
            if ((lowerlim * upperlim) > 0) {colorhist <- "red"; star="*"}
            else {colorhist <- NULL; ; star=""}
              
            uncentered <- FRBres$bootest$centered[vecposition,] + FRBres$bootest$vecest[vecposition]
            hist(uncentered, xlab=paste("Beta_",rinds[i],cinds[j],sep=""), col.main=colorhist, main=paste(resp[j],"~",expl[i],star))
            abline(v=FRBres$bootest$vecest[vecposition], col="blue", lwd=2, lty=3)
            abline(v=lowerlim, col="blue", lwd=2, lty=1)
            abline(v=upperlim, col="blue", lwd=2, lty=1)
        }
    }
    par(mfrow=(c(1,1)))
  }
  # or separate pages are used for each response variable
  else {
  
    if (nexpl<=1) rowcol <- c(1,1)
    else if (nexpl==2) rowcol <- c(2,1)
    else if (nexpl==3) rowcol <- c(2,2)
    else if (nexpl==4) rowcol <- c(2,2)
    else rowcol <- c(ceiling(nexpl/3),3)
    par(mfrow=rowcol)
    a <- try(hist(1:R), silent=TRUE);  devAskNewPage(ask = FALSE)
    if (class(a)=="try-error") {
      DoesNotFit <- TRUE
      giveWarning <- TRUE
    }
    else { 
      DoesNotFit <- FALSE
      giveWarning <- FALSE
    } 
    while (DoesNotFit & nexpl>0) {
      nexpl <- nexpl - 1
      if (nexpl<=1) rowcol <- c(1,1)
      else if (nexpl==2) rowcol <- c(2,1)
      else if (nexpl==3) rowcol <- c(2,2)
      else if (nexpl==4) rowcol <- c(2,2)
      else rowcol <- c(ceiling(nexpl/3),3)
      par(mfrow=rowcol)
      a <- try(hist(1:R), silent=TRUE)
      if (class(a)!="try-error") DoesNotFit <- FALSE
    }  
    if (nexpl==0) stop("Something is wrong: plot margins too small?")
  
    if (giveWarning) warning("Number of plots too large to fit on the pages, subset was selected: consider specifying (fewer) 
      variables in 'expl'; or enlarge graphics device")
  
    for(j in 1:nresp) {
      par(mfrow = rowcol)
      for (i in 1:nexpl) {
            vecposition <- (cinds[j] - 1)*p + rinds[i]
            if (confmethod=="basic") {
                lowerlim <- FRBres$bootest$CI.basic[vecposition,1]
                upperlim <- FRBres$bootest$CI.basic[vecposition,2]
            }
            else {
                lowerlim <- FRBres$bootest$CI.bca[vecposition,1]
                upperlim <- FRBres$bootest$CI.bca[vecposition,2]
            }
            if ((lowerlim * upperlim) > 0) colorhist <- "red2"
            else colorhist <- NULL
              
            uncentered <- FRBres$bootest$centered[vecposition,] + FRBres$bootest$vecest[vecposition]
            hist(uncentered, xlab=paste("Beta_",rinds[i],cinds[j],sep=""), col.main=colorhist, main=paste(resp[j],"~",expl[i]))
            abline(v=FRBres$bootest$vecest[vecposition], col="blue", lwd=2, lty=3)
            abline(v=lowerlim, col="blue", lwd=2, lty=1)
            abline(v=upperlim, col="blue", lwd=2, lty=1)
      }
      devAskNewPage(ask = TRUE) 
    }
  }
par(mfrow=(c(1,1)))
}
else warning("Could not plot confidence intervals; FRB was not performed") 

devAskNewPage(ask = currentAsk) 

}

