plotFRBangles <- function(x, pcs=1:min(12,length(x$eigval))) {

currentAsk <- devAskNewPage(ask = NULL)

FRBres <- x
q <- length(FRBres$eigval)
if (any(!(pcs %in% c(1:q)))) stop(paste("indices in 'pcs' should be in 1:",q,sep=""))

bootangles <- FRBres$angles
R <- ncol(bootangles)
nplots <- length(pcs)

if (nplots<=1) rowcol <- c(1,1)
else if (nplots==2) rowcol <- c(1,2)
else if (nplots==3) rowcol <- c(1,3)
else if (nplots==4) rowcol <- c(2,2)
else rowcol <- c(ceiling(nplots/3),3)
par(mfrow=rowcol)
a <- try(hist(1:R), silent=TRUE)
devAskNewPage(ask = FALSE)
if (class(a)=="try-error") {
   DoesNotFit <- TRUE
   giveWarning <- TRUE
}
else { 
   DoesNotFit <- FALSE
   giveWarning <- FALSE
} 
while (DoesNotFit & nplots>0) {
  nplots <- nplots - 1
  if (nplots<=1) rowcol <- c(1,1)
  else if (nplots==2) rowcol <- c(1,2)
  else if (nplots==3) rowcol <- c(1,3)
  else if (nplots==4) rowcol <- c(2,2)
  else rowcol <- c(ceiling(nplots/3),3)
  par(mfrow=rowcol)
  a <- try(hist(1:R), silent=TRUE)
  if (class(a)!="try-error") DoesNotFit <- FALSE
}
if (nplots==0) stop("Something is wrong: plot margins too small?")

if (giveWarning) warning("Number of plots too large to fit on the page, subset was selected: consider specifying (fewer) 
    variables in 'pcs'; or enlarge graphics device")

histbreaks <- (0:20)/20*pi/2

par(mfrow=rowcol)
for (i in 1:nplots) {
    hist(bootangles[pcs[i],], main=paste("PC",pcs[i],sep=""), xlab="angle(sample,resample)", breaks=histbreaks)
    abline(v=pi/2, lwd=2, col="red")
    abline(v=0, lwd=2, col="red")
}

devAskNewPage(ask = currentAsk)

}
