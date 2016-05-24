cont_flm_plot <-
function(smoothdata, matchresults, flmresults, xlim, ylim, ftest, 
nperm, lb, xat, legendx, legendy, L, xlab="Time", ylab="Activity"){

if(missing(smoothdata) || missing(matchresults) || missing(flmresults) || missing(xlim) || 
missing(ylim) || missing(ftest) || missing(nperm) || missing(lb) || missing(xat) || 
missing(legendx) || missing(legendy) || missing(L)) 
stop("Missing Arguments")

colort <- factor(matchresults$cov[,3])
ucont <- length(unique(colort))
covname <- names(matchresults$cov[3])
maintitle <- paste(ylab, "~", covname, sep="")
LofLegend <-ifelse(length(levels(colort)) <= 100, 10, 1)

minmedmax <- function(contvar){
contvar <- contvar
mincont <- signif(min(contvar), 3)
medcont <- signif(median(contvar), 3)
maxcont <- signif(max(contvar), 3)
contlab <- list(mincont, medcont, maxcont)
return(contlab)
}

contlabels <- minmedmax(matchresults$cov[,3])

par(mfrow=c(2,1), mar=c(4,4,3,1))
plot(0, 0, xlim=xlim, ylim=ylim, xaxt="n", xlab=xlab, 
ylab=ylab, type='n', main=maintitle)
for(i in 1:length(colort)) 
lines(predict(flmresults$freg$yhatfdobj, c(1:L))[,i], col=topo.colors(ucont)[colort[i]])
        
axis(1, at=xat, labels=lb)
colorsSamples <- topo.colors(ucont*LofLegend)
pnts <- cbind(x=c(legendx, (legendx + xlim[2]/20), (legendx + xlim[2]/20), legendx), y=c(legendy-ylim[2]/5, legendy, legendy, legendy-ylim[2]/5))
SDMTools::legend.gradient(pnts, colorsSamples, c(contlabels[[1]], contlabels[[3]]), paste(covname, "Value"))
        
geftFtestresults <- flm_ftest(smoothdata, "fourier", smoothdata$fd$fd$basis$nbasis, 4, ftest, nperm, xat, lb, 1.5)
        
return(geftFtestresults)
}
