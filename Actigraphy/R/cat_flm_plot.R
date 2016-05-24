cat_flm_plot <-
function(smoothdata, matchresults, flmresults, ftest, nperm, lb, xat, 
varname, col, ylim, L, xlab="Time", ylab="Activity"){

if(missing(smoothdata) || missing(matchresults) || missing(flmresults) || missing(ftest) || 
missing(nperm) || missing(lb) || missing(varname) || missing(col) || missing(ylim) || missing(L)) 
stop("Missing Arguments")

geft <- flmresults
maintitle<- paste(ylab, "~", varname, sep="")
leg <- names(matchresults$cov)[-1]
l <- length(geft$freg$betaestlist)
beta <- geft$freg$betaestlist
std <- geft$fregstd$betastderrlist

par(mfrow=c(2,1), mar=c(4,4,3,1))
plot(0, 0, xlim=c(0, L), ylim=ylim, xlab=xlab, ylab=ylab, type="n", main=maintitle, xaxt="n", cex.main=0.8)
lines(beta[[1]]$fd, col=col[1], lwd=2)
for(i in 2:l) 
lines(beta[[1]]$fd + beta[[i]]$fd, col=col[i], lwd=2)

axis(1, at=xat, labels=lb)
legend("topleft", leg, col=col, lty=rep(1, length(beta)), cex=0.8, lwd=2)

geftFtestresults <- flm_ftest(smoothdata, smoothdata$fd$fd$basis[[2]], smoothdata$fd$fd$basis$nbasis, 4, ftest, nperm, xat, lb, 1.5)

return(geftFtestresults)
}
