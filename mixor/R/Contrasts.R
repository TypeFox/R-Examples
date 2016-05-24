Contrasts <-
function(fit, contrast.matrix, digits=max(3, getOption("digits") - 2), signif.stars=TRUE, dig.tst = max(1, min(5, digits - 1))) {
        thresholds<-fit$Model[,"Estimate"]%*%contrast.matrix
        thresholds.se<-sqrt(diag(t(contrast.matrix)%*%fit$varcov%*%contrast.matrix))
        Z.thresh <- thresholds/thresholds.se
        P.thresh <- 2*(1-pnorm(abs(Z.thresh)))
        Thresholds<- cbind(Estimate=as.numeric(thresholds), SE=thresholds.se, Z=as.numeric(Z.thresh), P.value=as.numeric(P.thresh))
		dimnames(contrast.matrix)[[1]]<-dimnames(fit$Model)[[1]]
        dimnames(contrast.matrix)[[2]]<-1:dim(contrast.matrix)[2]
        dimnames(Thresholds)[[1]]<-dimnames(contrast.matrix)[[2]]
        dimnames(Thresholds)[[2]] <- c("Estimate", "Std. Error", "z value", "P(>|z|)")
print(contrast.matrix) 
cat("\n\n")
printCoefmat(Thresholds, digits=digits, signif.stars=signif.stars, dig.tst=dig.tst, cs.ind=1:2, tst.ind=3, Pvalues=TRUE, has.Pvalue=TRUE)
}
