ind.contrib <-
function(model,print.diff=FALSE,graph=TRUE,warning=25) {
  if ("lm" %in% class(model)) {
    coeffs <- model$coefficients
    coeffs.diff <- lm.influence(model)$coefficients
  } else if ("least.rect" %in% class(model)) {
    coeffs <- model$coefficients
    coeffs.mat <- matrix(ncol=2,nrow=nrow(model$model),dimnames=list(1:nrow(model$model),c("(Intercept)",colnames(model$model)[2])))
    x <- model$model[,2]
    y <- model$model[,1]
    for (i in 1:nrow(model$model)) {
	x2 <- x[-i]
	y2 <- y[-i]
	coeffs2 <- least.rect(y2~x2)$coefficients
	coeffs.mat[i,1] <- coeffs2[1]
	coeffs.mat[i,2] <- coeffs2[2]
    }
    coeffs.diff <- data.frame(coeffs[1]-coeffs.mat[,1],coeffs[2]-coeffs.mat[,2])
    colnames(coeffs.diff)[1] <- "(Intercept)"
    colnames(coeffs.diff)[2] <- colnames(model$model)[2]
  } else {
    stop("model not recognized")
  }
  coeffs.prop <- 100*coeffs.diff/coeffs
  if (graph) {
    plot(coeffs.prop[,1],ylim=c(1.1*min(coeffs.prop),1.1*max(coeffs.prop)),type="o",pch=16,cex=0.5,xlab="Individual",ylab="Difference in parameters (%)")
    abline(h=0,col="grey",lty=3)
    abline(h=-100,col="grey",lty=3)
    abline(h=warning,col="grey",lty=3)
    abline(h=-warning,col="grey",lty=3)
    abline(h=100,col="grey",lty=3)
    for (i in 2:ncol(coeffs.prop)) {
	lines(coeffs.prop[,i],col=i,type="o",pch=16,cex=0.5)
    }
    legend(0.75*nrow(model$model),1.05*max(coeffs.prop),colnames(coeffs.prop),col=1:ncol(coeffs.prop),lty=1)
    lignes <- which(abs(coeffs.prop)>warning)
    if (length(lignes)>0) {
	colonnes <- lignes%/%nrow(coeffs.prop)+1
	lignes <- lignes-(colonnes-1)*nrow(coeffs.prop)
	ecart <- abs(1.1*max(coeffs.prop)-1.1*min(coeffs.prop))*4/100
	for (i in 1:length(lignes)) {
	  text(lignes[i],coeffs.prop[lignes[i],colonnes[i]]+ecart*sign(coeffs.prop[lignes[i],colonnes[i]]),lignes[i],cex=0.5)
	}
    }
  }
  if (print.diff) {print(coeffs.prop,digits=5)}
  result <- list(coefficients=coeffs,coefficients.diff=coeffs.diff,coefficients.prop=coeffs.prop)
  return(result)
}
