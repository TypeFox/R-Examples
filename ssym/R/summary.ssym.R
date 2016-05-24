summary.ssym <-
function(object, ...){

if(object$family=='Normal'){
cat("\n     Family: ",object$family,"\n")}
else{
    if(object$family=='Contnormal' | object$family=='Sinh-t')
       cat("\n     Family: ",object$family,"(",object$xi[1],",",object$xi[2],")\n")
	else cat("\n     Family: ",object$family,"(",object$xi[1],")\n")
}

cat("Sample size: ",length(object$y),"\n")

if(object$censored==FALSE){
	cat(" Quantile of the Weights\n")
	temp <- round(quantile(object$weights),digits=2)
	cmat <- cbind(temp[1], temp[2], temp[3], temp[4], temp[5])
	colnames(cmat) <- c("0%", "25%", "50%", "75%", "100%")
	rownames(cmat) <- ""
	printCoefmat(cmat,digits=3)
}else cat("Censored %: ",round(100*mean(object$event),2),"\n")

cat("\n ************************* Median submodel ****************************\n")
if(object$orig=="linear") cat("link: ", attr(object$l1.mu,"link"),"\n")
if(object$p>0){
TAB		 <- cbind(Estimate <- object$theta.mu[1:object$p],
				  StdErr <- sqrt(diag(as.matrix(object$vcov.mu)))[1:object$p],
				  tval <- Estimate/StdErr,
				  p.value <- 2*pnorm(-abs(tval)))
colnames(TAB) <- c("Estimate", "Std.Err", "z-value", "Pr(>|z|)")
rownames(TAB) <- colnames(object$model.matrix.mu)[1:object$p]
  cat(" ******** Parametric component\n\n")		
printCoefmat(TAB, P.values=TRUE, has.Pvalue=TRUE,digits=5, signif.legend=FALSE, tst.ind=c(2,3))}

if(sum(object$qm) > 0){
TAB		 <- cbind(sp <- round(object$lambdas.mu,digits=5),
                  bd <- object$qm,
				  df <- round(object$gle.mu[(1+min(1,object$p)):(length(object$qm)+min(1,object$p))],digits=3),
				  Stat <- round(object$stes.mu[,1],digits=5),
				  p.value <- round(object$stes.mu[,2],digits=5))
colnames(TAB) <- c("Smooth.param", "Basis.dimen", "d.f.", "Statistic", "p-value")
rownames(TAB) <- colnames(object$model.matrix.mu)[(object$p+1):(length(object$qm)+object$p)]
  cat(" ******** Nonparametric component\n\n")		
printCoefmat(TAB, P.values=TRUE, has.Pvalue=TRUE,digits=4, signif.legend=FALSE,tst.ind=c(2,3))}

cat("\n **** Deviance: ",round(sum(object$deviance.mu),digits=2),"\n")

cat(" **************************** Dispersion submodel ***********************\n")
cat("link: ", attr(object$l1.phi,"link"),"\n")
if(object$l>0){
TAB		 <- cbind(Estimate <- object$theta.phi[1:object$l],
				  StdErr <- sqrt(diag(as.matrix(object$vcov.phi)))[1:object$l],
				  tval <- Estimate/StdErr,
				  p.value <- 2*pnorm(-abs(tval)))
colnames(TAB) <- c("Estimate", "Std.Err", "z-value", "Pr(>|z|)")
rownames(TAB) <- colnames(object$model.matrix.phi)[1:object$l]
  cat(" ******** Parametric component\n\n")		
printCoefmat(TAB, P.values=TRUE, has.Pvalue=TRUE,digits=5, signif.legend=FALSE,tst.ind=c(2,3))}

if(sum(object$q) > 0){
TAB		 <- cbind(sp <- round(object$lambdas.phi,digits=5),
                  bd <- object$q,
                  df <- round(object$gle.phi[(1+min(1,object$l)):(length(object$q)+min(1,object$l))],digits=3),
				  Stat <- round(object$stes.phi[,1],digits=5),
				  p.value <- round(object$stes.phi[,2],digits=5))
colnames(TAB) <- c("Smooth.param", "Basis.dimen", "d.f.", "Statistic", "p-value")
rownames(TAB) <- colnames(object$model.matrix.phi)[(object$l+1):(length(object$q)+object$l)]
  cat(" ******** Nonparametric component\n\n")		
printCoefmat(TAB, P.values=TRUE, has.Pvalue=TRUE,digits=4, signif.legend=FALSE,tst.ind=c(2,3))}

cat("\n **** Deviance: ",round(sum(object$deviance.phi),digits=2),"\n")

cat(" ************************************************************************\n")

if(object$censored==FALSE) temp2 <- qqnorm(qnorm(object$cdfz),plot.it=FALSE)
else{surv0 <- survfit(Surv(object$z_es,1-object$event)~1)
ids <- ifelse(surv0$n.event>0,TRUE,FALSE)
survs <- ifelse(1-surv0$surv[ids] < 1e-30,1 - 1e-30, 1- surv0$surv[ids])
survs <- ifelse(survs > 1 - 1e-15,1 - 1e-15, survs)
probs <- object$cdf(surv0$time[ids])
probs <- ifelse(probs < 1e-30,1e-30, probs)
probs <- ifelse(probs > 1 - 1e-15,1 - 1e-15, probs)
temp2 <- list(x=qnorm(survs),y=qnorm(probs))}
cat(" Overall goodness-of-fit statistic: ",round(mean(abs(sort(temp2$x)-sort(temp2$y))),digits=6),"\n")
cat("                 -2*log-likelihood: ",round(-2*sum(object$lpdf),digits=3),"\n")
cat("                               AIC: ",round(-2*sum(object$lpdf) + 2*(sum(object$gle.mu)+sum(object$gle.phi)),digits=3),"\n")
cat("                               BIC: ",round(-2*sum(object$lpdf) + log(length(object$y))*(sum(object$gle.mu)+sum(object$gle.phi)),digits=3),"\n")
}
