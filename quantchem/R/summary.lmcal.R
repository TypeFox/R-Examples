"summary.lmcal" <-
function (object,sort.models=FALSE,...) 
{

obj = object

durbin = function (obj) 
{
        mf <- model.frame(obj)
        y <- model.response(mf)
        X <- model.matrix(obj)
    n <- nrow(X)
    k <- ncol(X)
    resid <- residuals(obj)
    dw <- sum(diff(resid)^2)/sum(resid^2)
    ac <- (sum(resid[2:n] * resid[1:(n-1)]))/sum(resid^2)
    Q1 <- chol2inv(qr.R(qr(X)))
    AX <- matrix(as.vector(filter(X, c(-1, 2, -1))), ncol = k)
    AX[1, ] <- X[1, ] - X[2, ];  AX[n, ] <- X[n, ] - X[(n - 1), ]
    XAXQ <- t(X) %*% AX %*% Q1
    P <- 2 * (n - 1) - sum(diag(XAXQ))
    Q <- 2 * (3 * n - 4) - 2 * sum(diag(crossprod(AX) %*% 
          Q1)) + sum(diag(XAXQ %*% XAXQ))
    dmean <- P/(n - k)
    dvar <- 2/((n - k) * (n - k + 2)) * (Q - P * dmean)
    pval <-  pnorm(dw, mean = dmean, sd = sqrt(dvar))
    return(c(ac,dw,pval));
}

  inverse.power <- function(x, p) {
        if (p == 0) 
            exp(x)
        else (1 + p * x)^(1/p)
    }

coefficients = c()
residual = c()
fit = c()
sensitivity = c()
res = list()
variances = c();

xf = as.factor(obj$x);
xfp = (all(sapply(split(obj$x,xf),length) >= 2))

for (i in 1:length(obj$models)) {

if (inherits(obj$models[[i]],"lm")) {

coef = coefficients(summary(obj$models[[i]]));
rownames(coef)=c("x0","x1","x2","x3","x4")[1:dim(coef)[1]]
if (sort.models) { rownames(coef)=paste(names(obj$models)[i],rownames(coef)); }
else { rownames(coef)=paste(rownames(coef),names(obj$models)[i]); }
if (inherits(obj$models[[i]],"rlm")) { 
	coefficients = rbind(coefficients,cbind(coef,(1-pnorm(abs(coef[,3])))*2)); }
else { coefficients = rbind(coefficients,coef); }
}

if (
!is.null(attr(obj$models[[i]],"ly")) ||
!is.null(attr(obj$models[[i]],"bx")) ||
!is.null(attr(obj$models[[i]],"by"))
 ) { sens = NA; } else { sens = summary(obj$models[[i]])$sigma/derivative(obj$models[[i]],0); }

sha = shapiro.test(residuals(obj$models[[i]]));
residual = rbind(residual,c(
quantile(residuals(obj$models[[i]])),
sha$statistic,
sha$p.value));
if (is.null(summary(obj$models[[i]])$adj.r.squared)) { adj.r = NA; }
else { adj.r = summary(obj$models[[i]])$adj.r.squared; }
 
fit = rbind(fit,c(
summary(obj$models[[i]])$r.squared,
adj.r,
AIC(obj$models[[i]]),
summary(obj$models[[i]])$sigma
));
sensitivity = rbind(sensitivity,c(
sens,sens*3.3,sens*10,
durbin(obj$models[[i]])
));

if (xfp) {
b = bartlett.test(residuals(obj$models[[i]]) ~ xf)
variances=rbind(variances,c(
quantile(tapply(residuals(obj$models[[i]]),xf,var),c(0,0.5,1)),
b$statistic,b$p.value));

}
}

rownames(residual)=names(obj$models);
if (xfp) { variances=variances[c(1,9,12),];
rownames(variances)=c("Pure","Log","Box-Cox");
colnames(variances)=c("Min","Median","Max","K","Pr(>K)"); }
rownames(fit)=names(obj$models);
rownames(sensitivity)=names(obj$models);
colnames(fit)=c("R-Sq","Adj-R-Sq","AIC","Sigma");
colnames(sensitivity)=c("Sens","LOD","LOQ","Auto","DW","Pr(<DW)");

coefficients=coefficients[order(rownames(coefficients)),];

res$coefficients = coefficients;
res$residuals = residual;
res$variances = variances;
res$fit = cbind(fit,lof(obj));
res$sensitivity = sensitivity;
colnames(res$residuals)[7]="Pr(<W)"

cat(paste("\nOptimal Box-Cox power on x:",format(obj$px)));
cat(paste("\nOptimal Box-Cox power on y:",format(obj$py)));

cat(paste("\n\nOptimal power of weighting on x:",format(obj$wx)));
cat(paste("\nOptimal power of weighting on y:",format(obj$wy)));
cat(paste("\n\nMean relative error:\n",
		format(min(obj$weigh$mrex)),"on x\n",
		format(min(obj$weigh$mrey)),"on y\n",
		format(obj$weigh$mrex[41]),"unweighted"));
cat(paste("\n\nBetter weighting on y:",obj$yw,"\n"));

cat("\nCoefficients:\n");
printCoefmat(res$coefficients,tst.ind=1:3);
cat("\nResiduals:\n");
printCoefmat(res$residuals,tst.ind=1:6);
if (!is.null(res$variances)) { cat("\nVariances:\n");
printCoefmat(res$variances,tst.ind=1:3); }
cat("\nGoodness of fit:\n");
printCoefmat(res$fit,cs.ind=1:3,tst.ind=4:6);
cat("\nSensitivity and autocorrelation:\n");
printCoefmat(res$sensitivity,tst.ind=1:4);





invisible(res);

}

