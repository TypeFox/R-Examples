"print.lmcal" <-
function (x,...) 
{

obj = x

coefficients = c()
fit = c()

for (i in 1:length(obj$models)) {

if (inherits(obj$models[[i]],"lm")) {

coef = coefficients(obj$models[[i]]);
coef = c(coef,rep(NA,5-length(coef)));
coefficients = rbind(coefficients,coef); 
        if (is.null(summary(obj$models[[i]])$adj.r.squared)) { adj.r = NA; }
        else { adj.r = summary(obj$models[[i]])$adj.r.squared; }
fit = rbind(fit,c(
summary(obj$models[[i]])$r.squared,
adj.r,
AIC(obj$models[[i]]),
summary(obj$models[[i]])$sigma
));
}

}
rownames(coefficients)=names(obj$models);
colnames(coefficients)=c("x0","x1","x2","x3","x4");
rownames(fit)=names(obj$models);
colnames(fit)=c("R-Sq","Adj-R-Sq","AIC","Sigma");

#coefficients=coefficients[order(rownames(coefficients)),];

res=list();

res$coefficients = coefficients;
res$fit = fit;

cat(paste("\nOptimal Box-Cox power on x:",obj$px));
cat(paste("\nOptimal Box-Cox power on y:",obj$py));

cat(paste("\n\nOptimal power of weighting on x:",obj$wx));
cat(paste("\nOptimal power of weighting on y:",obj$wy));
cat(paste("\nBetter weighting on y:",obj$yw,"\n"));

cat("\nCoefficients:\n");
printCoefmat(res$coefficients);
cat("\nGoodness of fit:\n");
printCoefmat(res$fit);

}

