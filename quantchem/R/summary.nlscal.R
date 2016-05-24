"summary.nlscal" <-
function (object,sort.models=FALSE,...) 
{

obj = object

res=list()

durbin = function (obj) 
{
    resid <- residuals(obj)
    n <- length(resid)
    dw <- sum(diff(resid)^2)/sum(resid^2)
    ac <- (sum(resid[2:n] * resid[1:(n-1)]))/sum(resid^2)
    return(c(ac,dw));
}

xf = as.factor(obj$x);
xfp = (all(sapply(split(obj$x,xf),length) >= 2))

coefficients = NULL
variances = NULL
fit = NULL
residual = NULL
sens = NULL

for (i in 1:length(obj$models)) {

if (inherits(obj$models[[i]],"nls")) {

coef = coefficients(summary(obj$models[[i]]));
if (sort.models) { rownames(coef)=paste(names(obj$models)[i],rownames(coef)); }
else { rownames(coef)=paste(rownames(coef),names(obj$models)[i]); }
 coefficients = rbind(coefficients,coef);

        fit = rbind(fit,c(
                                AIC(obj$models[[i]]),
                                summary(obj$models[[i]])$sigma
                                ));
     }

if (inherits(obj$models[[i]],c("nls","loess"))) {

        sha = shapiro.test(residuals(obj$models[[i]]));
   residual = rbind(residual,c(
                                quantile(residuals(obj$models[[i]])),
                                sha$statistic,
                                sha$p.value
                                ));
rownames(residual)[nrow(residual)]=names(obj$models)[i]
}


}
coefficients=coefficients[order(rownames(coefficients)),];

        if (xfp) { 
                b = bartlett.test(residuals(obj$models[[1]]) ~ xf)
        variances=rbind(variances,c(
                quantile(tapply(residuals(obj$models[[1]]),xf,var),c(0,0.5,1)),
                b$statistic,b$p.value));

rownames(variances)="";
colnames(variances)=c("Min","Median","Max","K","Pr(>K)"); 
 }

res$residuals = residual;
colnames(res$residuals)[7]="Pr(<W)"


res$coefficients = coefficients;
res$variances = variances;
colnames(fit)=c("AIC","Sigma");
res$fit = cbind(fit,lof(obj));

if (inherits(obj$models$a1,"nls")) {
a=summary(obj$models$a1)$parameters[,1][1]
r=summary(obj$models$a1)$parameters[,1][2]
l=summary(obj$models$a1)$parameters[,1][3]
sen=summary(obj$models$a1)$sigma/(exp(l)*(a-r));
sens=rbind(sens,c(sen,3.3*sen,10*sen,durbin(obj$models$a1)));
rownames(sens)[nrow(sens)]="a1";
}

if (inherits(obj$models$a2,"nls")) {
a=summary(obj$models$a2)$parameters[,1][1]
l=summary(obj$models$a2)$parameters[,1][2]
sen=summary(obj$models$a2)$sigma/(a*exp(l));
sens=rbind(sens,c(sen,3.3*sen,10*sen,durbin(obj$models$a2)));
rownames(sens)[nrow(sens)]="a2";
}

if (inherits(obj$models$g1,"nls")) {
a=summary(obj$models$g1)$parameters[,1][1]
m=summary(obj$models$g1)$parameters[,1][2]
s=summary(obj$models$g1)$parameters[,1][3]
sen=summary(obj$models$g1)$sigma/
((a*exp(m/s))/(s*(exp(2*m/s)+2*exp(m/s)+1)));
sens=rbind(sens,c(sen,3.3*sen,10*sen,durbin(obj$models$g1)));
rownames(sens)[nrow(sens)]="g1";
}

if (inherits(obj$models$g2,"nls")) {
a=summary(obj$models$g2)$parameters[,1][1]
b=summary(obj$models$g2)$parameters[,1][2]
m=summary(obj$models$g2)$parameters[,1][3]
s=summary(obj$models$g2)$parameters[,1][4]
sen=summary(obj$models$g2)$sigma/
(exp(m/s)*(b-a)/(s*(exp(2*m/s)+2*exp(m/s)+1)));
sens=rbind(sens,c(sen,3.3*sen,10*sen,durbin(obj$models$g2)));
rownames(sens)[nrow(sens)]="g2";
}

if (inherits(obj$models$m1,"nls")) {
v=summary(obj$models$m1)$parameters[,1][1]
k=summary(obj$models$m1)$parameters[,1][2]
sen=summary(obj$models$m1)$sigma/(v/k);
sens=rbind(sens,c(sen,3.3*sen,10*sen,durbin(obj$models$m1)));
rownames(sens)[nrow(sens)]="m1";
}

colnames(sens) = c("Sens","LOD","LOQ","Auto","DW");

res$sens=sens;

cat("\nCoefficients:\n");
printCoefmat(res$coefficients,tst.ind=1:3);
cat("\nResiduals:\n");
printCoefmat(res$residuals,tst.ind=1:6);
if (!is.null(res$variances)) { cat("\nVariances:\n");
printCoefmat(res$variances); }
cat("\nGoodness of fit:\n");
printCoefmat(res$fit,tst.ind=2:4);
cat("\nSensitivity and autocorrelation:\n");
printCoefmat(res$sens,tst.ind=1:5);

invisible(res);

}

