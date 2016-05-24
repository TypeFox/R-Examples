"print.nlscal" <-
function (x,...) 
{

obj = x

res=list()
coefficients = c()
fit = c()

for (i in 1:length(obj$models)) {

if (inherits(obj$models[[i]],c("nls"))) {

coef = coefficients(obj$models[[i]]);
coef = c(coef,rep(NA,4-length(coef)));
coefficients = rbind(coefficients,coef); 
rownames(coefficients)[nrow(coefficients)]=paste(names(obj$models)[i],names(coef)[1],
names(coef)[2],names(coef)[3],names(coef)[4],sep="-");
fit = rbind(fit,c(
AIC(obj$models[[i]]),
summary(obj$models[[i]])$sigma
));
rownames(fit)[nrow(fit)]=names(obj$models)[i];
}

}

colnames(fit) = c("AIC","Sigma");
colnames(coefficients) = c("","","","");

res$coefficients = coefficients;
res$fit = fit;

cat("\nCoefficients:\n");
printCoefmat(res$coefficients);
cat("\nGoodness of fit:\n");
printCoefmat(res$fit);

}

