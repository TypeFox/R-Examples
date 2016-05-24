"derivative" <-
function(obj,x) {
der = 0;
coef = obj$coefficients[-1];
for (i in 1:length(coef)) { der = der + coef[i]*i*x^(i-1); }
return(der);
}

