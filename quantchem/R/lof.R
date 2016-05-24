"lof" <-
function(obj) 
{
res=c()
for (i in 1:length(obj$models)) {

if (inherits(obj$models[[i]],c("lm","nls"))) {

if (!is.null(attr(obj$models[[i]],"ly"))) {
lf = anova(obj$models[[i]],obj$f$l);
} else
if (!is.null(attr(obj$models[[i]],"by"))) {
lf = anova(obj$models[[i]],obj$f$b);
} else {
lf = anova(obj$models[[i]],obj$f$p);
}
res=rbind(res,c(lf[1,2],lf[2,2],lf[2,5],lf[2,6]));
rownames(res)[nrow(res)]=names(obj$models)[i];
}
}
colnames(res)=c("SSR","SSPE", "F","Pr(>F)");
return(res);
}

