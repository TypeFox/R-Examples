"confint.cal" <-
function (object,parm,level=0.95,sort.models=FALSE,...) 
{

obj=object

conftable=c()

for (i in 1:length(obj$models)) {

if (inherits(obj$models[[i]],c("lm","nls"))) {

coef = confint(obj$models[[i]],level=level);
if (inherits(obj$models[[i]],"lm")) rownames(coef)=c("x0","x1","x2","x3","x4")[1:dim(coef)[1]]
if (sort.models) { rownames(coef)=paste(names(obj$models)[i],rownames(coef)); }
else { rownames(coef)=paste(rownames(coef),names(obj$models)[i]); }

if (inherits(obj$models[[i]],"rlm")) {
coef[,1] = coef(summary(obj$models[[i]]))[,1] - pnorm((level+1)/2)*coef(summary(obj$models[[i]]))[,2]
coef[,2] = coef(summary(obj$models[[i]]))[,1] + pnorm((level+1)/2)*coef(summary(obj$models[[i]]))[,2]
}
conftable=rbind(conftable,coef);
}
}

conftable=conftable[order(rownames(conftable)),];

return(conftable);

}

