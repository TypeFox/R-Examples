"AIC.cal" <-
function(object,...,k = 2) 
{
obj=object;
res=c()
for (i in 1:length(obj$models)) {
if (inherits(obj$models[[i]],c("lm","nls"))) { res=rbind(res,AIC(obj$models[[i]],k=k)); }
else { res=rbind(res,NA); }
}
colnames(res)="AIC";
rownames(res)=names(obj$models);

return(res);

}

