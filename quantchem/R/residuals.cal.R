"residuals.cal" <-
function (object,...) 
{

obj=object

res = c()

for (i in 1:length(obj$models)) {

if (inherits(obj$models[[i]],c("lm","nls","loess")))
{
res = cbind(res,residuals(obj$models[[i]]));
colnames(res)[ncol(res)]=names(obj$models)[i];
}
}

res=as.data.frame(res);


return(res);

}

