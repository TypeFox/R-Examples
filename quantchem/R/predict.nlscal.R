"predict.nlscal" <-
function (object,dataset,...) 
{

obj=object

res=list()

if (inherits(obj$models$a1,"nls")) {
a=summary(obj$models$a1)$parameters[,1][1]
r=summary(obj$models$a1)$parameters[,1][2]
l=summary(obj$models$a1)$parameters[,1][3]
res$fit$a1=-exp(-l)*log((dataset-a)/(r-a));
}

if (inherits(obj$models$a2,"nls")) {
a=summary(obj$models$a2)$parameters[,1][1]
l=summary(obj$models$a2)$parameters[,1][2]
res$fit$a3=-exp(-l)*log((a-dataset)/a);
}

if (inherits(obj$models$g1,"nls")) {
a=summary(obj$models$g1)$parameters[,1][1]
m=summary(obj$models$g1)$parameters[,1][2]
s=summary(obj$models$g1)$parameters[,1][3]
res$fit$g1=s*log(dataset/(a-dataset))+m;
}

if (inherits(obj$models$g2,"nls")) {
a=summary(obj$models$g2)$parameters[,1][1]
b=summary(obj$models$g2)$parameters[,1][2]
m=summary(obj$models$g2)$parameters[,1][3]
s=summary(obj$models$g2)$parameters[,1][4]
res$fit$g2=s*log((dataset-a)/(b-dataset))+m;
}

if (inherits(obj$models$m1,"nls")) {
v=summary(obj$models$m1)$parameters[,1][1]
k=summary(obj$models$m1)$parameters[,1][2]
res$fit$m1=k*dataset/(v-dataset);
}


pred <- function(x,model,val) { predict(model,newdata=data.frame(x=x))-val }
low =c()
for (i in 1:length(dataset)) {
root = try(uniroot(pred,interval=range(obj$x),model=obj$models$s1,val=dataset[i])$root);
if (inherits(root,"try-error")) { low=c(low,NA); } else { low=c(low,root); }
}
res$fit$s1=low

res$fit=as.data.frame(res$fit);

return(res);
}

