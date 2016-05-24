"predict.lmcal" <-
function (object,dataset,conf.int = 0.95,...) 
{

obj = object

inverse.predict <- function(object,resp) {
	pre = c();
	for (i in 1:length(resp)) {
	coef = coefficients(object);
	coef[1]=coef[1]-resp[i];
	root = try(polyroot(coef));
	if (inherits(root,"try-error")) root = NA;
	if (!is.null(attr(object,"ly"))) xrange = range(log10(abs(obj$graph$grid)))
	else if (!is.null(attr(object,"bx"))) xrange = range(box.cox(obj$graph$grid,obj$px))
	else xrange = range(obj$x);
	root = root[which(Re(root)>xrange[1]&& Re(root)<xrange[2]&& abs(Im(root))<1e-10)];
	root = Re(root[1])
	pre=c(pre,root);
	}
	return(pre);
}

  inverse.power <- function(x, p) {
        if (p == 0) 
            exp(x)
        else (1 + p * x)^(1/p)
    }
box.cox <- function (x, p, start = 0) 
{
    if (p == 0) 
        log(x)
    else (x^p - 1)/p
}

res=list();
fit=c()
se=c()
upper=c();
lower=c();

for (i in 1:length(obj$models)) {

if (!is.null(attr(obj$models[[i]],"ly"))) datacurrent = log10(dataset) else
if (!is.null(attr(obj$models[[i]],"by"))) datacurrent = box.cox(dataset,obj$py)
else datacurrent = dataset;

f = inverse.predict(obj$models[[i]],datacurrent);

se = try(predict(obj$models[[i]],newdata=data.frame(x=f),se.fit=TRUE)$se/derivative(obj$models[[i]],f));

if (inherits(se,"try-error")) se=NA;

if (!is.na(df.residual(obj$models[[i]]))) {
up = f + qt((conf.int+1)/2,df.residual(obj$models[[i]]))*se;
lo = f - qt((conf.int+1)/2,df.residual(obj$models[[i]]))*se;
}
else{
up = f + qnorm((conf.int+1)/2)*se;
lo = f - qnorm((conf.int+1)/2)*se;
}


if (!is.null(attr(obj$models[[i]],"bx"))) { 
f=inverse.power(f,obj$px); 
up=inverse.power(up,obj$px);
lo=inverse.power(lo,obj$px);
}
if (!is.null(attr(obj$models[[i]],"ly"))) { f=10^f; up=10^up; lo=10^lo; }

fit = cbind(fit,f);
upper = cbind(upper,up);
lower = cbind(lower,lo);

}



res$fit=as.data.frame(fit);
res$upper=as.data.frame(upper);
res$lower=as.data.frame(lower);

names(res$fit)=names(obj$models);
names(res$upper)=names(obj$models);
names(res$lower)=names(obj$models);

return(res);

}

