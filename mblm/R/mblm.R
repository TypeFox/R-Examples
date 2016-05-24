"mblm" <-
function (formula,dataframe,repeated=TRUE) 
{
 if (missing(dataframe)) 
        dataframe <- environment(formula)

term<-as.character(attr(terms(formula),"variables")[-1]);
x=dataframe[[term[2]]];
y=dataframe[[term[1]]];

if(length(term) > 2) { stop("Only linear models are accepted"); }

xx = sort(x)
yy = y[order(x)]
n = length(xx)

slopes = c()
intercepts = c()
smedians = c()
imedians = c()

if (repeated) {

for (i in 1:n) {
	slopes = c()
	intercepts = c()
	for (j in 1:n) {
		if (xx[j] != xx[i]) { slopes = c(slopes,(yy[j]-yy[i])/(xx[j]-xx[i]));
					    intercepts = c(intercepts,(xx[j]*yy[i]-xx[i]*yy[j])/(xx[j]-xx[i])); }
	}
		smedians = c(smedians,median(slopes));
		imedians = c(imedians,median(intercepts));
	}

	slope = median(smedians);
	intercept = median(imedians);
	}

else	{

	for (i in 1:(n-1)) {
		for (j in i:n) {
			if (xx[j] != xx[i]) { slopes = c(slopes,(yy[j]-yy[i])/(xx[j]-xx[i])); }
		}
	}

	slope = median(slopes);
	intercepts = yy - slope*xx;
	intercept = median(intercepts);

}


res=list();

res$coefficients=c(intercept,slope);
names(res$coefficients)=c("(Intercept)",term[2]);

res$residuals=y-slope*x-intercept;
names(res$residuals)=as.character(1:length(res$residuals));

res$fitted.values=x*slope+intercept;
names(res$fitted.values)=as.character(1:length(res$fitted.values));

if (repeated) {
	res$slopes = smedians;
	res$intercepts = imedians;
	}
	else	{
	res$slopes = slopes;
	res$intercepts = intercepts;
	}

res$df.residual=n-2;
res$rank=2;
res$terms=terms(formula);
res$call=match.call();
res$model=data.frame(y,x);

res$assign=c(0,1);

  if (missing(dataframe)) {
    res$effects = lm(formula)$effects;
    res$qr = lm(formula)$qr;
  } else {
    res$effects = lm(formula, dataframe)$effects;
    res$qr = lm(formula, dataframe)$qr;
  }

res$effects[2]=sqrt(sum((res$fitted-mean(res$fitted))^2));
res$xlevels=list();

names(res$model)=term;
attr(res$model,"terms")=terms(formula);

class(res)=c("mblm","lm");

res
			
}

