"lmcal" <-
function (x,y,confint=0.95,gridratio=0.05) 
{

if ((min(x) <= 0) || (min(y) <= 0 )) 
stop ("Both variables must be positive!");

box.cox <- function (x, p, start = 0) 
{
    if (p == 0) 
        log(x)
    else (x^p - 1)/p
}

  inverse.power <- function(x, p) {
        if (p == 0) 
            exp(x)
        else (1 + p * x)^(1/p)
    }

  weigh <- function (x,y,gamma=seq(-4,4,by=0.1)) {

mrex=c(); mrey=c();

for (i in 1:length(gamma)) {

fx = lm(y~x,weights=x^gamma[i]);
fy = lm(y~x,weights=y^gamma[i]);
mrex = c(mrex,mean(abs(((y-fx$coefficients[1])/fx$coefficients[2]-x)/x)));
mrey = c(mrey,mean(abs(((y-fy$coefficients[1])/fy$coefficients[2]-x)/x)));

}

res=list(); res$gamma=gamma; res$mrex=mrex; res$mrey=mrey;
res = as.data.frame(res);
}


res=list()

res$models$p1 = lm(y~x)
res$models$p2 = lm(y~x+I(x^2))
res$models$p3 = lm(y~x+I(x^2)+I(x^3))
res$models$p4 = lm(y~x+I(x^2)+I(x^3)+I(x^4))

res$weigh = weigh(x,y);

res$wx=res$weigh$gamma[which.min(res$weigh$mrex)];
res$wy=res$weigh$gamma[which.min(res$weigh$mrey)];

res$yw = if (min(res$weigh$mrey) < min(res$weigh$mrex)) { TRUE } else {FALSE}

res$weights = if (res$yw) { y^res$wy } else { x^res$wx }

res$models$P1 = update(res$models$p1,weights=res$weights);
res$models$P2 = update(res$models$p2,weights=res$weights);
res$models$P3 = update(res$models$p3,weights=res$weights);
res$models$P4 = update(res$models$p4,weights=res$weights);

res$models$l1 = lm(I(log10(y)) ~ I(log10(x)));
attr(res$models$l1,"ly")=TRUE;
res$models$l2 = lm(I(log10(y)) ~ I(log10(x)) + I(log10(x)^2)); 
attr(res$models$l2,"ly")=TRUE;

res$by = boxcox(y~x,lambda=seq(-4,4,0.1),plotit=FALSE)
res$bx = boxcox(x~y,lambda=seq(-4,4,0.1),plotit=FALSE)

py = res$by$x[which.max(res$by$y)];
px = res$bx$x[which.max(res$bx$y)];
by = boxcox(y~x,lambda=seq(py-0.1,py+0.1,0.01),plotit=FALSE)
bx = boxcox(x~y,lambda=seq(px-0.1,px+0.1,0.01),plotit=FALSE)
py = by$x[which.max(by$y)];
px = bx$x[which.max(bx$y)];
by = boxcox(y~x,lambda=seq(py-0.01,py+0.01,0.001),plotit=FALSE)
bx = boxcox(x~y,lambda=seq(px-0.01,px+0.01,0.001),plotit=FALSE)
py = by$x[which.max(by$y)];
px = bx$x[which.max(bx$y)];
by = boxcox(y~x,lambda=seq(py-0.001,py+0.001,0.0001),plotit=FALSE)
bx = boxcox(x~y,lambda=seq(px-0.001,px+0.001,0.0001),plotit=FALSE)
res$py = by$x[which.max(by$y)];
res$px = bx$x[which.max(bx$y)];

res$models$bx = lm(y~I(box.cox(x,res$px))); attr(res$models$bx,"bx")=TRUE;
res$models$by = lm(I(box.cox(y,res$py))~x); attr(res$models$by,"by")=TRUE;

res$f$p = lm(y~factor(x))
res$f$l = lm(I(log10(y))~factor(I(log10(x))))
res$f$b = lm(I(box.cox(y,res$py))~factor(x))

res$graph$grid=seq(min(x)-diff(range(x))*gridratio,max(x)+diff(range(x))*gridratio,length=100);

res$models$r1 = rlm(y~x,method="MM")
res$models$r2 = rlm(y~x+I(x^2),method="MM")
res$models$r3 = rlm(y~x+I(x^2)+I(x^3),method="MM")
res$models$r4 = rlm(y~x+I(x^2)+I(x^3)+I(x^4),method="MM")

res$models$R1 = rlm(y~x,method="MM",weights=res$weights)
res$models$R2 = rlm(y~x+I(x^2),method="MM",weights=res$weights)
res$models$R3 = rlm(y~x+I(x^2)+I(x^3),method="MM",weights=res$weights)
res$models$R4 = rlm(y~x+I(x^2)+I(x^3)+I(x^4),method="MM",weights=res$weights)

res$fitted=c()
for (i in 1:length(res$models)) {
fitc = predict(res$models[[i]],interval="confidence",level=0.95,newdata=data.frame(x=res$graph$grid));
fitp = predict(res$models[[i]],interval="prediction",level=0.95,newdata=data.frame(x=res$graph$grid));
if (!is.null(attr(res$models[[i]],"ly"))) { fitc = 10^fitc; fitp = 10^fitp; }
if (!is.null(attr(res$models[[i]],"by"))) { 
fitc = inverse.power(fitc,py); 
fitp = inverse.power(fitp,py); 
}

res$graph$fitted=cbind(res$graph$fitted,fitc[,1]);
res$graph$lowerc=cbind(res$graph$lowerc,fitc[,2]);
res$graph$upperc=cbind(res$graph$upperc,fitc[,3]);
res$graph$lowerp=cbind(res$graph$lowerp,fitp[,2]);
res$graph$upperp=cbind(res$graph$upperp,fitp[,3]);

}

res$graph$fitted=as.data.frame(res$graph$fitted);
res$graph$lowerc=as.data.frame(res$graph$lowerc);
res$graph$upperc=as.data.frame(res$graph$upperc);
res$graph$lowerp=as.data.frame(res$graph$lowerp);
res$graph$upperp=as.data.frame(res$graph$upperp);

names(res$graph$fitted)=names(res$models);
names(res$graph$upperc)=names(res$models);
names(res$graph$lowerc)=names(res$models);
names(res$graph$upperp)=names(res$models);
names(res$graph$lowerp)=names(res$models);

res$x=x
res$y=y

class(res) = c("lmcal","cal");

return(res);

}

