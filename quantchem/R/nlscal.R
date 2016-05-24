"nlscal" <-
function (x,y,confint=0.95,gridratio=0.05) 
{

if ((min(x) <= 0) || (min(y) <= 0 )) 
stop ("Both variables must be positive!");

res=list()

res$models$a1 = try(nls(y~SSasymp(x,A,r0,lrc)));
res$models$a2 = try(nls(y~SSasympOrig(x,A,lrc)));
res$models$g1 = try(nls(y~SSlogis(x,A,xmid,scal)));
res$models$g2 = try(nls(y~SSfpl(x,A,B,xmid,scal)));
res$models$m1 = try(nls(y~SSmicmen(x,V,k0)));
res$models$s1 = try(loess(y~x,family="symmetric",surface="direct"));

res$graph$grid=seq(min(x)-diff(range(x))*gridratio,max(x)+diff(range(x))*gridratio,length=100);

res$graph$fitted$a1 = try(predict(res$models$a1,newdata=data.frame(x=res$graph$grid)));
res$graph$fitted$a2 = try(predict(res$models$a2,newdata=data.frame(x=res$graph$grid)));
res$graph$fitted$g1 = try(predict(res$models$g1,newdata=data.frame(x=res$graph$grid)));
res$graph$fitted$g2 = try(predict(res$models$g2,newdata=data.frame(x=res$graph$grid)));
res$graph$fitted$m1 = try(predict(res$models$m1,newdata=data.frame(x=res$graph$grid)));
res$graph$fitted$s1 = try(predict(res$models$s1,newdata=data.frame(x=res$graph$grid)));

res$x=x
res$y=y

res$f$p=lm(y~factor(x))

class(res) = c("nlscal","cal");

return(res);

}

