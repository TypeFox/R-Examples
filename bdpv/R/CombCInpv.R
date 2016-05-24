CombCInpv <-
function(x0, x1, p, conf.level=0.95,
 alternative=c("two.sided", "less", "greater"))
{
if(any(x0==0)||any(x1==0))
{return(CIlnpvak(x0=x0, x1=x1, p=p, conf.level=conf.level,
 alternative=alternative))}
else{return(CIlnpv(x0=x0, x1=x1, p=p, conf.level=conf.level,
 alternative=alternative))
}}

