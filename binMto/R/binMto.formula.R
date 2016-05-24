"binMto.formula" <-
function(formula, data, base=1, conf.level=0.95, alternative="two.sided",  method="Add4", adj="Dunnett", ...)
{

if(length(formula)!=3)
 {stop("formula must be something like 'response ~ treatment'")}
if(length(formula[[3]])!=1)
 {stop("There should be only ONE grouping variable in formula, e.g. 'response ~ treatment'")}

mf <- model.frame(formula=formula, data=data)

if(any (mf[,1] !=0 & mf[,1] !=1 ))
 {stop("response variable should have values 0,1 only")}

grouplist <- split(mf[,1], as.factor(mf[,2]), drop=TRUE)

names<- names(grouplist)

n <- as.numeric(lapply(grouplist, length))
x <- as.numeric(lapply(grouplist, sum))

return(binMto.default(x=x, n=n, names=names, base=base, conf.level=conf.level, alternative=alternative,  method=method, adj=adj))
}

