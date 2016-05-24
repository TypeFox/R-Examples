"t.test.ratio.formula" <-
function(formula, data, base=2,...)
{

if(length(formula)!=3)
 {stop("formula must be a two-sided formula with one numeric response variable and one one factor")}

if( all(c(1,2)!=base) )
 {stop("base must be one of 1,2")}

mf<-model.frame(formula,data=data)

if (!is.numeric(mf[,1]))
 {stop("response variable must be numeric")}

datalist=split(x=mf[,1], f=droplevels(mf[,2]), drop=TRUE)

groupnames<-names(datalist)

if(length(datalist)!=2)
 {stop("grouping variable must have exactly 2 levels")}

args<-list(...)

args$x <- datalist[[-base]]
args$y <- datalist[[base]]
args$namex <- groupnames[-base]
args$namey <- groupnames[base]



out <- do.call("t.test.ratio.default", args)

return(out)
}

