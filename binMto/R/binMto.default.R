"binMto.default" <-
function(x, n, names=NULL, base=1, conf.level=0.95, alternative="two.sided",  method="Add4", adj="Dunnett", ...)
{

alternative <- match.arg(arg=alternative, choices=c("two.sided", "less", "greater"))
method <- match.arg(arg=method, choices=c("Add4", "Add2", "NHS", "Wald"))
adj <- match.arg(arg=adj, choices=c("Dunnett", "Bonf", "Unadj"))

# include warnings:

if(length(x)!=length(n))
 {stop("number of successes x and number of trials n must be of same length")}

if(is.null(names)==FALSE && length(names)!=length(n)) 
 {stop("number of groups n and number of names must be the same")}

if(is.null(names)){names<-paste(1:length(n))}

if(length(base)!=1)
 {stop("specify only one control group")}

if(is.numeric(base)==FALSE)
 {stop("base must be numeric")}

if(base>length(n) || base<1)
 {stop("argument base mis-specified: base greater than number of groups")}

if(conf.level>1 || conf.level<0 || length(conf.level)>1)  
 {stop("conf.level must be a single number between 0 and 1, usually 0.95")}

xC <- x[base]
xT <- x[-base]

nC <- n[base]
nT <- n[-base]

namesC <- names[base]
namesT <- names[-base]


pC <- xC/nC
pT <- xT/nT

esti <- pT-pC

k <- length(nT)
# k= number of treatment groups


quant <- Mtoquant(nc=nC, nx=nT, pc=pC, px=pT, conf.level=conf.level, adj=adj, alternative=alternative)



upper <- numeric(length=k)
lower <- numeric(length=k)
compnames <- character(length=k)


if (method=="Add4")
{
 for (i in 1:k) 
  {
   temp <- Add4(nx=nT[i], X=xT[i], ny=nC, Y=xC, quantile=quant, alternative=alternative)
   lower[i] <- temp$conf.int[1]
   upper[i] <- temp$conf.int[2]
   compnames[i] <- paste(namesT[i], namesC, sep="-")
  }
}

if (method=="Add2")
{
 for (i in 1:k) 
  {
   temp <- Add2(nx=nT[i], X=xT[i], ny=nC, Y=xC, quantile=quant, alternative=alternative)
   lower[i] <- temp$conf.int[1]
   upper[i] <- temp$conf.int[2]
   compnames[i] <- paste(namesT[i], namesC, sep="-")
  }
}

if (method=="NHS")
{
 for (i in 1:k) 
  {
   temp <- NHS(nx=nT[i], X=xT[i], ny=nC, Y=xC, quantile=quant, alternative=alternative)
   lower[i] <- temp$conf.int[1]
   upper[i] <- temp$conf.int[2]
   compnames[i] <- paste(namesT[i], namesC, sep="-")
  }
}

if (method=="Wald")
{
 for (i in 1:k) 
  {
   temp <- Wald(nx=nT[i], X=xT[i], ny=nC, Y=xC, quantile=quant, alternative=alternative)
   lower[i] <- temp$conf.int[1]
   upper[i] <- temp$conf.int[2]
   compnames[i] <- paste(namesT[i], namesC, sep="-")
  }
}

conf.int <- cbind(esti, lower, upper)
colnames(conf.int) <- c( "estimate", "lower", "upper")
rownames(conf.int) <- compnames


out.list<- list(
conf.int=conf.int,
conf.level=conf.level,
quantile=quant,
alternative=alternative,
adj=adj,
method=method
)

class(out.list) <- "binMto" 

return(out.list)


}

