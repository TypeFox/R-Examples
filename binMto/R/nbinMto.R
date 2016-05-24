"nbinMto" <-
function(Ntotal=500, pH1 , ratio=1, alpha=0.05, power=0.8, alternative="two.sided", method="Add4", 
trace=FALSE)
{

pc <- pH1[1]
px <- pH1[-1]

k=length(px)

if(length(Ntotal)==1)
{
 Nstop<-Ntotal
 Nstart<-max(2,2*ratio)+2*k
  if(Nstop<Nstart) {stop("Specify larger Ntotal")}
}
else
{
 if(length(Ntotal)==2)
 {
 Nstop<-max(Ntotal)
 Nstart<-max(min(Ntotal),(max(2,2*ratio)+2*k) )
  if(Nstop<Nstart) {stop("Specify larger min of Ntotal")}
 }
  else{stop("Ntotal mis-specified: either the maximal total sample size or the range of allowed total sample size")}
}


if(alternative=="less" & any(px-pc > 0)) 
{cat("Arguments pH1 and alternative mis-specified","\n", "pi-p0 should be less than 0 for alternative='less'", "\n")}

if(alternative=="greater" & any(px-pc < 0))
{cat("Arguments pH1 and alternative mis-specified","\n", "pi-p0 should be greater than 0 for alternative='greater'", 
"\n")}


nvec<-numeric(length=k+1)



n <- max(round(Nstart/(ratio+k)), 2)

nvec <- c( max(round(ratio*n),2),rep(n, times=k) )

powit <- apprPower(n=nvec, pH1=pH1, alpha=alpha, alternative=alternative, method=method)

it.hist <- list()

i <- 0

while(powit<power && sum(nvec)<Nstop)
{
i <- i+1
n <- n+1
nvec <- c(max(round(ratio*n),2),rep(n, times=k) )
powit <- apprPower(n=nvec, pH1=pH1, alpha=alpha, alternative=alternative, method=method)
it.hist[[i]] <- c(nvec, sum(nvec), powit) 
}

if(trace==FALSE)
{
out <- it.hist[[i]]
names(out) <- c( paste("n",0:k,sep="") ,"Ntotal", "power" )
}

if(trace==TRUE)
{
print(it.hist)
out<-matrix(it.hist[[1]], nrow=1)
for(e in 1:length(it.hist))
 {
  out<-rbind(out,it.hist[[e]])
 }
colnames(out) <- c( paste("n",0:k,sep="") ,"Ntotal", "power" )
}

return(out)

}

