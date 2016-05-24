"ec.mto" <-
function(n,x, alternative="less")

{
if(length(x)!=length(n)){stop("length(n) != length(x)")}

xs <- sum(x)
ns <- sum(n)
 
event<-list()

for(i in 1:length(n))
 {event[[i]]<-0:min(xs,n[i])}

alltables <- expand.grid(event)

xs.alltables <- apply(alltables, MARGIN=1, sum)

cond.tables <- alltables[which(xs.alltables==xs),]

Prob.cond.tables<- apply(cond.tables, MARGIN=1, FUN=function(x)
  {prod(choose(n=n, k=x))/(choose(n=ns, k=xs))} )

teststat<-function(n,x)
{
nc <- n[1]
pc <- x[1]/nc

ntreat <- n[-1]
ptreat <- x[-1]/ntreat

tstat<-numeric(length=length(ptreat))
for(e in 1:length(ptreat))
 {
 tstat[e] <- (ptreat[e]-pc)/sqrt( pc*(1-pc)/nc + ptreat[e]*(1-ptreat[e])/ntreat[e] + 0.0001 )
 }
return(max(tstat))
}

tstatobs <- teststat(n=n, x=x)

#tstat.cond.tables <- apply(cond.tables, MARGIN=1, FUN=teststat(x=x))

tstat.cond.tables <- numeric(length=nrow(cond.tables))

for(r in 1:nrow(cond.tables))
 {
 tstat.cond.tables[r] <- teststat(n=n, x=as.numeric(cond.tables[r,]))
 }

if(alternative=="greater")
{ p.value <- sum(Prob.cond.tables[which(tstat.cond.tables>=tstatobs)]) }

if(alternative=="less")
{ p.value <- sum(Prob.cond.tables[which(tstat.cond.tables<=tstatobs)]) }

if(alternative=="two.sided")
{ p.value.less <- sum(Prob.cond.tables[which(tstat.cond.tables<=tstatobs)])
  p.value.greater <- sum(Prob.cond.tables[which(tstat.cond.tables>=tstatobs)])
  p.value<-min(2*p.value.less, 2*p.value.greater,1)
}

return(p.value)
}

