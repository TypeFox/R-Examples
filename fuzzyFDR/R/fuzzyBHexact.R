"fuzzyBHexact" <-
function(pvals,pprev,alpha=0.05,tol=10e-6,q.myuni=T,dp=20){

pvals <- round(pvals,d=dp)
pprev <- round(pprev,d=dp)
    
## get the intervals D_j
if(q.myuni) p.all <- sort(myUnique(c(pvals,pprev),tol=tol))
else p.all <- sort(unique(c(pvals,pprev)))

#print(p.all)

n.int <- length(p.all)-1
n.tests <- length(pvals)
p.minus <- p.all[-(n.int+1)]
p.plus <- p.all[-1]
print(paste("total no. intervals = ",n.int))

## get z.min and z.max for all i (for intervals D_j)
z.min <- rep(-1,n.tests)
z.max <- rep(-1,n.tests)
numalloc.poss <- 1
for(i in 1:n.tests){
  if(q.myuni){
    z.min[i] <- myMatch(pprev[i],p.minus,tol=tol)
    if(pvals[i]==0) z.max[i] <- 1
    else z.max[i] <- myMatch(pvals[i],p.plus,tol=tol)
  }
  else{
    z.min[i] <- match(pprev[i],p.minus)
    if(pvals[i]==0) z.max[i] <- 1
    else z.max[i] <- match(pvals[i],p.plus)
  }
  numalloc.poss <- numalloc.poss * (z.max[i]-z.min[i]+1)
}
print(paste("total no. possible alloc. = ",numalloc.poss))
#print(pvals)
#print(pprev)
#print(z.min)
#print(z.max)

## get global r.minus,r.plus (over all alloc)
r.minus <- rep(-1,n.int)
r.plus <- rep(-1,n.int)
for(j in 1:n.int){
  r.minus[j] <- sum(z.max<j)+1
  r.plus[j] <- sum(z.min<=j)
  #### this bit was wrong !!!
  ####dum <- z.max>=j & z.min<=j
  ####r.minus[j] <- min((1:n.tests)[dum])
  ####r.plus[j] <- max((1:n.tests)[dum])
}


## get sf,sc
a.minus <- r.minus*alpha/n.tests
a.plus <- r.plus*alpha/n.tests
dummy <- as.numeric(p.minus <= a.plus)
global.sf <- sum(sign(rev(cumsum(rev(dummy)))))
dummy <- as.numeric(p.plus <= a.minus)
global.sc <- sum(sign(rev(cumsum(rev(dummy)))))
print(paste("global sf = ",global.sf))
print(paste("global sc = ",global.sc))

## print for intervals D_j
leng <- r.plus-r.minus+1
#ytable <- cbind(p.minus,p.plus,r.minus,r.plus,leng,a.minus,a.plus)
ytable <- cbind(round(p.minus,d=4),round(p.plus,d=4),r.minus,r.plus,
                leng,round(a.minus,d=4),round(a.plus,d=4))
ytable <- data.frame(ytable)
names(ytable) <- c("p.minus","p.plus","r.minus","r.plus","leng","a.minus","a.plus")
print(ytable)

## redefine intervals rej,fuzzy,acc
new.n.int <- global.sf-global.sc + as.numeric(global.sc>0) + as.numeric(global.sf<n.int)
print(paste("reduced no. intervals = ",new.n.int))
new.p.minus <- rep(-1,new.n.int)
new.p.plus <- rep(-1,new.n.int)
if(global.sc>0){
 new.p.minus[1] <- p.minus[1]
 new.p.plus[1] <- p.plus[global.sc]
}
if(global.sf>global.sc){
 for(j in 1:(global.sf-global.sc)){
  new.p.minus[j+as.numeric(global.sc>0)] <- p.minus[global.sc+j]
  new.p.plus[j+as.numeric(global.sc>0)] <- p.plus[global.sc+j]
}}
if(global.sf<n.int){
 new.p.minus[new.n.int] <- p.minus[global.sf+1]
 new.p.plus[new.n.int] <- p.plus[n.int]
}

## get z.min and z.max for all i (for new intervals)
new.z.min <- rep(-1,n.tests)
new.z.max <- rep(-1,n.tests)
new.probsample <- array(-1,dim=c(n.tests,new.n.int))
new.numalloc.poss <- 1
for(i in 1:n.tests){
  if(z.min[i]<=global.sc) new.z.min[i] <- 1
  else if(z.min[i]<=global.sf & global.sc>0) new.z.min[i] <- z.min[i]-global.sc+1
  else if(z.min[i]<=global.sf & global.sc==0) new.z.min[i] <- z.min[i]
  else new.z.min[i] <- new.n.int
  if(z.max[i]<=global.sc) new.z.max[i] <- 1
  else if(z.max[i]<=global.sf & global.sc>0) new.z.max[i] <- z.max[i]-global.sc+1
  else if(z.max[i]<=global.sf & global.sc==0) new.z.max[i] <- z.max[i]
  else new.z.max[i] <- new.n.int
  for(j in new.z.min[i]:new.z.max[i]){
    new.probsample[i,j] <- max(0, (min(new.p.plus[j],pvals[i])-max(new.p.minus[j],pprev[i])) )/(pvals[i]-pprev[i])
  }
  #print(paste("new.probsample = ",new.probsample[i,]))
  new.numalloc.poss <- new.numalloc.poss * (new.z.max[i]-new.z.min[i]+1)
}
print(paste("reduced no. alloc. = ",new.numalloc.poss))

## get global r.minus,r.plus (over all alloc) for new intervals
new.r.minus <- rep(-1,new.n.int)
new.r.plus <- rep(-1,new.n.int)
for(j in 1:new.n.int){
  new.r.minus[j] <- sum(new.z.max<j)+1
  new.r.plus[j] <- sum(new.z.min<=j)
  #### this bit was wrong !!!
  ####dum <- new.z.max>=j & new.z.min<=j
  ####new.r.minus[j] <- min((1:n.tests)[dum])
  ####new.r.plus[j] <- max((1:n.tests)[dum])
}

## print for new intervals
#ztable <- cbind(new.p.minus,new.p.plus,new.r.minus,new.r.plus)
ztable <- cbind(round(new.p.minus,d=4),round(new.p.plus,d=4),new.r.minus,new.r.plus)
ztable <- data.frame(ztable)
names(ztable) <- c("new.p.minus","new.p.plus","new.r.minus","new.r.plus")
print(ztable)

## get the allocations matrices (for prob.rej.ties)
if(global.sc>0 & global.sf<n.int) len.max <- max(new.r.plus[-c(1,new.n.int)]-new.r.minus[-c(1,new.n.int)]+1)
else if(global.sc>0 & global.sf==n.int) len.max <- max(new.r.plus[-1]-new.r.minus[-1]+1)
else if(global.sc==0 & global.sf<n.int) len.max <- max(new.r.plus[-new.n.int]-new.r.minus[-new.n.int]+1)
else len.max <- max(new.r.plus-new.r.minus+1)
#print(paste("len.max = ",len.max))
if(len.max>2) getAlloc(len.max-1)


############################ enumerate allocations

lowerprod <- rep(1,n.tests)
lowerprod[1] <- 1
for(i in 2:n.tests){
  delt <- new.z.max[i-1]-new.z.min[i-1]+1
  lowerprod[i] <- lowerprod[i-1] * delt
}
#print("new total alloc = ")
#print(lowerprod[n.tests]*(new.z.max[n.tests]-new.z.min[n.tests]+1))

z <- array(-1,dim=c(new.numalloc.poss,n.tests))
for(i in 1:n.tests){
  delt <- new.z.max[i]-new.z.min[i]+1
  upperprod <- new.numalloc.poss/(lowerprod[i]*delt)
  ##zprime <- rep(-1,upperprod)
  zprime <- rep(-1,upperprod*delt)
  for(j in 1:delt){
    zprime[ (j-1)*upperprod + (1:upperprod) ] <- new.z.min[i]+j-1
  }
  z[,i] <- rep(zprime,lowerprod[i])
}
#print(z)

############################## loop over allocations

## initialise/declare vectors
len <- rep(-1,new.n.int)
r.minus <- rep(-1,new.n.int)
r.plus <- rep(-1,new.n.int)
tau <- rep(0,n.tests)
totalpalloc <- 0

print("starting loop over allocations")

for(d in 1:new.numalloc.poss){

  #print(d)
  #print(z[d,])
  #for(i in 1:n.tests){
  # print(new.probsample[i,z[d,i]])
  #}

  ## prob alloc
  logpalloc <- 0
  for(i in 1:n.tests){
    logpalloc <- logpalloc + log(new.probsample[i,z[d,i]])
  }
  palloc <- exp(logpalloc)

  #print(palloc)

  ## normalising constant
  totalpalloc <- totalpalloc + palloc

  ## ranks, length ties for each interval
  len[1] <- sum(z[d,]==1)
  r.minus[1] <- 1
  if(new.n.int>1){
    for(j in 2:new.n.int){
      len[j] <- sum(z[d,]==j)
      r.minus[j] <- r.minus[j-1] + len[j-1]
    }
  }
  r.plus <- c(r.minus[-1]-1,n.tests)

  #print(len)
  #print(r.minus)
  #print(r.plus)
  #print(new.p.minus)
  #print(new.p.plus)

  ##### CHANGED CALC. OF SF HERE
  ## get sf, sc (note def. of sc different from global one)
  dummy <- as.numeric(new.p.minus < r.plus*alpha/n.tests)
  if(global.sf<n.int) 
	sf <- sum(sign(rev(cumsum(rev(dummy[-new.n.int])))))
  else sf <- sum(sign(rev(cumsum(rev(dummy)))))
  dummy <- as.numeric(new.p.plus <= r.plus*alpha/n.tests)
  sc <- sum(sign(rev(cumsum(rev(dummy)))))
  #print(paste("sf = ",sf))
  #print(paste("sc = ",sc))

  ## get pi and eta for each interval

  probrej <- rep(-1,new.n.int)
  tzero <- rep(-1,new.n.int+1)
  tzero[new.n.int+1] <- 1

  if(sf<new.n.int){
   for(j in new.n.int:(sf+1)){
      probrej[j] <- 0
      tzero[j] <- 1
   }
  }
  if(sf>sc){
   for(j in sf:(sc+1)){
     if(len[j]==0){
       tzero[j] <- tzero[j+1]
     }
     else{
       dummy <- probRejTies(new.p.minus[j],new.p.plus[j],
                              r.minus=r.minus[j],len=len[j],
                              ntests=n.tests,alpha=alpha)
       probrej[j] <- 1 - tzero[j+1]*(1 - dummy$piprob)
       tzero[j] <- dummy$tzero * tzero[j+1]
     }
   }
  }
  if(sc>0){
   for(j in sc:1){
       probrej[j] <- 1
       tzero[j] <- 0
   }
  }

  ## get tau for each i
  for(i in 1:n.tests){
    tau[i] <- tau[i] + palloc*probrej[z[d,i]]  
  }
}

#print("total palloc = ")
#print(totalpalloc)

## normalise tau with totalprob 
## (since not summing over all possible alloc)
#tau <- tau/totalpalloc

print("")
print("Exact Method")
print(paste("alpha = ",alpha))
#xtable <- cbind(pvals,pprev,z.min,z.max,new.z.min,new.z.max,tau)
#print(xtable)

ind <- order(pvals,pprev)
xtable <- cbind(round(pvals,d=4),round(pprev,d=4),z.min,z.max,new.z.min,new.z.max,round(tau,d=4))
xtable <- data.frame(xtable)
names(xtable) <- c("pvals","pprev","z.min","z.max","new.z.min","new.z.max","tau")
print(xtable[ind,])

return(xtable)

}

