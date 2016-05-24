###These are tests of the basic function of xBalance under some
###restricted sitations and using functions that hew closely to the
###expressions in Hansen and Bowers 2008. For example, with only
###binary treatment. The idea here to show that the math in that
###article is equivalent to the output from xBalance.

require("RItools")

data(nuclearplants)

s2.fn<-function(x){sum((x-mean(x))^2)/(length(x)-1)} ##same as sd(x)
h.fn<-function(n,m){(m*(n-m))/n}

var1<-function(x,m){ ##var(d)
  h<-h.fn(n=length(x),m=m)
  (1/h)*s2.fn(x)
}

var2<-function(x,m){ ##var(Z'x) (i.e. var of the sum statistic)
  h<-h.fn(n=length(x),m=m)
  (h)*s2.fn(x)
}

#####First just looking at the unstratified calculations
xb1a<-xBalance(pr~ date+ t1 + t2 + cap + ne + ct + bw + cum.n,
               strata=list(nostrat=NULL),
               data=nuclearplants,
               report=c("adj.means","adj.mean.diffs","adj.mean.diffs.null.sd","std.diffs","z.scores","p.values"))

##print(xb1a,digits=4)

testxb1a<-t(sapply(nuclearplants[,dimnames(xb1a$results)$vars],function(thevar){
  myssn<-with(nuclearplants,sum((pr-mean(pr))*thevar))
  myadjdiff<-myssn/h.fn(m=sum(nuclearplants$pr),n=length(nuclearplants$pr))
  mynullvar1<-var1(x=thevar,m=sum(nuclearplants$pr))
  myz<-myadjdiff/sqrt(mynullvar1)
  return(c(adj.diff=myadjdiff,adj.diff.null.sd=sqrt(mynullvar1),z=myz))
}))

##print(testxb1a,digits=4)

all.equal(xb1a$results[,c("adj.diff","adj.diff.null.sd","z"),"nostrat"],testxb1a,check.attributes = FALSE)

###Now with strata.
xb2<-xBalance(pr~ date+ t1 + t2 + cap + ne + ct + bw + cum.n,
              strata=list(pt=~pt),
              data=nuclearplants,
              report=c("all"))


test2.fn<-function(zz,mm,ss){
  ##Some notes on xBalanceEngine
  ##dv = h/(n-1)
  ##unsplit(tapply(zz,ss,function(z){h.fn(m=sum(z),n=length(z))/(length(z)-1)}),ss) ##h/(n-1)
  ##tmat*tmat = squared mean deviations
  ##sapply(split(data.frame(mm),ss),function(x){sapply(x,function(var){(var-mean(var))^2})})
  ##so dv*tmat*tmat=(h/(n-1))*(x_i-\bar(x))^2=h*s^2

  myssn<-apply(mm,2,function(x){sum((zz-unsplit(tapply(zz,ss,mean),ss))*x)})

  hs<-tapply(zz,ss,function(z){h.fn(m=sum(z),n=length(z))})
  mywtsum<-sum(hs)

  myadjdiff<-myssn/mywtsum

  s2s<-sapply(data.frame(mm),function(x){sapply(split(x,ss),function(var){var(var)})})

  myssvar<-apply(s2s,2,function(s2){sum(hs*s2)})
  mynullsd2<-apply(s2s,2,function(s2){sqrt((1/(sum(hs)^2))*sum(hs*s2))}) ##from StatSci eq 6
  mynullsd1<-sqrt(myssvar*(1/mywtsum)^2) ##with (1/h)

  stopifnot(all.equal(mynullsd1,mynullsd2,check.attributes=FALSE))

  ##If numbers of treated and controls are the same for all blocks:
  if( length(unique(ss))>1 & all(diff(hs)==0) ){
    mynullsd3<-apply(s2s,2,function(s2){sqrt((1/(length(hs)))^2*sum((1/hs)*s2))}) ##(1/B^2)*\sum_{b=1}^B (1/h) s2
    stopifnot(all.equal(mynullsd3,mynullsd2))

    ##For fun, the version ignoring the stratification.
    m<-sum(zz)
    n<-length(zz)
    h<-(m*(n-m))/n
    mynullsd.nostrat<-sqrt((1/h)*apply(mm,2,var))

  }
  ##For pairs
  if((length(unique(ss))>1 & all(table(ss)==2)) & all(diff(hs)==0)){
    mys2s.pairs<-sapply(data.frame(mm),function(x){sapply(split(x,ss),function(var){sum((var-mean(var))^2)})})
    mynullsd4<-apply(s2s,2,function(s2){sqrt((2/(length(s2)^2))*sum(s2))}) ##(1/B^2)*\sum_{b=1}^B (1/h) s2 = (2/B^2)\sum_{b=1}^B \sum_{i=1}^2 s2
    stopifnot(all.equal(mynullsd4,mynullsd2))
  }
  myz2<-myadjdiff/mynullsd2
  myz3<-myadjdiff/mynullsd1
  myz1<-myssn/sqrt(myssvar)

  stopifnot(all.equal(myz1,myz2,myz3,check.attributes=FALSE))

  return(cbind(adj.diff=myadjdiff,adj.diff.null.sd=mynullsd2,z=myz2))
}

mymm<-model.matrix(pr~ date+ t1 + t2 + cap + ne + ct + bw + cum.n-1,data=nuclearplants)

test2.fn(zz=nuclearplants$pr,mm=mymm,ss=nuclearplants$pt)

all.equal(test2.fn(zz=nuclearplants$pr,mm=mymm,ss=nuclearplants$pt),
          xb2$results[,c("adj.diff","adj.diff.null.sd","z"),"pt"],check.attributes=FALSE)
