# contribution of block to log likelihood
blcontrib<-function(obj,p.patt)
{

  # only for data elements in partsList
  #    data elements have three components: $counts, $notnaidx, $s
  #    last element has only one component - covariate structure -  and has no name here


#    ## only non 0 observations
#    n<-obj$s[length(obj$s)]
#    idx<-(1:length(obj$counts))[obj$count>0]     # index which counts>0
#    idx2<-obj$s %in% idx
#    s2<-obj$s*idx2
#    # calculates probabilities for patterns by summation according to s
#    new.p.patt<-tapply(p.patt,s2,sum)
#    counts<-obj$counts[idx]
#    ll<-sum(counts*log(new.p.patt[-1]))


     ## all obs (slightly slower than the above)
     new.p.patt<-tapply(p.patt,obj$s,sum)
     ll<-sum(obj$counts*log(new.p.patt))

     #log likelihood for full (saturated) model
     p.cnts<-obj$counts/sum(obj$counts)
     fl<-sum(log(p.cnts[p.cnts>0])*obj$counts[obj$counts>0])

     list(ll=ll,fl=fl)
}
