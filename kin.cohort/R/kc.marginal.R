`kc.marginal` <-
function(t, delta, genes, r, knots, f, pw=rep(1,length(t)), set=NULL, B=1, maxit=1000, tol=1e-5, subset, logrank=TRUE, trace=FALSE){
#
# function for kin-cohort cumulative risk estimation.
# Piecewise exponential models are used for survial function.
#
# INPUT ARGUMENTS
#      t = time to event
#  delta = 0-1 indicator of failure
#  genes = list of carrier status of proband for each gene
#          may be a vector (one gene), matrix, or data.frame
#          ideally use labeled factors otherwise (1=non-carriers, 2=carriers) is assumed
#          it is important that non-carriers are coded first
#  knots = vector of knots for piecewise constant model of hazards
#          These should be internal cutoints (lower and upper limits will be added):
#          The first knot should not be smaller than the time of the first event
#          and the last knot should not be higher than the time of the last event
#      r = relative type (1=parents/child, 2=siblings)
#      f = vector of allele frequencies of mutation for each gene
#     pw = prior weights for relatives (assume 1 if no weights are needed)
#
#    set = family indicator (not used)
#      B = if >1 (bootstrap) don't perform genotype checks
#
#  maxit = max number of iterations of EM loop (defaults 1000)
#    tol = relative change in hazard estimates to stop the EM loop (defaults 1e-5)
#   orig = if FALSE, indexed version of the EM loop code (slower)
#
# OUTPUT ARGUMENTS
#   cumrisk = K by ngenes matrix containing estimates of cumulative risk probabilities
#             (K: number of knots, ngenes:combinations of genotypes)
#    hazard = K by ngenes matrix containing estimates of hazards
#     knots = copy of knots
#      conv = indicator of whether the EM algorithm converges (1=converged)
#     niter = numer of iterations needed for convergence
#

cl<-match.call()
ngenes <- length(f)			             # number of genes
ngeno.rel <- ngenes*2			          # number of relative genotypes
genes<-data.frame(genes)

if(B>1 & is.null(set)) stop("set variable needed for bootstrap estimates")


valid<- complete.cases(t, delta, genes, r, pw)
if (sum(valid) != length(t) | !missing(subset) | any(r==0))  {
  if (missing(subset)) subset<-rep(TRUE,length(t))
  valid<-valid & subset
  if (sum(valid) != length(t[subset])) warning("Some cases excluded because of missing values", call.=FALSE)
  if (any(r[subset]==0)) warning("Probands excluded from this analysis", call.=FALSE)
  valid<-valid & r!=0  # exlude probands

  t<- t[valid]
  delta<- delta[valid]
  r<- r[valid]
  pw <- pw[valid] 
  genes <-  genes[valid,,drop=FALSE]
  if (!is.null(set)) set<- set[valid]
}

# relatives
    rel.codes<-unique(r)
    if (any(rel.codes)>3 | any(rel.codes)<0)
        stop("Code relatives 0:proband (ignored here) 1:parents 2:siblings 3:offspring (recoded to 1)")

    r[r==3]<-1                    # pool parents & offspring
    r[r==4 | r==5]<-3             # pool 2nd degree relatives (grandparents & parent-sibs)


if (ngenes==1)
  ngeno.pro <- length(unique(genes[[1]]))    # number of proband genotypes
else if (ngenes==2)
  ngeno.pro <-  length(unique(genes[[1]])) * length(unique(genes[[2]])) 
else stop("2 genes max at the moment")


# genotype checks
if (ngenes==2 & dim(genes)[2]!=2) stop("Please provide genotypes in a data.frame or matrix")
for (i in 1:ngenes){
	if (!is.factor(genes[[i]])){
		warning("Genotypes better specified as factor. Levels assigned, check adequacy", call.=FALSE)
		genes[[i]]<-factor(genes[[i]])
		if (length(unique(genes[[i]]))==2) levels(genes[[i]]) <- c("Noncarrier","Carrier")
		if (length(unique(genes[[i]]))==3) levels(genes[[i]]) <- c("Noncarrier","Heterozygous","Homozygous")
	}
	fr <- table(genes[[i]])
	if (fr[2]>fr[1]) warning("Non carrriers should be coded first!", call.=FALSE)
}

# vectors of genotypes combinations for 2 genes
if (ngenes==2){
   lev<-list(levels(genes[[1]]), levels(genes[[2]]) )
   g0<-factor(paste(lev[[1]][genes[[1]]], lev[[2]][genes[[2]]]), levels = t(outer(lev[[1]], lev[[2]], paste)))
   genes[[1]][as.numeric(genes[[1]])==3]<-lev[[1]][2]
   genes[[2]][as.numeric(genes[[2]])==3]<-lev[[2]][2]   
   g0.rel<-as.numeric(factor(paste(lev[[1]][genes[[1]]], lev[[2]][genes[[2]]]), levels = t(outer(lev[[1]], lev[[2]], paste))))
} else {
   g0<-genes[[1]]
   lev<-list(levels(g0),NULL)
   genes[[1]][as.numeric(genes[[1]])==3]<-lev[[1]][2]
   g0.rel<-as.numeric(genes[[1]])
}


nk <- length(knots)
nobs<- length(t)
nknots<- c(0,knots,Inf)

# labels

g0.lab<- levels(g0)
g0<- as.numeric(g0)

if (ngenes==2) {
   n=matrix(NA,2,2)
   for (i in 1:2){
      nn=levels(genes[[i]])
      if(length(nn)==3)
         n[i,]<-c(nn[1],paste(nn[2:3],collapse="+"))
      else
         n[i,]<-nn
   }
   g0.lab = as.vector(t(outer(n[1,],n[2,], paste)))
} else {
      nn=levels(genes[[1]])
      if(length(nn)==3)
         g0.lab <-c(nn[1],paste(nn[2:3],collapse="+"))
}
sur.nam<-apply(cbind(rep("<",nk), knots),1,paste,collapse="")
haz.nam<- cbind(c(knots[1],knots[-nk]),c("",rep("-",nk-1)),c("",knots[-1]))
haz.nam<- apply( cbind(c("<",rep("[",nk-1)), apply(haz.nam,1,paste,collapse="") ,c("",rep(")",nk-1))),1,paste,collapse="")


# descriptive: events & person years
dpy <- pyear(t,delta,knots)
d<-dpy$d
py<-dpy$py
events<- apply(d*pw,2,sum)
p.years<- apply(py*pw,2,sum)
p.years<-c(p.years,NA,sum(p.years), round(sum(events)/sum(p.years)*100000,1)) # don't reorder !!
events<-c(events,NA,sum(events), NA)
epy<- cbind(events, p.years)

for (i in 1:ngeno.rel) {
   filter<-(g0.rel==i)
	events<-apply(d[filter,]*pw[filter],2,sum) 
	p.years<- apply(py[filter,]*pw[filter],2,sum)
   p.years<-c(p.years,NA,sum(p.years), round(sum(events)/sum(p.years)*100000,1))
   events<-c(events,NA,sum(events), NA)
	epy0<- cbind(events, p.years)
	colnames(epy0)[1]<-g0.lab[i]
	epy<- cbind(epy, epy0)
}
rownames(epy)<-c(haz.nam, paste(knots[nk],"+"),"","Total","x10e5 py")

# table of carrier probabilities for relative conditional on proband

tp <- mendelian.combine(f  ,c(length(lev[[1]]), length(lev[[2]])), c(2,2)) # dominant model for output


########################################################
# Logrank test
#
  if(logrank){
    lr<-list()
    lrtmp<- survdiff(Surv(t,delta)~genes[[1]])
    lr[[1]]<- 1-pchisq(lrtmp$chisq,length(lrtmp$obs))
    if (ngenes==2){
       lrtmp<- survdiff(Surv(t,delta)~genes[[2]])
       lr[[2]]<- 1-pchisq(lrtmp$chisq,length(lrtmp$obs))

       lrtmp<- survdiff(Surv(t,delta)~genes[[1]]+genes[[2]])
       lr[[3]]<- 1-pchisq(lrtmp$chisq,length(lrtmp$obs))

       lrtmp<- survdiff(Surv(t,delta)~genes[[1]]+strata(genes[[2]]))
       lr[[4]]<- 1-pchisq(lrtmp$chisq,length(lrtmp$obs))

       lrtmp<- survdiff(Surv(t,delta)~genes[[2]]+strata(genes[[1]]))
       lr[[5]]<- 1-pchisq(lrtmp$chisq,length(lrtmp$obs))
       names(lr)<-c("gen1", "gen2","gen1+gen2","gen1(gen2)","gen2(gen1)")
    }
    lr<-unlist(lr)

  } else lr<-NULL
  
########################################################
# Bootstrap starts here
#
marginal <- function(t, delta, g0, r, pw, tp, knots, nk, nknots, ngenes, ngeno.pro, ngeno.rel ) {


nobs<- length(t)
tc<-   cut(t,nknots,labels=FALSE,include.lowest=TRUE,right=FALSE) #bin2
ltc<-  tc-1

# pyear called just once here
dpy <- pyear(t,delta,knots)
d<-dpy$d
py<-dpy$py


pool  <- pwexp(d,py,pw,knots)  # start with all geno.rel equal

h<- h0<-   matrix(rep(pool$h,ngeno.rel),nk+1,ngeno.rel)
s<- s0<-   matrix(rep(pool$s,ngeno.rel),nk,ngeno.rel)

fv<- haz<- surv<- matrix(0,nobs,ngeno.rel) # individual contributions
w <- array(0,c(nobs,ngeno.rel,ngeno.pro))

rerror<-  1
niter<-   1
conv<-    0


#
# E-M loop
#
while  ((rerror > tol) & (niter < maxit)){

#
# indexed version, slower but may make ease the accomodation of multiple genotypes
#
   for (i in 1:ngeno.rel){
     # expected survival: individual contributions to L based on theta: h/s and t
     surv[tc>1,i] <- s[ltc[tc>1],i]*exp(-h[tc[tc>1],i]*(t[tc>1]-nknots[tc[tc>1]]+0.5))

     # s(0)=1
     if (any(tc==1)) surv[tc==1,i] <- exp(-h[tc[tc==1],i]*(t[tc==1]-nknots[tc[tc==1]]+0.5))

     haz[,i]<- h[tc,i]

     # q(y,theta)  theta: h/s at knots
     fv[, i]<- delta * haz[,i] * surv[, i] + (1 - delta) * surv[, i]
  }

  #cat("E-step\n")
  # E-step calculate mixing weights
  #
  # Pr(gR,gP) [depends on relative type] * q(y,theta)
  #       parent/child    siblings

   for (j in 1:ngeno.pro){        # for each proband genotype combination g0 = g1 x g2
     for (k in 1:3){             # cond carrier prob vary for relative type
        tpsum <- 0
        for (i in 1:ngeno.rel) { # calculate denominator: sum for each possible relative genotype
                 tpsum <- tpsum + tp[j,i,k]*fv[(r==k),i]
        }
        for (i in 1:ngeno.rel) { # calculate weight for relative genotype
            w[(r==k),i,j] <- tp[j,i,k]*fv[(r==k),i]/tpsum
        }
     }
   }
   w[is.na(w)]<-0

   # M-step
   #
	for (i in 1:ngeno.rel){ # fits a model for each combination of possible relative genotypes with weights
		tw <-0
		for (j in 1:ngeno.pro){
			tw = tw + w[,i,j]*(g0==j)  # selects weight according to proband genotype g0 (only 1 g0==1)
		}
		hat<-pwexp(d,py,tw*pw,knots)
		s[,i]<-hat$s
		h[,i]<-hat$h
	}
   #test convergence
   rerror = max( t(abs(h-h0)) / apply(abs(h0),2,max,0.1) )

   h0 = h
   s0 = s
   niter = niter + 1
}
#
# post processing
#
if  (rerror < tol) conv= 1

#
# average harzard ratio
#  weights proportional to cumulative number of events

logHR<- rep(NA,ngeno.rel-1)
for (i in 1:(ngeno.rel-1)){

 cd<-cumsum( apply(d*pw,2,sum) )[-ncol(d)]
 wcd<-(cd/sum(cd))
 lhr<-log(-log(s[,i+1]))-log(-log(s[,1]))
 logHR[i]<- sum( lhr[is.finite(lhr)]*wcd[is.finite(lhr)] )
}

# s-> cum risk
s<-1-s
for (i in 2:ngeno.rel){
   s<-cbind(s,s[,i]/s[,1])
}
# ignore last (open) interval for h
h<-h[-(nk+1),]
for (i in 2:ngeno.rel){
  h<-cbind(h,h[,i]/h[,1])
}

rownames(s)<-sur.nam
colnames(s)<-c(g0.lab,paste(rep("Cum.RR",ngeno.rel-1),g0.lab[-1]))
rownames(h)<- haz.nam
colnames(h)<-c(g0.lab,paste(rep("HR",ngeno.rel-1),g0.lab[-1]))
names(logHR)<-g0.lab[-1]

list(cumrisk=s, hazard=h, conv=conv, niter=niter, logHR=logHR)

} # marginal (bootstrap)
###################################################################

b.orig <- marginal(t, delta, g0, r, pw, tp, knots, nk, nknots, ngenes, ngeno.pro, ngeno.rel )

o<-list(cumrisk=b.orig$cumrisk, hazard=b.orig$hazard, knots=knots, conv=b.orig$conv, niter=b.orig$niter, ngeno.rel=ngeno.rel, events=epy, logHR=b.orig$logHR, f=f, call=cl, logrank=lr )
class(o)<-c("kin.cohort","chatterjee")

#### bootstrap code
if (B>1){
	index<- 1:nobs
	id<- unique(set)
	nfamily<- length(id)
	
   if( ngenes == 1 ){
      s1 <- s2 <- sr<- h1 <- h2 <- hr<- matrix(0,B,nk)
      #   sb <- hb <- array(0,c(B,nk,3))
   } else {
      #	s1 <- s2 <- s3 <- s4 <- sr2<- sr3<- sr4<- h1 <- h2 <- h3 <- h4 <- hr2<- hr3<- hr4<-  matrix(0,B,nk)
      sb <- hb <-array(0,c(nk,7,B))
   }
 
	logHR <- matrix(0,B,ngeno.rel-1)
	
	x<- matrix(NA,nfamily,max(table(set)))
	for (i in 1:nfamily){
	   ii<-index[set==id[i]]
      x[i,1:length(ii)] <- ii
	}
	
	i<-1
	while (i<=B){

      if (trace){
        if (i==1) cat("Sample: ")
        if (i %% 10==0) cat(i," ")
        if (i==B) cat("\n")
      }
	   s<- sample(1:nfamily, replace=TRUE) 
	   index_s <- NULL
	   for (j in 1:nfamily){
	      ii<-x[s[j],]
	      index_s<-c(index_s,ii[!is.na(ii)])
	   }
	
	   t_s=t[index_s]
	   delta_s=delta[index_s]
	   g0_s=g0[index_s]
	   r_s =r[index_s]
	   set_s=set[index_s]
	   pw_s=pw[index_s]
	
	   b<-marginal(t_s, delta_s, g0_s, r_s, pw_s, tp, knots, nk, nknots, ngenes, ngeno.pro, ngeno.rel )

      if( ngenes == 1 ){
         s1[i,]<- b$cumrisk[,1]
         s2[i,]<- b$cumrisk[,2]
         sr[i,]<- b$cumrisk[,3]
         h1[i,]<- b$hazard[,1]
         h2[i,]<- b$hazard[,2]
         hr[i,]<- b$hazard[,3]
      } else {
         sb[,,i] <- b$cumrisk
         hb[,,i] <- b$hazard
      }
      logHR[i,]<- b$logHR
      
	   i<-i+1
	}

   if( ngenes == 1 ){
   	colnames(s1)<-rownames(b$cumrisk)
   	colnames(s2)<-rownames(b$cumrisk)
   	colnames(sr)<-rownames(b$cumrisk)
   	colnames(h1)<-rownames(b$hazard)
   	colnames(h2)<-rownames(b$hazard)
   	colnames(hr)<-rownames(b$hazard)
   	sb<- list(Noncarrier=s1,Carrier=s2,"C.R.R."=sr)
   	hb<- list(Noncarrier=h1,Carrier=h2,"H.R."=hr)
	}


	colnames(logHR)<-names(b$logHR)
	
   o<-list(cumrisk=sb,hazard=hb, logHR=logHR ,estimate=o)
	class(o)<-c("kin.cohort.boot","chatterjee")

} # bootstrap

####

o
}

