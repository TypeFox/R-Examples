`kc.moments` <-
function(t, delta, genes, r, knots, f, pw=rep(1,length(t)), set=NULL, B=1, logrank=TRUE, subset, trace=FALSE){
#
# function for kin-cohort cumulative risk estimation.
# method by Wacholder et al. AJE 1998; 148: 632-9
#
# INPUT ARGUMENT (all vectors should be column vector)
#      t = time to event
#  delta = 0-1 indicator of failure
#     g0 = carrier status of proband (1=non-carriers, 2=carriers)
#  knots = vector of knots for point estimates
#      f = allele frequency of mutation
#     pw = prior weights for relatives (assume 1 if no weights are needed)
#
#      r = relative type (only used to filter relatives from data)
#    set = family indicator (not used)
#      B = bootstrap samples (not used)

cl<-match.call()
if(B>1 & is.null(set)) stop("set variable needed for bootstrap estimates")

# genotype checks
if (length(f)!=1 | dim(data.frame(genes))[2]!=1) stop("Only one gene supported for method of moments")

# relatives
    rel.codes<-unique(r)
    if (any(rel.codes)>3 | any(rel.codes)<0)
        stop("Code relatives 0:proband (ignored here) 1:parents 2:siblings 3:offspring (recoded to 1)")

    r[r==3]<-1   # pool parents & descendents
                 # could extend to 2nd-3rd degree relatives

    if (any(r==0)) warning("Probands excluded from this analysis", call.=FALSE)
    if (any(r>3))  warning("Second degree relatives excluded from this analysis", call.=FALSE)

valid<- complete.cases(t, delta, genes, r, pw)
  if (sum(valid) != length(t) | !missing(subset) | any(r==0) | any(r>3))  {
  if (sum(valid) != length(t)) warning("Some cases excluded because of missing values", call.=FALSE)
  if(!missing(subset)) valid<-valid & subset
  valid<-valid & r!=0 & r<=3 # exlude probands & 2nd D rel
  t<- t[valid]
  delta<- delta[valid]
  r<- r[valid]
  pw <- pw[valid]
  genes <-  genes[valid]
  if (!is.null(set)) set<- set[valid]
}

# genotype checks
	if (!is.factor(genes)){
		warning("Genotypes better specified as factor. Levels assigned, check adequacy", call.=FALSE)
		genes<-factor(genes)
		if (length(unique(genes))==3) {
 		    warning("3r genotype collapsed with 2nd", call.=FALSE)
 		    genes[genes==levels(genes)[3]]<-levels(genes)[2]
 		    }
		if (length(unique(genes))==2) levels(genes) <- c("Noncarrier","Carrier")
		else stop("something wrong with genotypes")
	} else if (length(unique(genes))==3) {
 		    warning("3r genotype collapsed with 2nd", call.=FALSE)
 		    genes[genes==levels(genes)[3]]<-levels(genes)[2]
 		    levels(genes)[2]<-paste(levels(genes)[2:3],collapse="+")
 		    genes<-factor(genes)
	}
	fr <- table(genes)
	if (fr[2]>fr[1]){
      warning("Non carriers should be coded first!, reordered, check labels", call.=FALSE)
      genes<-relevel(genes,2)
   }

nk <- length(knots)
nobs<- length(t)
nknots<- c(0,knots,Inf)

sur.nam<-apply(cbind(rep("<",nk), knots),1,paste,collapse="")
haz.nam<- cbind(c(knots[1],knots[-nk]),c("",rep("-",nk-1)),c("",knots[-1]))
haz.nam<- apply( cbind(c("<",rep("[",nk-1)), apply(haz.nam,1,paste,collapse="") ,c("",rep(")",nk-1))),1,paste,collapse="")

# descriptive: events & person years
dpy <- pyear(t,delta,knots)
d<-dpy$d
py<-dpy$py
events<- apply(d*pw,2,sum) 
events<-c(events,NA,sum(events))
p.years<- apply(py*pw,2,sum)
p.years<-c(p.years,NA,sum(p.years))
epy<- cbind(events, p.years)

for (i in 1:2) {
   filter<-(as.numeric(genes)==i)
	events<-apply(d[filter,]*pw[filter],2,sum) 
   events<-c(events,NA,sum(events))
	p.years<- apply(py[filter,]*pw[filter],2,sum)
   p.years<-c(p.years,NA,sum(p.years))
	epy0<- cbind(events, p.years)
	colnames(epy0)[1]<-levels(genes)[i]
	epy<- cbind(epy, epy0)
}
rownames(epy)<-c(haz.nam, paste(knots[nk],"+"),"","Total")


#####################################################

moments<-function(t,delta,genes,pw,knots,f, logrank){
  km<- survfit(Surv(t,delta)~genes,weights=pw)
  if(logrank){
    lr<- survdiff(Surv(t,delta)~genes) # no weights allowed here
    lr<- 1-pchisq(lr$chisq,length(lr$obs))
  }  else
    lr<-NA 
  #
  # Kaplan Meier estimates
  #
  r1<-1-summary(km[1], times=knots)$surv # non carriers
  r2<-1-summary(km[2], times=knots)$surv # carriers
  if (length(r1)!=length(knots)) r1<-c(r1,rep(r1[length(r1)],length(knots)-length(r1)))
  if (length(r2)!=length(knots)) r2<-c(r2,rep(r2[length(r2)],length(knots)-length(r2)))

  R<- cbind(r1,r2)
  rownames(R)<- knots
  colnames(R)<- levels(genes)
 #
 # kin-cohort cumulative risk 
 #
  r2s.plus<- c(-1,2)
  r2s.minus<- function(p)c((1+p)/(1-p),-2*p/(1-p))

  s<-cbind(R%*%r2s.minus(f),R%*%r2s.plus )

  s<-cbind(s,s[,2]/s[,1])
  colnames(s)<- c(levels(genes),"CumRisk Ratio")

  list(cumrisk=s, km=km, logrank=lr)
}

b<-moments(t,delta,genes,pw,knots,f, logrank=TRUE)
o<-list(cumrisk=b$cumrisk, knots=knots, km=b$km, logrank=b$logrank, events=epy, ngeno.rel=2, f=f, call=cl)
class(o)<-c("kin.cohort","wacholder")

#####################################################
# bootstrap code

if (B>1) {
	m<- length(t)
	index<- 1:m
	id<- unique(set)
	nfamily<- length(id)
	
	s1 <- s2 <- sr<- matrix(0,B,length(knots))
	
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
	   genes_s=genes[index_s]
	   set_s=set[index_s]
	   pw_s=pw[index_s]
	
	   b<-moments(t_s,delta_s,genes_s,pw_s,knots,f, logrank=FALSE)
	
	   s1[i,]<- b$cumrisk[,1]
	   s2[i,]<- b$cumrisk[,2]
	   sr[i,]<- b$cumrisk[,3]
	
	   i<-i+1
	}
	
	colnames(s1)<-rownames(b$cumrisk)
	colnames(s2)<-rownames(b$cumrisk)
	colnames(sr)<-rownames(b$cumrisk)
	
	o<-list(cumrisk=list(Noncarrier=s1,Carrier=s2,"C.R.R."=sr),estimate=o)
	class(o)<-c("kin.cohort.boot","wacholder")

} # bootstrap

o
}

