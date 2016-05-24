"bgtWidth" <-
function(n, s, p, conf.level=0.95, alternative="two.sided", method="CP")
{

 if( any(n<=3) )
  {stop("the number of groups n allowed in calculations must be integers greater than 1")}
 
 if( any(s<1) ){stop("group size s must be specified as integers > 0")}

 if( length(conf.level)!=1 || conf.level<0 || conf.level>1)
  {stop("conf.level must be a positive number between 0 and 1")}

 if( length(p)!=1 || p>1 || p<0)
  {stop("true proportion p must be specified as a single number between 0 and 1")}

  method<-match.arg(method, choices=c("CP","Blaker","AC","Score","Wald","SOC"))

  alternative<-match.arg(alternative, choices=c("two.sided","less","greater"))

# calculations:

 matnsp <- cbind(n,s,p)
 matnsp <- cbind("ns"=matnsp[,1]*matnsp[,2], matnsp)
 power <- numeric(length=nrow(matnsp))
 bias <- numeric(length=nrow(matnsp))

 expCIwidth<-numeric(length=nrow(matnsp))

 for (i in 1:length(expCIwidth))
  {
   expCIwidth[i]<-bgtWidthI(n=matnsp[[i,2]], s=matnsp[[i,3]], p=matnsp[[i,4]], conf.level=conf.level, alternative=alternative, method=method)$expCIWidth
  }

 return(as.matrix(cbind(matnsp,expCIwidth)))
}

"bgtWidthI" <-
function(n, s, p, conf.level=0.95, alternative="two.sided", method="CP")
{

# indicator function for the CI length at a special event
# in one sided case: length is defined as absolute difference between estimator and confidence bound

 L.Ind<-function(y, n, s, p, conf.level, alternative, method)

 {

  if(method=="Wald"){int=bgtWald(y=y, n=n, s=s, conf.level=conf.level, alternative=alternative)}
  if(method=="Wilson"){int=bgtWilson(y=y, n=n, s=s, conf.level=conf.level, alternative=alternative)}
  if(method=="AC"){int=bgtAC(y=y, n=n, s=s, conf.level=conf.level, alternative=alternative)}
  if(method=="SOC"){int=bgtSOC(y=y, n=n, s=s, conf.level=conf.level, alternative=alternative)}
  if(method=="CP"){int=bgtCP(y=y, n=n, s=s, conf.level=conf.level, alternative=alternative)}
  if(method=="Blaker"){int=bgtBlaker(y=y, n=n, s=s, conf.level=conf.level)}

  if(alternative=="less")
    {CIlength <- int[[2]]-p}

  if(alternative=="greater")
    {CIlength <- p-int[[1]]}

  if(alternative=="two.sided")
    {CIlength <- int[[2]]-int[[1]]}
  CIlength
 }

# Probability of a single event, the binomial group testing density:

 bgt.prob<-function(n,y,s,p.tr)
  {
  theta<-1-(1-p.tr)^s
  dbinom(x=y,size=n, prob=theta)
  }


#  calculate this for all possible events: 

yvec<-0:n

   Lvec<-numeric(length=length(yvec))   
   probvec<-numeric(length=length(yvec))

   for(i in 1:length(yvec))
    {Lvec[i]<-L.Ind(y=yvec[i], n=n, s=s, p=p, conf.level=conf.level, alternative=alternative, method=method)
     probvec[i]<-bgt.prob(y=yvec[i], n=n, s=s, p.tr=p)
    }
  expCILength=sum(Lvec * probvec)

# E(X)= sum(Xi * prob(Xi))

out<-list(expCIWidth=expCILength, alternative=alternative, n=n,s=s, p=p)

class(out)<-"binWidth"
return(out)
}

