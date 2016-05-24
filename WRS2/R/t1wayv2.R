t1wayv2 <-
function(x,tr=.2,grp=NA,MAT=FALSE,lev.col=1,var.col=2,nboot=100,SEED=TRUE,pr=TRUE,IV=NULL,loc.fun=median){
  #
  # Same a t1way, but computes explanatory power and related effect size
  #
  # For n1!=n2, this function calls t1way.effect.
  #
  #  A heteroscedastic one-way ANOVA for trimmed means
  #  using a generalization of Welch's method.
  #
  #  The data are assumed to be stored in $x$ in a matrix or in list mode.
  #
  # MAT=F, if x is a matrix, columns correspond to groups.
  # if MAT=T, assumes argument
  # lev.col
  # indicates which column of x denotes the groups. And
  #  var.col indicates the column where the data are stored.
  #
  #  IV, if specified, taken to be the independent variable
  #      That is, the group id values
  #      and x is assumed to be a vector containing all of the data
  #
  # if x has list mode:
  #  length(x) is assumed to correspond to the total number of groups.
  #  By default, the null hypothesis is that all groups have a common mean.
  #  To compare a subset of the groups, use grp to indicate which
  #  groups are to be compared. For example, if you type the
  #  command grp<-c(1,3,4), and then execute this function, groups
  #  1, 3, and 4 will be compared with the remaining groups ignored.
  #
  #  Missing values are automatically removed.
  #
  if(SEED)set.seed(2)
  if(MAT){
    if(!is.matrix(x))stop("With MAT=T, data must be stored in a matrix")
    if(length(lev.col)!=1)stop("Argument lev.col should have 1 value")
    temp=selby(x,lev.col,var.col)
    x=temp$x
    grp2=rank(temp$grpn)
    x=x[grp2]
  }
  if(!is.null(IV[1])){
    if(pr)print("Assuming x is a vector containing all of the data, the dependent variable")
    xi=elimna(cbind(x,IV))
    x=fac2list(xi[,1],xi[,2])
  }
  if(is.matrix(x))x<-listm(x)
  if(is.na(sum(grp[1])))grp<-c(1:length(x))
  if(!is.list(x))stop("Data are not stored in a matrix or in list mode.")
  J<-length(grp)
  h<-vector("numeric",J)
  w<-vector("numeric",J)
  xbar<-vector("numeric",J)
  pts=NULL
  nval=0
  for(j in 1:J)x[[j]]=elimna(x[[j]])
  for(j in 1:J){
    val<-x[[j]]
    val<-elimna(val)
    nval[j]=length(val)
    pts=c(pts,val)
    x[[j]]<-val # missing values have been removed
    h[j]<-length(x[[grp[j]]])-2*floor(tr*length(x[[grp[j]]]))
    # h is the number of observations in the jth group after trimming.
    w[j]<-h[j]*(h[j]-1)/((length(x[[grp[j]]])-1)*winvar(x[[grp[j]]],tr))
    xbar[j]<-mean(x[[grp[j]]],tr)
  }
  u<-sum(w)
  xtil<-sum(w*xbar)/u
  A<-sum(w*(xbar-xtil)^2)/(J-1)
  B<-2*(J-2)*sum((1-w/u)^2/(h-1))/(J^2-1)
  TEST<-A/(B+1)
  nu1<-J-1
  nu2<-1./(3*sum((1-w/u)^2/(h-1))/(J^2-1))
  sig<-1-pf(TEST,nu1,nu2)
  nv=lapply(x,length)
  #
  # Determine explanatory effect size
  #
  chkn=var(nval)
  if(chkn==0){
    top=var(xbar)
    bot=winvarN(pts,tr=tr)
    e.pow=top/bot
  }
  if(chkn!=0){
    vals=0
    N=min(nval)
    xdat=list()
    for(i in 1:nboot){
      for(j in 1:J){
        xdat[[j]]=sample(x[[j]],N)
        vals[i]=t1way.effect(xdat,tr=tr)$Var.Explained
      }}
    e.pow=loc.fun(vals,na.rm=TRUE)
  }
  list(TEST=TEST,nu1=nu1,nu2=nu2,n=nv,p.value=sig,Var.Explained=e.pow,
       Effect.Size=sqrt(e.pow))
}
