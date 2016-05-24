t1way.effect <-
function(x,tr=.2,grp=NA,MAT=FALSE,lev.col=1,var.col=2){
#
# Same as t1way, but computes explanatory power and related effect size
# Only use this function with = n's
# Use t1wayv2 in general, which calls this function when sample sizes differ.
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
if(MAT){
if(!is.matrix(x))stop("With MAT=T, data must be stored in a matrix")
if(length(lev.col)!=1)stop("Argument lev.col should have 1 value")
temp=selby(x,lev.col,var.col)
x=temp$x
grp2=rank(temp$grpn)
x=x[grp2]
}
if(is.matrix(x))x<-listm(x)
if(is.na(sum(grp[1])))grp<-c(1:length(x))
if(!is.list(x))stop("Data are not stored in a matrix or in list mode.")
J<-length(grp)
h<-vector("numeric",J)
w<-vector("numeric",J)
xbar<-vector("numeric",J)
pts=NULL
for(j in 1:J){
xx<-!is.na(x[[j]])
val<-x[[j]]
pts=c(pts,val)
x[[j]]<-val[xx]  # missing values have been removed
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
#
# Determine explanatory effect size
#
top=var(xbar)
bot=winvarN(pts,tr=tr)
e.pow=top/bot
list(TEST=TEST,nu1=nu1,nu2=nu2,siglevel=sig,Var.Explained=e.pow,
Effect.Size=sqrt(e.pow))
}
