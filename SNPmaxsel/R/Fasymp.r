###    Asymptotical distribution of the maximally selected chi-square statistic
###
### Copyright 2006-10 Anne-Laure Boulesteix 
###
### 
###
###
### This file is part of the `exactmaxstat' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
### 
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA


Fasymp<-function(t,a.vec,type=NULL,groups=NULL)
{
K<-length(a.vec)

if (t<0)
 stop("t must be positive")

if (round(sum(a.vec),8) != 1)
 stop("a.vec must sum to 1")

if (any(a.vec<0))
 stop("a.vec must have positive entries")

if (is.null(type)&is.null(groups))
 stop("type and groups can not be both null.")

if (!is.null(type)&!is.null(groups))
 stop("type and groups can not be both non-null.")

if (is.null(groups))
 {
 if (!is.element(type,c("all.pairs","all.partitions","inter.ord","inter.cat","ord","inter.ord.main")))
  {
  stop("type must be all.pairs, all.partitions, ord, inter.ord, inter.cat or main")
  } 
 groups<-groups(K=K,a.vec=a.vec,type=type)
 }     

if (!is.null(type))
 {
if (length(a.vec)!=9 && is.element(type,c("inter.ord","inter.cat","inter.ord.main")))
 stop("for type inter.ord and inter.cat, a.vec must be of length 9")
 }
   
A<-t(sapply(groups,FUN=makeA,a.vec=a.vec))
diagonal<-as.numeric(sapply(groups,FUN=makeB,a.vec=a.vec))
A<-diag(diagonal,nrow=length(diagonal),ncol=length(diagonal))%*%A

a.vec[a.vec==0]<-Inf
SigmaZ<-A%*%diag(1/a.vec)%*%t(A)
m<-nrow(A)
     
return(pmvnorm(lower=rep(-sqrt(t),m),upper=rep(sqrt(t),m),sigma=SigmaZ,mean=rep(0,m)) )

}

###############


makeA<-function(groups,a.vec)
{
group1<-groups$group1
group2<-groups$group2
p<-length(a.vec)

a<-numeric(p)
a[group1]<-a.vec[group1]/sum(a.vec[group1])
a[group2]<--a.vec[group2]/sum(a.vec[group2])

return(a)
}

##################

makeB<-function(groups,a.vec)
{
return((1/sum(a.vec[groups$group1])+1/sum(a.vec[groups$group2]))^(-0.5))

}


##################

groups<-function(K,a.vec,type="all.partitions")
{

if (K!=9&&(type=="inter.ord"|type=="inter.cat"))
 {
 stop("With K!=9, can not be inter.ord or inter.cat")
 }

if (type=="all.partitions")
 {
 groups_A<-c()
 groups<-list()
 for (i in 1:floor(K/2))
  {                      
  groups_A<-c(groups_A,combn(K,i,simplify=FALSE))
  }
  
 for (i in 1:length(groups_A))
  {
  groups[[i]]<-list()
  groups[[i]]$group1<-groups_A[[i]]
  groups[[i]]$group2<-setdiff(1:K,groups_A[[i]])
  }
 }
 
if (type=="all.pairs")
 {
 groups_matrix<-combn(K,2)
 groups<-list()
 for (i in 1:ncol(groups_matrix))
  {
  groups[[i]]<-list(group1=groups_matrix[1,i],group2=groups_matrix[2,i])
  }
 }

if (type=="ord")
 {
 groups<-list()
 for (k in 1:(K-1))
  {          
  groups[[k]]<-list(group1=1:k,group2=(k+1):K)
  }

 }

 
if (type=="inter.ord")
 {
 groups<-list()
 groups[[1]]<-list(group1=1,group2=setdiff(1:9,1))
 groups[[2]]<-list(group1=c(1,2),group2=setdiff(1:9,1:2))
 groups[[3]]<-list(group1=3,group2=setdiff(1:9,3))
 groups[[4]]<-list(group1=c(2,3),group2=setdiff(1:9,2:3))
 groups[[5]]<-list(group1=7,group2=setdiff(1:9,7))
 groups[[6]]<-list(group1=c(7,8),group2=setdiff(1:9,7:8))
 groups[[7]]<-list(group1=9,group2=setdiff(1:9,9))
 groups[[8]]<-list(group1=c(8,9),group2=setdiff(1:9,8:9))
 groups[[9]]<-list(group1=3,group2=setdiff(1:9,3))
 groups[[10]]<-list(group1=c(3,6),group2=setdiff(1:9,c(3,6)))
 groups[[11]]<-list(group1=c(1,4,2,5),group2=setdiff(1:9,c(1,4,2,5)))
 groups[[12]]<-list(group1=c(2,5,3,6),group2=setdiff(1:9,c(2,5,3,6)))
 groups[[13]]<-list(group1=c(4,7),group2=setdiff(1:9,c(4,7)))
 groups[[14]]<-list(group1=c(6,9),group2=setdiff(1:9,c(6,9)))
 groups[[15]]<-list(group1=c(4,7,5,8),group2=setdiff(1:9,c(4,7,5,8)))
 groups[[16]]<-list(group1=c(5,8,6,9),group2=setdiff(1:9,c(5,8,6,9)))
 }

if (type=="inter.cat")
 {
 groups<-list()
 groups[[1]]<-list(group1=c(1),group2=setdiff(1:9,1))
 groups[[2]]<-list(group1=c(1,2),group2=setdiff(1:9,c(1,2)))
 groups[[3]]<-list(group1=c(3),group2=setdiff(1:9,3))
 groups[[4]]<-list(group1=c(2,3),group2=c(1:9,c(2,3)))
 groups[[5]]<-list(group1=c(7),group2=setdiff(1:9,7))
 groups[[6]]<-list(group1=c(7,8),group2=setdiff(1:9,c(7,8)))
 groups[[7]]<-list(group1=c(8,9),group2=setdiff(1:9,c(8,9)))
 groups[[8]]<-list(group1=c(9),group2=setdiff(1:9,9))
 groups[[9]]<-list(group1=c(1,4),group2=setdiff(1:9,c(1,4)))
 groups[[10]]<-list(group1=c(3,6),group2=setdiff(1:9,c(3,6)))
 groups[[11]]<-list(group1=c(1,4,2,5),group2=setdiff(1:9,c(1,4,2,5)))
 groups[[12]]<-list(group1=c(2,5,3,6),group2=setdiff(1:9,c(2,5,3,6)))
 groups[[13]]<-list(group1=c(4,7),group2=setdiff(1:9,c(4,7)))
 groups[[14]]<-list(group1=c(6,9),group2=setdiff(1:9,c(6,9)))
 groups[[15]]<-list(group1=c(4,7,5,8),group2=setdiff(1:9,c(4,7,5,8)))
 groups[[16]]<-list(group1=c(5,8,6,9),group2=setdiff(1:9,c(5,8,6,9)))
 groups[[17]]<-list(group1=c(2),group2=setdiff(1:9,c(2)))
 groups[[18]]<-list(group1=c(2,5),group2=setdiff(1:9,c(2,5)))
 groups[[19]]<-list(group1=c(5,8),group2=setdiff(1:9,c(5,8)))
 groups[[20]]<-list(group1=c(8),group2=setdiff(1:9,c(8)))
 groups[[21]]<-list(group1=c(5),group2=setdiff(1:9,c(5)))
 groups[[22]]<-list(group1=c(4),group2=setdiff(1:9,c(4)))
 groups[[23]]<-list(group1=c(4,5),group2=setdiff(1:9,c(4,5)))
 groups[[24]]<-list(group1=c(5,6),group2=setdiff(1:9,c(5,6)))
 groups[[25]]<-list(group1=c(6),group2=setdiff(1:9,c(6)))
 }

if (type=="inter.ord.main")
 {
 groups<-list()
 groups[[1]]<-list(group1=1,group2=setdiff(1:9,1))
 groups[[2]]<-list(group1=c(1,2),group2=setdiff(1:9,1:2))
 groups[[3]]<-list(group1=3,group2=setdiff(1:9,3))
 groups[[4]]<-list(group1=c(2,3),group2=setdiff(1:9,2:3))
 groups[[5]]<-list(group1=7,group2=setdiff(1:9,7))
 groups[[6]]<-list(group1=c(7,8),group2=setdiff(1:9,7:8))
 groups[[7]]<-list(group1=9,group2=setdiff(1:9,9))
 groups[[8]]<-list(group1=c(8,9),group2=setdiff(1:9,8:9))
 groups[[9]]<-list(group1=3,group2=setdiff(1:9,3))
 groups[[10]]<-list(group1=c(3,6),group2=setdiff(1:9,c(3,6)))
 groups[[11]]<-list(group1=c(1,4,2,5),group2=setdiff(1:9,c(1,4,2,5)))
 groups[[12]]<-list(group1=c(2,5,3,6),group2=setdiff(1:9,c(2,5,3,6)))
 groups[[13]]<-list(group1=c(4,7),group2=setdiff(1:9,c(4,7)))
 groups[[14]]<-list(group1=c(6,9),group2=setdiff(1:9,c(6,9)))
 groups[[15]]<-list(group1=c(4,7,5,8),group2=setdiff(1:9,c(4,7,5,8)))
 groups[[16]]<-list(group1=c(5,8,6,9),group2=setdiff(1:9,c(5,8,6,9)))
 groups[[17]]<-list(group1=c(1:3),group2=c(4:9))
 groups[[18]]<-list(group1=c(1:6),group2=c(7:9))
 groups[[19]]<-list(group1=c(1,4,7),group2=c(2,3,5,6,8,9))
 groups[[20]]<-list(group1=c(3,6,9),group2=c(1,2,4,5,7,8))
 }

m<-length(groups)
if (m==1)
 {
 return(groups)
 }
duplicate<-c()

for (i in 1:m)
 {
 if (sum(a.vec[groups[[i]]$group1])==0|sum(a.vec[groups[[i]]$group2])==0)
  {
  duplicate<-c(duplicate,i)
  }
 groups[[i]]$group1<-groups[[i]]$group1[a.vec[groups[[i]]$group1]>0] 
 groups[[i]]$group2<-groups[[i]]$group2[a.vec[groups[[i]]$group2]>0]
 if (i>1)
  {
  for (j in 1:(i-1))
   {
   if ((setequal(groups[[i]]$group1,groups[[j]]$group1)&setequal(groups[[i]]$group2,groups[[j]]$group2))|(setequal(groups[[i]]$group2,groups[[j]]$group1)&setequal(groups[[i]]$group1,groups[[j]]$group2)))
    {
    duplicate<-c(duplicate,i)
    }
   }
  }
 }
duplicate<-union(duplicate,duplicate)
if (length(duplicate)==0)
 return(groups)
else
 return(groups[-duplicate])

}