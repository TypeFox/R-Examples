###    Computes maximally selected chi-square statistics
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



maxsel<-function(x1,x2=NULL,y,type="inter.ord")
{

if (length(setdiff(y,c(0,1)))>0)
 stop("The entries of y must be 0 or 1")

if (is.null(x2)&is.element(type,c("inter.ord","inter.cat","inter.ord.main")) )
 {
 stop("if x2 is null, type can not be inter.ord, inter.cat or main")
 }

if (!is.null(x2)&!is.element(type,c("inter.ord","inter.cat","inter.ord.main")) )
 {
 stop("if x2 is not null, type must be inter.ord, inter.cat or main")
 }

if (!is.null(x2))
 {
 x<-transf.inter(x1,x2)
 }
else
 {
 x<-x1
 }


dat<-cbind(x,y)
dat<- subset(dat, complete.cases(dat))
y<-dat[,2]
x<-dat[,1]


if (any(x<0))
 stop("x must be a vector of positive integers")

if (type=="all.partitions")
 {
 crit<-c()
 x<-as.factor(x)
 K<-nlevels(x)
 levels(x)<-1:K
 prop2<-as.numeric(tapply(x[y==1],x[y==1],length))/as.numeric(tapply(x,x,length))
 x<-factor(x,levels=order(prop2))
 x<-as.numeric(x)
 a<-list()
 for (k in 1:(K-1))
  {
  a[[k]]<-c(1:k)
  }
 crit<-as.numeric(sapply(a,FUN=chisq.test.a,x,y,partition=TRUE))
 maxselcrit<-max(crit)
 }
 
 
if (type=="all.pairs")
 {
 crit<-c()
 x<-as.factor(x)
 K<-nlevels(x)
 levels(x)<-1:K
 x<-as.numeric(x)
 a<-combn(K,2,simplify=FALSE)
 crit<-as.numeric(sapply(a,FUN=chisq.test.a,x,y,partition=FALSE))
 maxselcrit<-max(crit) 
 }
 
if (type=="ord")
 {        
 x<-as.factor(x)
 K<-nlevels(x)
 levels(x)<-1:K
 x<-as.numeric(x)
 a<-list()
 for (k in 1:(K-1))
  {
  a[[k]]<-c(1:k)
  }
 crit<-as.numeric(sapply(a,FUN=chisq.test.a,x,y,partition=TRUE))
 maxselcrit<-max(crit)
 }

if (type=="inter.ord")
 {
 a<-list(x1=1,x2=c(1,2),x3=3,x4=c(2,3),x5=7,x6=c(7,8),x7=c(8,9),x8=9,x9=c(1,2,4,5),
 x10=c(4,5,7,8),x11=c(5,6,8,9),x12=c(2,3,5,6),x13=c(1,4),x14=c(4,7),x15=c(3,6),
 x16=c(6,9))
 crit<-as.numeric(sapply(a,FUN=chisq.test.a,x,y,partition=TRUE))  
 maxselcrit<-max(crit)
 }
 
if (type=="inter.cat")
 {
 a<-list(x1=1,x2=c(1,2),x3=3,x4=c(2,3),x5=7,x6=c(7,8),x7=c(8,9),x8=9,x9=c(1,2,4,5),
 x10=c(4,5,7,8),x11=c(5,6,8,9),x12=c(2,3,5,6),x13=c(1,4),x14=c(4,7),x15=c(3,6),
 x16=c(6,9),x17=2,x18=c(2,5),x19=c(5,8),x20=8,x21=5,x22=4,x23=c(4,5),x24=c(5,6),x25=6)
 crit<-as.numeric(sapply(a,FUN=chisq.test.a,x,y,partition=TRUE))  
 maxselcrit<-max(crit)
 }

if (type=="inter.ord.main")
 {
 a<-list(x1=1,x2=c(1,2),x3=3,x4=c(2,3),x5=7,x6=c(7,8),x7=c(8,9),x8=9,x9=c(1,2,4,5),
 x10=c(4,5,7,8),x11=c(5,6,8,9),x12=c(2,3,5,6),x13=c(1,4),x14=c(4,7),x15=c(3,6),
 x16=c(6,9),x17=1:3,x18=7:9,x19=c(1,4,7),x20=c(3,6,9))
 crit<-as.numeric(sapply(a,FUN=chisq.test.a,x,y,partition=TRUE))  
 maxselcrit<-max(crit)
 }

 
return(maxselcrit)


}


#######

chisq.test.a<-function(a,x,y,partition)
{

if (partition==FALSE)
 {

 if (length(a)!=2)
  stop("stop")

 selected<-which(is.element(x,a))
 x<-x[selected]
 y<-y[selected]
 xa<-as.numeric(x==a[1])
 }
else
 {
 xa<-is.element(x,a)
 }
if (is.element(sum(xa),c(0,length(xa))))
 {
 return(0)
 }
else
 {
 return(chisq.test(xa,y,correct=FALSE)$statistic)
 }
}
