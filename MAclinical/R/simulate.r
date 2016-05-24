###    Simulating artificial data sets
###
### Copyright 2007-11 Anne-Laure Boulesteix 
###
### 
###
###
### This file is part of the `MAclinical' library for R and related languages.
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


simuldata_list<-function(niter=50,n=500,p=1000,psig=50,q=5,muX=0,muZ=0)
{

datalist<-list()
if (psig>p)
 stop("psig must be <p")

for (i in 1:niter)
 {
 set.seed(i)
 y<-numeric(n)
 y[sample(n,n/2)]<-1

 x<-matrix(rnorm(n*p),n,p)
 x[y==1,1:psig]<-x[y==1,1:psig]+muX

 z<-matrix(rnorm(n*q),n,q)
 z[y==1,]<-z[y==1,]+muZ
 datalist[[i]]<-list(y=y,x=x,z=z)
 }

return(datalist)

}


##################

simuldatacluster_list<-function(niter=50,n=500,p=1000,psig=50,q=5,muX=0,muZ=0)
{
if (psig>p)
 stop("psig must be <p")
datalist<-list()

for (i in 1:niter)
 {
 set.seed(i)
 y<-numeric(n)
 y[sample(n,n/2)]<-1
 yy<-y
 yy[y==1]<-sample(1:2,sum(y==1),prob=c(1/2,1/2),replace=TRUE)

 x<-matrix(rnorm(n*p),n,p)
 x[yy==2,1:psig]<-x[yy==2,1:psig]+muX

 z<-matrix(rnorm(n*q),n,q)
 z[yy==1,]<-z[yy==1,]+muZ
 datalist[[i]]<-list(y=y,x=x,z=z)
 }

return(datalist)

}


