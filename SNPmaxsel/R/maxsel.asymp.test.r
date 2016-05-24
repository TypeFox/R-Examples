###    Test of independence based on maximally selected statistics
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



maxsel.asymp.test<-function(x1,x2=NULL,y,type)
{

maxselcrit<-maxsel(x1=x1,x2=x2,y=y,type=type)

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
n<-length(x)
            
if (is.element(type,c("inter.ord","inter.cat","inter.ord.main")))
 {
 K<-9
 }
else 
 {
 K<-max(x)
 }
  
a.vec<-numeric(K)
                   
for (i in 1:K)
 {
 a.vec[i]<-sum(x==i)/n
 }
  
pval<-Fasymp(maxselcrit,a.vec=a.vec,type=type)
 
  
return(list(maxselstat=maxselcrit,value=pval))


}



