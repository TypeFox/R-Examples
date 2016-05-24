### Using an IPF-lasso model for prediction of new observations
###
### Copyright 2015-07 Anne-Laure Boulesteix 
###
### Using an IPF-lasso model for prediction of new observations
###
###
### This file is part of the `ipflasso' library for R and related languages.
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
 



ipflasso.predict<-function(object,Xtest)
{
coeff<-object$coeff[,object$ind.bestlambda]
if (object$family=="gaussian"|object$family=="binomial")
 {
 linpredtest<-coeff[1]+Xtest%*%coeff[-1]
 }

if (object$family=="cox")
 {
 linpredtest<-Xtest%*%coeff[-1]
 }


if (object$family=="gaussian"|object$family=="cox")
 {
 classtest<-NULL
 probabilitiestest<-NULL
 }

if (object$family=="binomial")
 {
 probabilitiestest<-as.numeric(plogis(linpredtest))
 classtest<-as.numeric(linpredtest>0)
 }

return(list(linpredtest=linpredtest,classtest=classtest,probabilitiestest=probabilitiestest))
}