###    Class prediction based on support vector machines
### 	using microarray data only
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


svm_x<-function(Xlearn,Zlearn=NULL,Ylearn,Xtest,Ztest=NULL,ordered=NULL,nbgene=NULL,...)
{

Ylearn<-as.numeric(factor(Ylearn))-1
nlearn<-length(Ylearn)

if (is.null(ordered))
 {
 ordered<-1:ncol(Xlearn)
 }

output.svm<-svm(x=Xlearn[,ordered],y=factor(Ylearn),kernel="linear",...)
prediction<-as.numeric(predict(object=output.svm,newdata=Xtest[,ordered]))-1

return(list(prediction=prediction))


}

