###    Class prediction based on random forests using clinical parameters only
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

rf_z<-function(Xlearn=NULL,Zlearn,Ylearn,Xtest=NULL,Ztest,...)
{
Ylearn<-as.numeric(factor(Ylearn))-1
nlearn<-length(Ylearn)
data.learn<-data.frame(Zlearn,y=factor(Ylearn))
data.test<-data.frame(Ztest)

output.forest<-cforest(formula=y~.,data=data.learn,controls=cforest_control(ntree=200,mincriterion=qnorm(0.5),mtry=floor(sqrt(ncol(data.learn)-1)),replace=FALSE))
importance<-varimp(output.forest)

prediction<-as.numeric(predict(object=output.forest,newdata=data.test))-1
importance<-varimp(output.forest)
OOB<-sum(predict(output.forest,OOB=TRUE)!=Ylearn)/nlearn
return(list(prediction=prediction,importance=importance,OOB=OOB))
}


