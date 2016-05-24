###    Class prediction based on logistic regression using clinical parameters only
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

logistic_z<-function(Xlearn=NULL,Zlearn,Ylearn,Xtest=NULL,Ztest,...)
{
Ylearn<-as.numeric(factor(Ylearn))-1
data.learn<-data.frame(Zlearn,y=Ylearn)
data.test<-data.frame(Ztest)

model<-glm(data=data.learn,formula=y~.,family=binomial)

prediction<-as.numeric(predict(object=model,newdata=data.test)>0)
return(list(prediction=prediction))

}