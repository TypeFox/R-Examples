###   Area under the curve (AUC)
###
### Copyright 2015-07 Anne-Laure Boulesteix 
###
### 
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
 


my.auc<-function(linpred,Y)
{
sub.auc<-function(x,y) {return(as.numeric(x>y)+0.5*as.numeric(x==y))}

return(sum(kronecker(linpred[Y==1],linpred[Y==0],FUN=sub.auc))/(sum(Y==1)*sum(Y==0)))
}


