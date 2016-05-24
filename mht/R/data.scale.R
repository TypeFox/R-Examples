# Author : F.Rohart,  Australian Institute for Bioengineering and Nanotechnology, The University of Queensland, Brisbane, QLD
# created: 18-09-2014
# last modification: 10-10-2014
# Copyright (C) 2014
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

## scaling the data and adding intercept if missing


#si la matrice de depart ne contient pas l'intercept (colonne de 1) on la rajoute
data.scale=function(data,warning)
{
    if(missing(warning)) warning=TRUE
    n=nrow(data)
    p=ncol(data)
    
    if(is.null(colnames(data)))
    {
        colnames(data)=paste("X", 1:p,sep="")
    }
    col.names=colnames(data)
    #if(((sum(data[,1])-sqrt(n))<1e-10) & ((mean(data[,1])-1/sqrt(n))<1e-10))
    if(sum(sum(data[,1]==1)==n))
    {
        temp=scale(data[,-1])
        data=cbind(data[,1],temp/sqrt(n-1))
        data[which(is.na(data))]=0
        intercept=TRUE
        means.data=c(0,attr(temp,"scaled:center")) #not center the intercept
        sigma.data=c(1,attr(temp,"scaled:scale")*sqrt(n-1)) #not scale the intercept
        colnames(data)=col.names
    }else{
        if(warning==TRUE)
        {
            message("intercept has been added")
        }
        temp=scale(data)
        data=temp/sqrt(n-1)
        data[which(is.na(data))]=0
        data=cbind(1,data)
        #data=cbind(1/sqrt(n),data)
        intercept=FALSE
        means.data=attr(temp,"scaled:center")
        sigma.data=attr(temp,"scaled:scale")*sqrt(n-1)
        colnames(data)=c("Intercept",col.names)
    }
    



    out=list(data=data,intercept=intercept,means.data=means.data,sigma.data=sigma.data)
}
