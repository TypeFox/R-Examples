# Author : F.Rohart,  Australian Institute for Bioengineering and Nanotechnology, The University of Queensland, Brisbane, QLD
# created: pre 01-01-2013
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

plot.mht <- plot.mht.order <-function(x,which.plot=1:4,...)
{
    
    alpha=as.numeric(colnames(x$coefficients))

    par(ask=TRUE)
    
    for(i in 1:length(alpha))
    {
        if(1%in%which.plot)
        {
            plot(x$fitted.values[,i],x$data$Y,xlab="Fitted values",ylab="Y",main=paste("Y vs fitted values, alpha =",alpha[i]))
            abline(c(0,1),lty=3)
        }
        
        if(2%in%which.plot)
        {
            plot(x$fitted.values[,i],x$residuals[,i],xlab="Fitted values",ylab="Residuals",main=paste("residuals vs fitted values, alpha =",alpha[i]))
        abline(h=0,lty=3)
        }
        
        if(3%in%which.plot)
        {
            qqnorm(x$residuals[,i]/sqrt(var(x$residuals[,i])),ylab="Standardised residuals",main=paste("Normal Q-Q plot, alpha =",alpha[i]))
        abline(c(0,1),lty=3)
        }
        
        if(4%in%which.plot)
        {
            barplot(x$coefficients[,i],las=2,ylim=range(x$coefficients[,i]),axes=TRUE,axis.lty=1,,main=paste("Coefficients, alpha =",alpha[i]),horiz=FALSE)
        }
        
        
    }
    par(ask=FALSE)

}


plot.bolasso <-function(x,...)
{
    col=topo.colors(ncol(x$data$X))
    ind=sort(x$mu,index.return=TRUE)
    mu=ind$x
    plot(mu,x$frequency[1,],xlim=range(mu)[2:1],type="n",col=col[1],ylim=c(0,1),xlab="Regularization parameter",ylab="Frequency of selection",main="Bolasso")
    for(i in 1:ncol(x$data$X))
    {
        lines(mu,x$frequency[i,][ind$ix],col=col[i])
    }
    
}

