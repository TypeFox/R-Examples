# Author : F.Rohart,  Australian Institute for Bioengineering and Nanotechnology, The University of Queensland, Brisbane, QLD
# created: 29-05-2014
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



#-------------------------------------------
#       plot the p.value from variable selection
#-------------------------------------------


plot.component.selection=function(x,pch,col,type,lty,...)
{
    if(class(x)!="component.selection") stop("problem")
    
    
    ncomp=ncol(x$object$nbr.var)
    alpha=x$alpha
    pval=x$pval
    
    if(missing(pch)) pch=1
    if(missing(col)) col=1
    if(missing(lty)) lty=1
    if(missing(type)) type="b"
    
    
    par(mar=c(5,5,3,2)+0.1)
    
    plot(1:ncomp,-log((pval),10),xlab=expression(bold("Components")),ylab=expression(bold(paste("-log"[10]," (p-value)"))),
            pch=pch,type=type,lty=lty,ylim=c(0,max(c(-log(pval,10),-log(alpha,10)),na.rm=TRUE)),col=col,xaxt="n",...)
            axis(side=1,at=1:(ncomp))
            
    abline(h=-log(alpha,10))
    
    if(x$opt>1)
    {
        text(x$opt,-log(pval[x$opt],10),paste(x$opt,"components"),pos=2,col=col)
    }else{
        text(x$opt,-log(pval[x$opt],10),paste(x$opt,"component"),pos=2,col=col)
    }

}


plot.variable.selection=function(x,pch,col,type,lty,...)
{
   if(class(x)!="variable.selection") stop("problem")
   

   
   ncomp=length(x$opt)
   limit=lapply(x$pval,length)
   frequency=x$object$frequency
   alpha=x$alpha
   
   if(missing(pch)) pch=1:ncomp
   if(missing(col)) col=1:ncomp
   if(missing(lty)) lty=1:ncomp
   if(missing(type)) type="b"
   
   
   par(mar=c(5,5,3,2)+0.1)
   
   for(compi in 1:ncomp)
   {
       variables=colnames(x$subsamplings[[compi]])
       ind=which(!is.na(x$pval[[compi]]))
       pval=x$pval[[compi]]
       
       if(compi==1)
       {
           plot(frequency[compi,variables[1:limit[[compi]]][ind]],-log((pval)[ind],10),xlab=expression(bold("Stability")),ylab=expression(bold(paste("-log"[10]," (p-value)"))),
           pch=pch[compi],type=type,lty=lty[compi],xlim=c(1,0),ylim=c(0,max(-log(unlist(x$pval),10),na.rm=TRUE)),col=col[compi],...)

       }else{
           
           points(frequency[compi,variables[1:limit[[compi]]][ind]],-log((pval)[ind],10),pch=pch[compi],type=type,col=col[compi],lty=lty[compi],...)

       }
       
       text(frequency[compi,variables[1:x$opt[compi]]][x$opt[compi]],-log(pval[x$opt[compi]],10),paste(x$opt[compi],"variables"),pos=3,col=col[compi])
       text(frequency[compi,variables[1:x$opt[compi]]][x$opt[compi]],-log(pval[x$opt[compi]],10),paste("comp",compi),pos=1,col=col[compi])

       
       
   }
   abline(h=-log(alpha,10))
   
}
