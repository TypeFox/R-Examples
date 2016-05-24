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
#       Figure Supp. Stability
#-------------------------------------------




plot.bootsPLS=function(x,light,pch,col,legend.position,title,...)
{
    
    if(missing(x)) stop("`x' is missing")
    H=ncol(x$nbr.var)
    
    if(missing(light)) light=0
    if(missing(pch)) pch=1:H
    if(missing(col)) col=1:H
    if(missing(legend.position)) legend.position="bottomright"
    if(missing(title)) title="Frequency of selection depending on the component"
    
    p=ncol(x$frequency)
    if(light>0) #design to remove the variable that are never selected from the plot
    {
        ind=apply(x$frequency,2,function(x){sum(x>light)})#genes that are selected more than light on at least one component
        ind=which(ind>0)
    }else{ind=1:p}
    
    ylab="Frequency"
    
    par(mar=c(10,4,4,2)+0.1)
    plot(ind,x$frequency[1,ind],ylab=ylab,xlab="",pch=pch[1],col=col[1],xaxt="n",ylim=c(0,1),main=title,...)
    axis(side=1,at=1:ncol(x$data$X),labels=substr(colnames(x$data$X),1,30),las=2)
    
    for(i in c(1:H,1)) points(ind,x$frequency[i,ind],col=col[i],pch=pch[i],...)

    if(legend.position!=FALSE)
    legend(legend.position,legend=paste("Comp",1:H),pch=pch,col=col)
    #legend(ncol(x$frequency)-2000,0.3,legend=paste("Comp",1:H),pch=pch,col=col)
    
   
    
}#end function plot.stability

# exple
# plot.stability(x$frequency)
