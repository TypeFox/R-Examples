######################################################################
#
# gfitpca.R
#
# copyright (c) 2014-07-10, Katalin Csillery, Louisiane Lemaire, 
# Olivier Francois and Michael GB Blum
#
#     This program is free software; you can redistribute it and/or
#     modify it under the terms of the GNU General Public License,
#     version 3, as published by the Free Software Foundation.
# 
#     This program is distributed in the hope that it will be useful,
#     but without any warranty; without even the implied warranty of
#     merchantability or fitness for a particular purpose.  See the GNU
#     General Public License, version 3, for more details.
# 
#     A copy of the GNU General Public License, version 3, is available
#     at http://www.r-project.org/Licenses/GPL-3
#
# Part of the R/abc package
#
######################################################################


gfitpca=function(target, sumstat, index, cprob=0.1, xlim=NULL, ylim=NULL, ...){

   loc2plot=function(x, y, cprob ,...)
   {
      fit=locfit(~lp(x, y, nn=.2), maxk=100, mint=100, maxit=100)
      lev=sort(fitted(fit))[floor(cprob*length(x))]
      plot(fit, lev=lev, m=100, drawlabels = FALSE, ...)
      return (list(fit=fit, lev=lev))
   }

   #when target is a vector 
   if (is.vector(target)){
      target=t(as.data.frame(target))
      if (is.data.frame(sumstat)){colnames(target)=names(sumstat)}
      if (is.matrix(sumstat)){colnames(target)=colnames(sumstat)}
   }

   #acp
   res.prcomp=prcomp(sumstat, scale=T, center=T)

   nmod=length(table(index))
   theindex=names(table(index))

   #plot   
   if (is.null(xlim)){xlim=ylim}
   if (is.null(ylim)){ylim=xlim}
   
   if (!is.null(xlim)){
      plot(0, type="n", xlim=xlim, ylim=ylim, xlab="PC1", ylab="PC2")
   }
   for (i in 1:nmod){
      ind=index==theindex[i]
      if ((i==1)&(is.null(xlim))){add=FALSE}else{add=TRUE}
      loc2plot(res.prcomp$x[ind,1], res.prcomp$x[ind,2], cprob, col=i, lty=1, lwd=2, add=add,...)
   }

   #observed data
   points(predict(res.prcomp, target)[1], predict(res.prcomp, target)[2], col=1, cex=2, pch=3, lwd=2)

   #legend
   legend("topright", legend=theindex, cex=1.5, col=c(1:nmod), lty=0, pch=15)


}
