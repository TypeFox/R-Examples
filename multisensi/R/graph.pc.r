# Multisensi R package ; file graph.pc.r (last modified: 2016-04-18) 
# Copyright INRA 2011-2015 
# Authors: C. Bidot, M. Lamboni, H. Monod
# MaIAGE, INRA, Univ. Paris-Saclay, 78350 Jouy-en-Josas, France
#
# More about multisensi in http://cran.r-project.org/web/packages/multisensi/
#
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use, 
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info". 
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability. 
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or 
# data to be ensured and,  more generally, to use and operate it in the 
# same conditions as regards security. 
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.
#
#===========================================================================
graph.pc <-function(x, nb.plot=15, nb.comp=NULL, xmax=NULL, beside=TRUE, cor.plot=FALSE, xtick=TRUE, type="l",...)
#===========================================================================
{
    ##x :          GSI objects
    ##nb.plot:       A number decribing the max number of factor bars to be ploted

    ## modified HM 1/3/2010 until the end of the function (cf. "toPlot" and "toWrite")

  if (is.null(nb.comp)){
    nbcomp <- ncol(x$L)
  } else {
    nbcomp <- min(nb.comp,ncol(x$L))
  }

  inertie <- rep(0,nbcomp)
  inertie[1] <- x$inertia[1]
  if(nbcomp>1){
    for(k in 2:nbcomp) {
      inertie[k] <- x$inertia[k]-x$inertia[k-1]
    }
  }
#  names.comp <- paste(paste("PC",1:nbcomp,sep=""),signif(inertie,3),sep=" (" )
#  names.comp <- paste(colnames(x$L)[1:nbcomp],signif(inertie,3),sep=" (" )
#  names.comp <- paste(names.comp," %)",sep="")
#  main.comp <- names.comp
  main.comp <- paste(colnames(x$L)[1:nbcomp],signif(inertie,3),sep=" (" )
  main.comp <- paste(main.comp," %)",sep="")


  if(cor.plot){ # on trace la sortie cor de multisensi = sdH*L
    if(x$normalized){
        toWrite <- "Correlation"
    }
    else{
        toWrite <- "Weighted loadings"
    }

    toPlot <- x$cor
    if(toPlot[floor(nrow(toPlot)/2),1]<0){ toPlot[,1] <- -toPlot[,1]}

    corrmin <- min(toPlot,na.rm=TRUE)
    corrmax <- max(toPlot,na.rm=TRUE)

    if(nbcomp>1){ if(toPlot[1,2]<0){toPlot[,2] <- -toPlot[,2]}}

    par(mfrow=c(2,nbcomp))

    for (k in 1:nbcomp){
        plot(toPlot[,k], ylim=c(corrmin, corrmax) , type=type, col="blue", main= main.comp[k],  ylab=toWrite, lwd=3 ,cex.axis=2 ,cex=4, xlab="", xaxt="n",...)#xlab="time",
        abline(h=0)
        if(xtick){
          axis(1,at=1:nrow(toPlot),labels=colnames(x$Y))
        }else{
          axis(1,at=NULL,labels=TRUE,tick=TRUE)
        }
    }

    for(k in 1:nbcomp){
        graph.bar(x ,k ,nb.plot, xmax=xmax, beside=beside, ...)
    }

  }else{ # on trace la variabilite des vecteurs de base L multiplies par leurs coeff H H*L
    par(mfrow=c(2,nbcomp))
    for(k in 1:nbcomp){
      quantH=quantile(x$H[,k],c(0.50,0,0.25,0.75,1,0.10,0.90))
      Lk=matrix(rep(x$L[, k],times=length(quantH)),ncol=length(quantH))
      qHLk=Lk%*%diag(quantH)

      matplot(1:nrow(qHLk),qHLk[,2:length(quantH)],type="l",main=main.comp[k],col=c("red","black","black","red","blue","blue"),lty=c(2,3,4,6,1,5),ylim=c(min(qHLk,na.rm=TRUE),max(qHLk,na.rm=TRUE)),xlab="",ylab="",xaxt="n")
      polygon(c(1:nrow(qHLk),seq(from=nrow(qHLk),to=1,by=-1)),c(qHLk[,3],qHLk[seq(from=nrow(qHLk),to=1,by=-1),4]),col="gray",lty=0)
      lines(1:nrow(qHLk),qHLk[, 1],lwd=2,col=1)
      if(xtick){
        ww=axTicks(1)
        ww[ww==0]=1
        axis(1,at=ww,labels=colnames(x$Y)[ww])
      }else{
        axis(1,at=NULL,labels=TRUE,tick=TRUE)
      }

#      if(k==1){
#        legend("topright",legend=c("median","min H", "1st quart.", "3rd quart.","max H"),col=c(1:length(quantH)),lty=c(1:length(quantH)),lwd=c(2,rep(1,length(quantH))))
#      }
    }

    for(k in 1:nbcomp){
        par(mar=c(5.1, 4.1, 0 ,2.1))
        graph.bar(x ,k ,nb.plot, xmax=xmax, beside=beside, xlab="SI",...)
    }

  }

}

