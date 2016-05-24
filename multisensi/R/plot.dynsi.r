# Multisensi R package ; file plot.dynsi.r (last modified: 2016-02-03) 
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
plot.dynsi <- function(x, normalized=FALSE, text.tuning=NULL, shade=FALSE, color=NULL, xtick=TRUE, total.plot=FALSE, ...)
#===========================================================================
{
    ## new argument "normalized" : FALSE => SI plotted within var(Y)
    ## changed by HM, 2013-05-03 : new arguments "shade" and "color"
    ## ENTREE
    ##               x:         dynsi  object

    ##SORTIE
    ##              Graphique
  VarY=apply(x$Y,2,var)
  jmax=which.max(VarY)

  if(x$call.info$analysis=="anova"){
    # dans le cadre de l'anova on traite les interactions
    filtre.main <- colSums(x$Att)==1 #apply(x$Att,2,sum)==1
    filtre.inter <- !filtre.main
    interaction <-  filtre.inter %*%as.matrix(x$SI)
#    interaction <- colSums(x$iSI)/2 # pourrait remplacer les 3 lignes precedentes dans les cas d'interaction 2a 2
    tab.indicesmi <- as.data.frame(rbind( as.matrix(x$mSI[order(x$mSI[,jmax]),]),interaction))
    rownames(tab.indicesmi)[nrow(tab.indicesmi)]<- "interaction"
  } else if(x$call.info$analysis=="sensitivity"){
    tab.indicesmi <- as.data.frame(x$mSI[order(x$mSI[,jmax]),])
    #rownames(tab.indicesmi)<- rownames( x$mSI)
  }

#  tab.indicesmi[1:nrow(x$mSI),]=tab.indicesmi[order(x$mSI[,jmax]),]
  indice.cum <-apply( tab.indicesmi,2,cumsum)

    ## ajout CB
    indice.cum=rbind(indice.cum,1)#100)
  if(x$call.info$analysis=="anova"){
    rownames(indice.cum)[nrow(indice.cum)]="residual"
  } else if(x$call.info$analysis=="sensitivity"){
    rownames(indice.cum)[nrow(indice.cum)]="interaction";
  }

    if(!normalized){
        #sdY=sqrt(apply(x$Y,2,var))
        for (j in 1:ncol(indice.cum)) {
          indice.cum[,j] <- indice.cum[,j]* VarY[j]#/100 
        }
    }


    ## recherche des positions pour afficher les noms des parametres

#    nb.outp <- ncol(tab.indicesmi)
#    names <- rownames(tab.indicesmi)
#    n.name <-1:(nrow(tab.indicesmi)+1)
    nb.outp <- ncol(indice.cum)
    names <- rownames(indice.cum)
    n.name <-1:nrow(indice.cum)

    if(is.null(text.tuning)){ text.tuning=floor(nb.outp/10)+1}
    if(text.tuning> floor(nb.outp/10)){text.tuning=floor(nb.outp/10)+1}

    for(k in 1:nrow(tab.indicesmi)){

        max.k <- max(tab.indicesmi[k,text.tuning:(nb.outp -text.tuning)],na.rm=TRUE)

        for(i in text.tuning:(nb.outp-text.tuning)){
            if(max.k==tab.indicesmi[k,i]){ n.name[k] <- i}
        }
    }

    min.inter <- min( indice.cum[nrow(indice.cum)-1,text.tuning:(nb.outp -text.tuning)],na.rm=TRUE )
    for(i in text.tuning:(nb.outp-text.tuning)){
        if(min.inter== indice.cum[nrow(indice.cum)-1,i]){ n.name[length(n.name)] <- i}
    }

    ## changed by HM, 2013-05-03
    if(shade & is.null(color)) densite <- c(7,1,5,10,15,19,23,4) 
    else densite <- rep(NULL, 8) 
    if(!is.null(color)) colset <- color(nrow(tab.indicesmi)) 
    else colset <- rep(NA, nrow(tab.indicesmi)) 

    #on ouvre un nouveau device si besoin
    if(dev.cur() == 1 | dev.interactive()) dev.new()
    old.par <- par(no.readonly = TRUE)

    if(normalized){
      # plot des outputs au-dessus des indices
      Y1=colMeans(x$Y)-sqrt(VarY)
      Y2=colMeans(x$Y)+sqrt(VarY)

      par(fig=c(0,1,0.7,1),mar=c(0,4.1, 4.1, 2.1))
      #pour limites et etiquettes
      ytck=pretty(c(min(x$Y,na.rm=TRUE),max(x$Y,na.rm=TRUE)))
      #max
      plot(1:length(VarY), apply(x$Y,2,max), type="l", xaxt="n", xlab="", ylab="Outputs", yaxp=c(min(ytck),max(ytck),1), lwd=1, ylim=c(min(ytck),max(ytck)), lty=4, cex.axis=.7,las=1)
      #moy +/- sdY
      if(shade || !is.null(color)){
        polygon(c(1:length(VarY),seq(from=length(VarY),to=1,by=-1)),c(Y1,Y2[seq(from=length(VarY),to=1,by=-1)]),col="gray",lty=0)
      }else{
        lines(1:length(VarY),Y1,lty=2)
        lines(1:length(VarY),Y2,lty=2)
      }
      #min
      lines(1:length(VarY),apply(x$Y,2,min),lty=4)
      #moy
      lines(1:length(VarY),colMeans(x$Y),lty=1,lwd=2)
      #emplacement plot des indices
      par(fig=c(0,1,0,0.7),mar=c(5.1, 4.1, 0 ,2.1),new=TRUE,las=1)
    }

    plot(c(1,nb.outp),
         c(0,0),
         type="l",
         ylab="Indices",
         ylim=c(0,max(indice.cum,na.rm=TRUE)),
         cex=3, xaxt="n", xlab="",
         ...)
    if(xtick){
      ww=axTicks(1)
      ww[ww==0]=1
      axis(1,at=ww,labels=colnames(x$Y)[ww])
    }else{
      axis(1,at=NULL,labels=TRUE,tick=TRUE)
    }
    polygon(c(1,1:nb.outp,nb.outp),
            c(0,indice.cum[1,],0),
            density=densite[1],
            angle=100,
            col=colset[1],
            cex=2)
    text(n.name[1],
         indice.cum[1,n.name[1]]/2,
         names[1],
         cex=1,
         col="black",
         font=2)

    for (k in 2:nrow(tab.indicesmi)){
#    for (k in 2:nrow(indice.cum)){

        j <- k-1
        ## changed by HM, 2013-05-03
        if(shade){
          quat <-floor(j/4)
          if(quat*4-j==0){ang <-0}
          if(-quat*4+j==1){ang <-45}
          if(-quat*4+j==2){ang <-90}
          if(-quat*4+j==3){ang <-135}
        }
        else{ ang <- NULL }

        polygon( c(1:nb.outp, rev(1:nb.outp) ) ,
                c(indice.cum[k-1,],rev(indice.cum[k,])),
                density=densite[ 2+(j-1)%%6 ],angle=ang,
                col=colset[k], border="blue")
        text(n.name[k] ,
             (indice.cum[k-1,n.name[k]]+indice.cum[k,n.name[k]])/2,names[k],
             cex=1,col="black",font=2)
      }

#Avant ajout de la ligne residual dans indice.cum :
#    polygon( c(1:nb.outp, rev(1:nb.outp) ) ,
#            c(indice.cum[nrow(indice.cum),],rep(100,nb.outp)),
#            density=densite[8], angle=30)
#    text(n.name[nrow(tab.indicesmi)+1],
#         indice.cum[nrow(indice.cum),n.name[nrow(tab.indicesmi)+1]]/2+50,
#         "residual",
#         cex=1,col="black",font=2)

#Apres ajout de la ligne residual dans indice.cum :
# cad k = nrow(indice.cum) = nrow(tab.indicesmi)+1
    polygon( c(1:nb.outp, rev(1:nb.outp) ) ,
            c(indice.cum[nrow(indice.cum)-1,],rev(indice.cum[nrow(indice.cum),])),
            density=densite[8], angle=30)
    text(n.name[nrow(tab.indicesmi)+1],
         (indice.cum[nrow(indice.cum)-1,n.name[nrow(indice.cum)]]+indice.cum[nrow(indice.cum),n.name[nrow(indice.cum)]])/2,
         names[nrow(indice.cum)],#"residual",
         cex=1,col="black",font=2)

    # une ligne renforcee pour le 1 ou variance
    lines(1:nb.outp,indice.cum[nrow(indice.cum),1:nb.outp],lwd=3,col="black")

    # le plot est fini on reinitialise les paramètres
    par(old.par)

  if(total.plot){
    YLabel="Total Indices"
#    if(x$call.info$analysis=="sensitivity"){
#      if(x$call.info$fct=="morris"){
#        YLabel="sigma"
#      }
#    }
    total.indices=x$tSI[order(x$tSI[,jmax],decreasing=TRUE),];
    if(!normalized){
      for (j in 1:ncol(total.indices)) {
            total.indices[,j] <- total.indices[,j]* VarY[j]#/100 
        }
    }
    if(dev.cur() == 1 | dev.interactive()) dev.new()
    old.par <- par(no.readonly = TRUE)
    plot(1:nb.outp,total.indices[1,], type="l", ylab=YLabel, ylim=c(0,max(total.indices,max(indice.cum,na.rm=TRUE),na.rm=TRUE)), cex=3, xaxt="n", xlab="")
    for(k in 2:nrow(total.indices)){
      lines(1:nb.outp, total.indices[k,], col=k, lty=k)
    }
    # une ligne renforcee pour le 1 ou variance
    lines(1:nb.outp,indice.cum[nrow(indice.cum),1:nb.outp],lwd=3,col="black")
    if(xtick){
      axis(1,at=1:nb.outp,labels=colnames(x$Y))
    }else{
      axis(1,at=NULL,labels=TRUE,tick=TRUE)
    }

    legend("right",legend=rownames(total.indices),col=1:nrow(total.indices),lty=1:nrow(total.indices))
    par(old.par)

  }

  if(x$call.info$analysis=="sensitivity"){
    if(x$call.info$fct=="morris"){
      if(dev.cur() == 1 | dev.interactive()) dev.new()
      old.par <- par(no.readonly = TRUE)
      par(mfrow=c(2,1))
      # mu.star
      plot(1:nb.outp,x$mSI[1,], type="l", ylab="mu.star", ylim=c(0,max(x$mSI,na.rm=TRUE)), cex=3)
      for(k in 2:nrow(x$mSI)){
        lines(1:nb.outp, x$mSI[k,], col=k, lty=k)
      }
      legend("right",legend=rownames(x$mSI),col=1:nrow(x$mSI),lty=1:nrow(x$mSI))
      # sigma
      plot(1:nb.outp,x$tSI[1,]^2, type="l", ylab="sigma^2", ylim=c(0,max(x$tSI^2,na.rm=TRUE)), cex=3)
      for(k in 2:nrow(x$tSI)){
        lines(1:nb.outp, x$tSI[k,]^2, col=k, lty=k)
      }
      par(old.par)
    }#endif morris
  }#endif sensitivity

}

