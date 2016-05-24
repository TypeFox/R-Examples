textvarsup <- function(resmca,var,sel=1:nlevels(var),axes=c(1,2),col='black',app=0,vname=NULL) {
   vs <- varsup(resmca,var)
   levs <- names(vs$weight) %in% levels(var)[sel]
   xy <- vs$coord[levs,axes]
   texte <- names(vs$weight)[levs]
   if(!(is.null(vname))) texte <- paste(vname,texte,sep='.')
   prop <- round(vs$weight/sum(vs$weight)*2+0.5,1)[levs]
   if(app==0) text(xy,texte,col=col,cex=1)
   if(app==1) text(xy,texte,col=col,cex=prop)
   if(app==2) {
      points(xy,pch=17,cex=prop)
      text(xy,texte,pos=3,col=col)
      }
   }
