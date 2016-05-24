### DPcaterpillar.R
### Caterpillar plots for random effects from a DPpackage object.
###
### Copyright: Alejandro Jara, 2006-2012.
###
### Last modification: 15-12-2006.
###
### This program is free software; you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation; either version 2 of the License, or (at
### your option) any later version.
###
### This program is distributed in the hope that it will be useful, but
### WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
### General Public License for more details.
###
### You should have received a copy of the GNU General Public License
### along with this program; if not, write to the Free Software
### Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
###
### The author's contact information:
###
###      Alejandro Jara
###      Department of Statistics
###      Facultad de Matematicas
###      Pontificia Universidad Catolica de Chile
###      Casilla 306, Correo 22 
###      Santiago
###      Chile
###      Voice: +56-2-3544506  URL  : http://www.mat.puc.cl/~ajara
###      Fax  : +56-2-3547729  Email: atjara@uc.cl
###

"DPcaterpillar"<-
function(object, midpoint="mean", hpd=TRUE , ask=TRUE, nfigr=1, nfigc=1, ...)
{
   hpdf<-function(x)
   {
        alpha<-0.05
        vec<-x
        n<-length(x)         
        alow<-rep(0,2)
        aupp<-rep(0,2)
        a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
                    alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
        return(c(a$alow[1],a$aupp[1]))
   }
   
   pdf<-function(x)
   {
        alpha<-0.05
        vec<-x
        n<-length(x)         
        alow<-rep(0,2)
        aupp<-rep(0,2)
        a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
                    alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
        return(c(a$alow[2],a$aupp[2]))
   }

   caterpillar.plot<-function (segs, main = "", xlab = "", vgrid = NA, mar = NA, 
          col = gray(0.6) , bg = "green", ...) 
   {
       if (missing(segs)) 
           stop("Usage: caterpillar.plot(segs,...)\n\twhere segs is a matrix of midpoints and limits")
       segdim <- dim(segs)
    
       left.labels <- colnames(segs)
       if (is.null(left.labels)){
          left.labels <- 1:segdim[2]
          colnames(segs)<-left.labels
       }   
       
       segs <- segs[, order(segs[1, ])]

       #oldpar <- par(no.readonly = TRUE)
       if (is.na(mar)) 
           mar = c(4, 6, 1 + 2 * (nchar(main) > 0), 5)
       par(mar = mar)
       plot(x = c(min(segs[2, ]), max(segs[3, ])), y = c(1, segdim[2]), 
           main = main, xlab = "", ylab = "", type = "n", axes = FALSE)
        
       box()
       if (!is.na(vgrid)) 
        abline(v = vgrid, lty = 2)
       axis(1)
    
       cexline<- 1- 0.0035*segdim[2]
       if(cexline<0.65)cexline<-0.65
    
       arrows(segs[2, ], 1:segdim[2], segs[3, ], 1:segdim[2], length = 0.005, 
           angle = 90, code = 3, col = gray(0.6) , lwd = 1.2 )
        
       points(segs[1, ], 1:segdim[2], pch = 21, col = gray(0.6) , bg = "green" , cex=cexline)

       cexleft<-1-0.009875*segdim[2]
       if(cexleft<0.0125)cexleft<-0.0125
       left.labels <- colnames(segs)
       plot.limits <- par("usr")
       mtext(left.labels, 2, line = 0.2, at = 1:segdim[2], adj = 1, 
             las = 1, cex = cexleft)
    
       if (nchar(xlab)) 
           mtext(xlab, 1, line = 2)
       #par(oldpar)
   }

   if(is(object, "DPrandom"))
   {

      nsubject <- object$nsubject
      nrandom <- object$nrandom
      type <- object$type
      centered <- object$centered
      nsave <- object$nsave
      randommat <- NULL

      idnames <- rownames(object$randomm)
      renames <- colnames(object$randomm)

      for(i in 1:nrandom)
      {
         for(j in 1:nsubject)
         {
              count<- (j-1)*nrandom + i 

              if(centered)
              {
                 vec<-matrix((object$randommat[,count]-object$thetamat[,i]),ncol=1)
              }
              else vec<-matrix(object$randommat[,count],ncol=1)
              colnames(vec)<- idnames[j]
              randommat<-cbind(randommat,vec)  
         }
      
      }

      if(hpd){             
         limm<-apply(randommat, 2, hpdf)
         coef.l<-limm[1,]
         coef.u<-limm[2,]
      }
      else
      {
         limm<-apply(randommat, 2, pdf)
         coef.l<-limm[1,]
         coef.u<-limm[2,]
      }
      
      if(midpoint=="mean")
      {
         coef.m <-apply(randommat, 2, mean)          
      }
      else
      {
         coef.m <-apply(randommat, 2, median)    
      }

      tables<-cbind(coef.m,coef.l,coef.u)

      start<-1
      end<-nsubject

      #oldpar <- par(no.readonly = TRUE)
      par(ask = ask)
      layout(matrix(seq(1,nfigr*nfigc,1), nrow=nfigr , ncol=nfigc ,byrow=TRUE))

      for(i in 1:nrandom)
      {
          work<-t(tables[start:end,])
          start<-end+1
          end<-start+nsubject-1
          caterpillar.plot(segs=work,xlab = renames[i]) 
      }
      #par(oldpar)     
   }
}

