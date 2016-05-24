### DProc.R                   
### Estimates the ROC curve using a Dirichlet process mixture of normals
### model for each group
###
### Copyright: Alejandro Jara, 2006-2012.
###
### Last modification: 22-12-2006.
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

"DProc"<-
function(x,y,fitx=NULL,fity=NULL,ngrid=1000,priorx,priory,
         mcmcx,mcmcy,statex,statey,
         statusx,statusy,data=sys.frame(sys.parent()),
         na.action=na.fail) 
UseMethod("DProc")

"DProc.default"<-
function(x,
         y,
         fitx=NULL,
         fity=NULL,
         ngrid=1000,
         priorx,
         priory,
         mcmcx,
         mcmcy,         
         statex,
         statey,
         statusx,
         statusy,         
         data=sys.frame(sys.parent()),
         na.action=na.fail)
{

         #########################################################################################
         # call parameters
         #########################################################################################
         cl <- match.call()
	  
         #########################################################################################
         # data structure
         #########################################################################################

         if(is.null(fitx))
         {
            x<-na.action(as.matrix(x))
            y<-na.action(as.matrix(y))

            n<-dim(x)[1]
            nvarx<-dim(x)[2]
            m<-dim(y)[1]
            nvary<-dim(y)[2]

            if(nvarx>1) stop("So far, this function is only for univariate ROC analysis")           
            if(nvary>1) stop("So far, this function is only for univariate ROC analysis")           

            varnamesx<-colnames(x)[1]
            if(is.null(varnamesx))
            {
               varnamesx<-all.vars(cl)[1]
            }
          
            varnamesy<-colnames(y)[1]
            if(is.null(varnamesy))
            {
               varnamesy<-all.vars(cl)[2]
            }
         }
         else
         {
            x<-fitx$y
            nvarx<-fitx$nvar
            n<-fitx$nrec
            varnamesx<-fitx$varnames

            y<-fity$y
            nvary<-fity$nvar
            m<-fity$nrec
            varnamesy<-fity$varnames

            if(nvarx>1) stop("So far, this function is only for univariate ROC analysis")           
            if(nvary>1) stop("So far, this function is only for univariate ROC analysis")           
         }
         
         #########################################################################################
         # calling the function
         #########################################################################################

         if(is.null(fitx))
         {
            cat("\n")
            cat(paste("Fitting the model for",varnamesx,sep=" "),"\n")
            fitx<-DPdensity(y=x,ngrid=ngrid,prior=priorx,mcmc=mcmcx,
                            state=statex,status=statusx)
     
            cat(paste("Fitting the model for",varnamesy,sep=" "),"\n")
            fity<-DPdensity(y=y,ngrid=ngrid,prior=priory,mcmc=mcmcy,
                            state=statey,status=statusy)
         }                   

         #########################################################################################
         # ROC analysis
         #########################################################################################
         cdf <- function(object,eval)
         {
            randommat <- NULL
            nrec<-object$nrec
            for(i in 1:2)
            {
               for(j in 1:(nrec+1))
               {
                    count<- (j-1)*2 + i 
                    vec<-matrix(object$save.state$randsave[,count],ncol=1)
                    randommat<-cbind(randommat,vec)  
               }
            }

            means <- randommat[,1:(nrec+1)]
            vars  <- randommat[,(nrec+2):(2*(nrec+1))]
            alpha <- object$save.state$thetasave[,5]
            nsave <- length(alpha)
            neval <- length(eval)
            F<-rep(0,neval)
            out2 <- cbind(matrix(rep(1/(alpha+nrec),nrec),ncol=nrec,byrow=FALSE),alpha/(alpha+nrec))
            for(i in 1:neval)
            {
                out1 <- pnorm(eval[i],mean=means,sd=sqrt(vars))
                F[i]<- mean(apply(out1*out2,1,sum))
            }
            return(F)
         }

         cdf2 <- function(object,eval)
         {
            nrec <- object$nrec
            samp <- object$save.state$randsave[,( 2*(nrec+1)+1 )] 
            neval <- length(eval)
            F<-rep(0,neval)             

            for(i in 1:neval)
            {
                F[i]<- mean(samp<=eval[i])
            }
            out<-list(F=F,samp=samp)
            return(out)
         }

         auc<-function(x,y)
         {
            n <- length(x)
            m <- length(y)
            z1 <- NULL
            z2 <- NULL

            x <- as.matrix(x)

            fun1 <- function(x)
            {
                 return(as.integer(x<y))
            }
            z1 <- apply(x,1,fun1)

            fun2 <- function(x)
            {
                 return(as.integer(x==y))
            }
            z2 <- apply(x,1,fun2)

            U1 <- sum(z1 + 0.5*z2)/(n*m)
            return(U1)
         }

         eval <- sort(c(x,y))
         eval <- sort(c(seq(min(eval),max(eval),length=ngrid),eval))
         eval <- as.double(names(table(eval)))
         
         anal1 <- cdf2(object=fitx,eval=eval)
         anal2 <- cdf2(object=fity,eval=eval)
         
         FX <- anal1$F 
         FY <- anal2$F 
         
         AUC <- auc(anal1$samp,anal2$samp)
         
         sens <- 1-FY
         spec <- FX
         P <- m / (m + n) 
         Q <- P * sens + (1-P) * (1-spec)
         k10 <- ( sens - Q ) / (1 - Q)
         k00 <- ( spec - (1-Q)) / (Q)
         kfin <- ( P*(1-Q)*k10 + (1-P)*Q*k00 ) / (P*(1-Q) + (1-P)*Q)
        
         #########################################################################################
         # save state
         #########################################################################################

         model.name<-"Bayesian semiparametric ROC curve estimation"		

         cuto<-c(eval[kfin==max(kfin)][1],1-FY[kfin==max(kfin)][1],FX[kfin==max(kfin)][1],AUC) 
         names(cuto)=c("Cutoff","Sensitivity","Specificity","AUC")

         z<-list(modelname=model.name,call=cl,n1=n,n2=m,
                 cuto=cuto,eval=eval,gridx=fitx$x1,gridy=fity$x1,
                 densx=fitx$dens,densy=fity$dens,cdfx=FX,cdfy=FY,
                 namex=varnamesx,namey=varnamesy,fitx=fitx,fity=fity)
         cat("\n")

         class(z)<-c("DProc")
	 return(z) 
}


###                    
### Tools
###
### Copyright: Alejandro Jara Vallejos, 2006
### Last modification: 22-12-2006.
###


"print.DProc"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\n",x$modelname,"\n\nCall:\n", sep = "")
    print(x$call)

    cat("\nPosterior Inference of Parameters:\n")
    print.default(format(x$cuto, digits = digits), print.gap = 2, 
            quote = FALSE)

    cat("\nNumber of Observations in x:",x$n1)
    cat("\nNumber of Observations in y:",x$n2)
    cat("\n\n")
    invisible(x)
}


"plot.DProc"<-function(x, ask=TRUE, param=NULL, hpd=TRUE, nfigr=1, nfigc=1, col="#bdfcc9", ...) 
{
   if(is(x, "DProc"))
   {

      par(ask = ask)
      layout(matrix(seq(1,nfigr*nfigc,1),nrow=nfigr,ncol=nfigc,byrow=TRUE))
      
      plot(x$gridx,x$densx,type="l",lwd=2,lty=1,
           xlim=c(min(c(x$gridx,x$gridy)),max(c(x$gridx,x$gridy))),
           ylim=c(0,max(c(x$densx,x$densy))),
           xlab="Values",ylab="density")
      lines(x$gridy,x$densy,lwd=2,lty=2)     
      abline(v=x$cuto[1],col="red",lwd=1)
      legend(x$eval[length(x$eval)*0.55], max(c(x$densx,x$densy)), lty=c(1,2),
             legend=c(x$namex,x$namey),bty="n",lwd=2)
      
      plot(x$eval,x$cdfx,type="l",lwd=2,lty=1,
           xlab=paste("Values in Group: ",x$namex),ylab="cdf",ylim=c(0,1))

      plot(x$eval,x$cdfy,type="l",lwd=2,lty=1,
           xlab=paste("Values in Group: ",x$namey),ylab="cdf",ylim=c(0,1))
      
      plot(1-x$cdfx,1-x$cdfy,typ="l",lwd=2,xlab="1-Specificity",ylab="Sensitivity",
           xlim=c(0,1),ylim=c(0,1))        
      points(1-x$cuto[3],x$cuto[2],col="red",pch=19,lwd=2) 
   }
}

