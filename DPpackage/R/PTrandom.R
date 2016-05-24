### PTrandom.R
### Extracts random effects from DPpackage objects PTlmm, PTolmm, and PTglmm.
###
### Copyright: Alejandro Jara, 2007-2012.
### Last modification: 04-05-2007.
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

"PTrandom"<-
function(object,centered=FALSE,predictive=FALSE,ngrid=1000,gridl=NULL)
UseMethod("PTrandom")

"PTrandom.default"<-
function(object,centered=FALSE,predictive=FALSE,ngrid=1000,gridl=NULL)
{
   comput<-0 
   if(is(object, "PTlmm"))
   {
      comput<-1
      random<-matrix(0,nrow=object$nsubject,ncol=object$nrandom)
      predp<-rep(0,object$nrandom)
      predsd<-rep(0,object$nrandom)
      predse<-rep(0,object$nrandom)
      predl<-rep(0,object$nrandom)
      predu<-rep(0,object$nrandom)
      predm<-rep(0,object$nrandom)
      
      randommat<-matrix(object$save.state$randsave,
                 nrow=object$mcmc$nsave,ncol=object$nrandom*(object$nsubject+1))

      dimnames(randommat)<-dimnames(object$save.state$randsave)
      
      thetamat<-matrix(object$save.state$thetasave,nrow=object$mcmc$nsave, 
                       ncol=object$dimen)
      
      counter<-0
      for(i in 1:object$nsubject){
          for(j in 1:object$nrandom){
              counter<-counter+1
              if(centered)
              {
                   random[i,j]<-mean(object$save.state$randsave[,counter]-
                                     object$save.state$thetasave[,j])
              }
              else
              {
                   random[i,j]<-mean(object$save.state$randsave[,counter])              
              }
          }
      }
      
      type<-1
      dimnames(random)<-list(object$namesre1,object$namesre2)
      z<-list(randomm=random,randommat=randommat,thetamat=thetamat,centered=centered,
              predictive=predictive,nsubject=object$nsubject,nrandom=object$nrandom,
              modelname=object$modelname,call=object$call,type=type,nsave=object$mcmc$nsave)
      
      if(predictive==TRUE)
      {
      	 for(i in 1:object$nrandom)
      	 { 		
      	     counter<-counter+1	
      	     if(centered)
      	     {
                predp[i]<-mean(object$save.state$randsave[,counter]-
                                object$save.state$thetasave[,i])      	     	

                predm[i]<-median(object$save.state$randsave[,counter]-
                                object$save.state$thetasave[,i])      	     	

                predsd[i]<-sqrt(var(object$save.state$randsave[,counter]-
                                object$save.state$thetasave[,i]))      	     	

                vec<-object$save.state$randsave[,counter]-object$save.state$thetasave[,i]
                
                n<-length(vec)
                
                alpha<-0.05
                
                alow<-rep(0,2)
                
                aupp<-rep(0,2)
                
       
                a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
                            alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
                predl[i]<-a$alow[1]            
                predu[i]<-a$aupp[1]
                
                predse[i]<-predsd[i]/sqrt(n)
      	     }
      	     else
      	     {
                predp[i]<-mean(object$save.state$randsave[,counter])      	     	

                predm[i]<-median(object$save.state$randsave[,counter])      	     	

                predsd[i]<-sqrt(var(object$save.state$randsave[,counter]))      	     	

                vec<-object$save.state$randsave[,counter]
                
                n<-length(vec)
                
                alpha<-0.05
                
                alow<-rep(0,2)
                
                aupp<-rep(0,2)
                
       
                a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
                            alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
                predl[i]<-a$alow[1]            
                predu[i]<-a$aupp[1]
                
                predse[i]<-predsd[i]/sqrt(n)
      	     }
      	 }
      	 
      	 predtable <- cbind(predp, predm, predsd, predse , predl , predu)
         dimnames(predtable) <- list(object$namesre2, c("Mean", "Median", "Std. Dev.", "Naive Std.Error",
                "95%HPD-Low","95%HPD-Upp"))
      	 
      	 z$prediction<-predtable
      }
   }

   if(is(object, "PTglmm"))
   {
      comput<-1   
      random<-matrix(0,nrow=object$nsubject,ncol=object$nrandom)
      predp<-rep(0,object$nrandom)
      predsd<-rep(0,object$nrandom)
      predse<-rep(0,object$nrandom)
      predl<-rep(0,object$nrandom)
      predu<-rep(0,object$nrandom)
      predm<-rep(0,object$nrandom)
      
      randommat<-matrix(object$save.state$randsave,
                 nrow=object$mcmc$nsave,ncol=object$nrandom*(object$nsubject+1))
      
      dimnames(randommat)<-dimnames(object$save.state$randsave)
      
      thetamat<-matrix(object$save.state$thetasave,nrow=object$mcmc$nsave, 
                       ncol=object$dimen)
      
      counter<-0
      for(i in 1:object$nsubject){
          for(j in 1:object$nrandom){
              counter<-counter+1
              if(centered)
              {
                   random[i,j]<-mean(object$save.state$randsave[,counter]-
                                     object$save.state$thetasave[,j])
              }
              else
              {
                   random[i,j]<-mean(object$save.state$randsave[,counter])              
              }
          }
      }
      
      type<-1
      dimnames(random)<-list(object$namesre1,object$namesre2)
      z<-list(randomm=random,randommat=randommat,thetamat=thetamat,centered=centered,
              predictive=predictive,nsubject=object$nsubject,nrandom=object$nrandom,
              modelname=object$modelname,call=object$call,type=type,nsave=object$mcmc$nsave)
      
      if(predictive==TRUE)
      {
      	 for(i in 1:object$nrandom)
      	 { 		
      	     counter<-counter+1	
      	     if(centered)
      	     {
                predp[i]<-mean(object$save.state$randsave[,counter]-
                                object$save.state$thetasave[,i])      	     	

                predm[i]<-median(object$save.state$randsave[,counter]-
                                object$save.state$thetasave[,i])      	     	

                predsd[i]<-sqrt(var(object$save.state$randsave[,counter]-
                                object$save.state$thetasave[,i]))      	     	

                vec<-object$save.state$randsave[,counter]-object$save.state$thetasave[,i]
                
                n<-length(vec)
                
                alpha<-0.05
                
                alow<-rep(0,2)
                
                aupp<-rep(0,2)
                
       
                a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
                            alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
                predl[i]<-a$alow[1]            
                predu[i]<-a$aupp[1]
                
                predse[i]<-predsd[i]/sqrt(n)
      	     }
      	     else
      	     {
                predp[i]<-mean(object$save.state$randsave[,counter])      	     	

                predm[i]<-median(object$save.state$randsave[,counter])      	     	

                predsd[i]<-sqrt(var(object$save.state$randsave[,counter]))      	     	

                vec<-object$save.state$randsave[,counter]
                
                n<-length(vec)
                
                alpha<-0.05
                
                alow<-rep(0,2)
                
                aupp<-rep(0,2)
                
       
                a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
                            alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
                predl[i]<-a$alow[1]            
                predu[i]<-a$aupp[1]
                
                predse[i]<-predsd[i]/sqrt(n)
      	     }
      	 }
      	 
      	 predtable <- cbind(predp, predm, predsd, predse , predl , predu)
         dimnames(predtable) <- list(object$namesre2, c("Mean", "Median", "Std. Dev.", "Naive Std.Error",
                "95%HPD-Low","95%HPD-Upp"))
      	 
      	 z$prediction<-predtable
      }
   }

   if(is(object, "PTolmm"))
   {
      comput<-1   
      random<-matrix(0,nrow=object$nsubject,ncol=object$nrandom)
      predp<-rep(0,object$nrandom)
      predsd<-rep(0,object$nrandom)
      predse<-rep(0,object$nrandom)
      predl<-rep(0,object$nrandom)
      predu<-rep(0,object$nrandom)
      predm<-rep(0,object$nrandom)
      
      randommat<-matrix(object$save.state$randsave,
                 nrow=object$mcmc$nsave,ncol=object$nrandom*(object$nsubject+1))
      
      dimnames(randommat)<-dimnames(object$save.state$randsave)
      
      thetamat<-matrix(object$save.state$thetasave,nrow=object$mcmc$nsave, 
                       ncol=object$dimen)
      
      counter<-0
      for(i in 1:object$nsubject){
          for(j in 1:object$nrandom){
              counter<-counter+1
              if(centered)
              {
                   random[i,j]<-mean(object$save.state$randsave[,counter]-
                                     object$save.state$thetasave[,j])
              }
              else
              {
                   random[i,j]<-mean(object$save.state$randsave[,counter])              
              }
          }
      }

      type<-1      
      dimnames(random)<-list(object$namesre1,object$namesre2)
      z<-list(randomm=random,randommat=randommat,thetamat=thetamat,centered=centered,
              predictive=predictive,nsubject=object$nsubject,nrandom=object$nrandom,
              modelname=object$modelname,call=object$call,type=type,nsave=object$mcmc$nsave)
      
      if(predictive==TRUE)
      {
      	 for(i in 1:object$nrandom)
      	 { 		
      	     counter<-counter+1	
      	     if(centered)
      	     {
                predp[i]<-mean(object$save.state$randsave[,counter]-
                                object$save.state$thetasave[,i])      	     	

                predm[i]<-median(object$save.state$randsave[,counter]-
                                object$save.state$thetasave[,i])      	     	

                predsd[i]<-sqrt(var(object$save.state$randsave[,counter]-
                                object$save.state$thetasave[,i]))      	     	

                vec<-object$save.state$randsave[,counter]-object$save.state$thetasave[,i]
                
                n<-length(vec)
                
                alpha<-0.05
                
                alow<-rep(0,2)
                
                aupp<-rep(0,2)
                
       
                a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
                            alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
                predl[i]<-a$alow[1]            
                predu[i]<-a$aupp[1]
                
                predse[i]<-predsd[i]/sqrt(n)
      	     }
      	     else
      	     {
                predp[i]<-mean(object$save.state$randsave[,counter])      	     	

                predm[i]<-median(object$save.state$randsave[,counter])      	     	

                predsd[i]<-sqrt(var(object$save.state$randsave[,counter]))      	     	

                vec<-object$save.state$randsave[,counter]
                
                n<-length(vec)
                
                alpha<-0.05
                
                alow<-rep(0,2)
                
                aupp<-rep(0,2)
                
       
                a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
                            alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
                predl[i]<-a$alow[1]            
                predu[i]<-a$aupp[1]
                
                predse[i]<-predsd[i]/sqrt(n)
      	     }
      	 }
      	 
      	 predtable <- cbind(predp, predm, predsd, predse , predl , predu)
         dimnames(predtable) <- list(object$namesre2, c("Mean", "Median", "Std. Dev.", "Naive Std.Error",
                "95%HPD-Low","95%HPD-Upp"))
      	 
      	 z$prediction<-predtable
      }
   }
   

   if(comput==1 && predictive==TRUE)
   {
   
      q<-object$nrandom
      frstlprob<-object$frstlprob

      if(q<=2)
      {
         m<-object$m
         nsubject<-object$nsubject 
         nfixed<-object$nfixed

         dimen1<-object$nrandom+object$nfixed
         dimen2<-0
         if(is(object, "PTlmm"))dimen2<-1
         if(is(object, "PTglmm"))dimen2<-object$dispp
         total<-dimen1+dimen2
         dimen3<-object$nrandom
         mumat<-object$save.state$thetasave[,(total+1):(total+dimen3)]
         total<-dimen1+dimen2+dimen3
         dimen4<-object$nrandom*(object$nrandom+1)/2
         sigmamat<-object$save.state$thetasave[,(total+1):(total+dimen4)]             
         cparvec<-object$save.state$thetasave[,(total+dimen4+1)]
         total<-total+dimen4+1
         dimen5<-object$nrandom*object$nrandom
         typepmat<-object$save.state$thetasave[,(total+1):(total+dimen5)]

         nsave<-length(cparvec)
         randsave<-object$save.state$randsave

         b<-matrix(0,nrow=nsubject,ncol=q)
         bz<-matrix(0,nrow=nsubject,ncol=q)
         iflagr<-rep(0,q) 
         linf<-rep(0,q)
         lsup<-rep(0,q)
         mu<-rep(0,q)
         parti<-rep(0,q)
         sigma<-matrix(0,nrow=q,ncol=q)
         sigmainv<-matrix(0,nrow=q,ncol=q)
         theta<-rep(0,q)
         thetaz<-rep(0,q)
         whicho<-rep(0,nsubject)
         whichn<-rep(0,nsubject)      
         workmhr<-rep(0,q*(q+1)/2) 
         workmr<-matrix(0,nrow=q,ncol=q) 
         workmr1<-matrix(0,nrow=q,ncol=q) 
         workmr2<-matrix(0,nrow=q,ncol=q)
         ortho<-matrix(0,nrow=q,ncol=q)  
         workvr<-rep(0,q)
         workvr1<-rep(0,q)

         if(q==1)
         {
            fs<-rep(0,ngrid)
            left<-min(z$randomm)-2*sqrt(var(z$randomm))
	    right<-max(z$randomm)+2*sqrt(var(z$randomm))
	    
	    if(is.null(gridl))
	    {
   	       grid<-seq(left,right,length=ngrid)  
   	    }
   	    else
   	    {

               gridl<-matrix(gridl,nrow=2)
               if(ncol(gridl)!=1)
               {
                  stop("Incorrect number of grid points limits")
               }   
   	       grid<-seq(gridl[1],gridl[2],length=ngrid)  
   	    }

            typep <- 1
            foo <- .Fortran("predictiveptu",
             	m          =as.integer(m),
   		nsubject   =as.integer(nsubject),
   		q          =as.integer(q),
   		nsave      =as.integer(nsave),
  		randsave   =as.double(randsave),
  		mumat      =as.double(mumat),
  		sigmamat   =as.double(sigmamat),
  		cparvec    =as.double(cparvec),
  		typep      =as.integer(typep),
  		ngrid      =as.integer(ngrid),
  		grid       =as.double(grid),
  		fs         =as.double(fs),	 
 	        iflagr     =as.integer(iflagr),
 	        parti      =as.integer(parti),
 	        whicho     =as.integer(whicho),
 	        whichn     =as.integer(whichn),
 	        b          =as.double(b),
 	        bz         =as.double(bz),
        	linf       =as.double(linf),
 	        lsup       =as.double(lsup),
 	        mu         =as.double(mu),
 	        sigma      =as.double(sigma),
 	        sigmainv   =as.double(sigmainv),
 	        theta      =as.double(theta),
 	        thetaz     =as.double(thetaz),
         	workmr     =as.double(workmr),
         	workmr1    =as.double(workmr1),
    		workmr2    =as.double(workmr2),
 	        workmhr    =as.double(workmhr),
     		workvr     =as.double(workvr),
     		workvr1    =as.double(workvr1),
     		fixed      =as.integer(frstlprob),
	        PACKAGE    ="DPpackage")
	
	    z$dens<-foo$fs
	    z$grid1<-foo$grid
         }
         
         if(q==2)
         {
            ngrid1<-as.integer(sqrt(ngrid)) 
            ngrid2<-ngrid1

	    if(is.null(gridl))
	    {
               left<-rep(0,2)
               right<-rep(0,2)
               left[1]<-min(z$randomm[,1])-2*sqrt(var(z$randomm[,1]))
               left[2]<-min(z$randomm[,2])-2*sqrt(var(z$randomm[,2]))
               right[1]<-max(z$randomm[,1])+2*sqrt(var(z$randomm[,1]))
               right[2]<-max(z$randomm[,2])+2*sqrt(var(z$randomm[,2]))
               grid1<-seq(left[1],right[1],length=ngrid1)  
               grid2<-seq(left[2],right[2],length=ngrid1)  
   	    }
   	    else
   	    {
               gridl<-matrix(gridl,nrow=2)
               if(ncol(gridl)!=2)
               {
                  stop("Incorrect number of grid points limits")
               }   

   	       grid1<-seq(gridl[1,1],gridl[2,1],length=ngrid1)  
   	       grid2<-seq(gridl[1,2],gridl[2,2],length=ngrid2)  
   	    }

            fs<-matrix(0,nrow=ngrid1,ncol=ngrid2)


            foo <- .Fortran("predictiveptb2",
             	m          =as.integer(m),
   		nsubject   =as.integer(nsubject),
   		q          =as.integer(q),
   		nsave      =as.integer(nsave),
  		randsave   =as.double(randsave),
  		mumat      =as.double(mumat),
  		sigmamat   =as.double(sigmamat),
  		cparvec    =as.double(cparvec),
  		typepmat   =as.double(typepmat),
  		ngrid1     =as.integer(ngrid1),
  		ngrid2     =as.integer(ngrid2),
  		grid1      =as.double(grid1),
  		grid2      =as.double(grid2),
  		fs         =as.double(fs),	 
 	        iflagr     =as.integer(iflagr),
 	        parti      =as.integer(parti),
 	        whicho     =as.integer(whicho),
 	        whichn     =as.integer(whichn),
 	        b          =as.double(b),
 	        bz         =as.double(bz),
        	linf       =as.double(linf),
 	        lsup       =as.double(lsup),
 	        mu         =as.double(mu),
 	        sigma      =as.double(sigma),
 	        sigmainv   =as.double(sigmainv),
 	        theta      =as.double(theta),
 	        thetaz     =as.double(thetaz),
         	workmr     =as.double(workmr),
         	ortho      =as.double(ortho),
 	        workmhr    =as.double(workmhr),
     		fixed      =as.integer(frstlprob),     		
	        PACKAGE    ="DPpackage")

            dens<-matrix(foo$fs,nrow=ngrid1,ncol=ngrid2)  
            dist2<-foo$grid1[2]-foo$grid1[1] 
            dist1<-foo$grid2[2]-foo$grid2[1] 
            dens1<- (dist1/2)*(dens[,1]+dens[,ngrid2]+2*apply(dens[,2:(ngrid2-1)],1,sum))
            dens2<- (dist2/2)*(dens[1,]+dens[ngrid1,]+2*apply(dens[2:(ngrid1-1),],2,sum))
             
	    z$dens<-matrix(foo$fs,nrow=ngrid1,ncol=ngrid2)
	    z$grid1<-foo$grid1
	    z$grid2<-foo$grid2
	    z$f1<-dens1
	    z$f2<-dens2
         }
      }
   }
   class(z)<-c("PTrandom") 
   return(z)
}



###
### Tools for PTrandom: print, plot
###
### Copyright: Alejandro Jara, 2007
### Last modification: 02-04-2007.


"print.PTrandom"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\n","Random effect information for the DP object:","\n\nCall:\n",sep = "")
    print(x$call)
    cat("\n")

    if(x$predictive)
    {
        cat("\nPredictive distribution:\n")
        print.default(format(x$prediction, digits = digits), print.gap = 2, 
            quote = FALSE)
    }
    else
    {
        cat("\nPosterior mean of subject-specific components:\n\n")
        print.default(format(x$randomm, digits = digits), print.gap = 2, 
            quote = FALSE)
    }        
    cat("\n\n")
    
    invisible(x)
}


"plot.PTrandom"<-function(x, ask=TRUE, hpd=TRUE, nfigr=2, nfigc=2, subject=NULL, col="#bdfcc9", ...) 
{

fancydensplot<-function(x, hpd=TRUE, npts=200, 
                        xlim=NULL,xlab="", ylab="", main="",col="#bdfcc9", ...)
# Author: AJV, 2006
#
{
	dens <- density(x,n=npts)
	densx <- dens$x
	densy <- dens$y

	meanvar <- mean(x)
	densx1 <- max(densx[densx<=meanvar])
	densx2 <- min(densx[densx>=meanvar])
        densy1 <- densy[densx==densx1]
        densy2 <- densy[densx==densx2]
        ymean <- densy1 + ((densy2-densy1)/(densx2-densx1))*(meanvar-densx1)
        

        if(hpd==TRUE)
	{
		alpha<-0.05
		alow<-rep(0,2)
        	aupp<-rep(0,2)
        	n<-length(x)
		a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(x),
		                     alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
		xlinf<-a$alow[1]            
		xlsup<-a$aupp[1]            
	}
	else
	{
		xlinf <- quantile(x,0.025)
		xlsup <- quantile(x,0.975)
	}

	densx1 <- max(densx[densx<=xlinf])
	densx2 <- min(densx[densx>=xlinf])
	densy1 <- densy[densx==densx1]
	densy2 <- densy[densx==densx2]
	ylinf <- densy1 + ((densy2-densy1)/(densx2-densx1))*(xlinf-densx1)

	densx1 <- max(densx[densx<=xlsup])
	densx2 <- min(densx[densx>=xlsup])
        densy1 <- densy[densx==densx1]
        densy2 <- densy[densx==densx2]
        ylsup <- densy1 + ((densy2-densy1)/(densx2-densx1))*(xlsup-densx1)

        if(is.null(xlim))
        {
           xlim<-c(min(densx), max(densx)) 
        }


        plot(0.,0.,xlim = xlim, ylim = c(min(densy), max(densy)),
                axes = F,type = "n" , xlab=xlab, ylab=ylab, main=main, cex=1.2)

        
        xpol<-c(xlinf,xlinf,densx[densx>=xlinf & densx <=xlsup],xlsup,xlsup)
        ypol<-c(0,ylinf,densy[densx>=xlinf & densx <=xlsup] ,ylsup,0)
             
        polygon(xpol, ypol, border = FALSE,col=col)
        
        lines(xlim,c(0,0),lwd=1.2)
        
        segments(xlim[1],0, xlim[1],max(densy),lwd=1.2)
        
        lines(densx,densy,lwd=1.2)
             
        segments(meanvar, 0, meanvar, ymean,lwd=1.2)
        segments(xlinf, 0, xlinf, ylinf,lwd=1.2)
        segments(xlsup, 0, xlsup, ylsup,lwd=1.2)

	axis(1., at = round(c(xlinf, meanvar,xlsup), 2.), labels = T,pos = 0.)
        axis(1., at = round(seq(xlim[1],xlim[2],length=15), 2.), labels = F,pos = 0.)
        axis(2., at = round(seq(0,max(densy),length=5), 2.), labels = T,pos =xlim[1])
}

   if(is(x, "PTrandom")){

       oldpar <- par(no.readonly = TRUE)
       par(ask = ask)

       if(x$predictive)
       {
          coef.p<-x$randomm
          n<-dim(coef.p)[2]
          pnames<-colnames(coef.p)
          start<-x$nsubject*n
          
          if(n==1)
          {
             nfigr<-1
             nfigc<-1          
             layout(matrix(seq(1,nfigr*nfigc,1), nrow=nfigr , ncol=nfigc ,byrow=TRUE))
          
             xx<-matrix(x$grid1,ncol=1)          
             dens<-x$dens
             title<-paste("Density of",pnames[1],sep=" ")
             plot(0.,0.,xlim = c(min(xx),max(xx)), ylim = c(0, (max(dens)+0.01)),
                axes = F,type = "n" , xlab="values", ylab="density", main=title, cex=1.2)
             lines(x$grid1,x$dens,lwd=2)
             lines(c(min(xx),max(xx)),c(0,0),lwd=1.2)                
             segments(min(xx),0, min(xx),max(dens)+0.01,lwd=1.2)
             axis(1., at = round(seq(min(xx),max(xx),length=15), 2.), labels = T,pos = 0.)
             axis(2., at = round(seq(0,max(dens)+0.01,length=5), 2.), labels = T,pos =min(xx))
             for(j in 1:x$nsubject)
             {
                 points(x$randomm[j,1],0,col="red",pch=20)              
             }
             
          }
          
          if(n==2)
           {
             layout(matrix(seq(1,nfigr*nfigc,1), nrow=nfigr , ncol=nfigc ,byrow=TRUE))                      

             xx<-matrix(x$grid1,ncol=1)
             yy<-matrix(x$grid2,ncol=1)
             z<-x$dens
             ngrid<-length(xx)
             
             dens1<-x$f1
             dens2<-x$f2
             
             colnames(xx)<-pnames[1]
             colnames(yy)<-pnames[2]

             title<-paste("Density of",pnames[1],sep=" ")
             plot(0.,0.,xlim = c(min(xx),max(xx)), ylim = c(0, (max(dens1)+0.01)),
                axes = F,type = "n" , xlab="values", ylab="density", main=title, cex=1.2)
             lines(xx,dens1,lwd=2)
             lines(c(min(xx),max(xx)),c(0,0),lwd=1.2)                
             segments(min(xx),0, min(xx),max(dens1)+0.01,lwd=1.2)
             axis(1., at = round(seq(min(xx),max(xx),length=15), 2.), labels = T,pos = 0.)
             axis(2., at = round(seq(0,max(dens1)+0.01,length=5), 2.), labels = T,pos =min(xx))
             
             for(j in 1:x$nsubject)
             {
                 points(x$randomm[j,1],0,col="red",pch=20)              
             }

             title<-paste("Density of",pnames[2],sep=" ")
             plot(0.,0.,xlim = c(min(yy),max(yy)), ylim = c(0, (max(dens2)+0.01)),
                axes = F,type = "n" , xlab="values", ylab="density", main=title, cex=1.2)
             lines(yy,dens2,lwd=2)             
             lines(c(min(yy),max(yy)),c(0,0),lwd=1.2)                
             segments(min(yy),0, min(yy),max(dens2)+0.01,lwd=1.2)
             axis(1., at = round(seq(min(yy),max(yy),length=15), 2.), labels = T,pos = 0.)
             axis(2., at = round(seq(0,max(dens2)+0.01,length=5), 2.), labels = T,pos =min(yy))
             
             for(j in 1:x$nsubject)
             {
                 points(x$randomm[j,2],0,col="red",pch=20)
             }    

             title<-paste("Density of",pnames[1],sep=" ")
             title<-paste(title,pnames[2],sep="-")
             contour(xx,yy,z,xlab=pnames[1],ylab=pnames[2],main=title)
             points(x$randomm[,1],x$randomm[,2],col="black",pch=20)              
             persp(xx,yy,z,xlab=pnames[1],ylab=pnames[2],zlab="density",theta=-30,phi=15,expand = 0.9, ltheta = 120,main=title)

          }
          if(n>2)
          {
             layout(matrix(seq(1,nfigr*nfigc,1), nrow=nfigr , ncol=nfigc ,byrow=TRUE))                    
             for(i in 1:n)
             {
                 title1<-paste("Trace of",pnames[i],sep=" ")
                 title2<-paste("Density of",pnames[i],sep=" ")
                 
                 if(x$centered)
                 {
                   vec<-(x$randommat[,start+i]-x$thetamat[,i])
                 }
                 else vec<-x$randommat[,start+i]

                 names(vec)<-pnames[i]
                
                 vectmp<-x$randomm[,i]
                 xlim<-c(min(vectmp)-6.5*sqrt(var(vectmp)),max(vectmp)+6.5*sqrt(var(vectmp)))

                 plot(vec,type='l',main=title1,xlab="MCMC scan",ylab=" ")
                 fancydensplot(vec,xlim=xlim,hpd=hpd,main=title2,xlab="values", ylab="density", col=col)
                 for(j in 1:x$nsubject)
                 {
                  points(x$randomm[j,i],0,col="red",pch=20)              
                 }
            }     
          }
       }
       
       else
       {
          layout(matrix(seq(1,nfigr*nfigc,1), nrow=nfigr , ncol=nfigc ,byrow=TRUE))
          coef.p<-x$randommat
          n<-dim(coef.p)[2]
          pnames<-colnames(coef.p)
          count<-0
          for(i in 1:x$nsubject)
          {
              for(j in 1:x$nrandom)
              {
                  count<-count+1
                  title1<-paste("Trace of",pnames[count],sep=" ")
                  title2<-paste("Density of",pnames[count],sep=" ")
                  if(x$centered)
                  {
                     vec<-(x$randommat[,count]-x$thetamat[,j])
                  }
                  else vec<-x$randommat[,count]
                  plot(vec,type='l',main=title1,xlab="MCMC scan",ylab=" ")
                  fancydensplot(vec,hpd=hpd,main=title2,xlab="values", ylab="density", col=col)
              }
          }
       }
       par(oldpar)  
   }
}


