### DPrandom.R
### Extracts random effects from a DPpackage object.
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

"DPrandom"<-
function(object,centered=FALSE,predictive=FALSE)
UseMethod("DPrandom")

"DPrandom.default"<-
function(object,centered=FALSE,predictive=FALSE)
{
   if(is(object, "DPlmm"))
   {
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
      class(z)<-c("DPrandom") 
      return(z)
     
   }

   if(is(object, "PTlmm"))
   {
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
      class(z)<-c("DPrandom") 
      return(z)
     
   }



   if(is(object, "DPglmm")){
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
      class(z)<-c("DPrandom") 
      return(z)
     
   }


   if(is(object, "PTglmm")){
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
      class(z)<-c("DPrandom") 
      return(z)
     
   }



   if(is(object, "DPolmm"))
   {
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
      class(z)<-c("DPrandom") 
      return(z)
     
   }


   if(is(object, "PTolmm"))
   {
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
      class(z)<-c("DPrandom") 
      return(z)
     
   }


   if(is(object, "DPdensity") || is(object, "DPMdencens"))
   {

       if (centered) 
       { 
	   stop("This option is not implemented for DPdensity or DPMdencens.\n")
       }


       counter<-0 
       dimen<-object$nvar+object$nvar*(object$nvar+1)/2
       random<-matrix(0,nrow=object$nrec,ncol=dimen)
       randommat=object$save.state$randsave

       for(i in 1:object$nrec)
       {
           for(j in 1:dimen)
           {
                counter<-counter+1
                random[i,j]<-mean(object$save.state$randsave[,counter])       
           }
       }
       
       if(is.null(object$varnames))
       {
          object$varnames<-all.vars(object$call)[1:object$nvar]
       }
       
       namesre<-NULL
       for(i in 1:object$nvar)
       {
           namesre<-c(namesre,paste("mu",object$varnames[i],sep="-"))
       }

       for(i in 1:object$nvar)
       {
           for(j in i:object$nvar)
           {
               if(i==j)namesre<-c(namesre,paste("var",object$varnames[i],sep="-"))
               if(i!=j)
               {
                   tempname<-paste(object$varnames[i],object$varnames[j],sep="-")
                   namesre<-c(namesre,paste("sigma",tempname,sep="-"))
               }    
           }
       }
       
       predtable=NULL
       
       if(predictive==TRUE)
       {       
           start<-object$nrec*(object$nvar+object$nvar*(object$nvar+1)/2)
           predp<-rep(0,dimen)
           predm<-rep(0,dimen)
           predsd<-rep(0,dimen)
           predse<-rep(0,dimen)
           predl<-rep(0,dimen)
           predu<-rep(0,dimen)
           for(i in 1:dimen)
           {
                predp[i]<-mean(object$save.state$randsave[,(start+i)])      	     	
                predm[i]<-median(object$save.state$randsave[,(start+i)])      	     	
                predsd[i]<-sqrt(var(object$save.state$randsave[,(start+i)]))      	     	
                vec<-object$save.state$randsave[,(start+i)]
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
           predtable <- cbind(predp, predm, predsd, predse , predl , predu)
           dimnames(predtable) <- list(namesre, c("Mean", "Median", "Std. Dev.", "Naive Std.Error",
                "95%HPD-Low","95%HPD-Upp"))
       }

       dimnames(random)<-list(seq(1,object$nrec),namesre)

       type<-2       
       z<-list(randomm=random,modelname=object$modelname,call=object$call,predictive=predictive,
               nsubject=object$nrec,nrandom=dimen,centered=centered,randommat=randommat,
               prediction=predtable,type=type,nsave=object$mcmc$nsave)
               
       class(z)<-c("DPrandom") 
       return(z)
   }


   if(is(object, "DPraschpoisson") || is(object, "DPMraschpoisson"))
   {

       if (centered) 
       { 
	   stop("This option is not implemented for DPraschpoisson.\n")
       }


       random<-matrix(0,nrow=object$nsubject,ncol=1)
       
       randommat<-matrix(object$save.state$randsave,
                  nrow=object$mcmc$nsave,ncol=object$nsubject+1)
      
       dimnames(randommat)<-dimnames(object$save.state$randsave)

       for(i in 1:object$nsubject){
           random[i]<-mean(object$save.state$randsave[,i])              
       }
      
       colnames(random)<-"theta"

       type<-2      
       z<-list(randomm=random,modelname=object$modelname,call=object$call,predictive=predictive,
               nsubject=object$nsubject,nrandom=1,centered=FALSE,randommat=randommat,
               type=type,nsave=object$mcmc$nsave)

       if(predictive==TRUE)
       {
          predp<-mean(object$save.state$randsave[,object$nsubject+1])      	     	

          predm<-median(object$save.state$randsave[,object$nsubject+1])      	     	

          predsd<-sqrt(var(object$save.state$randsave[,object$nsubject+1]))      	     	

          vec<-object$save.state$randsave[,object$nsubject+1]
          
          n<-length(vec)
          
          alpha<-0.05
          
          alow<-rep(0,2)
          
          aupp<-rep(0,2)
          
       
          a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
                      alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
          predl<-a$alow[1]            
          predu<-a$aupp[1]
          
          predse<-predsd/sqrt(n)
     	 
      	  predtable <- cbind(predp, predm, predsd, predse , predl , predu)
          dimnames(predtable) <- list("theta", c("Mean", "Median", "Std. Dev.", "Naive Std.Error",
                "95%HPD-Low","95%HPD-Upp"))
      	 
      	  z$prediction<-predtable
       }

       class(z)<-c("DPrandom") 
       return(z)
   }


   if(is(object, "FPTraschpoisson"))
   {

       if (centered) 
       { 
	   stop("This option is not implemented for FPTraschpoisson.\n")
       }

       random<-matrix(0,nrow=object$nsubject,ncol=1)
       
       randommat<-matrix(object$save.state$randsave,
                  nrow=object$mcmc$nsave,ncol=object$nsubject+1)
      
       dimnames(randommat)<-dimnames(object$save.state$randsave)

       for(i in 1:object$nsubject){
           random[i]<-mean(object$save.state$randsave[,i])              
       }
      
       colnames(random)<-"theta"

       type<-2      
       z<-list(randomm=random,modelname=object$modelname,call=object$call,predictive=predictive,
               nsubject=object$nsubject,nrandom=1,centered=FALSE,randommat=randommat,
               type=type,nsave=object$mcmc$nsave)

       if(predictive==TRUE)
       {
          predp<-mean(object$save.state$randsave[,object$nsubject+1])      	     	

          predm<-median(object$save.state$randsave[,object$nsubject+1])      	     	

          predsd<-sqrt(var(object$save.state$randsave[,object$nsubject+1]))      	     	

          vec<-object$save.state$randsave[,object$nsubject+1]
          
          n<-length(vec)
          
          alpha<-0.05
          
          alow<-rep(0,2)
          
          aupp<-rep(0,2)
          
       
          a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
                      alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
          predl<-a$alow[1]            
          predu<-a$aupp[1]
          
          predse<-predsd/sqrt(n)
     	 
      	  predtable <- cbind(predp, predm, predsd, predse , predl , predu)
          dimnames(predtable) <- list("theta", c("Mean", "Median", "Std. Dev.", "Naive Std.Error",
                "95%HPD-Low","95%HPD-Upp"))
      	 
      	  z$prediction<-predtable
       }

       class(z)<-c("DPrandom") 
       return(z)
   }


   if(is(object, "FPTrasch") || is(object, "DPMrasch"))
   {

       if (centered) 
       { 
	   stop("This option is not implemented for FPTrasch.\n")
       }

       random<-matrix(0,nrow=object$nsubject,ncol=1)
       
       randommat<-matrix(object$save.state$randsave,
                  nrow=object$mcmc$nsave,ncol=object$nsubject+1)
      
       dimnames(randommat)<-dimnames(object$save.state$randsave)

       for(i in 1:object$nsubject){
           random[i]<-mean(object$save.state$randsave[,i])              
       }
      
       colnames(random)<-"theta"

       type<-2      
       z<-list(randomm=random,modelname=object$modelname,call=object$call,predictive=predictive,
               nsubject=object$nsubject,nrandom=1,centered=FALSE,randommat=randommat,
               type=type,nsave=object$mcmc$nsave)

       if(predictive==TRUE)
       {
          predp<-mean(object$save.state$randsave[,object$nsubject+1])      	     	

          predm<-median(object$save.state$randsave[,object$nsubject+1])      	     	

          predsd<-sqrt(var(object$save.state$randsave[,object$nsubject+1]))      	     	

          vec<-object$save.state$randsave[,object$nsubject+1]
          
          n<-length(vec)
          
          alpha<-0.05
          
          alow<-rep(0,2)
          
          aupp<-rep(0,2)
          
       
          a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
                      alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
          predl<-a$alow[1]            
          predu<-a$aupp[1]
          
          predse<-predsd/sqrt(n)
     	 
      	  predtable <- cbind(predp, predm, predsd, predse , predl , predu)
          dimnames(predtable) <- list("theta", c("Mean", "Median", "Std. Dev.", "Naive Std.Error",
                "95%HPD-Low","95%HPD-Upp"))
      	 
      	  z$prediction<-predtable
       }

       class(z)<-c("DPrandom") 
       return(z)
   }


   if(is(object, "DPrasch"))
   {

       if (centered) 
       { 
	   stop("This option is not implemented for DPrasch.\n")
       }

       random<-matrix(0,nrow=object$nsubject,ncol=1)
       
       randommat<-matrix(object$save.state$randsave,
                  nrow=object$mcmc$nsave,ncol=object$nsubject+1)
      
       dimnames(randommat)<-dimnames(object$save.state$randsave)

       for(i in 1:object$nsubject){
           random[i]<-mean(object$save.state$randsave[,i])              
       }
      
       colnames(random)<-"theta"

       type<-2      
       z<-list(randomm=random,modelname=object$modelname,call=object$call,predictive=predictive,
               nsubject=object$nsubject,nrandom=1,centered=FALSE,randommat=randommat,
               type=type,nsave=object$mcmc$nsave)

       if(predictive==TRUE)
       {
          predp<-mean(object$save.state$randsave[,object$nsubject+1])      	     	

          predm<-median(object$save.state$randsave[,object$nsubject+1])      	     	

          predsd<-sqrt(var(object$save.state$randsave[,object$nsubject+1]))      	     	

          vec<-object$save.state$randsave[,object$nsubject+1]
          
          n<-length(vec)
          
          alpha<-0.05
          
          alow<-rep(0,2)
          
          aupp<-rep(0,2)
          
       
          a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
                      alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
          predl<-a$alow[1]            
          predu<-a$aupp[1]
          
          predse<-predsd/sqrt(n)
     	 
      	  predtable <- cbind(predp, predm, predsd, predse , predl , predu)
          dimnames(predtable) <- list("theta", c("Mean", "Median", "Std. Dev.", "Naive Std.Error",
                "95%HPD-Low","95%HPD-Upp"))
      	 
      	  z$prediction<-predtable
       }

       class(z)<-c("DPrandom") 
       return(z)
   }


   if(is(object, "DPmeta") || is(object, "DPMmeta")  || is(object, "PTmeta"))
   {

       if (centered) 
       { 
	   stop("This option is not implemented for DPmeta, DPMmeta or PTmeta.\n")
       }

       randommat<-matrix(object$save.state$randsave,
                  nrow=object$mcmc$nsave,ncol=object$nrec+1)
      
       dimnames(randommat)<-dimnames(object$save.state$randsave)
       
       random<-apply(object$save.state$randsave,2,mean)
       random<-random[-(object$nrec+1)]
       
       random<-as.matrix(random,ncol=1)
       colnames(random)<-names(object$coefficients)[1]

       type<-2      
       z<-list(randomm=random,modelname=object$modelname,call=object$call,predictive=predictive,
               nsubject=object$nrec,nrandom=1,centered=FALSE,randommat=randommat,
               type=type,nsave=object$mcmc$nsave)

       if(predictive==TRUE)
       {
          predp<-mean(object$save.state$randsave[,object$nrec+1])      	     	

          predm<-median(object$save.state$randsave[,object$nrec+1])      	     	

          predsd<-sqrt(var(object$save.state$randsave[,object$nrec+1]))      	     	

          vec<-object$save.state$randsave[,object$nrec+1]
          
          n<-length(vec)
          
          alpha<-0.05
          
          alow<-rep(0,2)
          
          aupp<-rep(0,2)
          
       
          a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
                      alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
          predl<-a$alow[1]            
          predu<-a$aupp[1]
          
          predse<-predsd/sqrt(n)
     	 
      	  predtable <- cbind(predp, predm, predsd, predse , predl , predu)
          dimnames(predtable) <- list("theta", c("Mean", "Median", "Std. Dev.", "Naive Std.Error",
                "95%HPD-Low","95%HPD-Upp"))
      	 
      	  z$prediction<-predtable
       }

       class(z)<-c("DPrandom") 
       return(z)
   }


   if(is(object, "LDDPrasch"))
   {

       if (centered) 
       { 
	   stop("This option is not implemented for LDDPrasch.\n")
       }

       random <- matrix(0,nrow=object$nsubject,ncol=1)
       
       randommat<-matrix(object$save.state$randsave,
                  nrow=object$mcmc$nsave,ncol=object$nsubject+object$q)
      
       dimnames(randommat)<-dimnames(object$save.state$randsave)

       for(i in 1:object$nsubject){
           random[i]<-mean(object$save.state$randsave[,i])              
       }
      
       colnames(random)<-"theta"

       type<-2      
       z<-list(randomm=random,modelname=object$modelname,call=object$call,predictive=predictive,
               nsubject=object$nsubject,nrandom=1,centered=FALSE,randommat=randommat,
               type=type,nsave=object$mcmc$nsave)

       if(predictive==TRUE)
       {
		  #predp<-mean(object$save.state$randsave[,object$nsubject+object$q])      	     	

          #predm<-median(object$save.state$randsave[,object$nsubject+])      	     	

          #predsd<-sqrt(var(object$save.state$randsave[,object$nsubject+1]))      	     	

          #vec<-object$save.state$randsave[,object$nsubject+1]
          
          #n<-length(vec)
          
          #alpha<-0.05
          
          #alow<-rep(0,2)
          
          #aupp<-rep(0,2)
          
       
          #a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
          #            alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
          #predl<-a$alow[1]            
          #predu<-a$aupp[1]
          
          #predse<-predsd/sqrt(n)
     	 
      	  #predtable <- cbind(predp, predm, predsd, predse , predl , predu)
          #dimnames(predtable) <- list("theta", c("Mean", "Median", "Std. Dev.", "Naive Std.Error",
          #      "95%HPD-Low","95%HPD-Upp"))
      	 
      	  #z$prediction<-predtable
       }

       class(z)<-c("DPrandom") 
       return(z)
   }

    
   if(is(object, "LDDPtwopl"))
   {

       if (centered) 
       { 
	   stop("This option is not implemented for LDDPtwopl.\n")
       }

       random <- matrix(0,nrow=object$nsubject,ncol=1)
       
       randommat<-matrix(object$save.state$randsave,
                  nrow=object$mcmc$nsave,ncol=object$nsubject+object$q)
      
       dimnames(randommat)<-dimnames(object$save.state$randsave)

       for(i in 1:object$nsubject){
           random[i]<-mean(object$save.state$randsave[,i])              
       }
      
       colnames(random)<-"theta"

       type<-2      
       z<-list(randomm=random,modelname=object$modelname,call=object$call,predictive=predictive,
               nsubject=object$nsubject,nrandom=1,centered=FALSE,randommat=randommat,
               type=type,nsave=object$mcmc$nsave)

       if(predictive==TRUE)
       {
		  #predp<-mean(object$save.state$randsave[,object$nsubject+object$q])      	     	

          #predm<-median(object$save.state$randsave[,object$nsubject+])      	     	

          #predsd<-sqrt(var(object$save.state$randsave[,object$nsubject+1]))      	     	

          #vec<-object$save.state$randsave[,object$nsubject+1]
          
          #n<-length(vec)
          
          #alpha<-0.05
          
          #alow<-rep(0,2)
          
          #aupp<-rep(0,2)
          
       
          #a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
          #            alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
          #predl<-a$alow[1]            
          #predu<-a$aupp[1]
          
          #predse<-predsd/sqrt(n)
     	 
      	  #predtable <- cbind(predp, predm, predsd, predse , predl , predu)
          #dimnames(predtable) <- list("theta", c("Mean", "Median", "Std. Dev.", "Naive Std.Error",
          #      "95%HPD-Low","95%HPD-Upp"))
      	 
      	  #z$prediction<-predtable
       }

       class(z)<-c("DPrandom") 
       return(z)
   }
    

}


###
### Tools for DPrandom: print, plot
###
### Copyright: Alejandro Jara Vallejos, 2006
### Last modification: 15-12-2006.


"print.DPrandom"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
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


"plot.DPrandom"<-function(x, ask=TRUE, hpd=TRUE, nfigr=2, nfigc=2, subject=NULL, col="#bdfcc9", ...) 
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

   if(is(x, "DPrandom")){

       oldpar <- par(no.readonly = TRUE)
       par(ask = ask)
       layout(matrix(seq(1,nfigr*nfigc,1), nrow=nfigr , ncol=nfigc ,byrow=TRUE))

       if(x$predictive)
       {
          coef.p<-x$randomm
          n<-dim(coef.p)[2]
          pnames<-colnames(coef.p)
          start<-x$nsubject*n
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
       
       else
       {
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


