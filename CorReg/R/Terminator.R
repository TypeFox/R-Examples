#' Destructing values to have missing ones
#' @export
#' @param wrath the ratio of missing values in the output
#' @param target the dataset (matrix or data.frame) in which missing values will be made
#' @param diag if =1 it creates a diagonal band of missing values (no complete line, no complete column, but not too much missing values)
#' @param Z adjacency matrix to coerce a maximum of 1 missing value per sub-regression for each individual
#' @return the matrix with missing values.


#' @examples
#'    \dontrun{
#'   rm(list=ls())#clean the workspace
#'   
#'    data<-mtcars
#'    require(CorReg)
#'   datamiss=Terminator(target = data,wrath=0.05)#5% of missing values
#'   showdata(datamiss)#plot positions of the missing values
#'   datamiss=Terminator(target = data,diag=1)#diag of missing values
#'   showdata(datamiss)#plot positions of the missing values (no full individuals, no full variable)

#'     }


Terminator<-function(target=NULL, wrath=0.1,diag=0,Z=NULL){
   if(is.null(target)){
      target="Sarah Connor"
   }
   if(is.null(dim(target))){
      if(target[1]=="Sarah Connor" ){
      print("I'll be back !")}
   }else if (is.null(dim(target))){
      if( target[1]=="bender"){
      Bender(a=wrath,b=diag,c=Z)}
   }else if(diag>0){
      n=nrow(target)
      p=ncol(target)
      quidiag=cbind(rep(1:n,length.out=max(n,p)),rep(1:p,length.out=max(n,p)))
      target[quidiag]=NA
      for(j in 2:diag){
         quidiag=cbind(rep(1:n,length.out=max(n,p)),rep(c(j:p,1:(j-1)),length.out=max(n,p)))
         quidiag=cbind(rep(1:n,length.out=max(n,p)),rep(1:p,length.out=max(n,p)))
         target[quidiag]=NA
      }
   }else if (wrath>0){
      target=as.matrix(target)
      n=nrow(target)
      p=ncol(target)
      quidiag=cbind(rep(1:n,length.out=max(n,p)),rep(1:p,length.out=max(n,p)))
      nbmank=floor(wrath*n*p)
      loc=target
      target=NA*target
      if(!is.null(Z)){
         Zc=colSums(Z)
         quiZ_vect=which(Z!=0,arr.ind=TRUE)#liste des impliques en ssreg
         quiZtot=unique(c(quiZ_vect))#liste des impliques en ssreg
         for (i in 1:n){
            quiZ=quiZtot#liste des impliques en ssreg
            quiblok=c()
            for(j in 1:length(quiZ)){#maxi pr manquants
               if(length(quiZ)>0){
                  mankloc=sample(quiZ,size=1)#on tue quelqu'un
                  quiZ=quiZ[quiZ!=mankloc]#le mort n'est plus candidat ni bloquable
                  if(Zc[mankloc]>0){#variable a gauche, on retire la regression
                     reste=which(Z[,mankloc]!=0)
                     quiblok=c(quiblok,reste)
                     quiZ=unique(c(reste,quiZ))[-c(1:length(reste))]#on enleve les bloques
                  }else{#variable a droite, on retire la gauche, les autres a droites, et les autres regressions touchees
                     impact=quiZ_vect[quiZ_vect[1,]==mankloc,2]
                     for(k in impact){
                        if(length(quiZ)>0){
                           ssreg=c(k,which(Z[,k]!=0))
                           ssreg=ssreg[ssreg!=mankloc]
                           quiblok=unique(c(quiblok,ssreg))
                           quiZ=unique(c(quiblok,quiZ))[-c(1:length(quiblok))]#on enleve les bloques
                        }else{
                           break
                        }
                     }
                  }                  
               }else{
                  break
               }
            }
            quidiag=rbind(quidiag,cbind(i,unique(quiblok)))
         }
      }
      target[quidiag]=loc[quidiag]
      mankmax=n*p-length(quidiag)
      
      if(nbmank<mankmax){
         candidat=which(is.na(target),arr.ind=TRUE)
         quimank=sample(1:nrow(candidat),size=nbmank)
         target=loc
         target[candidat[quimank,]]=NA
      }
   }
   return(target)
}