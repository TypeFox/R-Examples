PdataPed<-function(formula=NULL, data=NULL, id=data$id, sex=data$sex, offspring=data$offspring, timevar=data$timevar, USdam=FALSE, USsire=FALSE){

 if(is.null(timevar)){
   timevar<-rep(1,length(id))
 }

 if(length(data)!=0){
  object<-list(formula=formula, data=data, id=id, sex=sex, offspring=offspring, timevar=timevar, USdam=USdam, USsire=USsire)
 }else{
  object<-list(formula=NULL, data=NULL, id=NULL, sex=NULL, offspring=NULL, timevar=timevar, USdam=FALSE, USsire=FALSE)
 }

 class(object)<-c("PdataPed", "list")
 object

}

is.PdataPed<-function(x){ 
inherits(x, "PdataPed")
 }

GdataPed<-function(G=NULL, id=NULL, categories=NULL, perlocus=FALSE, marker.type="MSW"){

if(length(G)!=0){
   if(length(id)==0){
     if(length(G$id)==0){
       stop("id is empty and is not found in names(G)")
     }else{
       id<-G$id
       G<-G[,-which(names(G)=="id")]
     }
   }

   id<-as(id, "character")

   if(length(categories)==0){
     if(length(G$categories)==0){
       categories<-rep(1,length(id))
     }else{
       categories<-G$categories       
       G<-G[,-which(names(G)=="categories")]
     }
   }

   if(is.genotype(G[[1]])==FALSE & is.genotypeD(G[[1]])==FALSE){
     G<-genotype.list(G, marker.type=marker.type)
   }
 }  
 object<-list(G=G, id=id, categories=categories, perlocus=perlocus, marker.type=marker.type)
 class(object)<-c("GdataPed", "list")
 object

}

is.GdataPed<-function(x){ 
inherits(x, "GdataPed")
 }

startPed<-function(G=NULL,id=NULL, estG=TRUE, A=NULL, estA=TRUE, E1=NULL, estE1=TRUE, E2=NULL, estE2=TRUE, ped=NULL, estP=TRUE, beta=NULL, estbeta=TRUE,  USdam=NULL, estUSdam=TRUE, USsire=NULL, estUSsire=TRUE){

object<-list(G=G, id=id, estG=estG, A=A, estA=estA, E1=E1, estE1=estE1, E2=E2, estE2=estE2, ped=ped, estP=estP, beta=beta, estbeta=estbeta,  USdam=USdam, estUSdam=estUSdam, USsire=USsire, estUSsire=estUSsire)

 class(object)<-c("startPed", "list")
 object
}

is.startPed<-function(x){ 
inherits(x, "startPed")
 }

tunePed<-function(E1=NULL, E2=NULL, beta=NULL,  USdam=NULL, USsire=NULL){

  object<-list(E1=E1, E2=E2, beta=beta, USdam=USdam, USsire=USsire)
   class(object)<-c("tunePed", "list")
   object

}

is.tunePed<-function(x){ 
inherits(x, "tunePed")
}


priorPed<-function(E1=999, E2=999, beta=list(mu=999, sigma=999), USdam=list(mu=999, sigma=999),  USsire=list(mu=999, sigma=999)){

  if(length(E1)==2 & is.matrix(E1)==FALSE){
    E1<-t(as.matrix(E1))
  }
  if(length(E2)==2 & is.matrix(E2)==FALSE){
    E2<-t(as.matrix(E2))
  }


object<-list(E1=E1, E2=E2, beta=beta, USdam=USdam, USsire=USsire)
 class(object)<-c("priorPed", "list")
 object

}

is.priorPed<-function(x){
 inherits(x, "priorPed") 
}






