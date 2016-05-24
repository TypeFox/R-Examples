predict.dbplsr <-function(object,newdata,type.var="Z",...){
 
      # stop if the object is not a dblm object.
     if (!inherits(object, "dbplsr"))
        stop("use only with \"dbplsr\" objects")

     # controls 
     if (missing(newdata))
      stop("newdata matrix must be defined")
     if (type.var!="Z"&&type.var!="D2"&&type.var!="G")
      stop("type.var must be Z for explanatory values, D2 for square distance between individuals or G for centered Euclidean configuration")
     
     if (class(newdata)[1]=="dist"|| class(newdata)[1]=="dissimilarity")
      newdata<-disttoD2(newdata)
      
     # if new data have the explanatory values of the new indiviudals      
     if (type.var=="Z"){
         newdata<-as.data.frame(newdata)
         if (attr(object,"way")=="D2"||attr(object,"way")=="G")
           stop("If type.var=Z,the format of dblm call must be as dblm.default or dblm.formula format")
         z<-data.frame(attr(object,"zs"))
         if (ncol(z)!=ncol(newdata)){
            stop(gettextf("The number of columns of the newdata %d, must be the same that %d (number of Z's columns)", 
            ncol(newdata),ncol(z)))
         }
         
         # used metric in dblm
         metric<-attr(object,"metric")
         # names must be consistents  
         names(newdata)<-names(z)
         # row bind x and newdata
         zaux<-rbind(z,newdata)
         for (i in 1:ncol(z))
            class(zaux[,i])<-class(z[,i])
            
         # see if at least one explanatory variable is the class=factor
         if (metric!="gower")         
          zaux<-model.matrix(formula(paste("~",paste("zaux[,",1:ncol(zaux),
                    "]",sep="",collapse="+"),-1)))
        
         # distance between the new k individuals and the n originals from the
         # values of the variables x
         newdata<- as.matrix(daisy(zaux,metric=metric))[(nrow(z)+1):nrow(zaux),1:nrow(z)] 
                      
         # squared distance between the new and the originals individuals.              
         if (metric!="gower")newdata<-newdata^2  
    }    
    
    # recover y
    y<-object$y
    
    if (type.var=="D2"){
     if (attr(object,"way")=="Z"||attr(object,"way")=="G")
       stop("If type.var=D2,the format of dblm call must be as dblm.dist or dblm.D2 format") 
      if (!is.matrix(newdata))
       newdata <- t(as.matrix(newdata))
      if (length(y)!=ncol(newdata)){
       stop(gettextf("The number of columns of the newdata %d, must be the same that %d (Distance dimensions)", 
            ncol(newdata),length(y)))
       }
    }
    
    if (type.var=="G"){
     if (attr(object,"way")=="Z"||attr(object,"way")=="D2")
       stop("If type.var=G,the format of dblm call must be as dblm.Gram format") 
      if (!is.matrix(newdata))
       newdata <- t(as.matrix(newdata))
      if ((length(y)+1)!=ncol(newdata)){
       stop(gettextf("The number of columns of the newdata %d, must be the same that %d (Distance dimensions)", 
            ncol(newdata),length(y)+1))
       }
    }
 
    # weights
    weights <- object$weights
    weights <- weights/sum(weights)
    y0<- y - sum(weights*y)
     
    # Ghat i Gtit
    Gtit <- object$Gk
    # Ghat <- object$H%*%object$G0%*%object$H
    Ghat <- object$G0-Gtit
    G<-Ghat
     
    # if there is only one (k) new individual k<-1 
    if (is.null(nrow(newdata))) k<-1 
    else k <- nrow(newdata)
        
    # oensk, g diagonal of G inner products
    onesk <- rep(1,k)
    g <- object$gvec
    gk <- onesk %*% t(g)
     
    # newdata for format Gram
    if (type.var=="G"){
        onesn<-as.matrix(rep(1:length(y)))
          newdata<- as.matrix(newdata[,ncol(newdata)])%*%t(onesn)+
          gk -2*newdata[,1:length(y)]
        
       } 
     
    Dw <- diag(weights)
    Dsqw <- sqrt(Dw)
     
    F <- object$fk
    F <- as.matrix(as.data.frame(F))
    lamda2<-t(F)%*%F
    svdl2<-svd(lamda2)
    if (length(svdl2$d)==1)
     lamda<-svdl2$v%*%sqrt(svdl2$d)%*%svdl2$u
    else 
     lamda<-svdl2$v%*%diag(sqrt(svdl2$d))%*%svdl2$u
    
    lamba.svd<-svd(lamda)
    lambda.plus<- lamba.svd$u%*%diag(lamba.svd$d^-1)%*%lamba.svd$v
        
    F1 <- F%*%lambda.plus
    N1 <- t(F1)%*%Dw%*%Ghat%*%F1  
    Nplus <-  solve(N1)
    aux<-Dsqw%*%F1%*%Nplus%*%t(F1)%*%Dsqw 
    fit <-sum(weights*y) + .5*(gk-newdata)%*%(aux)%*%y0
    
    return(fit)
}

