

  ######################
  #### predict.dblm ####
  ######################

 ## Description: predict the response of the new data newdata. 
 ##
 ##     Inputs: type (Z or D)
 ##     Output: fit
 ## 

  predict.dblm<-function(object,newdata,type.var="Z",...){
     
     # stop if the object is not a dblm object.
     if (!inherits(object, "dblm"))
        stop("use only with \"dblm\" objects")

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
           stop("If type.var=Z,the format of dblm call must be as  dblm.formula format")
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
    
    # if there is only one (k) new individual k<-1 
    if (is.null(nrow(newdata))) k<-1 
    else k <- nrow(newdata)
    
    # oensk, g diagonal of G inner products
    onesk <- rep(1,k)
    g<-diag(attr(object,"G"))
    gk <- onesk %*% t(g)
    
    # y, weights, Dsqw, Fwplus and y0
    weights<- object$weights/sum(object$weights)
    Dsqw <-diag(sqrt(weights))
    Fwplus<-attr(object,"Fwplus")
    y0<- y - sum(weights*y)

    if (type.var=="G"){
     onesn<-as.matrix(rep(1:length(y)))
       newdata<- as.matrix(newdata[,ncol(newdata)])%*%t(onesn)+
       gk -2*newdata[,1:length(y)]
    
    } 
    #the predicted y values for the new data (fit)
    fit <- sum(weights*y) + .5 * (gk-newdata) %*% Dsqw %*% Fwplus %*% Dsqw %*% y0
    
    return(fit)

  }
