

 ########################
 #### predict.dbglm  ####
 ########################

 ## Description: predict the response of the new data newdata. 
 ##
 ##     Inputs: type,type_var
 ##     Output: fit
 ## 
 

predict.dbglm <- function(object,newdata,type.pred= c("link", "response"),
                  type.var="Z",...){
                  
     # stop if the object is not a dblm object.
     if (!inherits(object, "dbglm"))
        stop("use only with \"dbglm\" objects")
     
     type.pred <- match.arg(type.pred)
        
     # controls
     if (missing(newdata))
      stop("newdata matrix must be defined")
     if (type.var!="Z"&&type.var!="D2"&&type.var!="G")
      stop("type.var must be Z for explanatory values or D2 for square distance between individuals")

     if (class(newdata)[1]=="dist"|| class(newdata)[1]=="dissimilarity")
      newdata<-as.matrix(newdata)^2

     # if new data have the explanatory values of the new indiviudals
     if (type.var=="Z"){
         newdata<-as.data.frame(newdata)
         if (attr(object,"way")=="D2"||attr(object,"way")=="G")
           stop("If type.var=Z,the format of dbglm call must be as dbglm.default format")
         z<-data.frame(attr(object,"zs"))
         if (ncol(z)!=ncol(newdata)){
            stop(gettextf("The number of columns of the newdata %d, must be the same that %d (number of Z's columns)",
            ncol(newdata),ncol(z)))
         }
         
         # names must be consistents        
         names(newdata)<-names(z)
         # used metric in dblm
         metric<-attr(object,"metric")
         # row bind x and newdata
         zaux<-rbind(z,newdata)
         for (i in 1:ncol(z))
            class(zaux[,i])<-class(z[,i]) 
         
         if (metric!="gower")
           zaux<-model.matrix(formula(paste("~",paste("zaux[,",1:ncol(zaux),"]",
                    sep="",collapse="+"),-1)))
                    
         # distance between the new k individuals and the n originals from the 
         # values of the variables x 
         newdata<- as.matrix(daisy(zaux,metric=metric))[(nrow(z)+1):(nrow(z)+
                    nrow(newdata)),1:nrow(z)] 
         # squared distance between the new and the originals individuals. 
         if (metric!="gower") newdata<-newdata^2 
      }    
    
    if (!is.matrix(newdata))
        newdata<-t(as.matrix(newdata))
   
    # recover y
    y<-object$y
    
    if (type.var=="D2"){
     if (!is.matrix(newdata))
      newdata<-t(as.matrix(newdata))
     if (attr(object,"way")=="Z"||attr(object,"way")=="G")
       stop("If type.var=D2,the format of dbglm call must be as dbglm.dist or dbglm.D2 format")  
     if (length(y)!=ncol(as.matrix(newdata))){
            stop(gettextf("The number of columns of the newdata %d, must be the same that %d (Distance dimensions)",
            ncol(newdata),length(y)))
       }
    }                                          

    if (type.var=="G"){
     if (attr(object,"way")=="Z"||attr(object,"way")=="D2")
       stop("If type.var=G,the format of dblm call must be as dbglm.Gram format") 
      if (!is.matrix(newdata))
       newdata <- t(as.matrix(newdata))
      if ((length(y)+1)!=ncol(newdata)){
       stop(gettextf("The number of columns of the newdata %d, must be the same that %d (Distance dimensions)", 
            ncol(newdata),length(y)+1))
       }
       
    # if there is only one (k) new individual k<-1 
    if (is.null(nrow(newdata))) k<-1 
    else k <- nrow(newdata)
   
    # oensk, g diagonal of G inner products
    onesk <- rep(1,k)
    g<-diag(attr(object,"G"))
    gk <- onesk %*% t(g)
     
    onesn<-as.matrix(rep(1:length(y)))
    newdata<- as.matrix(newdata[,ncol(newdata)])%*%t(onesn)+
       gk -2*newdata[,1:length(y)]

    }

   # prediction for new individuals
   neweta <- predict(attr(object,"last_dblm"),newdata,type.var="D2")

   switch(type.pred, response = {
                  fit <- object$family$linkinv(neweta)
                 },
                 link ={
                  fit <- neweta
                 }
          )
   return (fit)
}
                 