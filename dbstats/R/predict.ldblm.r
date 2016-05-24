

  ######################
  #### predict.ldblm ####
  ######################

 ## Description: predict the response of the new data newdata1 
 ##             (newdata2 for kernel weights). 
 ##
 ##     Inputs: type and type_var (Z or D)
 ##     Output: fit and newS
 ## 




predict.ldblm<-function(object,newdata1,newdata2=newdata1,new.k.knn=3,
      type.var="Z",...){

     # stop if the object is not a dblm object.
     if (!inherits(object, "ldblm"))
        stop("use only with \"ldblm\" objects")

     y<-object$y
     # controls of newdata
     if (type.var!="Z"&&type.var!="D2"&&type.var!="G")
      stop("type.var must be Z for explanatory values, D2 for square distance between individuals or G for centered Euclidean configuration")
     if (missing(newdata1))
      stop("newdata1 matrix must be defined")  

     # transform an object of class dist to one of "D2"
     if (class(newdata1)[1]=="dist"|| class(newdata1)[1]=="dissimilarity")
      newdata1<-as.matrix(newdata1)^2
     
     if (class(newdata2)[1]=="dist"|| class(newdata2)[1]=="dissimilarity")
      newdata2<-as.matrix(newdata2)^2
          
     if (type.var=="Z"){
        if (attr(object,"way")=="G"||attr(object,"way")=="D2")
        stop("If type.var=Z,the format of dblm call must be as ldblm.dist format")
       newdata1<-as.data.frame(newdata1) 
       nr <- nrow(newdata1)
       nc <- ncol(newdata1) 
       z<-data.frame(attr(object,"zs"))
       if (ncol(z)!=nc){
            stop(gettextf("The number of columns of the newdata1 %d, must be the same that %d (number of Z's columns)", 
            ncol(as.matrix(newdata1)),ncol(z)))
       }
       
       # used metric in ldblm
       metric1<-attr(object,"metric1")
       metric2<-attr(object,"metric2")
       # names must be consistents  
       names(newdata1)<-names(z)
       # row bind x and newdata
       zaux<-rbind(z,newdata1)
       for (i in 1:ncol(z))
         class(zaux[,i])<-class(z[,i])    
      
       # distance between the new k individuals and the n originals from the 
       # values of the variables x  
       if (metric1!="gower"|| metric2!="gower")
         zaux2<-model.matrix(formula(paste("~",paste("zaux[,",1:ncol(zaux),"]",sep="",collapse="+"),-1)))
         
       if (metric1=="gower")                         
        newdata1<- as.matrix(daisy(zaux,metric=metric1))[(nrow(z)+1):(nrow(z)+nrow(newdata1)),1:nrow(z)] 
       else{
        newdata1<- as.matrix(daisy(zaux2,metric=metric1))[(nrow(z)+1):(nrow(z)+nrow(newdata1)),1:nrow(z)] 
        newdata1<-newdata1^2 
       }                         
       
       if (metric1==metric2)
        newdata2 <- newdata1
       else if (metric2=="gower")
        newdata2<- as.matrix(daisy(zaux,metric=metric2))[(nrow(z)+1):(nrow(z)+nrow(newdata1)),1:nrow(z)] 
       else{ 
        newdata2<- as.matrix(daisy(zaux2,metric=metric2))[(nrow(z)+1):(nrow(z)+nrow(newdata1)),1:nrow(z)]   
        newdata2<-newdata2^2 
       }
  
       if (!is.matrix(newdata1))
       newdata1<-t(as.matrix(newdata1))
          
       if (!is.matrix(newdata2))                 
        newdata2<-t(as.matrix(newdata2))
     } 
           
     if (type.var=="D2"){
       if (attr(object,"way")=="Z"||attr(object,"way")=="G")
        stop("If type.var=D2,the format of ldblm call must be as ldblm.dist format") 
      if (!is.matrix(newdata1))
         newdata1<-t(as.matrix(newdata1))
      if (!is.null(newdata2)){ 
       if (!is.matrix(newdata2))
         newdata2<-t(as.matrix(newdata2)) 
      }else
       newdata2<-newdata1   
    }
    
     if (type.var=="G"){
       if (attr(object,"way")=="Z"||attr(object,"way")=="D2")
        stop("If type.var=G,the format of ldblm call must be as ldblm.Gram format") 
      if (!is.matrix(newdata1))
         newdata1<-t(as.matrix(newdata1))
      if (!is.null(newdata2)){ 
       if (!is.matrix(newdata2))
         newdata1<-t(as.matrix(newdata2)) 
      }else
       newdata2<-newdata1  
      
     if (is.null(nrow(newdata1))) k<-1 
     else k <- nrow(newdata1)
   
     # oensk, g diagonal of G inner products
     onesk <- rep(1,k)
     gk1<- onesk %*%t(diag(attr(object,"G1")))
     gk2<- onesk %*%t(diag(attr(object,"G2")))
      
     onesn<-as.matrix(rep(1:length(y)))
      
     newdata1<- as.matrix(newdata1[,ncol(newdata1)])%*%t(onesn)+gk1
      -2*newdata1[,1:length(y)]

     newdata2<- as.matrix(newdata2[,ncol(newdata2)])%*%t(onesn)+gk2
       t(onesn%*%t(as.matrix(newdata2[,ncol(newdata2)]))) -2*newdata2[,1:length(y)]                          
    }
    
    nr <- nrow(newdata1)
    nc <- ncol(newdata1) 
     
    # controls: newdata1 dimension must be the same as newdata2.                
    if (nr!=nrow(newdata2))
     stop("newdata1 and newdata2 must have the same dimensions")
    if (nc!=ncol(newdata2))
     stop("number of newdata2 columns must be the same as in newdata1")
     
    dist2<-attr(object,"dist2")
    if (ncol(newdata2)!=ncol(dist2))
     stop("number of newdata2 columns must be the same as in dist2 of ldblm object")  

    # bandwidth h_opt  and kind of kernel of ldblm object
    h.opt<-object$h.opt
    kind.of.kernel<-attr(object,"kind.of.kernel")
     
    # max between the minimum bandwidth h.knn.new 
    # (3 nearest neighbors in the newdata1) and the using h.
    if (new.k.knn>1)
      new.h.knn<-h.knn.funct(newdata1^.5,k=new.k.knn) 
    else
     stop(" new.k.knn must be >1")    
    hi<-pmax(new.h.knn+1e-10,object$h.opt)
         
    # initialize n, y, fit, newS, newShat and recover S
    n<-length(object$fitted.values)
    fit <- rep(1,nr)
    newS <- matrix(0,nr,n)
    newShat <- matrix(0,nr,n)
    S<-object$S
    weights <-object$weights
    
    rel.gvar<-attr(object,"rel.gvar")
    eff.rank<-attr(object,"eff.rank")
    if (is.null(eff.rank)) method="rel.gvar"
    else method="eff.rank"   
        
    for (i in 1:nr){
      newS[i,] <- kernel.number(newdata1[i,]^.5/hi[i],j=kind.of.kernel)
       
       # find the observations without null weight
       iid<-which(newS[i,]>0)
       
       # auxiliar weights, auxiliar y and auxiliar dist2 only with the iid observations.
       weights_aux<-newS[i,iid]
       y_aux<-y[iid]
       dist2_aux<-dist2[iid,iid]
       newdataaux<-as.matrix(t(newdata2[i,][iid]))
       class(dist2_aux)<-"D2"
       if (!is.null(eff.rank)) eff.rank_aux=min(length(y_aux)-1,eff.rank)
       if (is.null(eff.rank))  eff.rank_aux=NULL
       
       # call dblm linear model to achieve the fitted values yhat and the Hat 
       # matrix Hhat for the i observation
       dblmaux <- dblm.D2(y=y_aux,D2=dist2_aux,weights=weights[iid]*weights_aux,method=method,#"OCV",
                    rel.gvar=rel.gvar,eff.rank=eff.rank_aux)#,full.search=T)
       fit[i]<-predict(dblmaux,newdata=newdataaux,type.var="D2")
    }

    ans<-list(fit=fit,newS=newS)
    class(ans)="predict.ldblm"
    return(ans)
}
                                         