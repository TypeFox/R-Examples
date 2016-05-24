predict.fregre.gsam<-function(object,newx=NULL,type="response",...){
 if (is.null(object)) stop("No fregre.gsam object entered")
 if (is.null(newx)) {
    yp=predict.gam(object,...)
    print("No newx entered")
    return(yp)
    }
 else {
 data=newx
 basis.x=object$basis.x
 basis.b=object$basis.b
 formula=object$formula.ini
 tf <- terms.formula(formula, specials = c("s", "te", "t2"))
 terms <- attr(tf, "term.labels")
 if (length(terms)==0) return(rep(object$coefficient,length=nrow(newx[[1]])) ) 
 special <- attr(tf, "specials")
 nt <- length(terms)
 specials<-rep(NULL,nt)
 if (!is.null(special$s)) specials[special$s-1]<-"s"
 if (!is.null(special$te)) specials[special$te-1]<-"te"
 if (!is.null(special$t2)) specials[special$t2-1]<-"t2"
 if (attr(tf, "response") > 0) {
        response <- as.character(attr(tf, "variables")[2])
        pf <- rf <- paste(response, "~", sep = "")
    } else pf <- rf <- "~"
vtab<-rownames(attr(tf,"factors"))
gp <- interpret.gam(formula)
#print(gp$smooth.spec)
len.smo<-length(gp$smooth.spec)
#if (len.smo==0) stop("No smoothing variables")
speci<-NULL
specials1<-specials2<-fnf1<-fnf2<-fnf<-bs.dim1<-bs.dim2<-vfunc2<-vfunc<-vnf<-NULL
nterms<-length(terms)
func<-nf<-sm<-rep(0,nterms)
 for (i in 1:nterms)
           if  (substr(terms[i],1,1)=="I")
           terms[i]<-substr(terms[i],3,nchar(terms[i])-1)
names(func)<-names(nf)<-names(sm)<-terms
if (len.smo==0) {warning("No smoothing variables");speci=NA}
else {
for (i in 1:nterms) if (!is.na(specials[i]))     speci<-c(speci,specials[i])
}
#if (object$nnf!=0)  ndata<-length(data)-1
#else ndata<-length(data)
#names.vfunc<-rep("",ndata)
   #if (object$nnf!=0) {  for (i in 1:ndata) names.vfunc[i]<-names(data)[i+1] }
#else  { for (i in 1:ndata) {    names.vfunc[i]<-names(data)[i]    }}
names.vfunc<-object$vfunc
ndata<-length(object$vfunc)+object$nnf
for (i in 1:nterms) {
         if (any(terms[i]==names(data$df))){
                     vnf<-c(vnf,terms[i] )
                     sm[i]<-nf[i]<-1
                     fnf1<-c(fnf1,0)
                     bs.dim1<-c(bs.dim1,0)
                     specials1<-c(specials1,"0")
                     }
         else {
         if (any(terms[i]==names.vfunc)){
#                     vfunc<-c(vfunc,terms[i] )
                     func[i]<-1;fnf2<-c(fnf2,0)
                     bs.dim2<-c(bs.dim2,0)
                     specials2<-c(specials2,"0")
    }     }
    }
if (len.smo!=0) for (i in 1:len.smo) {
if (speci[i]!="s") {
         if (any(gp$smooth.spec[[i]]$term==names(data$df))){
                     vnf<-c(vnf,gp$smooth.spec[[i]]$margin[[1]]$term)
                     bs.dim1<-c(bs.dim1,gp$smooth.spec[[i]]$margin[[1]]$bs.dim)
                     fnf<-c(fnf,1);fnf1<-c(fnf1,1)
                     specials1<-c(specials1,speci[i])
                     }
         else {
         if (any(gp$smooth.spec[[i]]$term==names.vfunc)){
#                     vfunc<-c(vfunc,gp$smooth.spec[[i]]$margin[[1]]$term)
                     bs.dim2<-c(bs.dim2,gp$smooth.spec[[i]]$margin[[1]]$bs.dim)
                     fnf<-c(fnf,2);fnf2<-c(fnf2,1)
                     specials2<-c(specials2,speci[i])
   }   }}
else{ if (any(gp$smooth.spec[[i]]$term==names(data$df))){
                     vnf<-c(vnf,gp$smooth.spec[[i]]$term)
                     bs.dim1<-c(bs.dim1,gp$smooth.spec[[i]]$bs.dim)
                     fnf<-c(fnf,1);fnf1<-c(fnf1,1)
                     specials1<-c(specials1,speci[i])
                     }
         else {
         if (any(gp$smooth.spec[[i]]$term==names.vfunc)){
#                     vfunc<-c(vfunc,gp$smooth.spec[[i]]$term)
                     bs.dim2<-c(bs.dim2,gp$smooth.spec[[i]]$bs.dim)
                     fnf<-c(fnf,2);fnf2<-c(fnf2,1)
                     specials2<-c(specials2,speci[i])
   }  } }
}
else speci=NA
vfunc<-object$vfunc
nfunc<-sum(func)
nnf<-length(vnf)
name.coef=nam=par.fregre=beta.l=list()
kterms=1
if (!is.null(vnf)) {
 first=FALSE
 XX=data[["df"]]
if   (attr(tf,"intercept")==0) {
     pf<- paste(pf,-1,sep="")
     }
 for ( i in 1:nnf) {
                                #      print(paste("Non functional covariate:",vnf[i]))
        if (fnf1[i]==1) sm1<-TRUE
        else sm1<-FALSE
        if (sm1) {
                pf <- paste(pf, "+",specials1[i],"(",vnf[i],",k=",bs.dim1[i],")",sep = "")
                }
              else   pf <- paste(pf, "+", vnf[i], sep = "")
           kterms <- kterms + 1
           }
}
else {  first=TRUE}
#                            print(paste("Functional covariate:",vfunc))
                            lenfunc<-length(vfunc)
if (lenfunc>0) {
k=1
 mean.list=vs.list=JJ=list()
 bsp1<-bsp2<-TRUE
 for (i in 1:lenfunc) {
	if(class(newx[[vfunc[i]]])[1]=="fdata"){
      tt<-data[[vfunc[i]]][["argvals"]]
      rtt<-data[[vfunc[i]]][["rangeval"]]
      fdataobj<-data[[vfunc[i]]]
      fdat<-data[[vfunc[i]]];      dat<-fdataobj$data
                if (nrow(dat)==1) rwn<-NULL
         else rwn<-rownames(dat)
      if (is.null(basis.x[[vfunc[i]]]))  basis.x[[vfunc[i]]]<-create.fdata.basis(fdat,l=1:7)
      else   if (basis.x[[vfunc[i]]]$type=="pc" | basis.x[[vfunc[i]]]$type=="pls") bsp1=FALSE
      if (is.null(basis.b[[vfunc[i]]])& bsp1)  basis.b[[vfunc[i]]]<-create.fdata.basis(fdat)
      else           if (basis.x[[vfunc[i]]]$type=="pc" | basis.x[[vfunc[i]]]$type=="pls") bsp2=FALSE
      if (bsp1 & bsp2) {
          if (is.null(rownames(dat)))    rownames(dat)<-1:nrow(dat)
          fdnames=list("time"=tt,"reps"=rwn,"values"="values")
          xcc<-fdata.cen(data[[vfunc[i]]],object$mean[[vfunc[i]]])
          mean.list[[vfunc[i]]]=xcc[[2]]

          if (!is.null( basis.x[[vfunc[i]]]$dropind)) {
              int<-setdiff(1:basis.x[[vfunc[i]]]$nbasis,basis.x[[vfunc[i]]]$dropind)
              basis.x[[vfunc[i]]]$nbasis<-length(int)
              basis.x[[vfunc[i]]]$dropind<-NULL
              basis.x[[vfunc[i]]]$names<-basis.x[[vfunc[i]]]$names[int]
              }
          if (!is.null( basis.b[[vfunc[i]]]$dropind)) {
              int<-setdiff(1:basis.b[[vfunc[i]]]$nbasis,basis.b[[vfunc[i]]]$dropind)
              basis.b[[vfunc[i]]]$nbasis<-length(int)
              basis.b[[vfunc[i]]]$dropind<-NULL
              basis.b[[vfunc[i]]]$names<-basis.b[[vfunc[i]]]$names[int]
              }
#x.fd = Data2fd(argvals = tt, y = t(xcc[[1]]$data),basisobj = basis.x[[vfunc[i]]])
#x.fd = Data2fd(argvals = tt, y = t(xcc[[1]]$data),basisobj = basis.x[[vfunc[i]]],fdnames=rwn)
      x.fd = Data2fd(argvals = tt, y = t(xcc[[1]]$data),basisobj = basis.x[[vfunc[i]]],fdnames=fdnames)
       r=x.fd[[2]][[3]]
       J=inprod(basis.x[[vfunc[i]]],basis.b[[vfunc[i]]])
          Z =t(x.fd$coefs) %*% J
          colnames(J)=colnames(Z) = name.coef[[vfunc[i]]]=paste(vfunc[i],".",basis.b[[vfunc[i]]]$names,sep="")
#          XX = cbind(XX,Z)
        if (fnf2[i]==1) sm2<-TRUE
        else sm2<-FALSE
          for ( j in 1:length(colnames(Z))){
              if (sm2) {
#             print( bs.dim2[i])
                pf <- paste(pf, "+",specials2[i],"(",name.coef[[vfunc[i]]][j],",k=",bs.dim2[i],")",sep = "")
                }
              else   pf <- paste(pf, "+", name.coef[[vfunc[i]]][j], sep = "")
           kterms <- kterms + 1
           }
        	JJ[[vfunc[i]]]<-J
				}
         else {
        l<-nrow(basis.x[[vfunc[i]]]$basis)
        vs <- t(basis.x[[vfunc[i]]]$basis$data)
        xcen<-fdata.cen(fdat,object$mean[[vfunc[i]]])[[1]]
        if ( basis.x[[vfunc[i]]]$type=="pls"){
            if (basis.x[[vfunc[i]]]$norm)  {  
              nn<-nrow(fdat)  
              sd.X <- sqrt(apply(basis.x[[vfunc[i]]]$fdataobj$data, 2, var))             
              xcen$data<- xcen$data/(rep(1, nn) %*% t(sd.X))
            } 
}
        Z<- inprod.fdata(xcen,object$vs.list[[vfunc[i]]])
        response = "y"
        colnames(Z) = name.coef[[vfunc[i]]]=paste(vfunc[i], ".",rownames(basis.x[[vfunc[i]]]$basis$data),sep ="")
#        XX = cbind(XX,Z)
        vs.list[[vfunc[i]]]=basis.x[[vfunc[i]]]$basis
        mean.list[[vfunc[i]]]=basis.x[[vfunc[i]]]$mean

        if (fnf2[i]==1) sm2<-TRUE
        else sm2<-FALSE
        for ( j in 1:length(colnames(Z))){
              if (sm2) {
                pf <- paste(pf, "+",specials2[i],"(",name.coef[[vfunc[i]]][j],",k=",bs.dim2[i],")",sep = "")
                }
              else   pf <- paste(pf, "+", name.coef[[vfunc[i]]][j], sep = "")
           kterms <- kterms + 1
       }
       }
#       print(colnames(Z))
        if (first) {
           XX=Z;       first=FALSE         }
        else {
         XX = cbind(XX,Z)}
  }
 	else {
 		if(class(data[[vfunc[i]]])[1]=="fd"){
      fdat<-data[[vfunc[i]]]
      if (is.null(basis.x[[vfunc[i]]]))  basis.x[[vfunc[i]]]<-fdat$basis
      else   if (class(basis.x[[vfunc[i]]])=="pca.fd") bsp1=FALSE
      if (is.null(basis.b[[vfunc[i]]])& bsp1)
         basis.b[[vfunc[i]]]<-create.fdata.basis(fdat,
         l=1:max(5,floor(basis.x[[vfunc[i]]]$nbasis/5)),type.basis=basis.x[[vfunc[i]]]$type,
         rangeval=fdat$basis$rangeval)
      else           if (class(basis.x[[vfunc[i]]])=="pca.fd") bsp2=FALSE
      if (bsp1 & bsp2) {
          r=fdat[[2]][[3]]
          if (!is.null( basis.x[[vfunc[i]]]$dropind)) {
              int<-setdiff(1:basis.x[[vfunc[i]]]$nbasis,basis.x[[vfunc[i]]]$dropind)
              basis.x[[vfunc[i]]]$nbasis<-length(int)
              basis.x[[vfunc[i]]]$dropind<-NULL
              basis.x[[vfunc[i]]]$names<-basis.x[[vfunc[i]]]$names[int]
              }
          if (!is.null( basis.b[[vfunc[i]]]$dropind)) {
              int<-setdiff(1:basis.b[[vfunc[i]]]$nbasis,basis.b[[vfunc[i]]]$dropind)
              basis.b[[vfunc[i]]]$nbasis<-length(int)
              basis.b[[vfunc[i]]]$dropind<-NULL
              basis.b[[vfunc[i]]]$names<-basis.b[[vfunc[i]]]$names[int]
              }
          J=inprod(basis.x[[vfunc[i]]],basis.b[[vfunc[i]]])
          mean.list[[vfunc[i]]]<-mean.fd(x.fd)
          x.fd<-center.fd(x.fd)
          Z =t(x.fd$coefs) %*% J
          colnames(J)=colnames(Z) = name.coef[[vfunc[i]]]=paste(vfunc[i],".",basis.b[[vfunc[i]]]$names,sep="")
#          XX = cbind(XX,Z)
          if (fnf2[i]==1) sm2<-TRUE
          else sm2<-FALSE
          for ( j in 1:length(colnames(Z))){
              if (sm2) {
                pf <- paste(pf, "+",specials2[i],"(",name.coef[[vfunc[i]]][j],",k=",bs.dim2[i],")",sep = "")
                }
              else   pf <- paste(pf, "+", name.coef[[vfunc[i]]][j], sep = "")
           kterms <- kterms + 1
           }
        	JJ[[vfunc[i]]]<-J
				}
      else {
        l<-ncol(basis.x[[vfunc[i]]]$scores)
        vs <- basis.x[[vfunc[i]]]$harmonics$coefs
        Z<-basis.x[[vfunc[i]]]$scores
        response = "y"
        colnames(Z) = name.coef[[vfunc[i]]]=paste(vfunc[i], ".",colnames(basis.x[[vfunc[i]]]$harmonics$coefs),sep ="")
#        XX = cbind(XX,Z)
        vs.list[[vfunc[i]]]=vs
        mean.list[[vfunc[i]]]=basis.x[[vfunc[i]]]$meanfd
        if (fnf2[i]==1) sm2<-TRUE
        else sm2<-FALSE
        for ( j in 1:length(colnames(Z))){
              if (sm2) {
#                print( bs.dim2[i])
                pf <- paste(pf, "+",specials2[i],"(",name.coef[[vfunc[i]]][j],",k=",bs.dim2[i],")",sep = "")
                }
              else   pf <- paste(pf, "+", name.coef[[vfunc[i]]][j], sep = "")
           kterms <- kterms + 1
    }       }
     if (first) {    XX=Z;              first=FALSE         }
     else XX = cbind(XX, Z)
    }
   else stop("Please, enter functional covariate")
   }
   }
        }
if (!is.data.frame(XX)) XX=data.frame(XX)
 yp=predict.gam(object=object,newdata=XX,type=type,...)
return(yp)
}
}
 
