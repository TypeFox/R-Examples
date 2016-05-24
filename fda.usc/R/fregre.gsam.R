fregre.gsam<-function (formula, family = gaussian(),data=list(),weights=NULL,basis.x=NULL,
basis.b=NULL,CV=FALSE,...){
# tf <- terms.formula(formula);print("tf")
 tf <- terms.formula(formula, specials = c("s", "te", "t2"))
 terms <- attr(tf, "term.labels")
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
#if (len.smo==0) { specials<-rep(NA,len=nterms)#stop("No smoothing variables")} ##############################
#}
nterms<-length(terms)
speci<-NULL
  
specials1<-specials2<-fnf1<-fnf2<-fnf<-bs.dim1<-bs.dim2<-vfunc2<-vfunc<-vnf<-NULL

func<-nf<-sm<-rep(0,nterms)
names(func)<-names(nf)<-names(sm)<-terms
ndata<-length(data)-1
if (ndata>0) {
names.vfunc<-rep("",ndata)
for (i in 1:ndata) names.vfunc[i]<-names(data)[i+1]
}
else names.vfunc<-NULL

if (len.smo==0) {
 specials<-rep(NA,nterms)
 speci<-rep("0",nterms)
 gp$smooth.spec[1:nterms]<-NULL 
 gp$smooth.spec[1:nterms]$term<-NULL 
 vnf<-terms 
 fnf1<-fnf<-nterms
  }
  else{
for (i in 1:nterms) if (!is.na(specials[i]))     speci<-c(speci,specials[i])
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
                     vfunc<-c(vfunc,terms[i] )
                     func[i]<-1;fnf2<-c(fnf2,0)
                     bs.dim2<-c(bs.dim2,0)
                     specials2<-c(specials2,"0")
    }     } }
for (i in 1:len.smo) {
if (speci[i]!="s") {
         if (any(gp$smooth.spec[[i]]$term==names(data$df))){
                     vnf<-c(vnf,gp$smooth.spec[[i]]$margin[[1]]$term)
                     bs.dim1<-c(bs.dim1,gp$smooth.spec[[i]]$margin[[1]]$bs.dim)
                     fnf<-c(fnf,1);fnf1<-c(fnf1,1)
                     specials1<-c(specials1,speci[i])
                     }
         else {
         if (any(gp$smooth.spec[[i]]$term==names.vfunc)){
                     vfunc<-c(vfunc,gp$smooth.spec[[i]]$margin[[1]]$term)
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
                     vfunc<-c(vfunc,gp$smooth.spec[[i]]$term)
                     bs.dim2<-c(bs.dim2,gp$smooth.spec[[i]]$bs.dim)
                     fnf<-c(fnf,2);fnf2<-c(fnf2,1)
                     specials2<-c(specials2,speci[i])
   }  } }
}
}
#print("specials1");print(specials1)
#print("specials2");print(specials2)
#print("bs dim 1");print(bs.dim1)
#print("bs dim 2");print(bs.dim2)
#print(vfunc);print(vnf);print("fnf");print(fnf)
#print("func");print(func);print("nf");print(nf)
#print("fnf1");print(fnf1) ;print("fnf2");print(fnf2)
nfunc<-sum(func)
nnf<-length(vnf)
name.coef=nam=par.fregre=beta.l=list()
kterms=1   
 if (nnf>0) {
 XX=data[["df"]]
if   (attr(tf,"intercept")==0) {
     pf<- paste(pf,-1,sep="")
     }
 for ( i in 1:nnf) {
#     print(paste("Non functional covariate:",vnf[i]))
        if (fnf1[i]==1 & len.smo!=0) sm1<-TRUE
        else sm1<-FALSE
        if (sm1) {
#             print( bs.dim1[i])
                pf <- paste(pf, "+",specials1[i],"(",vnf[i],",k=",bs.dim1[i],")",sep = "")
                }
              else   pf <- paste(pf, "+", vnf[i], sep = "")
           kterms <- kterms + 1
           }
}
else {
 XX=data.frame(data[["df"]][,response])
 names(XX)=response
}
#print(paste("Functional covariate:",vfunc))
 bsp1<-bsp2<-TRUE
 lenfunc<-length(vfunc)>0
if (lenfunc) {

k=1
 mean.list=vs.list=JJ=list()
 bsp1<-bsp2<-TRUE
 for (i in 1:length(vfunc)) {
	if(class(data[[vfunc[i]]])[1]=="fdata"){
      tt<-data[[vfunc[i]]][["argvals"]]
      rtt<-data[[vfunc[i]]][["rangeval"]]
      fdat<-data[[vfunc[i]]];      dat<-data[[vfunc[i]]]$data
      if (is.null(basis.x[[vfunc[i]]]))  basis.x[[vfunc[i]]]<-create.fdata.basis(fdat,l=1:7)
      else   if (basis.x[[vfunc[i]]]$type=="pc" | basis.x[[vfunc[i]]]$type=="pls") bsp1=FALSE
      if (is.null(basis.b[[vfunc[i]]])& bsp1)  basis.b[[vfunc[i]]]<-create.fdata.basis(fdat)
      else           if (class(basis.x[[vfunc[i]]])=="fdata.comp" | basis.x[[vfunc[i]]]$type=="pls") bsp2=FALSE
      if (bsp1 & bsp2) {
          if (is.null(rownames(dat)))    rownames(fdat$data)<-1:nrow(dat)
          fdnames=list("time"=tt,"reps"=rownames(fdat[["data"]]),"values"="values")
          xcc<-fdata.cen(data[[vfunc[i]]])
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
    	    x.fd = Data2fd(argvals = tt, y = t(xcc[[1]]$data),basisobj = basis.x[[vfunc[i]]],fdnames=fdnames)
          r=x.fd[[2]][[3]]
          J=inprod(basis.x[[vfunc[i]]],basis.b[[vfunc[i]]])
          Z =t(x.fd$coefs) %*% J
          colnames(J)=colnames(Z) = name.coef[[vfunc[i]]]=paste(vfunc[i],".",basis.b[[vfunc[i]]]$names,sep="")
          XX = cbind(XX,Z)
 #         print(vfunc[i]);          print(func);          print(func[vfunc[i]])
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
        l<-basis.x[[vfunc[i]]]$l
        lenl<-length(l)
        vs <- t(basis.x[[vfunc[i]]]$basis$data)
#        Z<-basis.x[[vfunc[i]]]$x[,1:l] ##
        Z<-basis.x[[vfunc[i]]]$x[,l,drop=FALSE]
        response = "y"
        colnames(Z) = name.coef[[vfunc[i]]]=paste(vfunc[i], ".",rownames(basis.x[[vfunc[i]]]$basis$data),sep ="")
        XX = cbind(XX,Z)
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
           }       }
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
          XX = cbind(XX,Z)
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
        XX = cbind(XX,Z)
        
        vs.list[[vfunc[i]]]=vs
        mean.list[[vfunc[i]]]=basis.x[[vfunc[i]]]$meanfd
        if (fnf2[i]==1) sm2<-TRUE
        else sm2<-FALSE
          for ( j in 1:length(colnames(Z))){
              if (sm2) {
                pf <- paste(pf, "+",specials2[i],"(",name.coef[[vfunc[i]]][j],",k=",bs.dim2[i],")",sep = "")
                }
              else   pf <- paste(pf, "+", name.coef[[vfunc[i]]][j], sep = "")
           kterms <- kterms + 1
       }           }
    }
   else stop("Please, enter functional covariate")
   }
   }   }
 if (!is.data.frame(XX)) XX=data.frame(XX)

    par.fregre$data=XX
    ndatos<-nrow(XX)
    yp<-rep(NA,ndatos)  
   if (length(vfunc)==0 & length(vnf)==0)      {
       pf<- as.formula(paste(pf,1,sep=""))
       z=gam(formula=pf,data=XX,family=family,weights=weights,...)
       z$XX=XX
       z$data<-data
       z$formula.ini<-pf
       z$formula<-pf           
       class(z)<-c(class(z),"fregre.gsam") 
#       print("ffffffffffffffffffff")
        return(z )
       }
      else 
    par.fregre$formula=as.formula(pf)     
    z=gam(formula=par.fregre$formula,data=XX,family=family,weights=weights,...)
#################### reconstruccion nombre variables sm
# names2<-NULL
# for (i in 1:nterms)
#  if  (sm[i]==0)  {
#            if (nf[i]==1)         names2 <- c(names2,terms[i])
#            else                  names2 <- c(names2,terms[i])
#            }
#  else {
#       n.smooth <- length(z$smooth)
#       if (n.smooth > 0) {
#        cat("\nEstimated degrees of freedom:\n")
#        for (i in 1:n.smooth) names2 <- c(names2,z[[1]][z$smooth[[i]]$first.para:z$smooth[[i]]$last.para])
#        }
#}
 if (lenfunc) {
  z$mean=mean.list
  z$basis.x=basis.x
  z$basis.b=basis.b
  z$JJ<-JJ
  z$vs.list=vs.list
  z$vfunc<-vfunc
}
 z$formula=pf
 z$formula.ini=formula
 z$data=z$data
 z$XX=XX
 z$nnf<-nnf
 class(z)<-c(class(z),"fregre.gsam")
z
}
