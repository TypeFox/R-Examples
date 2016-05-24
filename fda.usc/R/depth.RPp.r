 #hacer un dfunc,par.dfunc
depth.RPp<-function (lfdata,lfdataref=lfdata, nproj = 50, proj="vexponential",trim = 0.25,
dfunc="mdepth.TD",par.dfunc=list(scale=TRUE),draw = FALSE,ask=FALSE){  
verbose<-TRUE
if (class(lfdata)!="list") stop("lfdata input must be a list of fdata objects")
if (class(lfdataref)!="list") stop("lfdataref input must be a list of fdata objects")
 lenl<-length(lfdata)
 lenl2<-length(lfdataref)
 n<-nrow(lfdata[[1]])
 m<-nrow(lfdataref[[1]]) 
 mdist2<-matrix(0,n,n)
 amdist<-array(NA,dim=c(n,m,lenl))
 mdist<-list()
 nam1<-names(lfdata)
 nam2<-names(lfdataref) 
 if (is.null(nam1)) {names(lfdata)<-nam1<-paste("var",1:lenl,sep="")}
 if (is.null(nam2)) {names(lfdataref)<-nam2<-paste("var",1:lenl2,sep="")} 
# depth<-paste("mdepth.",dfunc,sep="")
# if (missing(dfunc))  dfunc<-depth.mode#dfunc<-rep("depth.mode",len=lenl)
# else if (length(dfunc)==1)   dfunc<-rep("depth.mode",len=lenl)
  if (is.fdata(lfdata[[nam1[1]]]))  {
    fdataobj<-lfdata[[nam1[1]]]
    fdataori<-lfdataref[[nam1[1]]]}
  else stop("No fdata object in the lfdata argument")
  data <- fdataobj[["data"]]
  data2 <- fdataori[["data"]]    
  names1 <- names2 <- names <- fdataobj[["names"]]
  names1$main <- "depth.RPD median"
  names2$main <- paste("RPD trim ", trim * 100, "%", sep = "")
  n <- nrow(data)
  m <- ncol(data)
  m2<-ncol(data2)
  n2<-nrow(data2) 
  if (is.null(m) && is.null(m2)) stop("ERROR IN THE DATA DIMENSIONS")
  modulo = function(z) {        sqrt(sum(z^2))    }
  if (is.null(n) || is.null(m))       stop("Input must be a matrix")
  tt = fdataobj[["argvals"]]
  rtt <- fdataobj[["rangeval"]]
  newfunc = array(NA, dim = c(n, m, lenl))
  newfunc2 = array(NA, dim = c(n2, m2, lenl))    
##################
  for (ider in 1:lenl) {
          newfunc[, , ider] = lfdata[[nam1[ider]]]$data     
          newfunc2[, , ider] = lfdataref[[nam1[ider]]]$data 
  }
  dep = rep(0, n)
  dep2 = rep(0, n2)    
  vproject = matrix(0, nrow = n, ncol = lenl)
  vproject2 = matrix(0, nrow = n2, ncol = lenl)   
    
#  z = matrix(rnorm(m * nproj), nrow = nproj, ncol = m)
#  modu = apply(z, 1, modulo)
#  z = z/modu
  if (is.fdata(proj)) {
     if (fdataobj$argvals!=proj$argvals || ncol(fdataobj)!=ncol(proj)) stop("Error en proj dimension")
     z<-proj
     nproj<-nrow(z)
    }
  else {	 z<-rproc2fdata(nproj,tt,sigma=proj,norm=TRUE)	}
  if (verbose) 
  pb = txtProgressBar(min = 0, max = nproj, style = 3)
  for (j in 1:nproj) {
      if (verbose)         setTxtProgressBar(pb, j - 0.5)
        for (ider in 1:lenl) {
            matriz = newfunc[, , ider]
            vproject[, ider] = matriz %*%z$data[j, ]
            vproject2[, ider] = newfunc2[, , ider] %*%z$data[j, ]              
        }
#        resul = dfunc2(vproject, ...)       
#        par.dfunc = list()
     if (substr(dfunc,1,1)=="m") {
        par.dfunc$x<- vproject       
        par.dfunc$xx <- vproject2
        }
     else{ 
        par.dfunc$fdataobj <- fdata(vproject)        
        par.dfunc$fdataori <- fdata(vproject2)
        }
#        par.dfunc$scale<-TRUE
        resul = do.call(dfunc, par.dfunc)
        if (dfunc=="depth.RP" | dfunc=="mdepth.HS" | dfunc=="mdepth.RP") par.dfunc$proj<-resul$proj        
        dep = dep + resul$dep
     if (verbose)         setTxtProgressBar(pb, j)
  }                                                               
  if (verbose)     close(pb)
    dep = dep/nproj
    k = which.max(dep)
    med = data[k,,drop=FALSE ]
    lista = which(dep >= quantile(dep, probs = trim,na.rm=TRUE))
    mtrim=apply(data[lista,,drop=FALSE],2,mean,na.rm=TRUE)
#else mtrim<-data     
    tr <- paste("RPD.tr", trim * 100, "%", sep = "")
    med <- fdata(med, tt, rtt, names1)
    mtrim <- fdata(mtrim, tt, rtt, names2)
    rownames(med$data) <- "RPD.med"
    rownames(mtrim$data) <- tr
if (draw){
  mf=5
  if (lenl>4) ask=TRUE
  if (ask) {par(mfrow = c(1, 1))
            dev.interactive()
            oask <- devAskNewPage(TRUE)
            on.exit(devAskNewPage(oask))}
   else{    mf<-switch(lenl,
   "1"={c(1,1)},
   "2"={c(1,2)},
   "3"={c(1,3)},
   "4"={c(2,2)})            
            par(mfrow =mf)                    }          
 for (idat in 1:lenl) {
   data<-lfdata[[idat]]$data
   med<-data[k,,drop=FALSE]      
  mtrim=apply(data[lista,,drop=FALSE],2,mean,na.rm=TRUE)
# mtrim=data[lista,]
   med<-fdata(med,tt,rtt,names1)
   mtrim<-fdata(mtrim,tt,rtt,names2)
   rownames(med$data)<-"RPd.med"
   rownames(mtrim$data)<-tr
   ind1<-!is.nan(dep)
   dep[is.nan(dep)]=NA
   cgray=1-(dep-min(dep,na.rm=TRUE))/(max(dep,na.rm=TRUE)-min(dep,na.rm=TRUE))
   plot(lfdata[[idat]][ind1, ], col =  gray(cgray[ind1]),lty=1, main = paste(nam1[idat],"curves with RPp Depth",sep=""))
   lines(mtrim,lwd=2,col="yellow")
   lines(med,col="red",lwd=2)
   legend("topleft",legend=c(tr,"Median"),lwd=2,col=c("yellow","red"),box.col=0)
 }
}
return(invisible(list(median = med, lmed = k, mtrim = mtrim,
        ltrim = lista, dep = dep,proj = z,dfunc=dfunc,par.dfunc=par.dfunc)))
}  
