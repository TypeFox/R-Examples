getWV <-
function(x, myx, kernel=kernel,  alpha=NULL, bw=NULL, bc="simple", getit=TRUE){
   ## to calculate the weights and dimensions
    ##This part is common to all iterations so should not be looped to save CPU time.
    if(is.character(kernel)) {
      kernel.char=kernel
      kernel=get.weight.function(kernel)
     }else{kernel.char=attr(kernel, "name")}
    nn<-length(myx)
    n1<-length(x)
    if(!"tempout"%in% dir()) dir.create("./tempout",showWarnings = FALSE)
    fname<-paste("./tempout/output.n_", length(x),"nn",  
      length(myx),"bw", ifelse(is.null(bw), "NA", bw),  "bc_", bc, "alpha", 
      ifelse(is.null(alpha), "NA",alpha), 
      "kernel",   kernel.char, ".RData",sep="_")
    if(getit|!file.exists(fname)){
        exnnn1 <- expand.grid(1:nn, 1:n1)
        #xs1=x: old sequence of x. ats=new seq of x. calculate the weight matrix. 
       app1 <- function(ij, xs1, ats, alpha1, ww1=kernel) {
	     at<-ats[ij[1]]
	    if(is.null(bw)) h=alpha2h(alpha=alpha, at=at, xseq=x) else  ##adaptive
       if(bc=="none") h=bw else 
       if(bc=="simple") {
         if(distance(at, min(myx))<bw) h=at-min(myx)+bw/10 else 
         if(distance(at, max(myx))<bw) h=max(myx)-at+bw/10 else
         h<-bw
       } else stop("Boundary correction must be one of none and simple");     
        return(ww1((at - xs1[ij[2]])/h) )
        }      
      fout <- apply(exnnn1, 1, FUN = function(ij){  
       at<-myx[ij[1]]  
       if(is.null(bw)) h=alpha2h(alpha=alpha, at=at, xseq=x) else 
       if(bc=="none") h=bw else 
       if(bc=="simple") {
        if(distance(at, min(myx))<bw) h=at-min(myx)+bw/10 else 
        if(distance(at, max(myx))<bw) h=max(myx)-at+bw/10 else
        h<-bw
       } else   stop("Boundary correction must be one of none and simple");     
     return(kernel((at - x[ij[2]])/h) )
      })    
    ####weighting function  nn by n1. 
      w1ij <- matrix(fout, nrow = nn, ncol = n1)   
       fout <- apply(exnnn1, 1, FUN = function(ij) {
       o<-w1ij[ij[1],ij[2]]*sum(w1ij[ij[1],]*(myx[ij[1]]-x)*(x[ij[2]]-x))
       return(o)
      })
     sx <- matrix(fout, nrow = nn, ncol = n1)
     sx <- sweep(sx, 1, rowSums(sx), FUN = "/")
     TX <- sweep(sx, 1,   sqrt(rowSums(sx^2)), FUN = "/")   
     k0<-sum( sqrt(sum((TX[1:(nn-1),]-TX[2:nn,])^2)) )
     exn1n1 <- expand.grid(1:n1, 1:n1)
     fout <- apply(exn1n1, 1, FUN = function(ij, t1 = x, 
          t2 = x, t3 = alpha) app1(ij, xs1=t1, ats=t2, alpha1=t3, ww1= kernel ))
     LL <- matrix(fout, nrow = n1, ncol = n1)
     LL <- sweep(LL, 1, rowSums(LL), FUN = "/")  
     LL <- sweep(LL, 1,sqrt(rowSums(LL)), FUN = "/")
     AA <- diag(rep(1, n1)) - LL
     delta <- sum((AA)^2)
     delta2<- sum((AA%*%AA)^2)    
     vv<-delta^2/delta2 
     o<-list(x=x, myx=myx, k0=k0, delta=delta, delta2=delta2,vv=vv, TX=TX, sx=sx,
       LL=LL, AA=AA, kernel=kernel, bc=bc)
     save(o,file= fname)
    } else  o<-get(load(fname))
    return(o)     
}
