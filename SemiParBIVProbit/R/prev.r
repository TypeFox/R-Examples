prev <- function(x, sw = NULL, type = "simultaneous", ind = NULL, delta = FALSE, n.sim = 100, prob.lev = 0.05, 
                      hd.plot = FALSE, main = "Histogram and Kernel Density of Simulated Prevalences", 
       xlab="Simulated Prevalences", ...){


if(x$Cont == "YES") stop("This function is not suitable for bivariate models with continuous margins.")
if(x$Cont == "NO" && x$VC$ccss == "yes" ) stop("This function is not suitable for bivariate selection models with continuous margin.")

lb <- wm <- ub <- qz <- sv <- Vv <- G <- X2sg <- 1
wms <- NA

if((x$Model=="B" || x$Model=="BPO") && x$triv == FALSE) stop("This function is suitable for sample selection models only.")

if(x$triv == TRUE && x$Model != "TSS") stop("This function is suitable for sample selection models only.")


if( !( type %in% c("naive","univariate","simultaneous") ) ) stop("Error in parameter type value. It should be one of: naive, univariate or simultaneous.")


if(x$Model == "BSS") Xsg <- x$X2s
if(x$Model == "TSS") Xsg <- x$X3s


if(type == "univariate"){

if(x$Model == "BSS"){
        etasg <- Xsg%*%coef(x$gam2)
        Vv    <- x$gam2$Vp
                    }

if(x$Model == "TSS"){
        etasg <- Xsg%*%coef(x$gam3)
        Vv    <- x$gam3$Vp
                    }

}  



if(type == "simultaneous" || type == "naive"){ # naive useful for sw below

if(x$Model == "BSS") etasg <- x$eta2
if(x$Model == "TSS") etasg <- x$eta3

Vv <- x$Vb

}




if(!is.null(sw)) { if( length(sw)!=length(etasg) ) stop("sw must have the same length as the number of observations used in fitting.")  }
if(is.null(sw)) sw <- rep(1,length(etasg)) 


if( !is.null(ind) ){ 

    if(is.logical(ind) == FALSE) stop("ind must be a logical variable.")
    if(length(table(ind))!=2 ) stop("ind must be a logical binary variable.")
    if( length(ind)!=length(etasg) ) stop("ind must have the same length as the number of observations used in fitting.")   

    etasg <- etasg[ind]
    Xsg <- Xsg[ind,]
    
    sw <- sw[ind]     

}


#######

wm <- weighted.mean(probm(etasg, x$margins[2])$pr, w=sw)

#######



if(type != "naive"){




if(delta == TRUE){

core <- colWeightedMeans( c( probm(etasg, x$margins[2], only.pr = FALSE)$d.n )*Xsg, w = sw, na.rm = FALSE) 


if(x$Model == "BSS" && type == "simultaneous"){

if( is.null(x$X3) ) zerod <- 0
if(!is.null(x$X3) ) zerod <- rep(0, x$X3.d2)

G <- c( rep(0,x$X1.d2), core, zerod) 

}


if(x$Model == "TSS" && type == "simultaneous"){

zerod <- c(0,0,0)
G <- c( rep(0, (x$X1.d2 + x$X2.d2)), core, zerod) 

}


if(type == "univariate")  G <- c( core )  

 
  sv <- sqrt( t(G)%*%Vv%*%G ) 
  
  qz <- qnorm(prob.lev/2, lower.tail = FALSE)

  lb <- wm - qz*sv 
  ub <- wm + qz*sv 

}







if(delta == FALSE){

  if(type == "simultaneous") coefm <- x$coefficients    
  
  
  if(type == "univariate"){ 
  
    if(x$Model == "BSS") coefm <- x$gam2$coefficients
    if(x$Model == "TSS") coefm <- x$gam3$coefficients
  
  }
  
  
   bs <- rMVN(n.sim, mean = coefm, sigma=Vv)
  
  
  if(type == "simultaneous"){
  
    if(x$Model == "BSS") bs <- bs[, x$X1.d2 + (1 : x$X2.d2) ]
    if(x$Model == "TSS") bs <- bs[, x$X1.d2 + x$X2.d2 + (1 : x$X3.d2) ]

  }
  
  
  
  
  
 
  ps  <- probm( Xsg%*%t(bs) , x$margins[2])$pr 
  wms <- colWeightedMeans( ps, w = sw, na.rm = FALSE)
  bb <- quantile(wms, probs = c(prob.lev/2,1-prob.lev/2), na.rm=TRUE )

  lb <- bb[1]
  ub <- bb[2] 
  
  if(hd.plot == TRUE){
  
  hist(wms*100, freq = FALSE, main=main, 
       xlab=xlab, 
       ylim=c(0,max(density(wms*100)$y,hist(wms*100, plot = FALSE)$density)), ...)
  lines(density(wms*100))

                     }

}






} # end type







if( type == "naive"){ 

if(x$Model == "BSS") inde <- x$inde
if(x$Model == "TSS") inde <- x$inde2
sw <- sw[inde]
if(x$Model == "BSS") resp <- x$y2
if(x$Model == "TSS") resp <- x$y3[inde] # because it is done a bit inefficiently at the moment

if( !is.null(ind) ){ 

if( length(ind) != length(resp) ) stop("ind must have the same length as the number of selected observations for the outcome.")  

resp <- resp[ind]
sw   <- sw[ind]
                    }


qz <- qnorm(prob.lev/2, lower.tail = FALSE)
wm <- weighted.mean(resp, w = sw)
sv <- sqrt( (wm*(1 - wm))/length(resp) )
lb <- wm - qz*sv 
ub <- wm + qz*sv 
  
}
  










  res <- c(lb, wm, ub)
  

  rm(lb,wm,ub,qz,sv,Vv,G,Xsg)

  out <- list(res=res, prob.lev=prob.lev, sim.prev = wms)
 
  class(out) <- "prev"

  out

}

























