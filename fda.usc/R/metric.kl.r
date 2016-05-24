metric.kl=function(fdata1, fdata2 = NULL,symm=TRUE, base=exp(1),eps=1e-10,...) 
{                                                                                
#verbose=TRUE
    C1 <- match.call()
    same <- FALSE
    if (is.fdata(fdata1)) {
        fdat <- TRUE
        tt <- tt1 <- fdata1[["argvals"]]
        rtt <- fdata1[["rangeval"]]
        nas1 <- apply(fdata1$data, 1, count.na)
        if (any(nas1)) 
            warning("fdata1 contain ", sum(nas1), " curves with some NA value \n")
        if (is.null(fdata2)) {
            fdata2 <- fdata1
            same <- TRUE
        }
        else if (!is.fdata(fdata2)) {
            fdata2 <- fdata(fdata2, tt1, rtt, fdata1$names)
        }
        nas2 <- apply(fdata2$data, 1, count.na)
        if (any(nas2)) {
            warning("fdata2 contain ", sum(nas2), " curves with some NA value \n")
        }
        DATA1 <- fdata1[["data"]]
        DATA2 <- fdata2[["data"]]
#print(DATA1[,1:3 ])        
#print(DATA2[,1:3 ])
        tt2 <- fdata2[["argvals"]]
        if (!same) {
            if (sum(tt1 != tt2) != 0) {
                stop("Error: different discretization points in the input data.\n")
            }
        }
    }
    else {
        fdat <- FALSE
        DATA1 <- fdata1
        if (is.null(fdata2)) {
            fdata2 <- fdata1
            same <- TRUE
        }
        DATA2 <- fdata2
    }
    numgr = nrow(DATA1)
    numgr2 = nrow(DATA2)
    np <- ncol(DATA1)   
#eps2<-eps
testfordim <- sum(dim(DATA1) == dim(DATA2)) == 2
	etiq1=rownames(DATA1)
	etiq2=rownames(DATA2)
twodatasets <- TRUE
if (testfordim)   twodatasets <- (all(DATA1== DATA2)) 


#if (verbose) {print(testfordim);print(twodatasets);print(symm)}

        eps2 <- as.double(.Machine[[1]] * 100)
        inf <- 1 - eps2                                      
        sup <- 1 + eps2                                      
        zz<-apply(DATA1,1,sum)
        if (all(zz > inf) & all(zz < sup)) {                
            densi=TRUE
        }                                                     
        else {
            densi = FALSE                                     
          warnings("Ponderation fdata1/sum(fdata1) is done")
          DATA1<-DATA1/apply(DATA1,1,sum)  #calcular integral
          }
        zz<-apply(DATA2,1,sum)          
         if (all(zz > inf) & all(zz < sup)) {   densi=TRUE        }                                                                       
         else {                                                                  
           densi = FALSE                                                       
           warnings("Ponderation fdata2/sum(fdata2) is done")       
           DATA2<-DATA2/apply(DATA2,1,sum)  #calcular integral                                     
           }                                                                        

     #& symm
#if (verbose) {print(testfordim)        ;print(twodatasets)}

   if (fdat) {
        dtt <- diff(tt)
        inf <- dtt - eps2
        sup <- dtt + eps2
        if (all(dtt > inf) & all(dtt < sup)) {
            equi = TRUE
        }
        else equi = FALSE
        mdist = array(0, dim = c(numgr, numgr2))
        predi <- TRUE
        if (testfordim) {
            if (twodatasets) {  
                predi <- FALSE
#if (verbose) {print(" no predi")  ;print("symm")  ;print(symm)  }
if (symm){              
                for (i in 1:numgr) {
                  ii = i + 1           
                  for (ii in i:numgr2) {
a4<-DATA1[i, ]  
a44<-DATA2[ii, ]
#      a4[a4<eps & a44<eps]<-min(eps2,min(a4[a4>0]))              
#                    a44[a4<eps & a44<eps]<-min(eps2,min(a44[a44>0]))  
if (any((a4<eps & a44>(eps+a4))|(a4>(eps+a44) & a44<eps)))   {
#          print("1")
          mdist[i, ii]<-Inf }
else  {
                    a4[a4<eps & a44<eps]<-min(eps2,min(a4[a4>0]))              
                    a44[a4<eps & a44<eps]<-min(eps2,min(a44[a44>0]))              
                    a4<-a4/sum(a4)                  
                    a44<-a44/sum(a44)
                    f = a4*log(a4/a44,base=base)
                    d1 = (int.simpson2(tt, f, equi))
                    f = a44*log(a44/a4,base=base)
                    d2= (int.simpson2(tt, f, equi))
                    mdist[i, ii] <-(d1+d2)/2    
  }                      
                  }
#print(                   mdist[i,] )
                }
              mdist = t(mdist) + mdist
                }
            else{
                for (i in 1:numgr) {
                  ii = i + 1
                  for (ii in 1:numgr2) {
a4<-DATA1[i, ]  
a44<-DATA2[ii, ] 
if (any((a4<eps & a44>(eps+a4))|(a4>(eps+a44) & a44<eps)))   {
#print("2")
mdist[i, ii]<-Inf }
else  { 
             a4[a4<eps & a44<eps]<-min(eps2,min(a4[a4>0]))              
             a44[a4<eps & a44<eps]<-min(eps2,min(a44[a44>0]))  
             a4<-a4/sum(a4)                  
             a44<-a44/sum(a44)                    
                    f <- a4*log(a4/a44,base=base)
                    mdist[i, ii]=(int.simpson2(tt, f, equi))                                     
#print( mdist[i, ii])                    
                  
                  }
#print(                   mdist[i,] )                  
                }
              }   }             

                
            }
        }
        if (predi) {
           if (symm){
            for (i in 1:numgr) {
                for (ii in 1:numgr2) {
a4<-DATA1[i, ]  
a44<-DATA2[ii, ]     
if (any((a4<eps & a44>(eps+a4))|(a4>(eps+a44) & a44<eps)))   {
#print(3)
mdist[i, ii]<-Inf }
else  { 
             a4[a4<eps & a44<eps]<-min(eps2,min(a4[a4>0]))              
             a44[a4<eps & a44<eps]<-min(eps2,min(a44[a44>0]))  
             a4<-a4/sum(a4)                  
             a44<-a44/sum(a44)                    
                    f = a4*log(a4/a44,base=base)
                    d1 = (int.simpson2(tt, f, equi))
                     f = a44*log(a44/a4,base=base)
                    d2= (int.simpson2(tt, f, equi))
                    mdist[i, ii] <-(d1+d2)/2
                }
            }
                        }}
           else{     
            for (i in 1:numgr) {
                for (ii in 1:numgr2) {
a4<-DATA1[i, ]  
a44<-DATA2[ii, ] 
if (any((a4<eps & a44>(eps+a4))|(a4>(eps+a44) & a44<eps)))   {
#print(4)
mdist[i, ii]<-Inf }
else  {    
                    a4[a4<eps & a44<eps]<-min(eps2,min(a4[a4>0]))              
                    a44[a4<eps & a44<eps]<-min(eps2,min(a44[a44>0]))  
 f = a4*log(a4/a44,base=base)
                    mdist[i, ii] = (int.simpson2(tt, f, equi))
}                     
                }
            }
          }
        }       
    }
    else {
   stop("No fdata class object")
     }
	rownames(mdist)<-etiq1
	colnames(mdist)<-etiq2     
    attr(mdist, "call") <- "metric.kl"
    attr(mdist, "par.metric") <- list( base=base,symm=symm)
    return(mdist)
}

