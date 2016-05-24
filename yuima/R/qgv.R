##:: function qgv
##:: Estimating the local Holder exponent of the path and the constant

qgv<- function(yuima, filter.type="Daubechies", order=2, a=NULL){

	call <- match.call()
	
	if(missing(yuima)){
		yuima.stop("yuima object is missing.")
	}
    
    if(class(yuima)!="yuima"){
        yuima.stop("an object of class yuima is needed.")
    }
    

    Ddiffx0 <- D(yuima@model@diffusion[[1]],yuima@model@state.variable) == 0
    Ddifft0 <- D(yuima@model@diffusion[[1]],yuima@model@time.variable) == 0
    
    bconst<-FALSE
    Hconst<-FALSE
    
    diffnoparam<-length(yuima@model@parameter@diffusion)==0
    
    if (Ddiffx0 & Ddifft0 & diffnoparam){
        
        if(eval(yuima@model@diffusion[[1]])==0){
            yuima.stop("the model has no fractional part.")
        }
        
        bconst<-TRUE
     }
    
    y<-yuima@model@drift
    dx <- D(y,yuima@model@state.variable)
    bt <- D(y,yuima@model@time.variable)==0
    bxx <- D(dx,yuima@model@state.variable)==0
    
    
    dbx <- D(yuima@model@diffusion[[1]],yuima@model@state.variable)==0
    dbt <- D(yuima@model@diffusion[[1]],yuima@model@time.variable)==0
    
    
    isfOU<-(bxx && bt) && (dbx && dbt)

    if (!is.na(yuima@model@hurst)){
        yuima.warn("No Hurst exponent will be estimated")
        Hconst<-TRUE
    }
    
    
    H <- yuima@model@hurst
	sdH<-NA
    
    if(Hconst & bconst){
        x<-c(H,eval(yuima@model@diffusion[[1]]))
        names(x)<-c("hurst","const")        
        sdx<-diag(c(NA,NA))
        colnames(sdx)<-names(x)
        rownames(sdx)<-names(x)
        sdx[2,1] <- sdx[1,2] <- NA
        obj <- list(coefficients=x,vcov=sdx,call=call)
        class(obj) <- "qgv"
        return(obj)
    }

    isregular<-yuima@sampling@regular 
    
    if (!isregular){
        yuima.stop("qgv method is only working for regular grid.")
    }
    
    if (!(order %in% 1:10)){
    yuima.warn("Classical filter implement are of order ranged in [1,10], order have been fixed to 2.")
    order=2
    }    

    if (missing(a)){

	if (filter.type=="Daubechies"){
		if (order==2){
			
		a<-1/sqrt(2)*c(.4829629131445341,
						  -.8365163037378077,
						   .2241438680420134,
						   .1294095225512603)
		}else{
			yuima.warn("There is no such order Daubechies' filter implemented, order have been fixed to 2.")	
				
		a<-1/sqrt(2)*c(.4829629131445341,
						  -.8365163037378077,
						   .2241438680420134,
						   .1294095225512603)
            order=2    
		
		
		} 
	}else if (filter.type=="Classical"){
		mesh<-0:order
		a=(-1)^(1-mesh)/2^order*choose(order,mesh)
	}else{
        yuima.warn("No such type of filter. Filter have been fixed to Daubechies' filter of order 2.")
        a<-1/sqrt(2)*c(.4829629131445341,
                            -.8365163037378077,
                            .2241438680420134,
                            .1294095225512603)
        order=2
    }


	}
	
	L<-length(a)
	a2<-rep(0,2*L)
	a2[seq(1,2*L,by=2)]<-a
	
	process<-yuima@data@zoo.data[[1]]
	N<-length(process)
	
	#Computation of the generalized quadratic variations
	
    V1<-sum(filter(process,a)^2,na.rm=TRUE)
	V2<-sum(filter(process,a2)^2,na.rm=TRUE)
    
    if(!Hconst){
        H<-1/2*log2(V2/V1)
    }    
    
    nconst <- "const"
    
    sdC<-NA
    C <- NA   
    if(isfOU){
        if( bconst||(H>=1)||(H<=0)){
            if(diffnoparam){
            C<-eval(yuima@model@diffusion[[1]])
            }    
        }else{
            #Compute the estimation of the constant C.
            delta<-yuima@sampling@delta
            constfilt<-sum(a%*%t(a)*abs(matrix(0:(L-1),L,L)-matrix(0:(L-1),L,L,byrow=TRUE))^(2*H))
            C<- sqrt(-2*V1/(N-L)/(constfilt*delta^(2*H)))
            nconst<-as.character(yuima@model@diffusion[[1]])  
        }
    } else {
      nconst<-as.character(yuima@model@diffusion[[1]]) 
    }
    
    
    if(isfOU){
        
        #Compute the standard error
        infty<-100
        
        C11<-rep(0,2*infty+1)
        C11bis<-rep(0,2*infty+1)
        C12<-rep(0,2*infty+1)
        C22<-rep(0,2*infty+1)
        l<-order+1
        
        for (i in (-infty:infty)){
            
            for (q in 0:l){
                for (r in 0:l){
                    C11[i+infty+1]<-C11[i+infty+1]+a[q+1]*a[r+1]*abs(q-r+i)^(2*H)	
                }	
            }
            
            for (q in 0:l){
                for (r in 0:l){
                    C12[i+infty+1]<-C12[i+infty+1]+a[q+1]*a[r+1]*abs(2*q-r+i)^(2*H)	
                }	
            }
            
            for (q in 0:l){
                for (r in 0:l){
                    C22[i+infty+1]<-C22[i+infty+1]+a[q+1]*a[r+1]*abs(2*q-2*r+i)^(2*H)	
                }	
            }
            
            for (q in 0:(2*l)){
                for (r in 0:(2*l)){
                    C11bis[i+infty+1]<-C11bis[i+infty+1]+a2[q+1]*a2[r+1]*abs(q-r+i)^(2*H)	
                }	
            }
            
        }    
        
        rho11<-1/2*sum((C11/C11[infty+1])^2)
        rho11dil<-1/2*sum((C11bis/C11bis[infty+1])^2)
        
        rho12<-1/2*sum((C12/2^H/C11[infty+1])^2)
        rho22<-1/2*sum((C22/2^H/2^H/C11[infty+1])^2)
        
        sigma1<-1/(log(2))^2*(rho11+rho11dil-2*rho12)  
        
        if((!Hconst)&(H>0)){
        sdH<-sqrt(sigma1)/sqrt(N)
        }
            
        if ((H>0)&(H<1)&(!bconst)){
            sigma2<-C^2/2*sigma1
            sdC<-sqrt(sigma2)/sqrt(N)*log(N)
        }

    }

x<-c(H,C)
names(x)<-c("hurst",nconst)        
sdx<-diag(c(sdH,sdC))
colnames(sdx)<-names(x)
rownames(sdx)<-names(x)
sdx[2,1] <- sdx[1,2] <- NA    
obj <- list(coefficients=x,vcov=sdx,call=call)
class(obj) <- "qgv"
return(obj)

}


print.qgv<-function(x,...){
    tmp <- rbind(x$coefficients, diag(x$vcov))
    rownames(tmp) <- c("Estimate", "Std. Error")
    cat("\nFractional OU estimation\n")
    print(tmp)
}

