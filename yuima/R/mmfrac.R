##estimate theta2 by LSE

mmfrac <- function(yuima,...){ 	

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
    }

    if (length(yuima@model@parameter@drift)!=1){
            yuima.stop("the drift is malformed.")
    }
    

    dx <- D(yuima@model@drift,yuima@model@state.variable)
    dbx <- D(yuima@model@diffusion[[1]],yuima@model@state.variable)==0
    dbt <- D(yuima@model@diffusion[[1]],yuima@model@time.variable)==0
    bxx <- D(dx,yuima@model@state.variable)==0
    bmix <- D(dx,yuima@model@time.variable)==0
    
    isfOU <- (bxx || bmix) && (dbx && dbt)

    if (!isfOU){yuima.stop("estimation not available for this model")}
    
	process<-yuima@data@zoo.data[[1]]
    

    T<-yuima@sampling@Terminal
	est<-qgv(yuima,...)

	H<-est$coeff[1]
	sigma<-est$coeff[2]
	
    if ((!is.na(H))&(!is.na(sigma))){
    
        theta<-(2*mean(process^2)/(sigma^2*gamma(2*H+1)))^(-1/2/H)
        sH<-sqrt((4*H-1)*(1+(gamma(1-4*H)*gamma(4*H-1))/(gamma(2-2*H)*gamma(2*H))))
        sdtheta<- sqrt(theta)*(sH/(2*H))/sqrt(T)
        
	}else{
        yuima.warn("Diffusion estimation not available, can not estimate the drift parameter")
    }
        
	x<-c(est$coeff,theta)
    names(x)<-c(names(est$coeff),yuima@model@parameter@drift)  
    sdx<-matrix(,3,3)
    diag(sdx)<-c(diag(est$vcov),sdtheta)
    colnames(sdx)<-names(x)
    rownames(sdx)<-names(x)
    
    obj <- list(coefficients=x,vcov=sdx,call=call)
    class(obj) <- "mmfrac"
    return(obj)

}

print.mmfrac<-function(x,...){
tmp <- rbind(x$coefficients, diag(x$vcov))
rownames(tmp) <- c("Estimate", "Std. Error")
cat("\nFractional OU estimation\n")
print(tmp)
}




