ivsplm <-function(formula,data=list(), index=NULL, endog = NULL, instruments= NULL, method = c("w2sls", "b2sls", "g2sls", "ec2sls"), lag = FALSE, listw = listw, effects = NULL, lag.instruments = FALSE){

# If the user do not make any choice in terms of method, when effects is Fixed the function calculates the w2sls. On the other hand, when effects is random the function calculates the ec2sls
if(length(method) !=1 && effects == "fixed") method <- "w2sls" 	
if(length(method) !=1 && effects == "random") method <- "ec2sls" 	
		
 if(!is.null(index)) {
    #require(plm)
    data <- plm.data(data, index)
    }
  
  index <- data[,1]
  tindex <- data[,2]

  names(index)<-row.names(data)
  ind <-index[which(names(index)%in%row.names(data))]
  tind<-tindex[which(names(index)%in%row.names(data))]
   spord <- order(tind, ind)
   data <-  data[spord,]


  ## record call
  cl <- match.call()

  ## check
  if(dim(data)[[1]]!=length(index)) stop("Non conformable arguments")
  
    mt <- terms(formula, data = data)
    mf <- lm(formula, data, na.action = na.fail, method = "model.frame")

    y <- model.extract(mf, "response")
    x <- model.matrix(mt, mf)
 
  N<-length(unique(ind))
  k<-dim(x)[[2]]
  T<-max(tapply(x[,1],ind,length))
  NT<-length(ind)



  balanced<-N*T==NT
if(!balanced) stop("Estimation method unavailable for unbalanced panels")

# print(listw)
#### creating the block diagonal matrix
if(lag){
I_T <- Diagonal(T)
Ws <- kronecker(I_T, listw)
}
# else Ws <- NULL

if(!lag){
	
if(is.null(endog)) stop("No endogenous variables specified. Please use plm instead of splm")

else {
	endog <- as.matrix(lm(endog, data, na.action = na.fail, method = "model.frame"))

if(!is.null(listw)){
I_T <- Diagonal(T)
Ws <- kronecker(I_T, listw)
	
}	
	}

if(is.null(instruments)) stop("No instruments specified")

else instruments <- as.matrix(lm(instruments, data, na.action = na.fail, method = "model.frame"))

}


else{

if(!is.null(endog)){

endog <- as.matrix(lm(endog, data, na.action = na.fail, method = "model.frame"))

if(is.null(instruments)) stop("No instruments specified for the additional variable")

else instruments <- as.matrix(lm(instruments, data, na.action = na.fail, method = "model.frame"))	

		
	}
	else instruments = NULL
	}



switch(method, 
w2sls = {
	result <- ivplm.w2sls(Y = y, X = x, H = instruments, endog = endog, lag = lag, listw = Ws, lag.instruments = lag.instruments, T = T, N = N, NT = NT)
	},
b2sls = {
	result <- ivplm.b2sls(Y = y,X =x, H = instruments, endog = endog, lag = lag, listw = Ws, lag.instruments = lag.instruments, T = T, N = N, NT = NT)
	},
ec2sls = {
	result <- ivplm.ec2sls(Y = y,X =x, H = instruments, endog = endog, lag = lag, listw = Ws, lag.instruments = lag.instruments, T = T, N = N, NT = NT)
	},
g2sls = {
	result <-ivplm.g2sls(Y = y,X =x, H = instruments, endog = endog, lag = lag, listw = Ws, lag.instruments = lag.instruments, T = T, N = N, NT = NT )
	},
stop("...\nUnknown method\n"))


    result$zero.policy <- FALSE
    result$robust <- FALSE
    result$legacy <- FALSE
    result$listw_style <- NULL
    result$call <- match.call()


class(result) <- "stsls"
result
}
