combat_np <-function(Y, b, covariates=NULL)
{
################################################################################

if(missing(Y)){stop("Argument 'Y' missing, with no default\n")}

if(missing(b)){stop("Argument 'b' missing, with no default\n")}

if(class(Y)!='matrix'){stop("'Y' must be of class 'matrix'\n")}

if(class(b)!='factor'){stop("'b' must be of class 'factor'\n")}

if(any(is.na(Y))){stop("NA values are not allowed in 'Y'\n")}

if(any(is.na(b))){stop("NA values are not allowed in 'b'\n")}

if(length(b)!=nrow(Y)){stop("length(b) is different from nrow(Y)\n")}

if(any(apply(Y,2,mode)!='numeric')){stop('Array expression columns contain non-numeric values!\n')}

################################################################################

if(!is.null(covariates)){
if(class(covariates)!="data.frame"){stop("'covariates' must be of class 'data.frame'\n")}

col.cov<-ncol(covariates)
for(i in 1:col.cov){if(class(covariates[,i])!="numeric" & class(covariates[,i])!="factor"){
stop("column ", i, " of 'covariates' must be of class 'factor' or 'numeric'\n")}}

if(any(is.na(covariates))){stop("NA values are not allowed in 'covariates'\n")}
}

################################################################################

    library(MASS)
    n   <- nrow(Y)
    g   <- ncol(Y)
    m1  <- nlevels(b)
    n.b <- summary(b)
    
    X            <- cbind(b,covariates)
    colnames(X)  <- paste("col", 1:ncol(X),sep = "")

    fmla         <- as.formula(paste("~-1+", paste(colnames(X),collapse = "+")))
    Xdes         <- model.matrix(fmla, X)                                       # (n x m1+m2)     m2 can be zero
    if(!is.null(covariates)){                                      
    XdesCov      <- as.matrix(Xdes[,(m1+1):ncol(Xdes)])                         # (n x m2)
    }                                    
    Xdes_batch   <- model.matrix(~-1 + b, b)                                    # (n x m1)
    n_batches    <- as.numeric(summary(b))                                      # (m1 x 1)
    mat_nbatches <- matrix(rep(n_batches,g),ncol=g)                             # (m1 x g)
    
    B.hat         <- ginv(Xdes)%*%Y                                             # (m1+m2 x g)     m2 can be zero        

    alpha.hat     <- t(n.b/n)%*%B.hat[1:m1,]                                    # (g x 1)   vector
    mat_alpha.hat <- matrix(rep(alpha.hat,n),byrow=T,nrow=n)                    # (n x g)
    if(!is.null(covariates)){
    Bcov.hat      <- B.hat[(m1+1):ncol(Xdes),]                                  # (m2 x g)
    }
        
    var.pooled <- matrix(rep(1/n,n),nrow=1)%*%((Y-Xdes%*%B.hat)^2)              # (g x 1)   vector
    mat_var    <- matrix(rep(var.pooled,n),byrow=T,nrow=n)                      # (n x g)
    
    if(!is.null(covariates)){
    Z_num  <- Y-mat_alpha.hat-XdesCov%*%Bcov.hat                                # (n x g)
    }else{
    Z_num  <- Y-mat_alpha.hat                                                   # (n x g)
    }
    Z      <- Z_num/sqrt(mat_var)                                               # (n x g)
    
#########################
    
    gamma_bg.hat <- (t(Xdes_batch)%*%Z)/mat_nbatches                            # (m1 x g)
        
    Zgamma         <- (Z-Xdes_batch%*%gamma_bg.hat)^2                           # (n x g)
    delta_bg.hat   <- (t(Xdes_batch)%*%Zgamma)/(mat_nbatches-1)                 # (m1 x g)
    
######################### 

#likelihood function used below is similar to this:
#L <-function(x,g.hat,d.hat){
#prod(dnorm(x,g.hat,d.hat))
#}

# Monte Carlo integration function to find the nonparametric adjustments
int.eprior <- function(sdat,g.hat,d.hat){
	g.star <- d.star <- NULL
	r <- nrow(sdat)
	for(i in 1:r){
		g <- g.hat[-i]
		d <- d.hat[-i]		
		x <- sdat[i,!is.na(sdat[i,])]
		n <- length(x)
		j <- numeric(n)+1
		dat <- matrix(as.numeric(x),length(g),n,byrow=T)
		resid2 <- (dat-g)^2
		sum2 <- resid2%*%j
		LH <- 1/(2*pi*d)^(n/2)*exp(-sum2/(2*d))
		LH[LH=="NaN"]=0
		g.star <- c(g.star,sum(g*LH)/sum(LH))
		d.star <- c(d.star,sum(d*LH)/sum(LH))
		#if(i%%1000==0){cat(i,'\n')}
		}
	adjust <- rbind(g.star,d.star)
	rownames(adjust) <- c("g.star","d.star")
	adjust	
	} 

#########################

#cat("Finding nonparametric adjustments\n")
      mat_gamma.star <- matrix(rep(0,m1*g),nrow=m1,ncol=g)
      mat_delta.star <- matrix(rep(0,m1*g),nrow=m1,ncol=g)
		  for (i in 1:m1){
		  ind_batch <- which(b==(levels(b)[i]))
			temp <- int.eprior(t(as.matrix(Z[ind_batch,])),gamma_bg.hat[i,],delta_bg.hat[i,])
			mat_gamma.star[i,] <- temp[1,]
			mat_delta.star[i,] <- temp[2,]
			} 

      mat_gamma.star <- as.matrix(mat_gamma.star)
      mat_delta.star <- as.matrix(mat_delta.star)
      
######################### 

    Zc_gamma.star                <- Z-Xdes_batch%*%mat_gamma.star
    Zc_gamma.star_var            <- sqrt(mat_var)*Zc_gamma.star
    Zc_gamma.star_var_delta.star <- Zc_gamma.star_var*(Xdes_batch%*%(1/sqrt(mat_delta.star))) 

    if(!is.null(covariates)){
    Y.star <- Zc_gamma.star_var_delta.star + mat_alpha.hat + XdesCov%*%Bcov.hat
    }else{
    Y.star <- Zc_gamma.star_var_delta.star + mat_alpha.hat 
    }

    return(Y.star)
}
