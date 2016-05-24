HubNetwork <-
function(p,sparsity,hubnumber,hubsparsity,type="Gaussian"){

if(sparsity<0 || sparsity>1) stop("sparsity has to take value between 0 and 1!")

if(hubsparsity<0 || hubsparsity>1) stop("hubsparsity has to take value between 0 and 1!")

if(type!="Gaussian" && type!="covariance" && type!="binary") stop("type can only takes argument Gaussian, covariance, or binary!")

if(hubnumber>p) stop("hubnumber cannot be larger than the number of features!")

if(hubnumber<0) stop("hubnumber cannot be negative!")

# Generate an Erdos Renyi type network with positive and negative entries  
  sparse <- rbinom(p*p,1,1-sparsity)*sample(c(-1,1),p*p,replace=TRUE)*runif(p*p,0.25,0.75)
  Theta <- matrix(data=sparse,nrow=p,ncol=p)
  Theta[lower.tri(Theta,diag=FALSE)] <- 0
  Theta <- Theta+t(Theta)
  
# Add in Hub Nodes and make the matrix symmetric  
  hubcol <- sample(1:p,hubnumber,replace=FALSE)
  Theta[,hubcol] <- rbinom(hubnumber*p,1,1-hubsparsity)*sample(c(-1,1),hubnumber*p,replace=TRUE)*runif(hubnumber*p,0.25,0.75)
  Theta <- (Theta+t(Theta))/2
  
  if(type=="binary"){ 
  	diag(Theta) <- sample(c(-1,1),p,replace=TRUE)*runif(p,0.25,0.75)
  	return(list(Theta=Theta,hubcol=hubcol))
  	}
  
# Make the matrix positive definite  
    diag(Theta) <- 0
    ee <- min(eigen(Theta,only.values=T)$values)
    diag(Theta) <- ifelse(ee < 0, -ee + 0.1, 0.1)

  if(type=="covariance"){
    Theta <- cov2cor(Theta)
  }

  return(list(Theta=Theta,hubcol=hubcol)) 
}
