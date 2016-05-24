EstCRMitem <-
function(data,
         max.item,
         min.item,
         max.EMCycle=500,
         converge=.01,
         type="Shojima",BFGS=TRUE) {  #Start Main function

n=ncol(data)  
N=nrow(data) 

if(is.data.frame(data)==FALSE) 
  stop("The input response data is not a data frame.
        Please use as.data.frame() and convert your response data to a data frame object before the analysis")

if(length(max.item)!=length(min.item)) 
  stop("The length of max.item vector is not equal to the length of min.item vector. Please check your inputs")

if(dim(data)[2]!=length(max.item)) 
  stop("The number of columns in the data is not equal to the length of max.item vector")

if(dim(data)[2]!=length(min.item)) 
  stop("The number of columns in the data is not equal to the length of min.item vector")

for(i in 1:n) {
  if(max(na.omit(data[,i]))> max.item[i]) 
    stop("The column ",i," has values higher than the maximum available score in the
          user specified max.item vector. Please check and clean your data.")
}

for(i in 1:n) {
  if(min(na.omit(data[,i]))< min.item[i]) 
    stop("The column ",i," has values smaller than the minimum available score in the
          user specified min.item vector. Please check and clean your data.")
}

if(is.numeric(data[,i])!=TRUE) 
  stop("The column vectors are not numeric. Please check your data")

desc <- as.data.frame(matrix(nrow=n,ncol=4))
colnames(desc) <- c("Mean","SD","Min","Max")
for(i in 1:n) {
  desc[i,1]=mean(data[,i],na.rm=TRUE)
  desc[i,2]=sd(data[,i],na.rm=TRUE)
  desc[i,3]=min(data[,i],na.rm=TRUE)
  desc[i,4]=max(data[,i],na.rm=TRUE)
}

for(i in 1:n){
  data[,i]= data[,i]-min.item[i]
}

max.item <- max.item-min.item
min.item <- rep(0,n)

for(i in 1:n) {
  if(length(which(data[,i]==max.item[i]))!=0) {
    data[which(data[,i]==max.item[i]),i]=max.item[i]-.01
  }
    if(length(which(data[,i]==0))!=0) {
      data[which(data[,i]==0),i]=.01
    }
}

data.original <- data

for(i in 1:n){data[,i]= log(data[,i]/(max.item[i]-data[,i]))}

desc$Mean.of.z <- NA
desc$SD.of.z <- NA
for(i in 1:n) {
desc[i,5]=mean(data[,i],na.rm=TRUE)
desc[i,6]=sd(data[,i],na.rm=TRUE)
}
rownames(desc) <- colnames(data.original) 


if(type=="Shojima") {

	loglikelihood <- function(data,ipar,mu,sigma) { 

		first.term <- N*sum(log(ipar[,1])+log(ipar[,3]))
		sec.term <- sum(rowSums(t(matrix(ipar[,1]^2,nrow=n,ncol=N))*
				(((t(matrix(ipar[,2],nrow=n,ncol=N))+(data*t(matrix(ipar[,3],nrow=n,ncol=N)))-matrix(rep(mu,n),nrow=N,ncol=n))^2)+
				sigma),na.rm=TRUE),na.rm=TRUE)/2
				first.term-sec.term - ((N*n/2)*log(2*pi))
	}

	estEM <- function(data,ipar) { #start internal function 2 

		sigma = 1/(sum(ipar[,1]^2)+1) 
		mu <-sigma*rowSums(t(matrix(ipar[,1]^2,ncol=N,nrow=n))*(t(matrix(ipar[,3],ncol=N,nrow=n))*data+t(matrix(ipar[,2],ncol=N,nrow=n))),na.rm=TRUE) # Equation 20 in Shojima's paper

		mumean <- mean(mu,na.rm=TRUE)
		muvar <- var(mu,na.rm=TRUE)
		zijmeanlist <- as.vector(colMeans(data,na.rm=TRUE))
		zijvarlist <- as.vector(diag(cov(data,use="pairwise.complete.obs")))
		zijmucovlist <- as.vector(cov(data,mu,use="pairwise.complete.obs"))

		gamma <- (muvar+sigma)/zijmucovlist  
		beta  <- mumean-gamma*zijmeanlist
		alpha <- 1/sqrt(gamma^2*zijvarlist-gamma*zijmucovlist)

		ipar <- cbind(alpha,beta,gamma)
		colnames(ipar) <- c("a","b","alpha")
		rownames(ipar) <- colnames(data)

		list(ipar,loglikelihood(data,ipar,mu,sigma),mu,sigma) 
	}

	
	l <- c()
	ipars <- vector("list",max.EMCycle)
      mus   <- vector("list",max.EMCycle)
	sig   <- c()

	d=1  
	iter <- 1  

	ipars[[iter]] <- t(matrix(nrow=3,ncol=n,c(1,0,1)))  
	ipars[[iter]][,2] <- -colMeans(data,na.rm=TRUE)
	mus[[iter]] <- rep(0,N)         
      sig[iter] <- 1
      l[iter]=loglikelihood(data,ipars[[iter]],mus[[iter]],sig[iter])

	while(abs(d)>converge && iter < max.EMCycle){

		est1 <- estEM(data,ipars[[iter]])
			ipars[[iter+1]] <- est1[[1]]
			l[iter+1] <- est1[[2]]
			d <- l[iter+1]-l[iter]
                  mus[[iter+1]] <- est1[[3]]
                  sig[iter+1] <- est1[[4]]
			iter <- iter+1
      }

	itempar <- vector("list",length(l))
		for(i in 1:length(l)) itempar[[i]]=ipars[[i]]
		for(i in 1:length(l)) itempar[[i]][,3]=1/ipars[[i]][,3]

	maximums <- vector("list",length(l))
	start <- t(matrix(nrow=3,ncol=n,c(1,0,1)))  
	start[,2] <- -colMeans(data,na.rm=TRUE)
	maximums[[1]]<- cbind(max(abs(itempar[[1]]-start)[,1]),max(abs(itempar[[1]]-start)[,2]),max(abs(itempar[[1]]-start)[,3]))
	for(i in 1:(length(l)-1)){
		maximums[[i+1]]=cbind(max(abs(itempar[[i+1]]-itempar[[i]])[,1]),
		max(abs(itempar[[i+1]]-itempar[[i]])[,2]),
		max(abs(itempar[[i+1]]-itempar[[i]])[,3]))
	}

name <- c()
name[1] <- paste0("EMCycle",1,"   Starting Parameters")
for(i in 2:length(l)) name[i]=paste("EMCycle",i,"   Largest Parameter Changes=",
	round(maximums[[i]],3)[1]," ",round(maximums[[i]],3)[2]," ",round(maximums[[i]],3)[3],
	sep="")
names(itempar) <- name

dif <- abs(l[length(l)]-l[length(l)-1])

	ipar.est <- itempar[[iter]]
	se.matrix <- matrix(nrow=n,ncol=3)

	for(j in 1:n) {

	 a   = ipar.est[j,1]
	 b   = ipar.est[j,2]
	 alp = ipar.est[j,3]

       s11 <- (N/a^2)+sum((data[,j]/alp+b-mus[[iter]])^2+sig[iter])
       s12 <-  2*a + sum(data[,j]/alp+b-mus[[iter]])
 	 s13 <-  2*a + sum(data[,j]*(data[,j]/alp+b-mus[[iter]]))
       s22 <-  N*a^2
	 s23 <- a^2*sum(data[,j])
       s33 <- N*alp + a^2*sum(data[,j]^2)


	 Hess <- matrix(c(s11,s12,s13,s12,s22,s23,s13,s23,s33),3,3,byrow=FALSE)
 	 se.matrix[j,] = sqrt(diag(solve(Hess)))
     }
}


if(type=="Wang&Zeng") {

  loglikelihood2 <- function(ipar,data,mu,sigma) { 

	ipar <- matrix(ipar,nrow=ncol(data),ncol=3)
	first.term <- N*sum(log(ipar[,1])+log(ipar[,3]))
	sec.term <- sum(rowSums(t(matrix(ipar[,1]^2,nrow=n,ncol=N))*
				(((t(matrix(ipar[,2],nrow=n,ncol=N))+(data*t(matrix(ipar[,3],nrow=n,ncol=N)))-matrix(rep(mu,n),nrow=N,ncol=n))^2)+sigma),
                            na.rm=TRUE),
                            na.rm=TRUE)/2
				first.term-sec.term - ((N*n/2)*log(2*pi))
	}

  grad <- function(ipar,data,mu,sigma) {
     
    ipar <- matrix(ipar,nrow=ncol(data),ncol=3)
    l0 = loglikelihood2(ipar,data,mu,sigma)
    g <- c()
    s=.00001
    for(i in 1:length(ipar)) {
      hold = ipar[i]
      h = s*hold
      ipar[i]=hold+h
      lj = loglikelihood2(ipar,data,mu,sigma)
      g[i]=(lj-l0)/h
      ipar[i]=hold
    }
    return(g)
  }
  
  hess <- function(ipar,data,mu,sigma){
   
    ipar <- matrix(ipar,nrow=ncol(data),ncol=3)
    
    s=.00001  
    g0 = grad(ipar,data,mu,sigma)
    fh = matrix(nrow=length(ipar),ncol=length(ipar))
    for(i in 1:length(ipar)){      
      hold = ipar[i]
      h = s*hold
      ipar[i]=hold+h
      gj=grad(ipar,data,mu,sigma)
      fh[,i]=(gj-g0)/h
      ipar[i]=hold
    }

    hh = (fh + t(fh))/2
    diag(hh)= diag(fh)
    return(hh)       
    
  }
  
  ipars    <- vector("list",max.EMCycle)
  gradient <- vector("list",max.EMCycle)
  H        <- vector("list",max.EMCycle)
  IH       <- vector("list",max.EMCycle)
  l        <- c()
  maxg     <- c()
  
    ipar <- t(matrix(nrow=3,ncol=n,c(1,0,1)))  
    ipar[,2] <- -colMeans(data,na.rm=TRUE)
    ipars[[1]] <- as.vector(ipar)

    # sigma1 = 1/(sum(ipars[[1]][1:4]^2)+1) 
    # mus <-sigma1*rowSums(t(matrix(ipars[[1]][1:4]^2,ncol=N,nrow=n))*((data*t(matrix(ipars[[1]][9:12],ncol=N,nrow=n)))+t(matrix(ipars[[1]][5:8],ncol=N,nrow=n))),na.rm=TRUE) # Equation 20 in Shojima's paper
    sigma1=1
    mus = rep(0,N)
    iter=1
    d=1
    l[iter]          <- loglikelihood2(ipars[[iter]],data,mu=mus,sigma=sigma1)
    gradient[[iter]] <- grad(ipars[[iter]],data,mu=mus,sigma=sigma1)
    H[[iter]]        <- hess(ipars[[iter]],data,mu=mus,sigma=sigma1)
    IH[[iter]]       <- solve(H[[iter]])
    maxg[iter]       <- abs(max(gradient[[iter]]))
  
      while(iter < max.EMCycle & abs(d)>converge) {   
    
        pos = seq(1,3*n,by=n)
        sigma1 = 1/(sum(ipars[[iter]][1:(pos[2]-1)]^2)+1) 
	  mus <-sigma1*rowSums(t(matrix(ipars[[iter]][1:(pos[2]-1)]^2,ncol=N,nrow=n))*((data*t(matrix(ipars[[iter]][(pos[3]):(3*n)],ncol=N,nrow=n)))+t(matrix(ipars[[iter]][(pos[2]):(pos[3]-1)],ncol=N,nrow=n))),na.rm=TRUE)

        iter <- iter+1
        step.size <- 1/(1:100)
        step=1
        ipars[[iter]]= ipars[[iter-1]]-(solve(H[[iter-1]])%*%as.matrix(gradient[[iter-1]]*step.size[step]))
        
        a     <- ipars[[iter]][1:(length(ipars[[iter]])/3)]
        b     <- ipars[[iter]][((length(ipars[[iter]])/3)+1):(2*(length(ipars[[iter]])/3))]
        alpha <- ipars[[iter]][(2*(length(ipars[[iter]])/3)+1):length(ipars[[iter]])]
        
        
          while((sum(a<0)!=0 | sum(alpha<0)!=0) & step <100) {
            step=step+1
            ipars[[iter]]= ipars[[iter-1]]-(solve(H[[iter-1]])%*%as.matrix(gradient[[iter-1]]*step.size[step])) 
             a     <- ipars[[iter]][1:(length(ipars[[iter]])/3)]
             b     <- ipars[[iter]][((length(ipars[[iter]])/3)+1):(2*(length(ipars[[iter]])/3))]
             alpha <- ipars[[iter]][(2*(length(ipars[[iter]])/3)+1):length(ipars[[iter]])]
          }

        l[iter] <- loglikelihood2(ipars[[iter]],data,mu=mus,sigma=sigma1)
        d <- l[iter]-l[iter-1]
        gradient[[iter]] <- grad(ipars[[iter]],data,mu=mus,sigma=sigma1)
        H[[iter]]        <- hess(ipars[[iter]],data,mu=mus,sigma=sigma1)
        maxg[iter]       <- max(abs(gradient[[iter]]))
        gradiff  = gradient[[iter]] - gradient[[iter-1]]
        pardiff  = ipars[[iter]]- ipars[[iter-1]]
       
        if(BFGS==FALSE) { 
           H[[iter]] <- hess(ipars[[iter]],data,mu=mus,sigma=sigma1) 
        } else {
         	Hstep  = ((H[[iter-1]]%*%pardiff%*%t(pardiff)%*%H[[iter-1]])/as.numeric((t(pardiff)%*%H[[iter-1]]%*%pardiff))) -
			    ((gradiff%*%t(gradiff))/as.numeric(t(pardiff)%*%gradiff))
            H[[iter]]= H[[iter-1]] - Hstep
        }
  }


	itempar <- vector("list",length(l))
		for(i in 1:length(l)) itempar[[i]]=matrix(ipars[[i]],nrow=ncol(data),ncol=3)
		for(i in 1:length(l)) itempar[[i]][,3]=1/itempar[[i]][,3]

	maximums <- vector("list",length(l))
	start <- t(matrix(nrow=3,ncol=n,c(1,0,1)))  
	start[,2] <- -colMeans(data,na.rm=TRUE)
	maximums[[1]]<- cbind(max(abs(itempar[[1]]-start)[,1]),max(abs(itempar[[1]]-start)[,2]),max(abs(itempar[[1]]-start)[,3]))
	for(i in 1:(length(l)-1)){
		maximums[[i+1]]=cbind(max(abs(itempar[[i+1]]-itempar[[i]])[,1]),
		max(abs(itempar[[i+1]]-itempar[[i]])[,2]),
		max(abs(itempar[[i+1]]-itempar[[i]])[,3]))
	}

name <- c()
name[1] <- paste0("EMCycle",1,"   Starting Parameters")
for(i in 2:length(l)) name[i]=paste("EMCycle",i,"   Largest Parameter Changes=",
	round(maximums[[i]],3)[1]," ",round(maximums[[i]],3)[2]," ",round(maximums[[i]],3)[3],
	sep="")
names(itempar) <- name

dif <- abs(l[length(l)]-l[length(l)-1])

	ipar.est <- itempar[[iter]]
	hessian <- hess(ipars[[iter]],data,mu=mus,sigma=sigma1) 
      se.matrix <- matrix(sqrt(diag(solve(-hessian))),nrow=n,ncol=3)
}

out <- list(data=data.original,
            descriptive=desc,
            param=itempar[[iter]],
            iterations=itempar,
            std.err = se.matrix,
            LL = l[length(l)],
            dif=dif)
class(out) <- "CRM"
return(out)
}

