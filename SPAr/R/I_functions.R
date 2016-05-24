print.SPA <- function(x,...){
	object <- x
	cat ('\ncall:	',paste(deparse(object$call),sep="\n"),'\n\n')
	if (object$type=="env"){
		cat('I2star=',object[[1]],'\tp_value_global=',object[[2]],'\tp_value_local=',object[[3]])
	}else{
	if (object$type=="I1")
		cat('I1=',object[[1]],'\t\tp_value=',object$pvalue)
	if (object$type=="I2")
		cat('I2=',object[[1]],'\t\tp_value=',object$pvalue)	
	if (object$type=="p*")
		cat('p*=',object[[1]][3],'\t\t\tp_value of p*=',object$pvalue[3],
			'\nI1=',object[[1]][1],'\tp_value of I1=',object$pvalue[3],
			'\nI1=',object[[1]][2],'\t\tp_value of I2=',object$pvalue[2])
	}
	cat('\n\n')
			
}

SPA.I <- function(x,y,nperm = 100,type="dichotomous",interaction=0){
# The function to compute I_1 
# @ x: a numeric genotype matrix with each row as an individual
#       and each column as a separate snp, genotype should be
#       coded as 0,1,2
# @ y: a vector of phenotype, the length of y should be the same
#       as the number of rows of x
# @nperm : the number of permutations
# @type : if the phenotype is dichotomous or continous
# @interaction: indicate if interaction is taken into account;
#               0 - default value. compute p*, it returns the value of I1, I2 and p* and their p-values
#				1 - only compute I1; 2 - only compute I2
 
 	call <- match.call()
	if (is.character(type))
		type <- f.get.type(type)
		else type <- NULL
	if (is.null(type))
		stop("Cannot recognize type!")
	if (!is.numeric(y))
		stop ("The phenotype must be a vector of numerical values!")
	if (!is.matrix(x))
		stop ("The genotype must be a matrix!")
	if (nrow(x)!=length(y))
		stop ("The length of the phenotype must match the number of rows of the genotype!")
	if (sum(x!=0 & x!=1 & x!=2)>0)
		stop ("The genotype can only take values 0,1,2!")
	if (type=="D" & length(setdiff(unique(y),c(0,1)))>0)
		stop("The phenotype can only take values 0 and 1 for dichotomous traits!")
	if (type == "C" & length(setdiff(unique(y),c(0,1)))==0)
		cat("The phenoype only has two unique values, so it is a dichotomous trait not continuous!\n")
	if (interaction!=1 & interaction!=2 & interaction!= 0)
		stop("Interaction must be 1, 2 or 3!")

	if (type == "D"){
		if (interaction==1){
#			cat('Compute I_1 for dichotomous trait...\n')
			type="I1"
			ans <- f.I1.dich(x,y,nperm)
			}
		else if (interaction==2){
			type="I2"
#			cat('Compute I_2 for dichotomous trait...\n')
			ans <- f.I2.dich(x,y,nperm)
			}
			else if (interaction==0){
				type="p*"
#				cat('Compute p* for dichotomous trait...\n')
				ans <- f.pstar.dich(x,y,nperm)
			}
		}
	if (type == "C"){
		if (interaction==1 ){
			type="I1"
#			cat('Compute I_1 for continuous trait...\n')
			ans <- f.I1.cont(x,y,nperm)
			}
		else if (interaction ==2){
			type="I2"
#			cat('Compute I_2 for continuous trait...\n')
			ans <- f.I2.cont(x,y,nperm)
		}
		else if (interaction == 0){
			type="p*"
#			cat('Compute p* for continuous trait...\n')
			ans <- f.pstar.cont(x,y,nperm)
		}
	}
	ans <- c(ans,list(nperm=nperm,call=call,type=type))
	class(ans) <- "SPA"
	invisible(ans)
}

#################################################################
f.get.type <- function(type){
	if (!is.na(charmatch(type,"dichotomous"))){
		return("D")
		}else
			if (!is.na(charmatch(type,"continuous"))){
				return("C")
				} else {
					return(NULL)
					}
}

f.I1.dich <- function(x,y,nperm){
	cat('Compute I_1 for dichotomous trait...\n')
	I1 <- f.I1.dich.single(x,y)
	I1.perm <- rep(NA,nperm)
	for (j in 1:nperm){
		yperm <- sample(y)
		I1.perm[j] <- f.I1.dich.single(x,y)
	}
	p <- sum(I1.perm>=rep(I1,nperm))/nperm
	return(list(I1=I1,pvalue=p))	
}

f.I1.dich.single <- function(x,y){
	xxx.temp <- as.matrix(x)
	myindex <- which(colSums(xxx.temp)!=0)
	xxx <- as.matrix(xxx.temp[,myindex])
	x.case <- as.matrix(xxx[y==1,])
	x.control <- as.matrix(xxx[y==0,])
	n.control <- sum(y==0)
	n.case <- sum(y==1)
	pd <- colSums(x.case)/(colSums(x.case)+colSums(x.control)) 
	pexp <- n.case/(n.case+n.control)
	ni <- colSums(xxx)
	wi <- ni/sum(ni)
	return( sum(wi*ni*(pd-pexp)^2))
}

f.I1.cont <- function(x,y,nperm){
	cat('Compute I_1 for continuous trait...\n')
	I1 <- f.I1.cont.single(x,y)
	I1.perm <- rep(NA,nperm)
	for (j in 1:nperm){
		yperm <- sample(y)
		I1.perm[j] <- f.I1.cont.single(x,yperm)
	}
	p <- sum(I1.perm>=rep(I1,nperm))/nperm
	return(list(I1=I1,pvalue=p))	
}

f.I1.cont.single <- function(x,y){
	xxx.temp <- as.matrix(x)
	myindex <- which(colSums(xxx.temp)!=0)
	xxx <- as.matrix(xxx.temp[,myindex])
	ybar <- mean(y)
	yibar <- apply(xxx,2,f.yibar,y=y)
	ni <- colSums(xxx)
	return(sum(ni^2*(yibar-ybar)^2))
}

f.yibar <- function(x,y){
	return (mean(y[x>0]))
}



#######################################################
f.I2.dich <- function(x,y,nperm){
  cat('Compute I_2 for dichotomous trait...\n')
  xxx.temp <- as.matrix(x)
  
  if (ncol(x)==1)
    stop("The number of columns of the genotype must 
    	be at least 2 to compute the interaciton score!")
  else {
	myindex <- which(colSums(xxx.temp)!=0)
  	xxx <- as.matrix(xxx.temp[,myindex])
  	if (ncol(xxx)==1)
  	  stop("The number of columns of the genotype must 
    		be at least 2 to compute the interaciton score!")
	mycombn <- combn(ncol(xxx),2)
	nij <- apply(mycombn,2,function(a) sum(xxx[,a]))
	wij <- nij/sum(nij)
	Iij <- apply(mycombn,2,function(a) f.Iij.dich(xxx[,a],y))
	I2 <- sum(wij*nij*Iij)
	
	I2.perm <- rep(NA,nperm)
	for (j in 1:nperm){
		yperm <- sample(y)
		Iij <- apply(mycombn,2,function(a) f.Iij.dich(xxx[,a],yperm))
		I2.perm[j] <- sum(wij*nij*Iij)
	}
	p <- sum(I2.perm>=rep(I2,nperm))/nperm
	return(list(I2=I2,pvalue=p))	
  }
}


f.Iij.dich <- function(x,y){
		n.case <- sum(y==1)
		n.control <- sum(y==0)
		x.case <- x[y==1,]
		x.control <- x[y==0,]
		n01d <- sum((x.case[,1]==0)*(x.case[,2]==1))+2*sum((x.case[,1]==0)*(x.case[,2]==2))
		n10d <- sum((x.case[,1]==1)*(x.case[,2]==0))+2*sum((x.case[,1]==2)*(x.case[,2]==0))
		n11d <- 2*sum((x.case[,1]==1)*(x.case[,2]==1))+3*sum((x.case[,1]==1)*(x.case[,2]==2))+3*sum((x.case[,1]==2)*(x.case[,2]==1))+4*sum((x.case[,1]==2)*(x.case[,2]==2))
		n01c <- sum((x.control[,1]==0)*(x.control[,2]==1))+2*sum((x.control[,1]==0)*(x.control[,2]==2))
		n10c <- sum((x.control[,1]==1)*(x.control[,2]==0))+2*sum((x.control[,1]==2)*(x.control[,2]==0))
		n11c <- 2*sum((x.control[,1]==1)*(x.control[,2]==1))+3*sum((x.control[,1]==1)*(x.control[,2]==2))+3*sum((x.control[,1]==2)*(x.control[,2]==1))+4*sum((x.control[,1]==2)*(x.control[,2]==2))
		p <- c()
		if ((n01d+n01c)!=0) 
			p <- c(p, n01d/(n01d+n01c))
			#(01,y=1)/(01)
		if ((n10d+n10c)!=0)
			p <- c(p, n10d/(n10d+n10c))
		if ((n11d+n11c)!=0)
			p <- c(p,n11d/(n11d+n11c))
		pexp <- n.case/(n.case+n.control)
		return(sum((p-pexp)^2))
}


f.I2.cont <- function(x,y,nperm){
  cat('Compute I_2 for continuous trait...\n')
  xxx.temp <- as.matrix(x)
  
  if (ncol(x)==1)
    stop("The number of columns of the genotype must 
    	be at least 2 to compute the interaciton score!")
  else {
	myindex <- which(colSums(xxx.temp)!=0)
  	xxx <- as.matrix(xxx.temp[,myindex])
  	if (ncol(xxx)==1)
  	  stop("The number of columns of the genotype must 
    		be at least 2 to compute the interaciton score!")
	mycombn <- combn(ncol(xxx),2)
	Iij <- apply(mycombn,2,function(a) f.Iij.cont(xxx[,a],y))
	I2 <- sum(Iij)
	
	I2.perm <- rep(NA,nperm)
	for (j in 1:nperm){
		yperm <- sample(y)
		Iij <- apply(mycombn,2,function(a) f.Iij.cont(xxx[,a],yperm))
		I2.perm[j] <- sum(Iij)
	}
	p <- sum(I2.perm>=rep(I2,nperm))/nperm
	return(list(I2=I2,pvalue=p))	
  }
}


f.Iij.cont <- function(x,yy){
		y00 <- mean(yy[x[,1]==0 & x[,2]==0])
		y01 <- mean(yy[x[,1]==0 & x[,2]>0])
		y10 <- mean(yy[x[,1]>0 & x[,2]==0])
		y11 <- mean(yy[x[,1]>0 & x[,2]>0])
		yave <- c(y01,y10,y11)
		nij <- sum((rowSums(x>0)>0))
		ybar <- mean(yy)
		return(nij^2*(sum((yave-ybar)^2,na.rm=T)))
}

#######################################################
f.pstar.dich <- function(x,y,nperm){
	cat('Compute p* for dichotomous trait...\n')
	  xxx.temp <- as.matrix(x)
  
  if (ncol(x)==1)
    stop("The number of columns of the genotype must 
    	be at least 2 to compute the interaciton score!")
  else {
	myindex <- which(colSums(xxx.temp)!=0)
  	xxx <- as.matrix(xxx.temp[,myindex])
  	if (ncol(xxx)==1)
  	  stop("The number of columns of the genotype must 
    		be at least 2 to compute the interaciton score!")
	mycombn <- combn(ncol(xxx),2)
	nij <- apply(mycombn,2,function(a) sum(xxx[,a]))
	wij <- nij/sum(nij)
	Iij <- apply(mycombn,2,function(a) f.Iij.dich(xxx[,a],y))
	I2 <- sum(wij*nij*Iij)
	
	I2.perm <- rep(NA,nperm)
	for (j in 1:nperm){
		yperm <- sample(y)
		Iij <- apply(mycombn,2,function(a) f.Iij.dich(xxx[,a],yperm))
		I2.perm[j] <- sum(wij*nij*Iij)
	}
	p2 <- sum(I2.perm>=rep(I2,nperm))/nperm


	I1 <- f.I1.dich.single(xxx,y)
	I1.perm <- rep(NA,nperm)
	for (j in 1:nperm){
		yperm <- sample(y)
		I1.perm[j] <- f.I1.dich.single(xxx,y)
	}
	p1 <- sum(I1.perm>=rep(I1,nperm))/nperm
	
	
	pstar <- min(p1,p2)	
	p1.perm <- rep(NA,nperm)
	p2.perm <- rep(NA,nperm)
	for (j in 1:nperm){
  		p1.perm[j] <- sum(I1.perm>=rep(I1.perm[j],nperm))/nperm
 	 	p2.perm[j] <- sum(I2.perm>=rep(I2.perm[j],nperm))/nperm
	}
	pstar.perm <- ifelse(p1.perm<p2.perm,p1.perm,p2.perm)
	p.pstar <-  sum(pstar.perm<=pstar)/nperm

	Is <- c(I1,I2,pstar)
	pvalues <- c(p1,p2,p.pstar)
	names(Is) <- c("I1","I2","pstar")
	names(pvalues) <- c("p_I1","p_I2","p_pstar")
	return(list(I=Is,pvalue=pvalues))
	}
}

f.pstar.cont <- function(x,y,nperm){
  cat('Compute p* for continuous trait...\n')
  xxx.temp <- as.matrix(x)
  
  if (ncol(x)==1)
    stop("The number of columns of the genotype must 
    	be at least 2 to compute the interaciton score!")
  else {
	myindex <- which(colSums(xxx.temp)!=0)
  	xxx <- as.matrix(xxx.temp[,myindex])
  	if (ncol(xxx)==1)
  	  stop("The number of columns of the genotype must 
    		be at least 2 to compute the interaciton score!")
	mycombn <- combn(ncol(xxx),2)
	Iij <- apply(mycombn,2,function(a) f.Iij.cont(xxx[,a],y))
	I2 <- sum(Iij)
	
	I2.perm <- rep(NA,nperm)
	for (j in 1:nperm){
		yperm <- sample(y)
		Iij <- apply(mycombn,2,function(a) f.Iij.cont(xxx[,a],yperm))
		I2.perm[j] <- sum(Iij)
	}
	p2 <- sum(I2.perm>=rep(I2,nperm))/nperm
	
	I1 <- f.I1.cont.single(x,y)
	I1.perm <- rep(NA,nperm)
	for (j in 1:nperm){
		yperm <- sample(y)
		I1.perm[j] <- f.I1.cont.single(x,yperm)
	}
	p1 <- sum(I1.perm>=rep(I1,nperm))/nperm
	
	pstar <- min(p1,p2)	
	p1.perm <- rep(NA,nperm)
	p2.perm <- rep(NA,nperm)
	for (j in 1:nperm){
  		p1.perm[j] <- sum(I1.perm>=rep(I1.perm[j],nperm))/nperm
 	 	p2.perm[j] <- sum(I2.perm>=rep(I2.perm[j],nperm))/nperm
	}
	pstar.perm <- ifelse(p1.perm<p2.perm,p1.perm,p2.perm)
	p.pstar <-  sum(pstar.perm<=pstar)/nperm

	Is <- c(I1,I2,pstar)
	pvalues <- c(p1,p2,p.pstar)
	names(Is) <- c("I1","I2","pstar")
	names(pvalues) <- c("p_I1","p_I2","p_pstar")
	return(list(I=Is,pvalue=pvalues))
	}
}



############################# Function I2star ##################
SPA.I.GE <- function(x,y,E,nperm=100){
# Function to compute I2star that deals with G*E interaction effect
# The phenotype is dichotomous
# returns the value of I2star and the p-value from both global and local permutation
	call <- match.call()

	if (length(setdiff(unique(y),c(0,1)))>0)
		stop ("The phenotype can only be dichotomous when testing with environmental factors!")
	
	cat('compute I2star with one environmental factor...\n')
	I2star <- f.I2star.single(x,y,E)	
	I2star.global.perm <- I2star.local.perm <- rep(NA,nperm)
	
	for (j in 1:nperm){
		I2star.global.perm[j] <- f.I2star.single(x,sample(y),E)
		I2star.local.perm[j] <- f.I2star.single(x,f.sampleE(y,E),E)
	}
	p.global <- sum(I2star.global.perm>=rep(I2star,nperm))/nperm
	p.local <- sum(I2star.local.perm>=rep(I2star,nperm))/nperm
	
	ans <- list(I2star=I2star,pvalue.global=p.global,pvalue.local=p.local,
				call=call,type="env")
	class(ans) <- "SPA"
	invisible(ans)
}


f.sampleE <- function(y,E){
# perform local permutation
	unE <- unique(E)
	lunE <- length(unE)
	ysample <- rep(NA,length(y))
	for (i in 1:lunE){
		ind <- (E==unE[i])
		ysample[ind] <- sample(y[ind])
	}
	return(ysample)	
}

f.I2star.single <- function(x,y,E){
	xxx.temp <- as.matrix(x)
	myindex <- which(colSums(xxx.temp)!=0)
	xxx <- as.matrix(xxx.temp[,myindex])
	ue <- unique(E)
	pexp <- mean(y)
	I <- 0
	for (i in ue){
		xe <- as.matrix(xxx[E==i,])
		ye <- y[E==i]
		xe.case <- as.matrix(xe[ye==1,])
		xe.control <- as.matrix(xe[ye==0,])
		ped <- colSums(xe.case)/(colSums(xe.case)+colSums(xe.control))
		nei <- colSums(xe)
		I <- I+sum(nei^2*(ped-pexp)^2,na.rm=TRUE)
	}
	return(I)
}


