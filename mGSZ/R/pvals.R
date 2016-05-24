emp <-
function(pos.data, neg.data){
    
    neg.dat.sz <- dim(neg.data)
    test_data <- rep(0, neg.dat.sz[2])
    for (k in 1:neg.dat.sz[2]){
        test_data[k] <- length(which(neg.data[,k] >= pos.data[k]))
    }
    test_data[test_data == 0] <- 0.5
    test_data <- (test_data)/(neg.dat.sz[1])
    #test_data <- -log10(test_data)
    test_data
    
    
}

############

emp.wrs <-
function(pos.data, neg.data){
    neg.dat.sz <- dim(neg.data)
    test_data <- rep(0, neg.dat.sz[2])
    for (k in 1:neg.dat.sz[2]){
        test_data[k] <- length(which(neg.data[,k] >= abs(pos.data[k])))
    }
    test_data[test_data == 0] <- 0.5
    test_data <- (test_data)/(neg.dat.sz[1])
    #test_data <- -log10(test_data)
    test_data
}

###############

KS.p.values <-
function(pos.data, perm.data){
	#library(MASS)
	tmp <- mGSZ.adj.mean.std(pos.data, perm.data)
	pos.data <- tmp$pos.out
	perm.data <- tmp$perm.out
	rm(tmp)
	col.ind <- pick.data.cols(perm.data)
    perm.data <- perm.data[,col.ind]
	pos.data <- pos.data[col.ind]
	emp.pval.class <- emp(pos.data, perm.data)
		
    output <- list(EMP.class=emp.pval.class,class.ind=col.ind)
    
    
    return(output)
}

###############

logEVcdf <-
function(x, mu, sigma) {
    
    x <- (mu - x)/sigma  # Note! Flip the sign here
    out <- rep(0, length(x))
    if(!(is.vector(x)) ){
        d.sz <- dim(x)
        out <- matrix(out, d.sz[1], d.sz[2])
    }
    po1 <- which(x < 5) # switch for appr. vs. norm. eq.
    out[po1] <- -log(1 - exp(-exp(x[po1])) )
    x <- x[-po1]
    out[-po1] <- -x + exp(x)/2 - exp(2*x)/24 + exp(4*x)/2880 #Polyn. Expansion
    out <- out/log(10)
    out <- 10^(-out)
    out
}

###############

logNORMcdf <-
function(x,mu,sigma){
	z <- (x-mu)/(sigma)
	
	#out <- rep(0,length(x))
	if(!(is.vector(x)) ){
        d.sz <- dim(x)
        out <- matrix(out, d.sz[1], d.sz[2])
    }
  	out <- pnorm(-z,lower.tail=T)
  	#out <- -log10(out)
    
  	out
}

###############

mAllez.p.values <-
function(pos.data, perm.data){
	#library(MASS)
	tmp <- mGSZ.adj.mean.std(pos.data, perm.data)
	pos.data <- tmp$pos.out
	perm.data <- tmp$perm.out
	rm(tmp)
	col.ind <- pick.data.cols(perm.data)
	perm.data <- perm.data[,col.ind]
	pos.data <- pos.data[col.ind]
	norm.param.all <- fitdistr(as.numeric(perm.data),'normal')$estimate
	norm.pval.class <- rep(0,length(pos.data))
	
	for (k in 1:length(pos.data)){
        norm.param.class <- fitdistr(as.vector(perm.data[,k]),'normal')$estimate
    	norm.pval.class[k] <- logNORMcdf(pos.data[k],norm.param.class[1],norm.param.class[2])

	}
	

    output <- list(NORM.class=norm.pval.class,class.ind=col.ind)
    output
	
}

################

mGSA.p.values <-
function(pos.data, perm.data) {
    #library(ismev)
    #library(MASS)
    tmp <- mGSZ.adj.mean.std(pos.data, perm.data)
    scaled.pos.data <- tmp$pos.out
    scaled.perm.data <- tmp$perm.out
    rm(tmp)
    get.perm.p.vals = FALSE
    perm.dat.sz <- dim(scaled.perm.data)
    col.ind <- pick.data.cols(scaled.perm.data) # discard cols with null var

    # First fit the parameters to the whole data
    scaled.perm.data <- scaled.perm.data[,col.ind]
    scaled.pos.data <- scaled.pos.data[col.ind]
    ev.p.val.class  <- rep(0, length(scaled.pos.data))
        
    for (k in 1:length(scaled.pos.data)){
        
        # Following fits parameters to perms of each class
        ev.param.class  <- gum.fit(as.vector(scaled.perm.data[,k]),show=FALSE)$mle
        ev.p.val.class[k]  <- logEVcdf(scaled.pos.data[k], ev.param.class[1], ev.param.class[2]  )
        
        }
        output <- list(EV.class = ev.p.val.class,class.ind=col.ind)
        
        output
}

###############

mGSZ.p.values <-
function(pos.data, perm.data) {
    #library(ismev)
    #library(MASS)
    tmp <- mGSZ.adj.mean.std(pos.data, perm.data)
    pos.data <- tmp$pos.out
    perm.data <- tmp$perm.out
    rm(tmp)
    get.perm.p.vals = FALSE
    perm.dat.sz <- dim(perm.data)
    col.ind <- pick.data.cols(perm.data)
    perm.data <- perm.data[,col.ind]
    pos.data <- pos.data[col.ind]
    ev.p.val.class  <- rep(0,length(pos.data))
        
    
    for (k in 1:length(pos.data)){

        # Following fits parameters to perms of each class
        ev.param.class  <- gum.fit(as.vector(perm.data[,k]),show=FALSE)$mle
        ev.p.val.class[k]  <- logEVcdf(pos.data[k], ev.param.class[1], ev.param.class[2]  )
        }

        output <- list(EV.class = ev.p.val.class,class.ind=col.ind)
        
    
    output
}

###############

SS.p.values <-
function(pos.data, perm.data){
	#library(MASS)
    tmp <- mGSZ.adj.mean.std(pos.data, perm.data)
    
	scaled.pos.data <- tmp$pos.out
	scaled.perm.data <- tmp$perm.out
	rm(tmp)
	col.ind <- pick.data.cols(scaled.perm.data)
	div.length <- length(scaled.perm.data)
	scaled.perm.data <- scaled.perm.data[,col.ind]
	scaled.pos.data <- scaled.pos.data[col.ind]
    emp.pval.class <- emp(scaled.pos.data, scaled.perm.data)
    output <- list(EMP.class = emp.pval.class,class.ind=col.ind )
    output
	
}

##############

WRS.p.values <-
function(pos.data, perm.data){
	#library(MASS)
	tmp <- mGSZ.adj.mean.std(pos.data, perm.data)
	pos.data <- tmp$pos.out
	perm.data <- tmp$perm.out
	rm(tmp)
	col.ind <- pick.data.cols(perm.data)
    perm.data <- perm.data[,col.ind]
	pos.data <- pos.data[col.ind]
    emp.pval.class <- emp.wrs(pos.data, perm.data)
    norm.pval.class <- rep(0,length(pos.data))

    output <- list(EMP.class=emp.pval.class,class.ind=col.ind)
    
    output
	
}

################

mGSZ.adj.mean.std <-
function(pos.data, perm.data) {
   perm.sz <- dim(perm.data)
   mean.prof <- apply(perm.data, 2, mean)
   std.prof    <- apply(perm.data, 2, sd)
   null.ind <- which(std.prof == 0)  
   std.prof[null.ind] <- 0.1
   perm.out <- matrix(0, perm.sz[1], perm.sz[2])
   pos.out <- rep(0, perm.sz[2])
   for( k in 1:perm.sz[2]){
     perm.out[,k] <- ( perm.data[,k] - mean.prof[k]) /std.prof[k]
     pos.out[k]    <- ( pos.data[k] - mean.prof[k] )/std.prof[k]   
    }
   out  <- list(pos.out = pos.out, perm.out = perm.out)
   out
}

##############

pick.data.cols <-
function(data) {

    # Discard cols where var == 0
    # These disturb distribution fitting
    
    tmp <- apply(data, 2, var)
    out <- which(tmp > 0)
    out
}

#############
