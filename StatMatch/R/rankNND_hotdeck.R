`rankNND.hotdeck` <-
function (data.rec, data.don, var.rec, var.don=var.rec, don.class=NULL, weight.rec=NULL, weight.don=NULL, constrained=FALSE, constr.alg="Hungarian", keep.t=FALSE)
{
    if(is.na(match(var.rec, colnames(data.rec)))) stop("The variable var.rec is not available in data.rec")
    if(is.na(match(var.don, colnames(data.don)))) stop("The variable var.don is not available in data.don")
    if(!is.null(weight.rec) &&  is.na(match(weight.rec, colnames(data.rec)))) {
        stop("The variable weight.rec is not available in data.rec")
    }
    if(!is.null(weight.don) &&  is.na(match(weight.don, colnames(data.don)))) {
        stop("The variable weight.don is not available in data.don")
    }
    
 #   if(constrained && (constr.alg=="Hungarian" || constr.alg=="hungarian")) require(clue)
 #   if(constrained && (constr.alg=="lpSolve" || constr.alg=="lpsolve")) require(lpSolve)
    
    if(!is.null(dim(data.rec))){
		nr <- nrow(data.rec)
		r.lab <- row.names(data.rec)
	}
	else{
		nr <- length(data.rec)
		r.lab <- names(data.rec)
	}
	if(!is.null(dim(data.don))){
		nd <- nrow(data.don)
		d.lab <- row.names(data.don)
	}
	else{
		nd <- length(data.don)
		d.lab <- names(data.don)
	}
	
	if(is.null(r.lab)) r.lab <- paste("rec", 1:nr, sep="=")
	else r.lab <- paste("rec", r.lab, sep="=")
	row.names(data.rec) <- r.lab
	
	if(is.null(d.lab)) d.lab <- paste("don", 1:nd, sep="=")
	else d.lab <- paste("don", d.lab, sep="=")
	row.names(data.don) <- d.lab

    if(is.null(weight.rec)) {
        weight.rec <- "weight.rec"
        data.rec <- data.frame(data.rec, weight.rec=1)
    }
    if(is.null(weight.don)) {
        weight.don <- "weight.don"
        data.don <- data.frame(data.don, weight.don=1)
    }

############################################################
rankNND.hd <- function (x.rec, x.don, w.rec=NULL, w.don=NULL, constr=FALSE, c.alg=NULL)
{ 
    nr <- length(x.rec)
	nd <- length(x.don)
	if(nr>nd) cat("Warning: the number of donors is less \n than the number of recipients", fill=TRUE)

    if(is.null(w.rec)) w.rec <- rep(1,nr)
    if(is.null(w.don)) w.don <- rep(1,nd)

# computes empirical cumulative ...
    cdf.rec <- numeric(nr)
    for(i in 1:nr) cdf.rec[i] <- sum(w.rec[x.rec<=x.rec[i]])/sum(w.rec)

    cdf.don <- numeric(nd)
    for(i in 1:nd) cdf.don[i] <- sum(w.don[x.don<=x.don[i]])/sum(w.don)

	mdist <- abs(outer(X=cdf.rec, Y=cdf.don, FUN="-"))
 	dimnames(mdist) <- list(1:nr, 1:nd)

# unconstrained matching    
    if(!constr){
        dist.rd <- numeric(nr)
	    nad <- rep(NA, nr)
	    don.pos <- numeric(nr)
	    for(i in 1:nr){
			vd <- mdist[i,]
            min.d <- min(vd) # smallest distance recipient-donor
			dist.rd[i] <- min.d
			appo <- c(1:nd)[vd==min.d]
			nad[i] <- length(appo) # number of availabe donors
            if(length(appo)==1) don.pos[i] <- appo
			else don.pos[i] <- sample(appo, 1)
	    }
    }

# CONSTRAINED nearest neighbour matching.
# the functions in library lpSolve are used
    
    if(constr && (c.alg=="lpSolve" || c.alg=="lpsolve")){
        if(nr==nd) appo <- lp.assign(cost.mat=mdist)
        else if(nr<nd){
            r.sig <- rep("==", nr)
            r.rhs <- rep(1, nr)
            c.sig <- rep("<=", nd)
            c.rhs <- rep(1, nd)
            appo <- lp.transport(cost.mat=mdist, row.signs=r.sig, row.rhs=r.rhs, col.signs=c.sig, col.rhs=c.rhs)
        }   
        else if(nr > nd){
            warning("There more recipients than donors!")
            cat("some donors will be used more than once", fill=TRUE)
            r.sig <- rep("==", nr)
            r.rhs <- rep(1, nr)
            c.sig <- rep(">=", nd)
            c.rhs <- rep(1, nd)
            appo <- lp.transport(cost.mat=mdist, row.signs=r.sig, row.rhs=r.rhs, col.signs=c.sig, col.rhs=c.rhs)
        }
        sol <- appo$solution
        ss <- c(t(sol))
        cc <- c(t(col(sol)))
        dist.rd <- mdist[cbind(1:nr, cc[as.logical(ss)] )]
        don.pos <- (1:nd)[c(cc[as.logical(ss)])]
        nad=NULL
    }
    
    # the function solve_LSAP() in package clue is used
    if(constr && (c.alg=="Hungarian" || c.alg=="hungarian")){
        if(nr > nd) stop("It is required that no. of donors is greater \n or equal than the no. of recipients")
        sol <- solve_LSAP(x=mdist, maximum=FALSE)
        don.pos <- as.integer(sol)
        dist.rd <- mdist[cbind(1:nr, don.pos)]
        nad=NULL
    }
    
    
# output
    list(don.pos=don.pos, dist.rd=dist.rd, noad=nad, call=match.call())
}

################ rankNND.hd ends here #############################

	if(is.null(don.class)){ 
		out <- rankNND.hd(x.rec=data.rec[,var.rec], x.don=data.don[,var.don],
                            w.rec=data.rec[,weight.rec], w.don=data.don[,weight.don],
		                    constr=constrained, c.alg=constr.alg)
        mtc.lab <- cbind(r.lab, d.lab[out$don.pos])
		mmm <- substring(mtc.lab, 5)
		if(is.null(rownames(data.rec)) && is.null(rownames(data.don)))  mtc.ids <- matrix(as.numeric(mmm), ncol=2, byrow=TRUE)
		else mtc.ids <- mmm
		dimnames(mtc.ids) <- list(NULL, c("rec.id", "don.id"))
		dist.rd <- out$dist.rd
		noad <- out$noad
	}
	else{
		if(length(don.class)==1){
			l.rec <- split(data.rec[ ,c(var.rec, weight.rec)], f=data.rec[ ,don.class], drop=TRUE)
			l.don <- split(data.don[ ,c(var.don, weight.don)], f=data.don[ ,don.class], drop=TRUE)
		}
		else{
			l.rec <- split(data.rec[ ,c(var.rec, weight.rec)], f=as.list(data.rec[ ,don.class]), drop=TRUE)
			l.don <- split(data.don[ ,c(var.don, weight.don)], f=as.list(data.don[ ,don.class]), drop=TRUE)
		}
#		if(length(l.rec)!=length(l.don)){
#			cat("The no. of donation classes in recipient data is not equal to \n the no. of donation classes in donor data", fill=TRUE)
#			stop("Possible reason: the variables used to classify units are not \n defined as factors or are factors with different levels")
#		}
#		if(!identical(names(l.rec), names(l.don)))
#			cat("Warning: the donation classes seem built using \n different factors with differnt levels")

		tst <- which( !(names(l.rec) %in% names(l.don)) )
		if(length(tst)) {
		    cat("The following donation classes:", fill=TRUE) 
		    cat(names(l.rec)[tst], fill=TRUE)
		    cat("in data.don are empty, i.e. there are no donors", fill=TRUE)
		    stop()
		}
		
        nn.r <- unlist(lapply(l.rec, nrow))
        nn.d <- unlist(lapply(l.don, nrow))

#        l.rec <- l.rec[nn.r>0]
#        l.don <- l.don[nn.r>0]
#        nn.r <- nn.r[nn.r>0]
#        nn.d <- nn.d[nn.r>0]
#        if(any(nn.d==0)) {
#			stop("For some donation classes there are NO donors available. \n Please modify the definition of the donation classes")
#		}	
		H <- length(l.rec)
		mtc.ids <- as.list(numeric(H))
		dist.rd <- as.list(numeric(H))
	    noad <- as.list(numeric(H))
		for(h in 1:H){
		    lab.h <- names(l.rec)[h]
		    if(keep.t) cat("Selecting donors for donation class: ", lab.h, fill=TRUE)
    		out <- rankNND.hd(x.rec=l.rec[[lab.h]][,var.rec], x.don=l.don[[lab.h]][,var.don],
                              w.rec=l.rec[[lab.h]][,weight.rec], w.don=l.don[[lab.h]][,weight.don],
    		                  constr=constrained, c.alg=constr.alg)
            mtc.ids[[h]] <- cbind(rownames(l.rec[[lab.h]]), rownames(l.don[[lab.h]])[out$don.pos])
			dist.rd[[h]] <- out$dist.rd
			noad[[h]] <- out$noad
		}

        mmm <- unlist(lapply(mtc.ids, t))
		mmm <- substring(mmm, 5)
		mtc.ids <- matrix(mmm, ncol=2, byrow=TRUE)
		if(is.null(rownames(data.rec)) && is.null(rownames(data.don)))  mtc.ids <- matrix(as.numeric(mmm), ncol=2, byrow=TRUE)
		dimnames(mtc.ids) <- list(NULL, c("rec.id", "don.id"))
		dist.rd <- unlist(dist.rd)
		if(!constrained) noad <- unlist(noad)
        else noad <- NULL
	}
	end.out <- list(mtc.ids=mtc.ids, dist.rd=dist.rd, noad=noad, call=match.call())
	end.out
}
