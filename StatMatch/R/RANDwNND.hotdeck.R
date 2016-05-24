`RANDwNND.hotdeck` <-
function (data.rec, data.don, match.vars=NULL, don.class=NULL, dist.fun="Manhattan", cut.don="rot", k=NULL, weight.don=NULL, keep.t=FALSE,  ...)
{
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
	
    p <- length(match.vars)
    if(!is.null(match.vars)){
        if(dist.fun=="Euclidean" || dist.fun=="euclidean" || dist.fun=="Manhattan" || dist.fun=="manhattan" || dist.fun=="Mahalanobis" || dist.fun=="mahalanobis" || dist.fun=="minimax" || dist.fun=="MiniMax" || dist.fun=="Minimax"){
            cat("Warning: The ", dist.fun, " distance is being used", fill=TRUE)
            cat("All the categorical matching variables in rec and don data.frames, \n if present, are recoded into dummies", fill=TRUE)
        }
        if(dist.fun=="exact" || dist.fun=="exact matching"){
            cat("Warning: the exact matching distance is being used", fill=TRUE)
            cat("all the matching variables in rec and don are converted to \n character variables and are treated as categorical nominal", fill=TRUE)
        }
        if(dist.fun=="difference" || dist.fun=="diff"){
            if(p>1) cat("When dist.fun='difference' just a single matching 
                        \n variable should be used", fill=TRUE)
            if(!(cut.don %in% c("lt","<","le","<=","ge",">=","gt",">"))) {
                cat("When dist.fun='difference' the argument cut.don should be \n
                        a comparison operator", fill=TRUE)
            }
        }
    }
################
RANDwNND.hd <- function (rec, don, dfun="Manhattan", cut.don="rot", k=NULL, w.don=NULL, ...)
{ 
    x.rec <- rec
    x.don <- don
    p <- ncol(rec)
  	nr <- nrow(rec)
	nd <- nrow(don)
	if(nr>nd) cat("Warning: the number of donors is less than the number of recipients", fill=TRUE)

    r.lab <- rownames(x.rec)
	if(is.null(r.lab)) r.lab <- paste("rec", 1:nr, sep="=")
	d.lab <- rownames(x.don)
	if(is.null(d.lab)) d.lab <-  paste("don", 1:nr, sep="=")
    if(is.null(w.don)) ww <- rep(1,nd)
    else ww <- w.don
# compute matrix of distances between obs. in x.don and obs. in x.rec
# function dist() in package "proxy" is used! 
    if(cut.don=="rot"){
        k <- ceiling(sqrt(nd))
        if(k==0) stop("k=sqrt(no. of dons) is equal to 0")
    }
    else if(cut.don=="span"){ 
        if(k==0 || k>1) stop("When cut.don='span' then  0 < k <= 1")
        k <- ceiling(nd*k)
        if(k>nd) k <- nd
        if(k==0) stop("k=round(N. of dons * k) is equal to 0")
    }
    else if(cut.don=="exact"){
        if(k==0 || k>nd) stop("When cut.don=exact, k should be such that  1 < k <= no. of dons") 
    }
    
    if(dfun=="Euclidean" || dfun=="Manhattan"){
   #     require(proxy)
        x.rec <- fact2dummy(x.rec, all=FALSE)
        x.don <- fact2dummy(x.don, all=FALSE)
        mdist <- dist(x=x.rec, y=x.don, method=dfun, ...)
        dimnames(mdist) <- list(r.lab, d.lab)
    }
	else if(dfun=="Mahalanobis" || dfun=="mahalanobis"){
        if(is.data.frame(x.rec)) x.rec <- fact2dummy(x.rec, all=FALSE)
        if(is.data.frame(x.don)) x.don <- fact2dummy(x.don, all=FALSE)
        mdist <- mahalanobis.dist(data.x=x.rec, data.y=x.don, ...)
        dimnames(mdist) <- list(r.lab, d.lab)
	}
	else if(dfun=="minimax" || dfun=="MiniMax" || dfun=="Minimax"){
        x.rec <- fact2dummy(x.rec, all=FALSE)
        x.don <- fact2dummy(x.don, all=FALSE)
        mdist <- maximum.dist(data.x=x.rec, data.y=x.don, ...)
        dimnames(mdist) <- list(r.lab, d.lab)
	}
    else if(dfun=="exact" || dfun=="exact matching"){
        dxr <- dim(x.rec)
        x.rec <- as.character(as.matrix(x.rec))
        dim(x.rec) <- dxr
        dxd <- dim(x.don)
        x.don <- as.character(as.matrix(x.don))
        dim(x.don) <- dxd
        xx <- data.frame(rbind(x.rec, x.don))
		x.rec <- xx[1:nr,]
		x.don <- xx[-(1:nr),]
		mdist <- gower.dist(data.x=x.rec, data.y=x.don)
        dimnames(mdist) <- list(r.lab, d.lab)
    }
    else if(dfun=="Gower" || dfun=="gower"){
        # if(p==1 && is.factor(x.rec)) x.rec <- list(x.rec)
        # if(p==1 && is.factor(x.don)) x.don <- list(x.don)
        mdist <- gower.dist(data.x=x.rec, data.y=x.don, ...)
        mdist[is.nan(mdist)] <- 1 # NaN can occur when p=1 and x.rec and x.don is of type logical
        mdist[is.na(mdist)] <- 1 # NA can occur when p=1 and x.rec and x.don is of type logical
        dimnames(mdist) <- list(r.lab, d.lab)
    }
    else if(dfun=="RANN" || dfun=="ANN"){
 #       require(RANN)
        if(cut.don=="min") k0 <- 10
        else if (cut.don=="k.dist") stop("When dist.fun='RANN' it is not possible to to set \n cut.don = 'k.dist' ")
        else k0 <- k
        dd <- nn2(data=x.don, query=x.rec, k=k0, ...)
        mdist <- dd$nn.dists
    }
    else if(dfun=="difference" || dfun=="diff"){
        
        xA <- data.matrix(x.rec)
        xB <- data.matrix(x.don)
        mdist <- outer(as.numeric(xA), as.numeric(xB), FUN="-")
    }

    else {
  #      require(proxy)
        mdist <- dist(x=x.rec, y=x.don, method=dfun, ...)
        dimnames(mdist) <- list(r.lab, d.lab)
    }

    min.d <- numeric(nr)
    max.d <- numeric(nr)
    sd.d <- numeric(nr)
    cut.d <- numeric(nr)
    dist.rd <- numeric(nr)
    nad <- rep(NA, nr)
    don.lab <- numeric(nr)
    
	for(i in 1:nr){
        vd <- mdist[i,]
        if(sum(is.na(vd))==nd) stop("The missing values on the mtc. vars determine missing distances")
#        cat("rec obs: ", i, fill=TRUE) #add check
#        cat("number of donors: ", nd, fill=TRUE) #add check
#        cat("number of non-miss distances: ", sum(!is.na(vd)), fill=TRUE) #add check
#        cat("distances:", summary(vd), fill=TRUE) #add check
        min.dist <- min(vd, na.rm=TRUE) # smallest distance recipient-donor
        min.d[i] <- min.dist
        if(dfun=="RANN" || dfun=="ANN") {
            max.d[i] <- NA
            sd.d[i] <- NA
        }
        else {
            max.d[i] <- max(vd, na.rm=TRUE)
		    sd.d[i] <- sd(vd, na.rm=TRUE)
        }
        if(cut.don=="min"){
            tst <- (vd==min.dist) & !is.na(vd)
            short.vd <- vd[tst]
            if(dfun=="RANN" || dfun=="ANN"){
                idx <- dd$nn.idx[i,]
                appo <- d.lab[idx]
                short.ww <- ww[idx]
            }
            else {
                appo <- d.lab[tst]
                short.ww <- ww[tst]
            }
            dist.rd[i] <- min.dist
            cut.d[i] <- min.dist
		}	
        
        else if(cut.don=="lt" || cut.don=="<"){
            tst <- vd<0
            appo <- d.lab[tst]
            short.vd <- vd[tst]
            short.ww <- ww[tst]
            cut.d[i] <- NA
        }
        else if(cut.don=="le" || cut.don=="<="){
            tst <- vd<=0
            appo <- d.lab[tst]
            short.vd <- vd[tst]
            short.ww <- ww[tst]
            cut.d[i] <- NA
        }
        else if(cut.don=="ge" || cut.don==">="){
            tst <- vd>=0
            appo <- d.lab[tst]
            short.vd <- vd[tst]
            short.ww <- ww[tst]
            cut.d[i] <- NA
        }
        else if(cut.don=="gt" || cut.don==">"){
            tst <- vd>0
            appo <- d.lab[tst]
            short.vd <- vd[tst]
            short.ww <- ww[tst]
            cut.d[i] <- NA
        }

        else if(cut.don=="k.dist"){
			if(k<min.dist) {
			    cat("Warning: the value of k,", k, fill=TRUE)
			    cat("is smaller than the minimum distance:", min.d, fill=TRUE)
			}			
            tst <- (vd<=k) & !is.na(vd)
		    appo <- d.lab[tst]
			short.vd <- vd[tst]
            short.ww <- ww[tst]
            cut.d[i] <- k
        }
        else {
            if(dfun=="RANN" || dfun=="ANN"){
                pos <- dd$nn.idx[i,]
                appo <- d.lab[pos]
                short.vd <- dd$nn.dists[i,]
                short.ww <- ww[pos]
                cut.d[i] <- dd$nn.dists[i,ncol(dd$nn.dists)]
            }
            else{
                pos <- order(vd, na.last=NA)
			    if(length(vd)<k) kk <- length(appo)
			    else kk <- k
                pos <- pos[1:kk]
#            cat("k:", k, fill=TRUE) #add check
#            cat("closest donors", pos, fill=TRUE) #add check
                appo <- d.lab[pos]
			    short.vd <- vd[pos]
			    short.ww <- ww[pos]
                cut.d[i] <- short.vd[kk]
            }
        }

        nad[i] <- length(appo) # number of availabe donors
        if(length(appo)==0) {
            nad[i] <- 0
			don.lab[i] <- NA
		    dist.rd[i] <- NA
            cat("Warning: there are no available donors for the", r.lab[i], "recipient!", fill=TRUE)
        }       
        else if(length(appo)==1){
			 nad[i] <- 1
			 don.lab[i] <- appo 
			 dist.rd[i] <- short.vd
		}	
        else{
            nn.dd <- length(appo)
            nad[i] <- nn.dd
            choi <- sample(1:nn.dd, 1, replace=TRUE, prob=short.ww)
            don.lab[i] <- appo[choi]
            dist.rd[i] <- short.vd[choi]
		}
    }       
    rec.lab <- r.lab

# output
    mtc.ids <- cbind(rec.id=rec.lab, don.id=don.lab)
	sum.dist <- cbind(min=min.d, max=max.d, sd=sd.d, cut=cut.d, dist.rd=dist.rd)
	list(mtc.ids=mtc.ids, sum.dist=sum.dist, noad=nad, call=match.call())
}
############################################
# RANDwNND.hd ends here
############################################	

##############################
# function pps.draw

pps.draw <- function(n, w){
    tot <- sum(w)
    udraw <- runif(n)
    pos <- cut(udraw*tot, breaks=c(0,cumsum(w)))
    as.integer(pos)
}
# pps.draw(n,w) ends here
############################

	if(is.null(don.class)){ 
		if(is.null(match.vars)){
            if(is.null(weight.don)) don.lab <- sample(d.lab, nr, replace=TRUE)
			else don.lab <- d.lab[pps.draw(n=nr, w=data.don[,weight.don])]
            mtc.ids <- cbind(rec.id=r.lab, don.id=don.lab)
			noad <- rep(nd, nr)
			sss.dist <- NULL
		}
		else{
            REC <- data.rec[, match.vars, drop=FALSE]
			DON <- data.don[, match.vars, drop=FALSE]
            if(is.null(weight.don)) out <- RANDwNND.hd(rec=REC, don=DON, dfun=dist.fun, cut.don=cut.don, k=k, w.don=NULL, ...)
            else out <- RANDwNND.hd(rec=REC, don=DON, dfun=dist.fun, cut.don=cut.don, k=k, w.don=data.don[,weight.don], ...)
			mtc.ids <- out$mtc.ids
			noad <- out$noad
			sss.dist <- out$sum.dist
		}
		mmm <- substring(c(mtc.ids), 5)
		mtc.ids <- matrix(mmm, ncol=2)
		if(is.null(rownames(data.rec)) && is.null(rownames(data.don)))  mtc.ids <- matrix(as.numeric(mmm), ncol=2)
	}
	else{
		if(length(don.class)==1){
			l.r.lab <- split(r.lab, f=data.rec[ ,don.class], drop=TRUE)
			l.rec <- split(data.rec, f=data.rec[ ,don.class], drop=TRUE)
			l.d.lab <- split(d.lab, f=data.don[ ,don.class], drop=TRUE)
			l.don <- split(data.don, f=data.don[ ,don.class], drop=TRUE)
		}
		else{
			l.r.lab <- split(r.lab, f=as.list(data.rec[ ,don.class]), drop=TRUE)
			l.rec <- split(data.rec, f=as.list(data.rec[ ,don.class]), drop=TRUE)
			l.d.lab <- split(d.lab, f=as.list(data.don[ ,don.class]), drop=TRUE)
			l.don <- split(data.don, f=as.list(data.don[ ,don.class]), drop=TRUE)
		}
		tst <- which( !(names(l.rec) %in% names(l.don)) )
		if(length(tst)) {
		    cat("The following donation classes:", fill=TRUE) 
		    cat(names(l.rec)[tst], fill=TRUE)
		    cat("in data.don are empty, i.e. there are no donors", fill=TRUE)
		    stop()
		}
		
#        if(length(l.rec)!=length(l.don)){
#			cat("The no. of donation classes in recipient data is not equal \n to the no. of donation classes in donor data", fill=TRUE)
#			stop("Possible reason: the variables used to classify units are not \n defined as factors or are factors with different levels")
#		}
#		if(!identical(names(l.rec), names(l.don)))
#			cat("Warning: the donation classes seem built \n using different factors with differnt levels")
			
		nn.r <- unlist(lapply(l.r.lab, length))
		nn.d <- unlist(lapply(l.d.lab, length))
				
#		if(sum((nn.r>0) & (nn.d==0))>0) {
#			stop("In some donation classes there are NO available donors. \n Please modify the definition of the donation classes")
#		}	
#		l.r.lab <- l.r.lab[nn.r>0]
#		l.d.lab <- l.d.lab[nn.r>0]
#		l.rec <- l.rec[nn.r>0]
#		l.don <- l.don[nn.r>0]
				
		H <- length(l.rec)
		mtc.ids <- as.list(numeric(H))
		sum.dist <- as.list(numeric(H))
		noad <- as.list(numeric(H))
		
		for(h in 1:H){
            lab.h <- names(l.rec)[h]
            if(keep.t) cat("Selecting donors for donation class: ", lab.h, fill=TRUE)
			if(is.null(match.vars)){
                if(is.null(weight.don)){
				    don.lab <- sample(l.d.lab[[lab.h]], nn.r[[lab.h]], replace=TRUE)
                }
                else{
                    pos <- pps.draw(n=nn.r[[lab.h]], w=l.don[[lab.h]][,weight.don])
                    don.lab <- l.d.lab[[lab.h]][pos]
                }
			    mtc.ids[[h]] <- cbind(rec.id=l.r.lab[[lab.h]], don.id=don.lab)
    		    sum.dist[[h]] <- NA
			    noad[[h]] <- rep(nn.d[[lab.h]], nn.r[[lab.h]])
			}
			else{
    			REC <- l.rec[[lab.h]][,match.vars, drop=FALSE]
				DON <- l.don[[lab.h]][,match.vars, drop=FALSE]
				if(is.null(weight.don)) out <- RANDwNND.hd(rec=REC, don=DON, dfun=dist.fun, w.don=NULL, cut.don=cut.don, k=k, ...)
				else out <- RANDwNND.hd(rec=REC, don=DON, dfun=dist.fun, w.don=l.don[[lab.h]][,weight.don], cut.don=cut.don, k=k, ...)
				mtc.ids[[h]] <- out$mtc.ids
				sum.dist[[h]] <- out$sum.dist
				noad[[h]] <- out$noad
			}	
		}
		
		mmm <- unlist(lapply(mtc.ids, t))
		mmm <- substring(mmm, 5)
		mtc.ids <- matrix(mmm, ncol=2, byrow=TRUE)
		if(is.null(rownames(data.rec)) && is.null(rownames(data.don)))  mtc.ids <- matrix(as.numeric(mmm), ncol=2, byrow=TRUE)
		sss.dist <- NA
		if(!is.null(match.vars)){
			
			sss <- unlist(lapply(sum.dist, t))
			sss.dist <- matrix(sss, ncol=5, byrow=TRUE)
			colnames(sss.dist)<- colnames(out$sum.dist)
		}
		noad <- unlist(noad)
	}
	dimnames(mtc.ids) <- list(NULL, c("rec.id", "don.id"))
	list(mtc.ids=mtc.ids, sum.dist=sss.dist, noad=noad, call=match.call())
}

