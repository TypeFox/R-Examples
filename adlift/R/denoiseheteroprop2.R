`denoiseheteroprop2` <-function (x, f, pred=1, neigh=1, int=TRUE, clo=FALSE, keep=2, rule = "median", gamvec=rep(1,times=length(x)), verbose=TRUE,index=1:length(x),returnall=FALSE) 
{

newcoeff <- NULL
ndetlist <- list()
tclist <- NULL
nonornewcoeff <- NULL
norcoeff <- NULL

n<-length(x)
out <- fwtnp2(x, f, pred, neigh, int, clo, keep,do.W=TRUE,varonly=FALSE)
w <- out$Wnew

sigmamat<-diag(gamvec^2)

varmat <- w %*% (sigmamat)
varmat <- (varmat) %*% t(w)
indsd <- sqrt(diag(varmat))
norcoeff <- out$coeff/indsd
nonorcoeff <- out$coeff
po<-out$pointsin
lr <- out$lengthsremove
lca<-out$lca
rem <- lca[,1]
al <- artlev(lr, rem)
levno <- length(al)

for (i in 1:levno) {
	ndetlist[[i]] <- norcoeff[al[[i]]]
}

sd <- mad(ndetlist[[1]])
sd1 <- mad(ndetlist[[1]], center = 0)
sd2 <- mad(norcoeff[rem])

for (i in 1:levno) {
	tclist <- ebayesthresh(ndetlist[[i]], prior = "cauchy", 
            a = NA, sdev = sd, threshrule = rule)
	newcoeff[al[[i]]] <- tclist * indsd[al[[i]]]
}

newcoeff[po]<-out$coeff[po]
    
neighbrs<-list()
for(kk in 1:length(rem)){
	neighbrs[[kk]]<-.C("nbrsfromlca",as.double(t(lca)),as.integer(ncol(lca)),as.integer(kk),
                    n=as.integer(rep(0,times=lca[kk,2])),PACKAGE="adlift")$n
}
                    
adds <- findadds(rem, neighbrs, po, index)
add <- max(adds)
if (verbose) {
	cat("doing ", add, " steps...\n", sep = "")
}
                            
fhat <- invtnp2(x, newcoeff, out$lengths, lr, po,lca, add)

if(returnall){                            
	return(list(fhat=fhat,w = w, indsd = indsd, al = al, sd = sd))
}
else{
	return(fhat$coeff)
}

}

