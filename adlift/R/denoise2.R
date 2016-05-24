`denoise2` <-function (x, f, pred=1, neigh=1, int=TRUE, clo=FALSE, keep=2, rule = "median", index = 1:length(x), verbose = FALSE, returnall = FALSE) 
{

d <- length(x) - keep
newcoeff <- indsd <- sd <- NULL
ndetlist <- al <- list()
tclist <- NULL


temp <- fwtnp2(x, f, pred, neigh, int, clo, keep,do.W=FALSE,varonly=TRUE)
out <- temp
newcoeff <- out$coeff
lr <- out$lengthsremove
lca <- out$lca
rem <- lca[, 1]

if(d>0){
	indsd<-sqrt(temp$v)
        norcoeff <- out$coeff/indsd
        al <- artlev(lr, rem)
        levno <- length(al)
        for (i in 1:levno) {
            ndetlist[[i]] <- norcoeff[al[[i]]]
        }
        sd <- mad(ndetlist[[1]])
        if (sd == 0) {
            sd <- NA
        }
        for (i in 1:levno) {
		tclist <- ebayesthresh(ndetlist[[i]], prior = "cauchy", a = NA, sdev = sd, threshrule = rule)
		newcoeff[al[[i]]] <- tclist * indsd[al[[i]]]
        }
        newcoeff[out$pointsin] <- out$coeff[out$pointsin]
}
neighbrs<-list()

for(kk in 1:length(rem)){
	neighbrs[[kk]]<-.C("nbrsfromlca",as.double(t(lca)),as.integer(ncol(lca)),as.integer(kk),
			n=as.integer(rep(0,times=lca[kk,2])),PACKAGE="adlift")$n
}

adds <- findadds(rem, neighbrs, out$po, index)
add <- max(adds)
if (verbose) {
        cat("doing ", add, " steps...\n", sep = "")
}
fhat <- invtnp2(x, newcoeff, out$lengths, lr, out$pointsin, lca, add)

if (returnall == FALSE) {
        return(fhat$coeff)
}
else {
        return(list(fhat=fhat,indsd = indsd, al = al, sd = sd))
}

}

