`corsup.fnc` <-
function(corres, bycol = TRUE, supp, plot = TRUE, font = 3, labels="", cex = 1.0) {

	if (!is(corres, "corres")) 
		stop("argument should be a correspondence object")

	xtab = corres@data$input
	
	if (bycol) {

		csupp = supp

		if (length(dim(csupp)) == 0) supName = "SUPP"
		else supName = colnames(csupp)

        # Projections of supplementary columns.  One or more 
        # supplementary columns can be used.  
        # Example of use on gobelets data (25 x 6):
        # suppl <- caSuppCol(gobelets[,-3], gobelets[,3])
        # This uses all except column 3 as principal; and 
        # column 3 as supplementary.
        # FM, 2003/12.

        tot <- sum(xtab)
        fIJ <- xtab/tot
        fI <- apply(fIJ, 1, sum)
        fJ <- apply(fIJ, 2, sum)
        fJsupI <- sweep(fIJ, 1, fI, FUN="/")
        fIsupJ <- sweep(fIJ, 2, fJ, FUN="/")
        s <- as.matrix(t(fJsupI)) %*% as.matrix(fIJ)
        s1 <- sweep(s, 1, sqrt(fJ), FUN="/")
        s2 <- sweep(s1, 2, sqrt(fJ), FUN="/")
        sres <- eigen(s2)
        sres$values[sres$values < 1.0e-8] <- 0.0
        # cat("Eigenvalues follow (trivial first eigenvalue removed).\n")
        # cat(sres$values[-1], "\n")
        # cat("Eigenvalue rate, in thousandths.\n")
        # tot <- sum(sres$values[-1])
        # cat(1000*sres$values[-1]/tot,"\n")
        # Eigenvectors divided rowwise by sqrt(fJ):
        evectors <- sweep(sres$vectors, 1, sqrt(fJ), FUN="/")
        rproj <- as.matrix(fJsupI) %*% evectors

        # Note: we must coerce csupp to matrix type, which 
        # propagates to csuppIJ
        csuppIJ <- as.matrix(csupp)/tot
        if (ncol(csuppIJ) > 1) csuppJ <- apply(csuppIJ, 2, sum)
        if (ncol(csuppIJ) == 1) csuppJ <- sum(csuppIJ)
        csuppproj <- t(csuppIJ) %*% rproj
        temp <- csuppproj
        # Divide rows by mass; and then cols. by sqrt of evals.
        csuppproj <- sweep ( sweep(temp,1,csuppJ,FUN="/"),2,
             sqrt(sres$values),FUN="/")

        # Value of this function on return: table of projections,
        # rows = set of supplementary columns; columns = set of factors.
        # (More than 1 supplementary column => labels will be retained.)
        # Adjust for trivial factor.
		if (plot) {
			res = csuppproj[,-1]
			if (length(supName) == 1) {
				text(res[1], res[2], supName, font=font,cex=cex)
			} else {
				if (length(labels) > 1) 
					supName = labels
				text(res[,1], res[,2], supName, font=font,cex=cex)
			}
		} else {
        	return(csuppproj[,-1])
		}
	} else { 

		rsupp = supp

		if (length(dim(rsupp)) == 0) supName = "SUPP"
		else supName = rownames(rsupp)

       # Projections of supplementary rows.  One or more supplementary rows
       # can be used.  Example of use:
       # x <- read.table("c:/mandible77s.dat")                 # Read data
       # xca <- ca(x[1:77,])                                   # Corr. analysis
       # xcar <- caSuppRow(x[1:77,], x[78:86,])                # Suppl. rows
       # plot(c(xca$rproj[,1],xca$cproj[,1]),                  # Prepare plot
       #       c(xca$rproj[,2],xca$cproj[,2]),
       #       type="n", xlab="Factor 1 (50.0% of inertia)",
       #       ylab="Factor 2 (21.0% of inertia)")
       # text(xca$rproj[,1],xca$rproj[,2],dimnames(x)[[1]])    # Plot prin. rows
       # text(xca$cproj[,1],xca$cproj[,2],dimnames(x)[[2]],font=4) # Plot cols.
       # text(xcar[,1],xcar[,2],dimnames(xcar)[[1]],font=3)    # Plot supp. rows
       # title("77 mandibles, 9 supplementary (groups), crossed by 9 attributes")
       # FM, 2003/12.

       tot <- sum(xtab)
       fIJ <- xtab/tot
       fI <- apply(fIJ, 1, sum)
       fJ <- apply(fIJ, 2, sum)
       fJsupI <- sweep(fIJ, 1, fI, FUN="/")
       fIsupJ <- sweep(fIJ, 2, fJ, FUN="/")
       s <- as.matrix(t(fJsupI)) %*% as.matrix(fIJ)
       s1 <- sweep(s, 1, sqrt(fJ), FUN="/")
       s2 <- sweep(s1, 2, sqrt(fJ), FUN="/")
       sres <- eigen(s2)
       sres$values[sres$values < 1.0e-8] <- 0.0
       # cat("Eigenvalues follow (trivial first eigenvalue removed).\n")
       # cat(sres$values[-1], "\n")
       # cat("Eigenvalue rate, in thousandths.\n")
       # tot <- sum(sres$values[-1])
       # cat(1000*sres$values[-1]/tot,"\n")
       # Eigenvectors divided rowwise by sqrt(fJ):
       evectors <- sweep(sres$vectors, 1, sqrt(fJ), FUN="/")
       # rproj <- as.matrix(fJsupI) %*% evectors
       temp  <- as.matrix(s2) %*% sres$vectors
       # Following divides rowwise by sqrt(fJ) and 
       # columnwise by sqrt(eigenvalues):
       # Note: first column of cproj is trivially 1-valued.
       cproj <- sweep ( sweep(temp,1,sqrt(fJ),FUN="/"), 2,
                        sqrt(sres$values),FUN="/")

       # Note: we must coerce rsupp to matrix type, which 
       # propagates to rsuppIJ
       rsuppIJ <- as.matrix(rsupp)/tot
       if (nrow(rsuppIJ) > 1) rsuppI <- apply(rsuppIJ, 1, sum)
       if (nrow(rsuppIJ) == 1) rsuppI <- sum(rsuppIJ)
       
       rsuppproj <- rsuppIJ %*% cproj 
       temp <- rsuppproj
       # Divide cols. by mass; and then rows. by sqrt of evals.
       rsuppproj <- sweep ( sweep(temp,1,rsuppI,FUN="/"),2,
                     sqrt(sres$values),FUN="/")

       # Value of this function on return: table of projections,
       # rows = set of supplementary rows; columns = set of factors.
       # (More than 1 supplementary rows => labels will be retained.)
       # Adjust for trivial factor.
	   if (plot) {
			res = rsuppproj[,-1]
			if (length(supName) == 1) {
				text(res[1], res[2], supName, font=font,cex=cex)
			} else {
				if (length(labels) > 1) 
					supName = labels
				text(res[,1], res[,2], supName, font=font,cex=cex)
			}
			# text(rsupproj[1,-1], rsupproj[2,-1], supName, font=3)
	   } else {
            return(rsuppproj[,-1])
	   }
	}
}

