#diffAmova.R
#package: msap
#Author: Andrés Pérez-Figueroa (anpefi@uvigo.es)

#Uses pegas' AMOVA and report differentiation

diffAmova <- function(DM, groups, nDec, pairwise){
  
  amova2<- function(formula, data = NULL, nperm = 1000, is.squared = FALSE)
  {
    ## THIS IS A MDOFIED VERSION OF amova {pegas} for evluating formula in the parent environment, as the original amova in pegas needs the formula memeber to be in the Global Environment (and this caused some restrictions in the CRAN)
    y.nms <- as.character(as.expression(formula[[2]]))
    rhs <- formula[[3]]
    gr.nms <- as.character(as.expression(rhs))
    ## keep the highest level first:
    if (length(rhs) > 1) gr.nms <- unlist(strsplit(gr.nms, "/"))
    
    y <- get(y.nms, envir=parent.env(environment()))
    if (!is.squared) y <- y^2 # square the distances
    if (class(y) == "dist") y <- as.matrix(y)
    if (!is.matrix(y))
      stop("the lhs of the formula must be either a matrix or an object of class 'dist'.")
    n <- dim(y)[1] # number of individuals
    
    gr <-
      if (is.null(data)) as.data.frame(sapply(gr.nms, get, envir = parent.env(environment())))
    else data[gr.nms]
    Nlv <- length(gr) # number of levels
    
    ### 5 local functions
    ## a simplified version of tapply(X, INDEX, FUN = sum):
    foo <- function(x, index) unlist(lapply(split(x, index), sum))
    
    getSSD <- function(y, gr, Nlv, N, n) {
      SSD <- numeric(Nlv + 2)
      SSD[Nlv + 2] <- sum(y/(2 * n)) # total SSD
      ## calculate SSD *within* each level:
      for (i in 1:Nlv) {
        p <- gr[, i] # extract the grouping at this level
        SSD[i + 1] <- sum((y/(2 * N[[i]])[p])[outer(p, p, "==")])
      }
      ## now differentiate to get the SSDs:
      if (Nlv > 1)
        for (i in Nlv:2) SSD[i] <- SSD[i] - SSD[i + 1]
      SSD[1] <- SSD[Nlv + 2] - sum(SSD[-(Nlv + 2)])
      SSD
    }
    
    getDF <- function(gr, Nlv, N, n) {
      df <- numeric(Nlv + 2)
      df[1:Nlv] <- unlist(lapply(N, length)) # number of groups
      ## get the df from the lowest level to the highest one:
      df[Nlv + 1] <- n
      for (i in (Nlv + 1):2) df[i] <-  df[i] - df[i - 1]
      df[1] <- df[1] - 1
      df[Nlv + 2] <- n - 1 # total df
      df
    }
    
    getNcoefficient <- function(gr, Nlv, N, n) {
      Nig <- N[[Nlv]] # numbers from the lowest level (ie, pop)
      Nig2 <- Nig^2
      npop <- length(Nig)
      if (Nlv == 1) # do classic variance components estimation
        ncoef <- (n - sum(Nig2)/n)/(npop - 1)
      else {
        if (Nlv == 2) { # use eq.9a-c in Excoffier et al. 1992
          ncoef <- numeric(3) # n, n', and n'', respectively
          G <- nlevels(gr[, 1]) # the number of groups
          ## use the fact that match() returns the 1st match of the
          ## 1st argument into the 2nd one, so this returns a factor
          ## indicating to which group each pop belongs to
          g <- gr[, 1][match(1:npop, as.integer(gr[, 2]))]
          npopBYgr <- tabulate(g) # number of pops in each group
          A <- sum(foo(Nig2, g)/foo(Nig, g))
          ncoef[1] <- (n - A)/sum(npopBYgr - 1) # corrects eq.9a in Excoffier et al. 1992
          ncoef[2] <- (A - sum(Nig2)/n)/(G - 1)
          ncoef[3] <- (n - sum(foo(Nig, g)^2/n))/(G - 1)
        } else { # use approximate formulae
          ncoef <- numeric(Nlv + 1)
          ## the coefficients are indexed from the highest to the lowest
          ## level to ease computation of the variance components
          ncoef[Nlv] <- (n - sum(Nig2)/n)/(npop - 1) # same than above
          ncoef[Nlv + 1] <- 1 # for the within pop level
          for (i in 1:(Nlv - 1)) {
            group <- gr[, i]
            g <- group[match(1:npop, as.integer(gr[, i + 1]))]
            A <- sum(foo(Nig, g)^2)/sum(foo(Nig, g))
            ncoef[i] <- (n - A)/(nlevels(group) - 1)
          }
        }
      }
      names(ncoef) <- letters[1:length(ncoef)] # suggestion by Rodney Dyer
      ncoef
    }
    
    getVarComp <- function(MSD, Nlv, ncoef) {
      ## get the variance components bottom-up
      if (Nlv == 1)
        sigma2 <- c((MSD[1] - MSD[2])/ncoef, MSD[2]) # 'n' changed to 'ncoef' (fix by Qixin He, 2010-07-15)
      else {
        sigma2 <- numeric(Nlv + 1)
        if (Nlv == 2) {
          sigma2[3] <- MSD[3]
          sigma2[2] <- (MSD[2] - sigma2[3])/ncoef[1]
          sigma2[1] <- (MSD[1] - MSD[3] - ncoef[2]*sigma2[2])/ncoef[3]
        } else {
          sigma2[Nlv + 1] <- MSD[Nlv + 1]
          for (i in Nlv:1) {
            sel <- i:(Nlv + 1)
            sigma2[i] <- (MSD[i] - sum(ncoef[sel]*sigma2[sel]))/ncoef[i]
          }
        }
      }
      names(sigma2) <- c(names(gr), "Error")
      sigma2
    }
    
    ## fit the AMOVA model to the data:
    N <- lapply(gr, tabulate) # number of individuals in each group
    SSD <- getSSD(y, gr, Nlv, N, n)
    df <- getDF(gr, Nlv, N, n)
    MSD <- SSD/df
    ncoef <- getNcoefficient(gr, Nlv, N, n)
    sigma2 <- getVarComp(MSD, Nlv, ncoef)
    
    ## output the results:
    res <- list(tab = data.frame(SSD = SSD, MSD = MSD, df = df,
                                 row.names = c(names(gr), "Error", "Total")),
                varcoef = ncoef, varcomp = sigma2, call = match.call())
    class(res) <- "amova"
    
    if (nperm) {
      rSigma2 <- matrix(0, nperm, length(sigma2))
      ## First shuffle all individuals among all pops to test
      ## the within pop var comp; if there is only one level,
      ## this is used to test both var components.
      j <- if (Nlv == 1) 1:2 else Nlv + 1
      for (i in 1:nperm) {
        rY <- perm.rowscols(y, n)
        ## the hierarchical structure is not changed so just
        ## need to recalculate the SSD and var comp
        rSSD <- getSSD(rY, gr, Nlv, N, n)
        rSigma2[i, j] <- getVarComp(rSSD/df, Nlv, ncoef)[j]
      }
      if (Nlv > 1) {
        ## for the lowest level, we just permute individuals within the
        ## level just above, so the hierarchical structure is unchanged
        j <- Nlv
        ## 'L' contains the indices of the individuals for each level just above
        L <- lapply(levels(gr[, j - 1]), function(x) which(gr[, j - 1] == x))
        for (i in 1:nperm) {
          rind <- unlist(lapply(L, sample))
          rY <- y[rind, rind]
          rSSD <- getSSD(rY, gr, Nlv, N, n)
          rSigma2[i, j] <- getVarComp(rSSD/df, Nlv, ncoef)[j]
        }
        if (Nlv > 2) {
          for (j in (Nlv - 1):2) {
            above <- gr[, j - 1]
            L <- lapply(levels(above), function(x) which(above == x))
            for (i in 1:nperm) {
              rind <- integer(0)
              for (k in L) rind <- c(rind, sample(k))
              rind <- unlist(lapply(L, sample))
              rY <- y[rind, rind]
              rGR <- gr[rind, ]
              rN <- lapply(rGR, tabulate)
              rSSD <- getSSD(rY, rGR, Nlv, rN, n)
              rDF <- getDF(rGR, Nlv, rN, n)
              rNcoef <- getNcoefficient(rGR, Nlv, rN, n)
              rSigma2[i, j] <- getVarComp(rSSD/rDF, Nlv, rNcoef)[j]
            }
          }
        }
        ## for the highest level permute the level below
        L <- lapply(levels(gr[, 1]), function(x) which(gr[, 1] == x))
        for (i in 1:nperm) {
          rind <- unlist(sample(L))
          rY <- y[rind, rind]
          rGR <- gr[rind, ]
          rN <- lapply(rGR, tabulate)
          rSSD <- getSSD(rY, rGR, Nlv, rN, n)
          rDF <- getDF(rGR, Nlv, rN, n)
          rNcoef <- getNcoefficient(rGR, Nlv, rN, n)
          rSigma2[i, 1] <- getVarComp(rSSD/rDF, Nlv, rNcoef)[1]
        }
      }
      P <- numeric(Nlv + 1)
      for (j in 1:(Nlv + 1))
        P[j] <- sum(rSigma2[, j] >= sigma2[j])/nperm
      P[Nlv + 1] <- NA
      res$varcomp <- data.frame(sigma2 = res$varcomp, P.value = P)
    }
    res
  }
  
	
	ntt <- length(levels(groups))

	#assign("DM", DM, envir=globalenv())
	#assign("groups", groups, envir=globalenv())
	
  amv<-amova2(DM ~ groups, data=data.frame(groups), nperm=10000) #from pegas

	cat("AMOVA TABLE \td.f. \tSSD \t\tMSD \t\tVariance\n")
	phiST <- as.numeric(amv$varcomp[1,1]/(amv$varcomp[1,1]+amv$varcomp[2,1]))
	cat("among groups\t",amv$tab[1,3],"\t",format(amv$tab[1,1], digits=nDec),"\t",format(amv$tab[1,2], digits=nDec),"\t",format(amv$varcomp[1,1], digits=nDec),"\n")
	cat("within groups\t",amv$tab[2,3],"\t",format(amv$tab[2,1], digits=nDec),"\t",format(amv$tab[2,2], digits=nDec),"\t",format(amv$varcomp[2,1], digits=nDec),"\n")
	cat("Total        \t",amv$tab[3,3],"\t",format(amv$tab[3,1], digits=nDec),"\t",format(amv$tab[3,2], digits=nDec), " \n")
	cat("\n")
	pval<- amv$varcomp[1,2]
	pval<- ifelse(pval<0.0001, "(P<0.0001)", paste("(P=",format(pval,digits=nDec),")"))
	cat("Phi_ST = ", format(phiST, digits=nDec), " ", pval,"\n")
	
	if(ntt<3) pairwise=FALSE
	if(pairwise){
		cat("\nPairwise Phi_ST\n------------------------------------\n")

		for (i in 1:(ntt-1) ){
			for(j in seq(i+1,ntt ) ) {
				p1 <- which(groups==levels(groups)[i])
				p2 <- which(groups==levels(groups)[j])
				ttos <- groups[c(p1,p2)]
				M <- subset(DM,c(p1,p2))
				a <- amova2(M ~ ttos, nperm=10000)
				pval<-as.numeric(a$varcomp$P.value[1])
	
				pval<-ifelse(pval<0.0001, "(P<0.0001)", paste("(P=",format(pval, digits=nDec),")"))
				phiST <- as.numeric(a$varcomp[1,1]/(a$varcomp[1,1]+a$varcomp[2,1]))
				
				cat(levels(groups)[i], " - ",levels(groups)[j],": ",format(phiST,trim=T,digits=nDec),"\t ",pval," \n")
			}
		}
}#end if pairwise
  
  
  ###############################################
  
  

	
  
  #############################################
  
  
  
  
  
}




