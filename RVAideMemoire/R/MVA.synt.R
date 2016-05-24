# pls: R2

#  - Total variance:
#     * PCA (dudi.pca[ade4],prcomp[stats],princomp[stats],pca[mixOmics],pca[labdsv],rda[vegan])
#     * sPCA (spca[mixOmics])
#     * PCoA (dudi.pco[ade4],pcoa[ape],pco[labdsv],cmdscale[stats]°,wcmdscale[vegan]°,capscale[vegan])
#		° if computed with 'eig=TRUE'
#	  Does not take into account imaginary axes
#
#  - Intergroup variance:
#     * LDA (lda[MASS],discrimin[ade4])
#     * PLS-DA (plsda[mixOmics])
#     * CDA (discrimin[ade4],discrimin.coa[ade4])
#
#  - X variance + Y variance:
#	* CPPLS (mvr[pls])
#	* PLSR (mvr[pls],plsR[plsRglm] !que Y!, !!! pls[mixOmics])
#	* PCR (mvr[pls])
#
#  - Kurtosis:
#	* IPCA (ipca[mixOmics])
#	* sIPCA (sipca[mixOmics])
#
#  - Stress:
#	* nMDS (isoMDS[MASS],monoMDS[vegan],metaMDS[vegan],nmds[labdsv])


print.MVA.synt <- function(x,...) {
  for (i in 1:(length(x)-1)) {
    x.i <- x[[i]]
    cat(paste0("Criterion: ",x.i$crit,"\n"))
    if ("tab" %in% names(x.i)) {
	if (x$rows>nrow(x.i$tab)) {
	  rows <- nrow(x.i$tab)
	} else {
	  rows <- x$rows
	}
	tab <- x.i$tab
	if ("Proportion" %in% colnames(tab)) {tab[,"Proportion"] <- round(tab[,"Proportion"],2)}
	if ("Cumulative" %in% colnames(tab)) {tab[,"Cumulative"] <- round(tab[,"Cumulative"],2)}
	if ("Kurtosis" %in% colnames(tab)) {tab[,"Kurtosis"] <- round(tab[,"Kurtosis"],2)}
	print(tab[1:rows,],row.names=FALSE)
    }
    if ("stress" %in% names(x.i)) {
	cat(paste0("Stress: ",round(x.i$stress,4),"\n"))
    }
    cat("\n")
  }
}

MVA.synt <- function(x,rows=5) {
  x <- MVA.ident(x)
  res <- MVA.get.synt(x)
  res$rows <- rows
  class(res) <- "MVA.synt"
  return(res)
}

MVA.get.synt <- function(x,...) {
  UseMethod("MVA.get.synt")
}

MVA.get.synt.default <- function(x,...) {
  stop("unknown multivariate analysis")
}

MVA.get.synt.PCA.ade4 <- MVA.get.synt.PCoA.ade4 <- function(x,...) {
  res <- list()
  vars <- x$eig
  vars <- vars/sum(vars)
  vars.each <- 100*vars
  vars.cum <- 100*cumsum(vars)
  tab <- data.frame(Axis=1:length(vars),Proportion=vars.each,Cumulative=vars.cum)
  res[[1]] <- list(crit="total variance (%)",tab=tab)
  return(res)
}

MVA.get.synt.PCA.prcomp <- function(x,...) {
  res <- list()
  vars <- summary(x)$importance
  vars.each <- 100*vars["Proportion of Variance",]
  vars.cum <- 100*vars["Cumulative Proportion",]
  tab <- data.frame(Axis=1:ncol(vars),Proportion=vars.each,Cumulative=vars.cum)
  res[[1]] <- list(crit="total variance (%)",tab=tab)
  return(res)
}

MVA.get.synt.PCA.princomp <- function(x,...) {
  res <- list()
  vars <- summary(x)$sdev^2
  vars <- vars/sum(vars)
  vars.each <- 100*vars
  vars.cum <- 100*cumsum(vars)
  tab <- data.frame(Axis=1:length(vars),Proportion=vars.each,Cumulative=vars.cum)
  res[[1]] <- list(crit="total variance (%)",tab=tab)
  return(res)
}

MVA.get.synt.PCA.mixOmics <- function(x,...) {
  res <- list()
  comp <- x$ncomp
  comp.max <- min(ncol(x$X),nrow(x$X)-1)
  x <- update(x,ncomp=comp.max)
  vars <- x$sdev
  vars <- vars/sum(vars)
  vars.each <- 100*vars
  vars.cum <- 100*cumsum(vars)
  tab <- data.frame(Axis=1:length(vars),Proportion=vars.each,Cumulative=vars.cum)
  tab <- tab[1:comp,]
  res[[1]] <- list(crit="total variance (%)",tab=tab)
  return(res)
}

MVA.get.synt.PCA.labdsv <- function(x,...) {
  res <- list()
  vars <- x$sdev^2
  vars <- vars/x$totdev
  vars.each <- 100*vars
  vars.cum <- 100*cumsum(vars)
  tab <- data.frame(Axis=1:length(vars),Proportion=vars.each,Cumulative=vars.cum)
  res[[1]] <- list(crit="total variance (%)",tab=tab)
  return(res)
}

MVA.get.synt.PCA.vegan <- function(x,...) {
  res <- list()
  vars <- x$CA$eig
  vars <- vars/sum(vars)
  vars.each <- 100*vars
  vars.cum <- 100*cumsum(vars)
  tab <- data.frame(Axis=1:length(vars),Proportion=vars.each,Cumulative=vars.cum)
  res[[1]] <- list(crit="total variance (%)",tab=tab)
  return(res)
}

MVA.get.synt.sPCA.mixOmics <- function(x,...) {
  res <- list()
  vars.each <- 100*c(x$varX[1],diff(x$varX))
  vars.cum <- 100*x$varX
  tab <- data.frame(Axis=1:length(vars.each),Proportion=vars.each,Cumulative=vars.cum)
  res[[1]] <- list(crit="total variance (%)",tab=tab)
  return(res)
}

MVA.get.synt.IPCA.mixOmics <- MVA.get.synt.sIPCA.mixOmics <- function(x,...) {
  res <- list()
  kurt <- x$kurtosis
  tab <- data.frame(Axis=1:length(kurt),Kurtosis=kurt)
  res[[1]] <- list(crit="kurtosis",tab=tab)
  return(res)
}

MVA.get.synt.PCoA.ape <- function(x,...) {
  res <- list()
  val <- if ("Cum_corr_eig" %in% colnames(x$values)) {
    x$values$Cum_corr_eig
  } else {
    x$values$Cumul_eig
  }
  vars.each <- 100*c(val[1],diff(val))
  vars.cum <- 100*val
  tab <- data.frame(Axis=1:length(vars.each),Proportion=vars.each,Cumulative=vars.cum)
  res[[1]] <- list(crit="total variance (%)",tab=tab)
  return(res)
}

MVA.get.synt.PCoA.labdsv <- function(x,...) {
  res <- list()
  vars <- x$eig
  vars[vars<0] <- 0
  vars <- vars/sum(vars)
  vars.each <- 100*vars
  vars.cum <- 100*cumsum(vars)
  tab <- data.frame(Axis=1:length(vars),Proportion=vars.each,Cumulative=vars.cum)
  res[[1]] <- list(crit="total variance (%)",tab=tab)
  return(res)
}

MVA.get.synt.PCoA.stats <- function(x,...) {
  if (!"eig" %in% names(x) || is.null(x$eig)) {
    stop("no eigenvalues available, compute the analysis with 'eig=TRUE'")
  }
  res <- list()
  vars <- x$eig
  vars[vars<0] <- 0
  vars <- vars/sum(vars)
  vars.each <- 100*vars
  vars.cum <- 100*cumsum(vars)
  tab <- data.frame(Axis=1:length(vars),Proportion=vars.each,Cumulative=vars.cum)
  res[[1]] <- list(crit="total variance (%)",tab=tab)
  return(res)
}

MVA.get.synt.PCoA.vegan <- function(x,...) {
  res <- list()
  if ("CA" %in% names(x)) {
    vars <- x$CA$eig

  } else {
    if (!"eig" %in% names(x) || is.null(x$eig)) {
	stop("no eigenvalues available, compute the analysis with 'eig=TRUE'")
    }
    vars <- x$eig
    vars[vars<0] <- 0
  }
  vars <- vars/sum(vars)
  vars.each <- 100*vars
  vars.cum <- 100*cumsum(vars)
  tab <- data.frame(Axis=1:length(vars),Proportion=vars.each,Cumulative=vars.cum)
  res[[1]] <- list(crit="total variance (%)",tab=tab)
  return(res)
}

MVA.get.synt.nMDS.mono.vegan <- function(x,...) {
  res <- list()
  res[[1]] <- list(crit="stress",stress=x$stress)
  return(res)
}

MVA.get.synt.nMDS.iso.vegan <- MVA.get.synt.nMDS.MASS <-
MVA.get.synt.nMDS.labdsv <- function(x,...) {
  res <- list()
  res[[1]] <- list(crit="stress",stress=x$stress/100)
  return(res)
}

MVA.get.synt.LDA.MASS <- function(x,...) {
  res <- list()
  form <- LDA.format(x)
  sco <- form$li
#  tot.inert <- inertia(sco)
#  tot.ax <- apply(sco,2,function(y) inertia(as.data.frame(y)))
#  tot.vars <- tot.ax/tot.inert
#  tot.vars.each <- 100*tot.vars
#  tot.vars.cum <- 100*cumsum(tot.vars)
#  tot.tab <- data.frame(Axis=1:length(tot.vars),Proportion=tot.vars.each,Cumulative=tot.vars.cum)
#  res[[1]] <- list(crit="total variance (%)",tab=tot.tab)
  Y <- form$grouping
  gp.n <- table(Y)
  gp.means <- as.data.frame(aggregate(sco~Y,FUN=mean)[,-1])
  gp.inert <- inertia(gp.means,w=gp.n)
  gp.ax <- apply(gp.means,2,function(y) inertia(as.data.frame(y),w=gp.n))
  gp.vars <- gp.ax/gp.inert
  gp.vars.each <- 100*gp.vars
  gp.vars.cum <- 100*cumsum(gp.vars)
  gp.tab <- data.frame(Axis=1:length(gp.vars),Proportion=gp.vars.each,Cumulative=gp.vars.cum)
  res[[1]] <- list(crit="intergroup variance (%)",tab=gp.tab)
  return(res)
}

MVA.get.synt.LDA.ade4 <- MVA.get.synt.CDA.ade4 <- function(x,...) {
  res <- list()
  Y <- eval(x$call$fac)
  x <- update(x,scannf=FALSE,nf=nlevels(Y)-1)
  sco <- x$li
#  tot.inert <- inertia(sco)
#  tot.ax <- apply(sco,2,function(y) inertia(as.data.frame(y)))
#  tot.vars <- tot.ax/tot.inert
#  tot.vars.each <- 100*tot.vars
#  tot.vars.cum <- 100*cumsum(tot.vars)
#  tot.tab <- data.frame(Axis=1:length(tot.vars),Proportion=tot.vars.each,Cumulative=tot.vars.cum)
#  res[[1]] <- list(crit="total variance (%)",tab=tot.tab)
  gp.n <- table(Y)
  gp.means <- as.data.frame(aggregate(as.matrix(sco)~Y,FUN=mean)[,-1])
  gp.inert <- inertia(gp.means,w=gp.n)
  gp.ax <- apply(gp.means,2,function(y) inertia(as.data.frame(y),w=gp.n))
  gp.vars <- gp.ax/gp.inert
  gp.vars.each <- 100*gp.vars
  gp.vars.cum <- 100*cumsum(gp.vars)
  gp.tab <- data.frame(Axis=1:length(gp.vars),Proportion=gp.vars.each,Cumulative=gp.vars.cum)
  res[[1]] <- list(crit="intergroup variance (%)",tab=gp.tab)
  return(res)
}

MVA.get.synt.PLSDA.mixOmics <- function(x,...) {
  res <- list()
  nco <- ncol(x$X)
  while(inherits(try(update(x,ncomp=nco),silent=TRUE),"try-error")) {
    nco <- nco-1
  }
  x <- update(x,ncomp=nco)
  sco <- x$variates$X
#  tot.inert <- inertia(sco)
#  tot.ax <- apply(sco,2,function(y) inertia(as.data.frame(y)))
#  tot.vars <- tot.ax/tot.inert
#  tot.vars.each <- 100*tot.vars
#  tot.vars.cum <- 100*cumsum(tot.vars)
#  tot.tab <- data.frame(Axis=1:length(tot.vars),Proportion=tot.vars.each,Cumulative=tot.vars.cum)
#  res[[1]] <- list(crit="total variance (%)",tab=tot.tab)
  Y <- eval(x$call$Y)
  gp.n <- table(Y)
  gp.means <- as.data.frame(aggregate(as.matrix(sco)~Y,FUN=mean)[,-1])
  gp.inert <- inertia(gp.means,w=gp.n)
  gp.ax <- apply(gp.means,2,function(y) inertia(as.data.frame(y),w=gp.n))
  gp.vars <- gp.ax/gp.inert
  gp.vars.each <- 100*gp.vars
  gp.vars.cum <- 100*cumsum(gp.vars)
  gp.tab <- data.frame(Axis=1:length(gp.vars),Proportion=gp.vars.each,Cumulative=gp.vars.cum)
  res[[1]] <- list(crit="intergroup variance (%)",tab=gp.tab)
  return(res)
}

MVA.get.synt.CPPLS.pls <- MVA.get.synt.PLSR.pls <-
MVA.get.synt.PCR.pls <- function(x,...) {
  res <- list()
  Y.vars <- 100*drop(pls::R2(x,intercept=FALSE,estimate="train")$val)
  if (is.matrix(Y.vars)) {
    crit <- "Y cumulative total variance (%)"
    Y.vars <- round(as.data.frame(t(Y.vars)),2)
    Y.tab <- as.data.frame(cbind(Axis=1:nrow(Y.vars),Y.vars))
  } else {
    crit <- "Y total variance (%)"
    Y.vars.each <- c(Y.vars[1],diff(Y.vars))
    Y.vars.cum <- Y.vars
    Y.tab <- data.frame(Axis=1:length(Y.vars),Proportion=Y.vars.each,Cumulative=Y.vars.cum)
  }
  res[[1]] <- list(crit=crit,tab=Y.tab)
  X.vars <- x$Xvar/x$Xtotvar
  X.vars.each <- 100*X.vars
  X.vars.cum <- 100*cumsum(X.vars)
  X.tab <- data.frame(Axis=1:length(X.vars),Proportion=X.vars.each,Cumulative=X.vars.cum)
  res[[2]] <- list(crit="X total variance (%)",tab=X.tab)
  return(res)
}

MVA.get.synt.PLSR.plsRglm <- function(x,...) {		# !!!! manque la variance de X
  res <- list()
  Y.vars <- as.vector(cor(x$dataY,x$tt,use="pairwise")^2)
  Y.vars.each <- 100*Y.vars
  Y.vars.cum <- 100*cumsum(Y.vars)
  Y.tab <- data.frame(Axis=1:length(Y.vars),Proportion=Y.vars.each,Cumulative=Y.vars.cum)
  res[[1]] <- list(crit="Y total variance (%)",tab=Y.tab)
#  X.vars <- x$Xvar/x$Xtotvar
#  X.vars.each <- 100*X.vars
#  X.vars.cum <- 100*cumsum(X.vars)
#  X.tab <- data.frame(Axis=1:length(X.vars),Proportion=X.vars.each,Cumulative=X.vars.cum)
#  res[[2]] <- list(crit="X total variance (%)",tab=X.tab)
  return(res)
}



