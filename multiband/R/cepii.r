#' Light curve data from two bands
#' 
#' A dataset containing band I and V data from OGLE-LMC-T2CEP-009.
#' There are two data frames, iband and vband, each of which has three columns. (time,magnitude,magnitude error).
#' 
#' 
#' @docType data
#' @keywords datasets
#' @format Two data frames with 3 variables, and 554 (iband) and 65 (vband) rows.
#' @name cepii
#' @examples
#' ## Load I and V bands
#' iband <- cepii[[1]]
#' vband <- cepii[[2]]
#' 
#' ## Try finer grid sizes as well, e.g. 'nOmega <- 500'
#' nOmega <- 10
#' omega_seq <- seq(3.5,3.65,length.out=nOmega)
#' sol_ls <- vector(mode="list",length=nOmega)
#' RSS_ls <- double(nOmega)
#' 
#' ## Drastically subsample the data and see if the methods can find the period.
#' subi <- seq(1,nrow(iband),by=20)
#' 
#' ## 1. Lomb Scargle on I band and V band separately
#' for (i in 1:nOmega) {
#'   sol_ls[[i]] <- lomb_scargle(iband[subi,1],iband[subi,2],iband[subi,3],omega_seq[i])
#'   RSS_ls[i] <- sol_ls[[i]]$RSS
#' }
#' plot(omega_seq,RSS_ls,xlab=expression(omega),ylab='RSS',main='I band',pch=16)
#' 
#' subv <- seq(1,nrow(vband),by=4)
#' for (i in 1:nOmega) {
#'   sol_ls[[i]] <- lomb_scargle(vband[subv,1],vband[subv,2],vband[subv,3],omega_seq[i])
#'   RSS_ls[i] <- sol_ls[[i]]$RSS
#' }
#' plot(omega_seq,RSS_ls,xlab=expression(omega),ylab='RSS',main='V band',pch=16)
#' 
#' ## 2. Naive pooled Lomb Scargle versus Fused
#' tms <- vector(mode="list",length=2)
#' tms[[1]] <- iband[subi,,drop=FALSE]
#' tms[[2]] <- vband[subv,,drop=FALSE]
#' B <- length(tms)
#' t <- c(); m <- c(); sigma <- c()
#' for (b in 1:B) {
#'   t <- c(t,tms[[b]][,1])
#'   m <- c(m,tms[[b]][,2])
#'   sigma <- c(sigma,tms[[b]][,3])
#' }
#' 
#' sol_ls <- vector(mode="list",length=nOmega)
#' sol_bcd <- vector(mode="list",length=nOmega)
#' RSS_seq <- double(nOmega)
#' RSS_ls_seq <- double(nOmega)
#' gamma1 <- 1000
#' gamma2 <- 10
#' for (i in 1:nOmega) {
#'   sol_ls[[i]] <- lomb_scargle(t,m,sigma,omega_seq[i])
#'   RSS_ls_seq[i] <- sol_ls[[i]]$RSS
#'   beta0_ls <- rep(sol_ls[[i]]$beta0,B)
#'   A_ls <- rep(sol_ls[[i]]$A,B)
#'   rho_ls <- rep(sol_ls[[i]]$rho,B)
#'   sol_bcd[[i]] <- bcd_inexact(tms,beta0_ls,A_ls,rep(1,B),rho_ls,omega_seq[i],gamma1,gamma2,
#'     max_iter=1e4,tol=1e-10)
#'   RSS_seq[i] <- pnll(tms,sol_bcd[[i]]$beta0,sol_bcd[[i]]$A,rep(1,B),
#'     sol_bcd[[i]]$rho,omega_seq[i],0,0)
#'   print(paste0("Completed ", i))
#' }
#' 
#' plot(omega_seq,RSS_seq,xlab=expression(omega),ylab='RSS',main='I & V band fusion',pch=16)
#' plot(omega_seq,RSS_ls_seq,xlab=expression(omega),ylab='RSS',main='naive Lomb-Scargle',pch=16)
#' ix_min <- which(RSS_seq==min(RSS_seq))
#' sol_bcd_final <- sol_bcd[[ix_min]]
NULL
