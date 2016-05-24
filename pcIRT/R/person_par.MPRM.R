#' @param set0 if \code{set0=TRUE} for those raw scores patterns with 0
#' observations (except in the reference category) the person parameter value
#' is set minimal. With this procedure it is possible to estimate at least the
#' remaining person parameters of these raw score pattern.  Note: only relevant
#' for person parameter estimation of MPRM. The person parameters for each raw score vector are constrained to sum zero

#'@rdname perspar
#'@method person_par MPRM
#'@export


person_par.MPRM <-
function(object, ..., set0 = FALSE){

call <- match.call()

kateg.zahl <- length(table(object$data))
item.zahl <- ncol(object$data)

row.table1 <- apply(object$data+1, 1, function(s) tabulate(s,nbins=kateg.zahl))
row.table <- apply(object$data+1,1,function(x) sprintf("%04d",(tabulate(x,nbins=kateg.zahl))))
pat   <- apply(row.table,2, function(n) paste0(n, collapse=""))
patt  <- table(pat)

ind <- lapply(names(patt), function(hpat) which(pat %in% hpat))
rto <- row.table1[,unlist(lapply(ind, function(el) el[[1]]))]

col.table <- apply(object$data+1, 2, function(s) tabulate(s,nbins=kateg.zahl))

if(set0 == FALSE){
      out <- which(apply(rto,2, function(nu) {any(nu==0) | nu[kateg.zahl] == item.zahl}))
      rto2 <- rto[,-out]
      patt2 <- patt[-out]
      startv <- rep(0, nrow(rto2[-kateg.zahl,])*ncol(rto2))
      fixP  <- NULL
    } else {
      out <- which(apply(rto,2, function(nu) {nu[kateg.zahl]==0 | nu[kateg.zahl] == item.zahl}))
      rto2    <- rto[,-out]
      patt2   <- patt[-out]
      fixP    <- which(as.vector(rto2[-kateg.zahl,])==0)
      startv  <- rep(0, length(as.vector(rto2[-kateg.zahl,])[-fixP]))
    }

estpar <- object$itempar[-kateg.zahl, -item.zahl]*(-1)

pl <- function(perspv=startv, rto=rto2, patt=patt2, estpar=estpar, col.table=col.table, fixP=fixP){

newp <- numeric(length(perspv)+length(fixP))
newp[fixP] <- -1000
newp[newp != -1000] <- perspv

perspar <- cbind(matrix(newp, ncol=nrow(rto)-1, nrow=ncol(rto), byrow=TRUE),rep(0,ncol(rto)))

fir <- sum(rowSums(perspar*t(rto))*patt)

itmat <- cbind(rbind(matrix(estpar, nrow=nrow(rto)-1), rep(0, item.zahl-1)), rep(0, nrow(rto)))

sec <- sum(itmat*col.table)

dif <- log(colSums(exp(mapply(function(pp,it) {perspar[pp,]+itmat[,it]}, pp=rep(seq_len(ncol(rto)), each=item.zahl), it=rep(seq_len(ncol(itmat)), ncol(rto))))))

lauf <- seq(1,length(dif), by=ncol(itmat))
sdif <- sum(dif[1:(ncol(itmat))])
for (i in lauf[-1]){
  sdif <- c(sdif,sum(dif[i:(i+ncol(itmat)-1)]))
}
thir <- sum(as.vector(sdif*patt))

fir + sec - thir
}

der1pl <- function(perspv=startv, rto=rto2, patt=patt2, estpar=estpar, col.table=col.table, fixP=fixP){

  newp <- numeric(length(perspv)+length(fixP))
  newp[fixP] <- -1000
  newp[newp != -1000] <- perspv

  perspar <- cbind(matrix(newp, ncol=nrow(rto)-1, nrow=ncol(rto), byrow=TRUE),rep(0,ncol(rto)))

  itmat <- cbind(rbind(matrix(estpar, nrow=nrow(rto)-1), rep(0, item.zahl-1)), rep(0, nrow(rto)))

  difz <- exp(mapply(function(pp,it) {perspar[pp,]+itmat[,it]}, pp=rep(seq_len(ncol(rto)), each=item.zahl), it=rep(seq_len(ncol(itmat)), ncol(rto))))
  lauf <- seq(1,ncol(difz), by=ncol(itmat))
  difz2 <- colSums(difz)
  difz3 <- t(apply(difz, 1, function(ze) {ze/difz2}))
  za <- rowSums(difz3[,1:(ncol(itmat))])

  for (i in lauf[-1]){
    za <- cbind(za,rowSums(difz3[,i:(i+ncol(itmat)-1)]))
  }

  if(is.null(fixP)){
    as.vector(rto[-kateg.zahl,]) - as.vector(za[-kateg.zahl,])
  } else {
    as.vector(rto[-kateg.zahl,])[-fixP] - as.vector(za[-kateg.zahl,])[-fixP]
  }
}


persest <- optim(startv, pl, gr=der1pl, rto=rto2, patt=patt2, estpar=estpar, col.table=col.table, fixP=fixP, method="BFGS", control=list(maxit=100, fnscale=-1), hessian=TRUE)


param <- numeric(length(persest$par)+length(fixP))
ppse <- diag(solve(persest$hessian)*(-1))
ppse <- numeric(length(param))
param[fixP] <- -1000
#ppse[fixP] <- 0
ppse[param != -1000] <-  diag(solve(persest$hessian)*(-1))
param[param != -1000] <- persest$par

ptablePP <- cbind(matrix(param, ncol=nrow(rto2)-1, nrow=ncol(rto2), byrow=TRUE),rep(0,ncol(rto2)))
ptableSE <- matrix(ppse, ncol=nrow(rto2)-1, nrow=ncol(rto2), byrow=TRUE)

ptable <- cbind(ptablePP, ptableSE)

row.names(ptable) <- apply(rto2,2, function(n) paste0(n, collapse="|"))
colnames(ptable) <- c(paste(rep("pers.par.cat",kateg.zahl), 1:kateg.zahl, sep=""),paste(rep("SE cat", kateg.zahl-1),1:(kateg.zahl-1)))

pvec <- apply(row.table1,2, function(n2) paste0(n2, collapse="|"))

ptable2 <- t(apply(ptable, 1, function(p1) scale(p1[1:kateg.zahl], scale=FALSE)))
ptable3 <- cbind(ptable2, ptable[,(kateg.zahl+1):length(colnames(ptable))])
colnames(ptable3) <- colnames(ptable)

pparList <- ptable3[match(pvec,rownames(ptable3)),]

res_par <- list(ptable=ptable3, pparList=pparList,fun_calls=persest$counts, call=call)

class(res_par) <- "person_par"

res_par

}
