#'@rdname perspar
#'@method person_par CRSM
#'@export

person_par.CRSM <-
function(object, ...){

call <- match.call()

scores <- unique(rowSums(object$data_p))
scores <- scores[scores != 0 & scores != ncol(object$data)] # extreme scores out

S0n <- function(t, paraI, itp)       {exp(t*(paraI - itp) + t*(1-t)*object$disppar)}
S1n <- function(t, paraI, itp) {t*   exp(t*(paraI - itp) + t*(1-t)*object$disppar)}
S2n <- function(t, paraI, itp) {t^2* exp(t*(paraI - itp) + t*(1-t)*object$disppar)}

#starting values
para1      <- rep(0, length(scores))
iter       <- 0

while( !exists("para") || max(abs(para1-para)) > 0.0001 ){
  para <- para1
  iter <- iter+1

  Sf <- sapply(para, function(pa){
    s0 <- sapply(object$itempar, function(it) integrate(S0n, paraI=pa, itp=it, lower=0, upper=1, stop.on.error=F)$value)
    s1 <- sapply(object$itempar, function(it) integrate(S1n, paraI=pa, itp=it, lower=0, upper=1, stop.on.error=F)$value)
    s2 <- sapply(object$itempar, function(it) integrate(S2n, paraI=pa, itp=it, lower=0, upper=1, stop.on.error=F)$value)
    su1 <- sum(s1/s0)
    su2 <- sum(s2/s0 - (s1/s0)^2)*(-1)
    list(su1=su1,su2=su2)
  })

  para1    <- as.vector(para - (scores - unlist(Sf[1,]))*(unlist(Sf[2,])^(-1)))
}

ppse <- sqrt(unlist(Sf[2,])^(-1)*(-1))

ptable <- cbind(scores,para1,ppse)
colnames(ptable) <- c("raw score", "person par", "SE")
ptableO <- ptable[order(ptable[,1]),]

pparList <- ptable[match(rowSums(object$data_p), scores),]

res_par <- list(ptable=ptableO, pparList=pparList, fun_calls=iter, call=call)

class(res_par) <- "person_par"

res_par
}
