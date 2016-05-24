################################################################################
# summary.ca(): Summarizing ca objects (ca package 0.70)
################################################################################
summary.ca <- function(object, scree = TRUE, rows = TRUE, columns = TRUE, ...){
  obj <- object
  nd  <- obj$nd
  if (is.na(nd)){
    nd <- 2
    } else {
    if (nd > length(obj$sv)){
      nd <- length(obj$sv)
      }
    }
 # principal coordinates:
  K   <- nd
  I   <- dim(obj$rowcoord)[1] ; J <- dim(obj$colcoord)[1]
  svF <- matrix(rep(obj$sv[1:K], I), I, K, byrow = TRUE)
  svG <- matrix(rep(obj$sv[1:K], J), J, K, byrow = TRUE)
  rpc <- obj$rowcoord[,1:K] * svF
  cpc <- obj$colcoord[,1:K] * svG
 # rows:
  strnascii <- function(x){
    foo1 <- unlist(strsplit(x, ""))
    foo2 <- grep('\\w', foo1)
    foo  <- paste(foo1[foo2], collapse = "")
    return(foo)
    }
  if(!rows) {
    r.out <- NULL
    } else {
    rnames.temp <- unlist(lapply(obj$rownames, strnascii))
    r.names     <- abbreviate(rnames.temp, 4)
    sr          <- obj$rowsup
    if (!is.na(sr[1])){
      r.names[sr] <- paste("(*)", r.names[sr], sep = "")
      }
    r.mass <- obj$rowmass
    if (length(obj$rowsup)>0){
      i0 <- obj$rowsup
      r.mass[i0] <- NA
      }
    r.inr  <- obj$rowinertia / sum(obj$rowinertia, na.rm = TRUE)
    r.ccc  <- matrix(NA, nrow = length(r.names), ncol = nd * 3)
    for (i in 1:nd){
      r.ccc[,3 * (i - 1) + 1] <- rpc[,i]
      r.ccc[,3 * (i - 1) + 2] <- rpc[,i]^2 / obj$rowdist^2
      r.ccc[,3 * (i - 1) + 3] <- obj$rowmass * rpc[,i]^2 /obj$sv[i]^2
      }
  if (nd > 1) {
    r.qlt <- apply(r.ccc[,((1:nd-1) * 3 + 2)], 1, sum) } 
    else {
    r.qlt <- r.ccc[,((1:nd-1) * 3 + 2)] 
    }
  r1              <- paste(" k=", 1:nd, sep = "")
  r2              <- rep("cor", nd)
  r3              <- rep("ctr", nd)
  rcclab          <- as.vector(rbind(r1, r2, r3))
  dimnames(r.ccc) <- list(r.names, rcclab)
  r.out           <- data.frame(r.names, round(1000 * r.mass, 0), round(1000 * r.qlt, 0),
                                round(1000 * r.inr, 0), round(1000 * r.ccc, 0))
  dimnames(r.out) <- list(as.character(1:length(r.names)), c("name", "mass", " qlt", " inr", rcclab))
  }
 # columns:
  if(!columns) {
    c.out <- NULL
  } else {
    cnames.temp <- unlist(lapply(obj$colnames, strnascii))
    c.names     <- abbreviate(cnames.temp, 4)
    sc          <- obj$colsup
    if (!is.na(sc[1])){
      c.names[sc] <- paste("(*)", c.names[sc], sep = "")
      }
    c.mass <- obj$colmass
    if (length(obj$colsup) > 0){
      i0 <- obj$colsup
      c.mass[i0] <- NA
      }
    c.inr   <- obj$colinertia / sum(obj$colinertia, na.rm = TRUE)
    c.ccc   <- matrix(NA, nrow = length(c.names), ncol = nd * 3)
    for (i in 1:nd){
      c.ccc[,3 * (i - 1) + 1] <- cpc[,i]
      c.ccc[,3 * (i - 1) + 2] <- cpc[,i]^2 / obj$coldist^2
      c.ccc[,3 * (i - 1) + 3] <- obj$colmass * cpc[,i]^2 /  obj$sv[i]^2
      }
    if (nd > 1) { 
      c.qlt <- apply(c.ccc[,((1:nd - 1) * 3 + 2)], 1, sum) } 
      else {
      c.qlt <- c.ccc[,((1:nd - 1) * 3 + 2)] 
      }
    c1     <- paste(" k=", 1:nd, sep = "")
    c2     <- rep("cor", nd)
    c3     <- rep("ctr", nd)
    ccclab <- as.vector(rbind(c1, c2, c3))
    dimnames(c.ccc) <- list(c.names, ccclab)
    c.out  <- data.frame(c.names, round(1000 * c.mass, 0), round(1000 * c.qlt, 0),
                         round(1000 * c.inr, 0), round(1000 * c.ccc, 0))
    dimnames(c.out) <- list(as.character(1:length(c.names)), c("name", "mass", " qlt", " inr", ccclab))
    }
 # scree plot:
  if (scree) {
    values     <- obj$sv^2
    values2    <- 100*(obj$sv^2)/sum(obj$sv^2)
    values3    <- cumsum(100*(obj$sv^2)/sum(obj$sv^2))
    scree.out  <- cbind(1:length(obj$sv), values, values2, values3)
    } else {
    scree.out <- NULL
    }
 # output:
  out <- list(scree   = scree.out,
              rows    = r.out,
              columns = c.out)
  class(out) <- "summary.ca"
  return(out)
  }
################################################################################
