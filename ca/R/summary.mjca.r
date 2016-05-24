################################################################################
# summary.mjca(): Summarizing mjca objects (ca package 0.70)
################################################################################
summary.mjca <- function(object, scree = TRUE, rows = FALSE, columns = TRUE, ...){
  obj <- object
  nd  <- obj$nd
  if (is.na(nd)){
    nd <- 2
    } else {
    if (nd > length(obj$sv)){
      nd <- length(obj$sv)
      }
    }
  if (obj$lambda != "JCA"){
    K   <- obj$nd.max
    } else {
    K   <- obj$nd
    }
  I   <- length(obj$rownames)
  J   <- dim(obj$colcoord)[1]
  cpc <- obj$colpcoord
 # row profiles:
  if (rows){
    rpc <- obj$rowpcoord
    r.names <- abbreviate(obj$rownames, 3)
    sr      <- obj$rowsup
    r.mass  <- obj$rowmass
    r.inr   <- obj$rowinertia / sum(obj$rowinertia, na.rm = TRUE)
    r.ccc   <- matrix(NA, nrow = length(r.names), ncol = nd * 3)
    for (i in 1:nd)  {
      r.ccc[,3 * (i - 1) + 1] <- rpc[,i]
      r.ccc[,3 * (i - 1) + 2] <- obj$rowmass * rpc[,i]^2 / obj$rowinertia
      r.ccc[,3 * (i - 1) + 3] <- obj$rowmass * rpc[,i]^2 / obj$sv[i]
      if (obj$lambda == "indicator"){
        r.ccc[,3 * (i - 1) + 3] <- obj$rowmass * rpc[,i]^2 / sqrt(obj$sv[i])
        }
      }
    if (nd > 1) {
      r.qlt <- apply(r.ccc[,((1:nd-1) * 3 + 2)], 1, sum) 
      } else {
      r.qlt <- r.ccc[,((1:nd-1) * 3 + 2)] 
      }
    r1     <- paste(" k=", 1:nd, sep = "")
    r2     <- rep("cor", nd)
    r3     <- rep("ctr", nd)
    rcclab <- as.vector(rbind(r1, r2, r3))
    dimnames(r.ccc) <- list(r.names, rcclab)
    r.out  <- data.frame(r.names, round(1000 * r.mass, 0), round(1000 * r.qlt, 0),
                         round(1000 * r.inr, 0), round(1000 * r.ccc, 0))
    dimnames(r.out) <- list(as.character(1:length(r.names)), c("name", "mass", " qlt", " inr", rcclab))
    } else {
    r.out <- NULL
    } 
 ### COLUMNS:
  if (columns){
    c.names  <- obj$levelnames    
    sc       <- obj$colsup
  if (!is.na(sc[1])){
    c.names[sc] <- paste("(*)", c.names[sc], sep = "")
    }
  c.mass   <- obj$colmass
  c.inr    <- obj$colinertia  / sum(obj$colinertia, na.rm = TRUE)
    if (obj$lambda != "JCA"){
      c.ccc    <- matrix(NA, nrow = length(c.names), ncol = nd * 3)
      for (i in 1:nd){
        c.ccc[,3 * (i - 1) + 1] <- cpc[,i]
        c.ccc[,3 * (i - 1) + 2] <- obj$colcor[,i]
        c.ccc[,3 * (i - 1) + 3] <- obj$colctr[,i]
        }
      if (nd > 1) { 
        c.qlt <- apply(c.ccc[,((1:nd - 1) * 3 + 2)], 1, sum) 
        } else {
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
      } else { #JCA CASE BELOW:
      c.ccc  <- cpc[,1:nd]
      ccclab <- paste(" k=", 1:nd, sep = "")
      dimnames(c.ccc) <- list(c.names, ccclab)
      c.out  <- data.frame(c.names, round(1000 * c.mass, 0), round(1000 * c.inr, 0), 
                           round(1000 * c.ccc, 0),round(1000 * obj$colcor, 0), round(1000 * obj$colctr, 0) )
      dimnames(c.out) <- list(as.character(1:length(c.names)), c("name", "mass", " inr", ccclab, "cor", "ctr"))
      }
    } else {
      c.out <- NULL
    } # END COLUMNS
 ### SCREE PLOT:
  sev.0 <- round(100*obj$inertia.et, 1)
  if (scree) {
    values     <- obj$sv^2
    values2    <- round(100*values/sum(values), 1)
    scree.out  <- cbind(1:length(obj$sv), round(values, 6), values2, 
                        round(cumsum(100*values/sum(values)), 1))
    if (obj$lambda == "adjusted"){
      values     <- round(obj$sv^2, 6)
      values2    <- round(100*obj$inertia.e, 1)
      values3    <- round(cumsum(100*obj$inertia.e), 1)
      scree.out  <- cbind(1:length(obj$sv), round(values, 6), values2, values3)
      }
    if (obj$lambda == "JCA"){
      values     <- round(obj$sv^2, 6)
      values2    <- rep(NA, length(values))
      values3    <- rep(NA, length(values))
      scree.out  <- cbind(1:length(obj$sv), round(values, 6), values2, values3)
      }
    } else {
    scree.out <- NULL
    }
 ### OUTPUT:
  out <- list(scree   = scree.out,
              rows    = r.out,
              columns = c.out, 
              sev     = sev.0, 
              JCA     = obj$JCA.iter,
              tin     = obj$inertia.t, 
              JCA.nd  = obj$nd, 
              JCA.ind = sum(diag(obj$subinertia)),
              JCA.nit = obj$JCA.iter[1],
              JCA.eps = obj$JCA.iter[2])
  class(out) <- "summary.mjca"
  return(out)
  }

