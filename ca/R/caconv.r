################################################################################
# caconv(): Converting between CA/MCA data types (ca 0.64)
# Arguments
#   - x   : A matrix or data frame.
#   - from: The type of input data in 'x': a frequency table ("freq"), or a 
#           response pattern matrix ("rpm"), or an indicator matrix ("ind"), or 
#           a Burt matrix ("Burt").
#   - to  : The data type into which 'x' should be converted. (Note: Conversion 
#           from "Burt" to "ind" or "rpm" is not supported).
#   - nlev: A vector containing the number of levels for each categorical 
#           variable (for 'from=ind' or 'from="Burt"'). If NA, 'nlev' is 
#           computed from the data.
#   - vars: A vector of length 2 specifying the index of the variables to use 
#           for converting to "freq" (i.e. to a regular two-way frequency table). 
#   -  ...: Further arguments (ignored).
# Value
#   - A matrix or data frame containing the converted data (with the type 
#     specified in 'to'). 
################################################################################
caconv <- function(x, 
            from = c("freq", "rpm", "ind", "Burt"), 
            to   = c("rpm", "ind", "Burt", "freq"), 
            nlev = NA, 
            vars = c(1,2),
            ...){
##### Input check:
  from <- match.arg(from)
  to   <- match.arg(to)
  out  <- NA
 # Reset factor levels for 'from="rpm"':
  if (from == "rpm"){
    x <- data.frame(lapply(x, factor))
    }
 # Fix row-/columnames:
  if (is.null(rownames(x))){
    rownames(x) <- 1:nrow(x)
    }
  if (is.null(colnames(x))){
    colnames(x) <- 1:ncol(x)
    }
 # Get nlev for 'nlev=NA & from="ind"|"Burt"':
  if (is.na(nlev)[1] & (from == "ind" | from == "Burt")){
    if (from == "ind"){
      B0 <- x %*% t(x)
      } else {
      B0 <- x
      }
    diag(B0) <- 0
    n        <- length(diag(B0))
    ind.lo   <- 1
    nlev     <- numeric(0)
    for (i in 2:n){
      B1 <- B0[ind.lo:i, ind.lo:i]
      if ((sum(B1 == 0) != (i-ind.lo+1)^2) | (i == n)){
        nlev   <- c(nlev, i - ind.lo + ifelse(i == n, 1, 0))
        ind.lo <- i
        }
      }
    }
##### Conversion:
### Source: Frequency matrix
  if (from == "freq"){
   # "freq" => "freq":
    if (to == "freq"){
      out <- x
      } else {
      var1 <- factor(rep(rownames(x), apply(x, 1, sum)), levels = rownames(x))
      var2 <- factor(rep(rep(colnames(x), nrow(x)), as.vector(t(as.matrix(x)))), levels = colnames(x))
      foo  <- data.frame(V1 = var1, V2 = var2)
     # "freq" => "rpm":
      if (to == "rpm"){
        out  <- foo
        } else {
        n.0      <- unlist(lapply(foo, nlevels))
        I.0      <- nrow(foo)
        out0     <- matrix(0, nrow = I.0, ncol = sum(n.0))
        foo1     <- lapply(foo, as.numeric)
        offset.b <- c(0, cumsum(n.0))
        offset   <- c(0, n.0[-length(n.0)])
        for (i in 1:ncol(foo)){
          out0[1:I.0 + (I.0 * (offset[i] + foo1[[i]] - 1))] <- 1
          }
        rownames(out0) <- rownames(foo)
        colnames(out0) <- paste(rep(colnames(foo), n.0), unlist(lapply(foo, levels)), sep = ".")
      # "freq" => "ind":
        if (to == "ind"){
          out <- out0
          } else {
        # "freq" => "Burt": 
          if (to == "Burt"){
            out <- t(out0)%*%out0
            }
          }
        }
      }
    } else {
### Source: Response pattern matrix
  if (from == "rpm"){
   # "rpm" => "rpm":
    if (to == "rpm"){ 
      out <- x
      } else {
      foo      <- x
      n.0      <- unlist(lapply(foo, nlevels))
      I.0      <- nrow(foo)
      out0     <- matrix(0, nrow = I.0, ncol = sum(n.0))
      foo1     <- lapply(foo, as.numeric)
      offset.b <- c(0, cumsum(n.0))
      offset   <- c(0, n.0[-length(n.0)])
      for (i in 1:ncol(foo)){
        out0[1:I.0 + (I.0 * (offset[i] + foo1[[i]] - 1))] <- 1
        }
      rownames(out0) <- rownames(foo)
      colnames(out0) <- paste(rep(colnames(foo), n.0), unlist(lapply(foo, levels)), sep = ".")
     # "rpm" => "ind":
      if (to == "ind"){ 
        out <- out0
        } else {
        out1 <- t(out0)%*%out0
       # "rpm" => "Burt":
        if (to == "Burt"){
          out <- out1
          } else {
         # "rpm" => "freq":
          if (to == "freq"){
            lo  <- cumsum(c(1,n.0))
            hi  <- cumsum(c(n.0,0))
            lut <- rbind(lo,hi)
            out <- out1[lut[,vars[1]][1]:lut[,vars[1]][2],lut[,vars[2]][1]:lut[,vars[2]][2]]
            }
          }
        }
      }
    } else {
### Source: Indicator matrix
  if (from == "ind"){
   # "ind" => "ind":
    if (to == "ind"){ 
      out <- x
      } else {
     # "ind" => "rpm":
      if (to == "rpm"){ 
        nn  <- length(nlev)
        mat <- matrix(NA, nrow = sum(nlev), ncol = nn)
        for (i in 1:nn){
          lo      <- rep(0, sum(nlev[(1:nn)[(1:nn) < i]]))
          mid     <- 1:nlev[i]
          up      <- rep(0, sum(nlev[(1:nn)[(1:nn) > i]]))
          mat[,i] <- c(lo,mid,up)
          }
        out0 <- data.frame(x %*% mat)
        out0 <- data.frame(lapply(out0, as.factor))
        foo0 <- matrix(unlist(strsplit(colnames(x), ".", fixed = TRUE)), ncol = 2, byrow = TRUE)
        lut0 <- cumsum(c(1,nlev))
        for (i in 1:ncol(out0)){
          colnames(out0)[i] <- foo0[lut0[i],1]
          levels(out0[,i])  <- foo0[lut0[i]:(lut0[i+1]-1),2]
          }
        out <- out0
        } else {
       # "ind" => "Burt":
        out0 <- x %*% t(x)
        if (to == "Burt"){ 
          out <- out0
          } else {
         # "ind" => "freq":
          if (to == "freq"){
            lo  <- cumsum(c(1,nlev))
            hi  <- cumsum(c(nlev,0))
            lut <- rbind(lo,hi)
            out <- out0[lut[,vars[1]][1]:lut[,vars[1]][2],lut[,vars[2]][1]:lut[,vars[2]][2]]
            }
          }
        }
      }
    } else {
### Source: Burt matrix
  if (from == "Burt"){
   # "Burt" => "Burt":
    if (to == "Burt"){ 
      out <- x
      } else {
     # "Burt" => "freq":
      if (to == "freq"){
        lo  <- cumsum(c(1,nlev))
        hi  <- cumsum(c(nlev,0))
        lut <- rbind(lo,hi)
        out <- x[lut[,vars[1]][1]:lut[,vars[1]][2],lut[,vars[2]][1]:lut[,vars[2]][2]]
        } else {
       # "Burt" => "ind":
        if (to == "ind"){
          stop("Option 'Burt'=>'ind' not implemented.")
          } else {
         # "Burt" => "rpm":
          if (to == "rpm"){ 
            stop("Option 'Burt'=>'rpm' not implemented.")
            } 
          }
        }
    }}}}}
  return(out)
  }
################################################################################

