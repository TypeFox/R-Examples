#  Create an object of class 'conect' containing the atoms IDs defining the connectivity of a molecular system.

conect <- function(...)
  UseMethod("conect")

conect.default <- function(eleid.1, eleid.2, ...)
{
  if(is.null(eleid.1) & is.null(eleid.2))
    return(NULL)
  eleid.1 <- as.integer(eleid.1)
  eleid.2 <- as.integer(eleid.2)
  con <- data.frame(eleid.1, eleid.2)
#   con[con$eleid.1 > con$eleid.2,] <- rev(con[con$eleid.1 > con$eleid.2,])
  con <- con[order(con$eleid.2),]
  con <- con[order(con$eleid.1),]
  rownames(con) <- NULL
  class(con) <- c("conect","data.frame")
  if(nrow(con) == 0)
    con <- NULL
  return(con)
}

conect.coords <- function(x, radii = 0.75, safety = 1.2, by.block = FALSE, ...) {
  if(!is.coords(x))
    stop("'x' must be an object of class 'coords'")

  radii <- radii*safety
  data <- cbind(x, radii)

  findCon <- function(data) {
    nat <- nrow(data)
    if(nat==0) return(NULL)
    r <- sqrt(
        outer(data$x1, data$x1, "-")^2 +
        outer(data$x2, data$x2, "-")^2 +
        outer(data$x3, data$x3, "-")^2
    )
    bond.dist <- outer(data$radii, data$radii, "+")
    M <- lower.tri(r) & (r < bond.dist)
    if(all(!M)) return(NULL)
    eleid <- matrix(rownames(data), nrow = nat, ncol = nat)   
    eleid.1 <- as.integer(t(eleid)[M])
    eleid.2 <- as.integer(  eleid [M]) 
    return(conect.default(eleid.1, eleid.2))
  }

  if(!by.block) {
    con <- findCon(data)
  }
  else {
    get.con <- function(shift = c(0,0,0), x, radii, step) {
      coords.range <- range(x)
      coords.range <- t(t(coords.range) - shift)

      x.cuts <- seq(coords.range["min","x"], coords.range["max","x"] + step, step)
      y.cuts <- seq(coords.range["min","y"], coords.range["max","y"] + step, step)
      z.cuts <- seq(coords.range["min","z"], coords.range["max","z"] + step, step)
    
      x.index <- cut(x$x1, x.cuts, include.lowest = TRUE)
      y.index <- cut(x$x2, y.cuts, include.lowest = TRUE)
      z.index <- cut(x$x3, z.cuts, include.lowest = TRUE)
      
      f <- interaction(x.index, y.index, z.index)
      
      data <- cbind(x, radii)
      
      con <- split(data, f)
      con <- lapply(con, findCon)
      con <- do.call(rbind, con)

      return(con)
    }

    width <-10
    shift <- expand.grid(0:1,0:1,0:1)*width/2

    con <- apply(shift, 1, get.con, x, radii, width)
    con <- unique(do.call(rbind, con))
    conect(con$eleid.1, con$eleid.2)
    rownames(con) <- NULL
  }
  
  return(con)
}

conect.pdb <- function(x, safety = 1.2, by.block = FALSE, ...) {
  symb <- toSymbols(x$atom$elename)
  symb[is.na(symb)] <- "Xx"
  rcov <- elements[match(symb, elements[,"symb"]), "rcov"]
  con <- conect(coords(x), rcov, safety, by.block)
  con <- conect(x$atoms$eleid[con$eleid.1], x$atoms$eleid[con$eleid.2])
  return(con)
}

is.conect <- function(x)
{
  to.return <- any(attr(x,which="class") == "conect")
  return(to.return)
}
