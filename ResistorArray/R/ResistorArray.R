"array.resistance" <-
function (x.offset, y.offset, rows.of.resistors, cols.of.resistors, 
    give.pots = FALSE) 
{
    total.n.resistors <- rows.of.resistors * cols.of.resistors
    A <- makefullmatrix(rows.of.resistors, cols.of.resistors)
    current.in.x <- cols.of.resistors%/%2
    current.in.y <- rows.of.resistors%/%2
    current.out.x <- current.in.x + x.offset
    current.out.y <- current.in.y + y.offset
    current.in.1d <- .get.1d.index(current.in.y, current.in.x, 
        rows.of.resistors, cols.of.resistors)
    current.out.1d <- .get.1d.index(current.out.y, current.out.x, 
        rows.of.resistors, cols.of.resistors)
    if (give.pots) {
        return(matrix(resistance(A, current.in.1d, current.out.1d, 
            give.pots = TRUE), rows.of.resistors, cols.of.resistors))
    }
    else {
        return(resistance(A, current.in.1d, current.out.1d, give.pots = FALSE))
    }
}
"circuit" <-
function (L, v, currents = 0, use.inverse = FALSE, give.internal = FALSE) 
{
    free.v <- is.na(v)
    fixed.v <- !free.v
    A <- L[fixed.v, fixed.v, drop = FALSE]
    B <- L[fixed.v, free.v, drop = FALSE]
    D <- L[free.v, free.v, drop = FALSE]
    v.known <- v[fixed.v, drop = FALSE]
    I <- v
    I[] <- currents
    if (any(free.v)) {
      if(use.inverse){
          v[free.v] <- crossprod(solve(D), I[free.v] - crossprod(B, v.known))
      } else {
          v[free.v] <- solve(D, I[free.v] - crossprod(B, v.known))
      }
    }
    I[fixed.v] <- crossprod(A, v.known) + B %*% v[free.v]
    out <- list(potentials = v, currents = I)
    if (give.internal) {
        jj <- L
        jj[] <- v
        pot.diffs <- jj - t(jj)
        return(list(potentials = v, currents = I, internal.currents = L * 
            pot.diffs, power = -L * pot.diffs^2, total.power = -sum(L * 
            pot.diffs^2)/2))
    }
    else {
        return(out)
    }
}
"cube" <-
function (x = 1) 
{
    out <- matrix(0, 8, 8)
    out[platonic("cube")] <- -1/x
    out <- out + t(out)
    diag(out) <- -apply(out, 2, sum)
    out
}
"currents" <-
function (L, earth.node, input.node) 
{
    edges <- which(L != 0 & row(L) < col(L), arr.ind = TRUE)
    rownames(edges) <- NULL
    colnames(edges) <- c("from", "to")
    volts <- resistance(L, earth.node = earth.node, input.node = input.node, 
        give.pots = TRUE)
    cbind(edges, (volts[edges[, 1]] - volts[edges[, 2]]) * L[edges])
}
"currents.matrix" <-
function (L, earth.node, input.node) 
{
    out <- L
    out[] <- resistance(L, earth.node, input.node, give.pots = TRUE)
    out <- L * (out - t(out))
    out
}
"dodecahedron" <-
function (x = 1) 
{
    out <- matrix(0, 20, 20)
    out[platonic("dodecahedron")] <- -1/x
    out <- out + t(out)
    diag(out) <- -apply(out, 2, sum)
    return(out)
}

".get.1d.index" <- function(rownum, colnum, R, C)
{
    (colnum - 1) * R + rownum
}
"icosahedron" <-
function (x = 1) 
{
    out <- matrix(0, 12, 12)
    out[platonic("icosahedron")] <- -1/x
    out <- out + t(out)
    diag(out) <- -apply(out, 2, sum)
    return(out)
}
"ladder" <-
function (n, x = 1, y = 1, z = NULL) 
{
    out <- matrix(0, n, n)
    jj.series <- rbind(cbind(1, 2:n), cbind(2:(n - 1), 3:n))
    if(is.null(z)){
      out[jj.series[1:(n - 1), ]] <- -1/x
      out[jj.series[n:(2 * n - 3), ]] <- -1/y
    } else {
      out[jj.series] <- -1/z
    }
    out <- out + t(out)
    diag(out) <- -apply(out, 2, sum)
    out
}

"makefullmatrix" <-
function (R, C)
{
    RC <- R*C
    out <- diag(4, nrow = RC)
    diagdist <- row(out) - col(out)
    
    out[diagdist == +1] <- -1
    out[diagdist == -1] <- -1
    out[diagdist == +R] <- -1
    out[diagdist == -R] <- -1
    out[diagdist == +RC - R] <- -1  
    out[diagdist == -RC + R] <- -1
    out[1 , RC] <- -1
    out[RC,  1] <- -1
    
    return(out)
  
}

"makefullmatrix_strict" <- function(R,C, toroidal)
{
    RC <- R*C
    jj <- as.vector(matrix(seq_len(RC),R,C)[-R,])
    iv <- cbind(jj,jj+1)         # interior vertical lines
      
    jj <- as.vector(matrix(seq_len(RC),R,C)[,-C])
    ih <- cbind(jj,jj+R)         # interior horizontal lines
    index <- rbind(iv,ih)
    if(toroidal){
      jj <- seq(from=1 , by=R , len=C)
      wv <- cbind(jj,jj+R-1)         # wrapped vertical lines

      jj <- seq_len(R)
      wh <- cbind(jj,jj+RC-R)         # wrapped horizontal lines

      index <- rbind(index , wv , wh)
    }
    out <- matrix(0,RC,RC)
    out[rbind(index,index[,2:1])] <- -1
    diag(out) <- -rowSums(out)
    return(out)
}
    

      
"octahedron" <-
function (x = 1) 
{
    out <- matrix(0, 6, 6)
    out[platonic("octahedron")] <- -1/x
    out <- out + t(out)
    diag(out) <- -apply(out, 2, sum)
    return(out)
}
"platonic" <-
function (a) 
{
    switch(a, tetrahedron = matrix(c(1, 1, 1, 2, 2, 3, 2, 3, 
        4, 3, 4, 4), ncol = 2), cube = matrix(c(1, 1, 1, 2, 2, 
        3, 3, 4, 5, 5, 6, 7, 2, 4, 5, 3, 6, 4, 7, 8, 6, 8, 7, 
        8), ncol = 2), octahedron = matrix(c(1, 1, 1, 1, 2, 2, 
        2, 3, 3, 4, 4, 5, 2, 3, 4, 5, 3, 5, 6, 4, 6, 5, 6, 6), 
        ncol = 2), dodecahedron = matrix(c(1, 1, 1, 2, 2, 3, 
        3, 4, 4, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 12, 13, 
        14, 15, 16, 16, 17, 18, 19, 2, 5, 6, 3, 7, 4, 8, 5, 9, 
        10, 11, 12, 12, 13, 13, 14, 14, 15, 11, 15, 16, 17, 18, 
        19, 20, 17, 20, 18, 19, 20), ncol = 2), icosahedron = matrix(c(1, 
        1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 
        6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 2, 3, 4, 5, 8, 3, 
        5, 6, 9, 4, 6, 7, 7, 8, 10, 8, 9, 11, 7, 9, 12, 10, 12, 
        10, 11, 11, 12, 11, 12, 12), ncol = 2))
}
"resistance" <-
function (A, earth.node, input.node, current.input.vector = NULL,
           give.pots = FALSE) 
{
    potentials <- rep(0, nrow(A))
    A <- A[-earth.node, -earth.node]
    if (is.null(current.input.vector)) {
        current.input.vector <- potentials
        current.input.vector[input.node] <- 1
    }
    else {
        give.pots <- TRUE
    }
    potentials[-earth.node] <- solve(A, current.input.vector[-earth.node])
    if (give.pots) {
        return(potentials)
    }
    else {
        return(potentials[input.node])
    }
}
"series" <-
function (x) 
{
    n <- length(x)
    out <- matrix(0, n + 1, n + 1)
    jj.series <- cbind(1:n, 2:(n + 1))
    out[jj.series] <- -1/x
    out <- out + t(out)
    diag(out) <- -apply(out, 2, sum)
    out
}
"tetrahedron" <-
function (x = 1) 
{
    out <- matrix(0, 4, 4)
    out[platonic("tetrahedron")] <- -1/x
    out <- out + t(out)
    diag(out) <- -apply(out, 2, sum)
    out
}
"Wu" <-
function (L) 
{
    jj <- eigen(L)
    n <- nrow(L)
    evals <- jj$values
    evecs <- t(jj$vectors)
    evals[n] <- Inf
    evecs <- evecs/sqrt(evals)
    out <- apply(evecs, 2, function(y) {
        apply(evecs, 2, function(x) {
            sum((x - y)^2)
        })
    })
    return(out)
}

"hypercube" <- function(n){

  jj <- as.matrix(expand.grid(rep(list(0:1),n))) 
  
  o <-
    -apply(jj, 1, function(y) {
      apply(jj, 1, function(x) {
        sum(x!=y)==1
      })
    })

  jj.names <- apply(jj,1,paste,collapse="")
  rownames(o) <- jj.names
  colnames(o) <- jj.names
  diag(o) <- -apply(o,1,sum,na.rm=TRUE)
  return(o)
}
