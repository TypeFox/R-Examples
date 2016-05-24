"adiag" <-
function (..., pad = as.integer(0), do.dimnames = TRUE) 
{
    args <- list(...)
    if (length(args) == 1) {
        return(args[[1]])
    }
    if (length(args) > 2) {
        jj <- do.call("Recall", c(args[-1], list(pad = pad)))
        return(do.call("Recall", c(list(args[[1]]), list(jj), 
            list(pad = pad))))
    }
    a <- args[[1]]
    b <- args[[2]]
    if (is.null(b)) {
        return(a)
    }
    if (is.null(dim(a)) & is.null(dim(b))) {
        dim(a) <- rep(1, 2)
        dim(b) <- rep(1, 2)
    }
    if (is.null(dim(a)) & length(a) == 1) {
        dim(a) <- rep(1, length(dim(b)))
    }
    if (is.null(dim(b)) & length(b) == 1) {
        dim(b) <- rep(1, length(dim(a)))
    }
    if (length(dim.a <- dim(a)) != length(dim.b <- dim(b))) {
        stop("a and b must have identical number of dimensions")
    }
    s <- array(pad, dim.a + dim.b)
    s <- do.call("[<-", c(list(s), lapply(dim.a, seq_len), list(a)))
    ind <- lapply(seq(dim.b), function(i) seq_len(dim.b[[i]]) + 
        dim.a[[i]])
    out <- do.call("[<-", c(list(s), ind, list(b)))
    n.a <- dimnames(a)
    n.b <- dimnames(b)
    if (do.dimnames & !is.null(n.a) & !is.null(n.b)) {
        dimnames(out) <- mapply(c, n.a, n.b, SIMPLIFY = FALSE)
        names(dimnames(out)) <- names(n.a)
    }
    return(out)
}

"allsubhypercubes" <-
function (a) 
{
    if (!minmax(dim(a))) {
        stop("only cubes of equal dimensions allowed")
    }
    n <- dim(a)[1]
    d <- length(dim(a))
    tri <- c("", "i", "n-i+1")
    q <- expand.grid(sapply(1:d, function(x) {
        tri
    }, simplify = FALSE))
    jj <- apply(apply(q, 2, paste), 1, paste, collapse = ",")
    wanted <- grep("i.*i", jj)
    jj <- jj[wanted]
    number.of.subhypercubes <- length(jj)
    f <- function(i, a, string) {
        n <- dim(a)[1]
        execute.string <- paste("jj <- a[", string, "]", collapse = "")
        eval(parse(text = execute.string))
        d <- round(log(length(jj))/log(n))
        if(d > 0.5){
          return(array(jj, rep(n, d)))
        } else {
          return(jj)
        }
    }
    dummy <- function(p) {
        x <- sapply(1:n, f, a = a, string = jj[p], simplify = FALSE)
        along.dim <- 1 + sum(dim(x[[1]]) > 1)
        return(do.call("abind", c(x, along = along.dim)))
    }
    out <- lapply(1:number.of.subhypercubes, dummy)
    names(out) <- jj
    return(out)
}

"allsums" <-
function (m, func = NULL, ...) 
{
    n <- nrow(m)
    if(is.null(func)){
      rowsums <- rowSums(m)
      colsums <- colSums(m)
      func <- sum
    } else {
      rowsums <- apply(m, 1, FUN=func, ...)
      colsums <- apply(m, 2, FUN=func, ...)
    }
    f1 <- function(i) {
        func(diag.off(m, i, nw.se = TRUE), ...)
    }
    f2 <- function(i) {
        func(diag.off(m, i, nw.se = FALSE), ...)
    }
    majors <- sapply(0:(n - 1), FUN=f1)
    minors <- sapply(0:(n - 1), FUN=f2)
    return(list(rowsums = rowsums, colsums = colsums, majors = majors, 
        minors = minors))
}

"aplus" <-
function(...){
  args <- list(...)
  if (length(args) == 1) {
    return(args[[1]])
  }
  if (length(args) > 2) {
    jj <- do.call("Recall", c(args[-1]))
    return(do.call("Recall", c(list(args[[1]]), list(jj))))
  }
  
  a <- args[[1]]
  b <- args[[2]]
  
  dima <- dim(a)
  dimb <- dim(b)
  
  stopifnot(length(dima)==length(dimb))
  
  out <- array(0,pmax(dima,dimb))
  return(
         do.call("[<-",c(list(out),lapply(dima,seq_len),list(a)))+
         do.call("[<-",c(list(out),lapply(dimb,seq_len),list(b)))
         )
}
  
"arev" <-
function(a, swap=TRUE)
{
    if(is.vector(a)){return(rev(a))}
    d <- dim(a)
    n <- length(d)
    N <- seq_len(n)
    if(is.logical(swap)){
      if(length(swap)==1){swap <- rep(swap,n)}
    } else {
      swap <- N %in% swap
    }
    f <- function(i){
      if(d[i]>0){
        return(swap[i]*rev(seq_len(d[i])) + (!swap[i])*seq_len(d[i]))
      } else {
        return(0)
      }
    }
    do.call("[", c(list(a), sapply(N, f, simplify=FALSE),drop=FALSE))
}

"arot" <- 
function (a, rights=1, pair=1:2) 
{
    d <- dim(a)
    n <- length(d)
    jj <- 1:n
    jj[pair] <- shift(jj[pair],1)
    rights <- rights%%4
    if(rights==0){
      return(a)
    } else if (rights==1){
      return(aperm(arev(a,pair[2]),jj))
    } else if (rights==2){
      return(arev(a,pair))
    } else if (rights==3){
      return(aperm(arev(a,pair[1]),jj))
    } else {
      stop("rights must be one of 0,1,2,3")
    }
}

"ashift" <-
function (a, v=rep(1,length(dim(a))))
{
    if (is.vector(a)) {
        return(shift(a, v))
    }
    v <- c(v,rep(0,length(dim(a))-length(v)))
    f <- function(i) {
        shift(seq_len(dim(a)[i]), v[i])
    }
    do.call("[", c(list(a), sapply(seq_along(dim(a)), f, simplify = FALSE)))
}

"as.standard" <-
function (a, toroidal=FALSE, one_minus=FALSE) 
{
    if(one_minus){
      a1 <- as.standard(a,          toroidal=toroidal,one_minus=FALSE)
      a2 <- as.standard(1L+max(a)-a,toroidal=toroidal,one_minus=FALSE)
      if(a1 %lt% a2){
        return(a1)
      } else {
        return(a2)
      }
    }  
    a <- drop(a)
    d.a <- dim(a)
    
    if(any(d.a) < 1){ return(a) }
    
    if(!toroidal &  (max(d.a) <= 1 )){ return(a) }

    d <- length(d.a)

    if(toroidal){
      jj <- which(a==min(a),arr.ind=TRUE)
      if(nrow(jj)==1){
        a <- ashift(a,1-jj)    # move the "1" to top-left.
      } else {
        stop("minimum not unique")
      }

      # now pivot so a[2,1,1] < a[d[1],1,1] etc:
      f <- function(a){cbind(c(1,a-1))}
      ind <- matrix(a[1+do.call("adiag",sapply(d.a, f, simplify=FALSE))],nrow=2)
      jj <- ind[1,] > ind[2,]
      a <- ashift(arev(a,jj),jj+0)
    } else {  # not toroidal
      corners <- as.matrix(do.call("expand.grid", lapply(1:d, function(i) c(1, d.a[i]))))
      pos.of.min <- corners[which.min(a[corners]), ]
      d.a[pos.of.min > 1] <- -1    
      a <- arev(a, d.a<0)
    }

    # now aperm so adjacent elements are in the correct order:
    return(aperm(a, order(-a[1 + diag(nrow = d)])))
}

"circulant" <-
function (vec)
{
   if(length(vec)==1){vec <- seq(length=vec)}
   n <- length(vec)
   a <- matrix(0,n,n)
   out <- process(1-row(a)+col(a),n)
   out[] <- vec[out]
   return(out)
}

latin <- circulant

"diag.off" <-
function (a, offset = 0, nw.se = TRUE) 
{
    n <- dim(a)[1]
    if (nw.se == TRUE) {
        indices <- cbind(1:n, 1:n)
    }
    else {
        indices <- cbind(1:n, n:1)
    }
    jj <- process(sweep(indices, 2, c(0, offset), "+"), n)
    return(a[jj])
}
"arow" <-
function (a, i) 
{
    p <- 1:prod(dim(a))
    n <- length(dim(a))
    d <- dim(a)[i]
    permute <- c(i, (1:n)[-i])
    a <- aperm(a, permute)
    a[] <- p
    permute[permute] <- 1:n
    return(force.integer(aperm(process(a, d), permute)))
}

"force.integer" <-
function (x) 
{
    out <- as.integer(x)
    attributes(out) <- attributes(x)
    return(out)
}

"hudson" <-
function (n = NULL, a = NULL, b = NULL) 
{
    if (is.null(n)) {
        n <- length(a)
    }
    if (is.null(a)) {
        a <- c(n - 1, 0:(n - 2))
    }
    if (is.null(b)) {
        b <- c(2:(n - 1), n, 1)
    }
    perm <- c(n - 1, n, 1:(n - 2))
    f <- function(i) {
        recurse(perm=perm, i=i, start = a)
    }
    g <- function(i) {
        recurse(perm = perm, i=i, start = b)
    }
    jj <- 0:(n - 1)
    aa <- t(sapply(jj, f))
    bb <- t(sapply(-jj, g))
    return(n * aa + bb)
}

"is.2x2.correct" <-
function (m, give.answers = FALSE) 
{
    window <- c(2, 2)
    sums <- subsums(m, window)
    answer <- minmax(sums)
    if (give.answers == FALSE) {
        return(answer)
    }
    else {
        return(list(answer = answer, tbt.sums = sums))
    }
}

"is.associative" <-
function (m) 
{
    if(is.list(m)){
      return(sapply(m,match.fun(sys.call()[[1]])))
    }
    is.magic(m) & minmax(c(m + arev(m)))
}

"is.square.palindromic" <-
function (m, base=10, give.answers=FALSE)
{
    n <- nrow(m)
    S <- function(i){ashift(diag(n),c(i,0))}
    f.maj <- function(i){
      is.persymmetric(m %*% S(i) %*% t(m))
    }

    f.min <- function(i){
      is.persymmetric(t(m) %*% S(i) %*% m)
    }

    row.sufficient <- is.persymmetric(t(m) %*% m)
    col.sufficient <- is.persymmetric(m %*% t(m))
    major.diag.sufficient <- all(sapply(1:nrow(m),f.maj))
    minor.diag.sufficient <- all(sapply(1:nrow(m),f.min))
    
    sufficient <- row.sufficient & col.sufficient & major.diag.sufficient & minor.diag.sufficient
      
    b <- base^(0:(n-1))
    R <- diag(n)[n:1,]

    is.necessary <- function(mat,tol=1e-8){
      as.vector(abs(
      (crossprod(b, R %*% mat %*% R) %*% b)/
        (crossprod(b, mat) %*% b)-1)<tol)
    }
    
    row.necessary <- is.necessary(crossprod(m,m))
    col.necessary <- is.necessary(m %*% t(m))

    f1 <- function(i) {
      diag.off(m, i, nw.se = TRUE)
    }
    f2 <- function(i) {
      diag.off(m, i, nw.se = FALSE)
    }
    m.tilde.major <- sapply(0:(n-1),f1)
    m.tilde.minor <- sapply(0:(n-1),f2)

    major.diag.necessary <-  is.necessary(crossprod(m.tilde.major %*% t(m.tilde.major)))
    minor.diag.necessary <-  is.necessary(crossprod(m.tilde.minor %*% t(m.tilde.minor)))

    necessary=row.necessary & col.necessary & major.diag.necessary & minor.diag.necessary

    if(give.answers){
      return(list(
                  necessary = necessary,
                  sufficient = sufficient,
                  row.necessary = row.necessary,
                  col.necessary = col.necessary,
                  major.diag.necessary = major.diag.necessary,
                  minor.diag.necessary = minor.diag.necessary,
                  row.sufficient = row.sufficient,
                  col.sufficient = col.sufficient,
                  major.diag.sufficient = major.diag.sufficient,
                  minor.diag.sufficient = minor.diag.sufficient
                  ))
    } else {
      return( sufficient )
    }
}

"is.bree.correct" <-
function (m, give.answers = FALSE) 
{
    diag.dist <- nrow(m)%/%2
    offsets <- matrix(c(0, diag.dist), 2, 2)
    diag.sums <- subsums(m, offsets)
    answer <- minmax(diag.sums)
    if (give.answers == FALSE) {
        return(answer)
    }
    else {
        return(list(answer = answer, diag.sums = diag.sums))
    }
}

"is.centrosymmetric" <-
function (m)
{
    if(is.list(m)){
      return(sapply(m,match.fun(sys.call()[[1]])))
    }
    all(m==arev(m))
}

"is.circulant" <-
function(m,dir=rep(1,length(dim(m))))
{
  return(all(m == ashift(m,dir)))
}

"is.diagonally.correct" <-
function (a, give.answers = FALSE, func = sum, boolean = FALSE, ...) 
{
    if (!minmax(dim(a))) {
        stop("only cubes of equal dimensions allowed")
    }
    n <- dim(a)[1]
    d <- length(dim(a))
    b <- c(1, -1)
    f <- function(dir) {
        (dir > 0) * (1:n) + (dir < 0) * (n:1)
    }
    g <- function(jj) {
        func(a[sapply(jj, f)], ...)
    }
    ans <- expand.grid(rep(list(b),d))
    diag.sums <- apply(ans, 1, g)
    dim(diag.sums) <- c(length(diag.sums)/(2^d), rep(2, d))
    if (boolean) {
        answer <- all(diag.sums)
    }
    else {
        answer <- minmax(diag.sums)
    }
    if (give.answers) {
        return(list(answer = answer, diag.sums = drop(diag.sums)))
    }
    else {
        return(answer)
    }
}

"is.latin" <-
function (m, give.answers = FALSE) 
{
    is.latinhypercube(a = m, give.answers = give.answers)
}

"is.latinhypercube" <-
function (a, give.answers = FALSE) 
{
    f <- function(x) {
        minmax(c(1, diff(sort(x))))
    }
    is.consecutive <- is.semimagichypercube(a, func = f, give.answers = TRUE)$rook.sums
    answer <- all(is.consecutive)
    if (give.answers) {
        return(list(answer = answer, is.consecutive))
    }
    else {
        return(answer)
    }
}

"is.magic" <-
function (m, give.answers = FALSE, func = sum,  boolean = FALSE)
{
  if(is.list(m)){
    out <- lapply(m,match.fun(sys.call()[[1]]),
                  give.answers = give.answers,
                  func         = func,
                  boolean      = boolean
                  )
    if(give.answers){
      return(out)
    } else {
      return(unlist(out))
    }
  }
  sums <- allsums(m, func = func)
  jj <- c(sums$rowsums, sums$colsums, sums$majors[1], sums$minors[1])
  if (boolean) {
    answer <- all(jj)
  }
  else {
    answer <- minmax(jj)
  }
  if (give.answers) {
    return(c(answer = answer, sums))
  }
  else {
    return(answer)
  }
}

"is.magichypercube" <-
function (a, give.answers = FALSE, func = sum, boolean = FALSE, ...) 
{
    diag.sums <- is.diagonally.correct(a, give.answers = TRUE, 
        func = func, boolean = boolean, ...)$diag.sums
    jj.semi <- is.semimagichypercube(a, give.answers = TRUE, 
        func = func, boolean = boolean, ...)
    answer <- minmax(diag.sums) & jj.semi$answer
    if (give.answers) {
        return(list(answer = answer, rook.sums = jj.semi$rook.sums, 
            diag.sums = diag.sums))
    }
    else {
        return(answer)
    }
}

"is.mostperfect" <-
function (m, give.answers = FALSE) 
{
    if (give.answers) {
        ibc <- is.bree.correct(m, give.answers = TRUE)
        i2c <- is.2x2.correct(m, give.answers = TRUE)
        ipd <- is.panmagic(m, give.answers = TRUE)
        return(list(answer = ibc$answer & i2c$answer, rowsums = ipd$rowsums, 
            colsums = ipd$colsums, majors = ipd$majors, minors = ipd$minors, 
            diag.sums = ibc$diag.sums, tbt.sums = i2c$tbt.sums))
    }
    else {
        return(is.bree.correct(m) & is.2x2.correct(m))
    }
}

"is.normal" <-
function (m) 
{
    if(is.list(m)){
      return(sapply(m,match.fun(sys.call()[[1]])))
    }
    minmax(c(1, diff(sort(m))))
}

"is.ok" <-
function (vec, n = length(vec), d = 2) 
{
    return(sum(vec) == magic.constant(n, d = d))
}

"is.panmagic" <-
function (m, give.answers = FALSE, func = sum, boolean = FALSE) 
{
    sums <- allsums(m, func = func)
    jj <- c(sums$rowsums, sums$colsums, sums$majors, sums$minors)
    if (boolean) {
        answer <- all(jj)
    }
    else {
        answer <- minmax(jj)
    }
    if (give.answers) {
        return(c(answer = answer, sums))
    }
    else {
        return(answer)
    }
}

"is.pandiagonal" <- is.panmagic

"is.perfect" <-
function (a, give.answers = FALSE, func = sum, boolean = FALSE) 
{
    d <- length(dim(a))
    putative.magic.constant <- func(do.call("[", c(list(a), alist(a = )$a, 
        rep(1, d - 1))))
    jj.is.ok <- function(jj, jj.give) {
        if (length(dim(jj)) == 1) {
            if (boolean) {
                return(func(jj))
            }
            else {
                if (jj.give) {
                  return(func(jj))
                }
                else {
                  return(func(jj) == putative.magic.constant)
                }
            }
        }
        else {
            return(is.semimagichypercube(jj, func = func, boolean = boolean, 
                give.answers = jj.give))
        }
    }
    semi.stuff <- is.semimagichypercube(a, give.answers = TRUE, 
        func = func, boolean = boolean)
    diag.stuff <- unlist(lapply(allsubhypercubes(a), jj.is.ok, 
        jj.give = FALSE))
    answer <- semi.stuff$answer & all(diag.stuff)
    if (give.answers) {
        diag.sums <- lapply(allsubhypercubes(a), jj.is.ok, jj.give = TRUE)
        return(list(answer = answer, rook.sums = semi.stuff$rook.sums, 
            diag.sums = unlist(diag.sums, recursive = FALSE)))
    }
    else {
        return(answer)
    }
}

"is.persymmetric" <-
function (m)
{
    jj <- m[,nrow(m):1]
    all(jj==t(jj))
}

"is.semimagic" <-
function (m, give.answers = FALSE, func = sum, boolean = FALSE) 
{
    sums <- allsums(m, func = func)
    jj <- c(sums$rowsums, sums$colsums)
    if (boolean) {
        answer <- all(jj)
    }
    else {
        answer <- minmax(jj)
    }
    if (give.answers) {
        return(c(answer = answer, sums))
    }
    else {
        return(answer)
    }
}

"is.semimagic.default" <-
function(m)
{
    minmax(c(rowSums(m),colSums(m)))
}
  
"is.semimagichypercube" <-
function (a, give.answers = FALSE, func = sum, boolean = FALSE, ...) 
{
    d <- length(dim(a))
    f <- function(i) {
        apply(a, (1:d)[-i], FUN=func, ...)
    }
    jj <- sapply(1:d, f)
    if (minmax(dim(a))) {
        dim(jj) <- c(dim(a)[-1], d)
        if (boolean) {
            answer <- all(jj)
        }
        else {
            answer <- minmax(jj)
        }
    }
    else {
        if (boolean) {
            answer <- all(unlist(jj))
        }
        else {
            answer <- minmax(unlist(jj))
        }
    }
    if (give.answers) {
        return(list(answer = answer, rook.sums = jj))
    }
    else {
        return(answer)
    }
}

"is.standard" <-
function (a,toroidal=FALSE,one_minus=FALSE) 
{
    if(is.list(a)){
      return(sapply(a,match.fun(sys.call()[[1]])))
    }

    if(one_minus){
      ans1 <- is.standard(a,toroidal=toroidal)
      ans2 <- a %le% as.standard(max(a)+1L-a,toroidal=toroidal)
      return(ans1 & ans2)
    }
    
    if(toroidal){return(is.standard.toroidal(a))}
    a <- drop(a)
    d.a <- dim(a)
    if(any(d.a==0)){return(TRUE)}
    d <- length(d.a)
    corners <- as.matrix(do.call("expand.grid", lapply(1:d, function(i) c(1, 
        d.a[i]))))
    corners.correct <- a[1] <= min(a[corners])

    jj <- 1 + diag(nrow = d)
    adjacent.correct <- all(diff(a[jj])<0)
    
    return(corners.correct & adjacent.correct)
}

"is.standard.toroidal" <- function(a){
    if(is.list(a)){
      return(sapply(a,match.fun(sys.call()[[1]])))
    }

    first.element.correct <- identical(which(a==min(a)) , 1L)

    jj <- 1 + diag(nrow = length(dim(a)))
    adjacent.correct <- all(diff(a[jj])<0)

    f <- function(a){cbind(c(1,a-1))}
    ind <- matrix(a[1+do.call("adiag",sapply(dim(a), f, simplify=FALSE))],nrow=2)
    octahedron.correct <- all(ind[1,] < ind[2,])

    return(first.element.correct & adjacent.correct & octahedron.correct)
}

"lozenge" <-
function (m) 
{
    if(length(m)>1){
      return(sapply(m,match.fun(sys.call()[[1]])))
    }
    n <- 2 * m + 1
    out <- matrix(NA, n, n)
    jj <- cbind(m:-m, 0:(2 * m)) + 1
    odd.a <- jj[1:(1 + m), ]
    odd.b <- odd.a
    odd.b[, 2] <- odd.b[, 2] + 1
    odd.b <- odd.b[-(m + 1), ]
    odd.coords <- rbind(odd.a, odd.b)
    even.a <- jj[(m + 2):(2 * m + 1), ]
    even.b <- jj[(m + 1):(2 * m + 1), ]
    even.b[, 2] <- even.b[, 2] + 1
    even.coords <- rbind(even.a, even.b)
    f <- function(a, x) {
        x + a
    }
    all.odd.coords <- do.call("rbind", sapply(0:m, f, x = odd.coords, 
        simplify = FALSE))
    all.even.coords <- do.call("rbind", sapply(0:m, f, x = even.coords, 
        simplify = FALSE))
    all.even.coords <- process(all.even.coords, n)
    diam.odd <- 1:(1 + 2 * m * (1 + m))
    out[all.odd.coords[diam.odd, ]] <- 2 * diam.odd - 1
    diam.even <- 1:(2 * m * (1 + m))
    out[all.even.coords[diam.even, ]] <- 2 * diam.even
    return(force.integer(out))
}

"magic" <-
function (n) 
{
    if(length(n)>1){
      return(sapply(n,match.fun(sys.call()[[1]])))
    }
    n <- round(n)
    if (n == 2) {
        stop("Normal magic squares of order 2 do not exist")
    }
    if (n%%2 == 1) {
        return(as.standard(magic.2np1(floor(n/2))))
    }
    if (n%%4 == 0) {
        return(as.standard(magic.4n(round(n/4))))
    }
    if (n%%4 == 2) {
        return(as.standard(magic.4np2(round((n - 2)/4))))
    }
    stop("This cannot happen")
}

"magic.2np1" <-
function (m, ord.vec = c(-1, 1), break.vec = c(1, 0), start.point = NULL) 
{
     if(length(m)>1){
       return(sapply(m,match.fun(sys.call()[[1]]),
                     ord.vec = ord.vec,
                     break.vec = break.vec,
                     start.point = start.point
                     ))
     }
     n <- 2 * m + 1
    if (is.null(start.point)) {
        start.row <- 0
        start.col <- n + 1
    }
    else {
        start.row <- start.point[1] - 1
        start.col <- m + start.point[2] + 1
    }
    f <- function(n) {
        ord.row <- seq(from = start.row, by = ord.vec[1], length = n)
        ord.col <- seq(from = start.col, by = ord.vec[2], length = n)
        out <- cbind(rep(ord.row, n) - (n - 1), rep(ord.col, 
            n) + m)
        break.row <- ord.vec[1] - break.vec[1]
        break.col <- ord.vec[2] - break.vec[2]
        adjust <- cbind(rep(seq(from = 0, by = break.row, len = n), 
            each = n), rep(seq(from = 0, by = break.col, len = n), 
            each = n))
        return(process(out - adjust, n))
    }
    a <- matrix(NA, n, n)
    a[f(n)] <- 1:(n * n)
    return(a)
}

"magic.4n" <-
function (m) 
{
    if(length(m)>1){
      return(sapply(m,match.fun(sys.call()[[1]])))
    }
    n <- 4 * m
    a <- matrix(1:(n^2), n, n)
    jj.1 <- kronecker(diag(2), matrix(1, 2, 2))
    jj <- as.logical(kronecker(matrix(1, m + 1, m + 1), jj.1)[2:(n + 
        1), 2:(n + 1)])
    a[jj] <- rev(a[jj])
    return(force.integer(a))
}

"magic.4np2" <-
function (m) 
{
    if(length(m)>1){
      return(sapply(m,match.fun(sys.call()[[1]])))
    }
    n <- 4 * m + 2
    f <- function(n) {
        if (n == 1) {
            return(matrix(c(4, 2, 1, 3), 2, 2))
        }
        if (n == 2) {
            return(matrix(c(1, 2, 4, 3), 2, 2))
        }
        if (n == 3) {
            return(matrix(c(1, 3, 4, 2), 2, 2))
        }
        return(NULL)
    }
    lux.n <- function(m) {
        lux <- matrix(1, 2 * m + 1, 2 * m + 1)
        lux[(m + 2), ] <- 2
        if (m > 1) {
            lux[(m + 3):(2 * m + 1), ] <- 3
        }
        lux[m + 1, m + 1] <- 2
        lux[m + 2, m + 1] <- 1
        return(lux)
    }
    i <- function(a, r) {
        jj <- which(a == r, arr.ind = TRUE)
        indices <- (cbind(jj[, 1] + (jj[, 3] - 1) * 2, jj[, 2] + 
            (jj[, 4] - 1) * 2))
        o <- order(indices[, 1] * nrow(jj) + indices[, 2])
        return(indices[o, ])
    }
    a <- apply(lux.n(m), 1:2, FUN = f)
    dim(a) <- c(2, 2, 2 * m + 1, 2 * m + 1)
    out <- matrix(NA, n, n)
    sequ <- as.vector(t(magic.2np1(m))) * 4 - 4
    out[i(a, 1)] <- sequ + 1
    out[i(a, 2)] <- sequ + 2
    out[i(a, 3)] <- sequ + 3
    out[i(a, 4)] <- sequ + 4
    return(force.integer(out))
}

"magic.8" <-
function (...) 
{
    f <- function(...) {
        0:1
    }
    j <- array(t(expand.grid(rep(list(0:1),16))),
        c(4, 4, 65536))
    all.rowsums.eq.2 <- apply(apply(j, c(1, 3), sum) == 2, 2, 
        all)
    all.colsums.eq.2 <- apply(apply(j, c(2, 3), sum) == 2, 2, 
        all)
    both.sums.eq.2 <- all.rowsums.eq.2 & all.colsums.eq.2
    j <- j[c(1:4, 4:1), c(1:4, 4:1), both.sums.eq.2] > 0
    n <- dim(j)[3]
    magics <- array(1:64, c(8, 8, n))
    ref <- function(magics, j) {
        magics[j] <- rev(magics[j])
        return(magics)
    }
    fun <- function(i) {
        ref(magics[, , i], j[, , i])
    }
    return(array(sapply(1:n, fun), c(8, 8, n)))
}

"magic.constant" <-
function (n, d = 2, start = 1) 
{
    if (is.array(n)) {
        return(Recall(n = dim(n)[1], d = length(dim(n))))
    }
    n * (start + (n^d - 1)/2)
}

"magiccube.2np1" <-
function (m) 
{

    if(length(m)>1){
      return(sapply(m,match.fun(sys.call()[[1]])))
    }

    n <- 2 * m + 1
    jj <- array(1:n, rep(n, 3))
    x <- arow(jj, 1)
    y <- arow(jj, 2)
    z <- arow(jj, 3)
    return(force.integer(((x - y + z - 1) - n * floor((x - y + z - 1)/n)) * 
                         n * n + ((x - y - z) - n * floor((x - y - z)/n)) * n + 
                         ((x + y + z - 2) - n * floor((x + y + z - 2)/n)) + 1))
}

"magichypercube.4n" <-
function (m, d = 3) 
{
    n <- 4 * m
    a <- array(0, rep(2, d))
    jj.f <- function(i) {
        arow(a, i)
    }
    x <- apply(sapply(1:d, jj.f, simplify = TRUE), 1, sum)
    dim(x) <- rep(2, d)
    a[x%%2 == 1] <- 1
    i <- kronecker(array(1, rep(m + 1, d)), kronecker(a, array(1, 
        rep(2, d)))) == 1
    i <- do.call("[", c(list(i), lapply(1:d, function(jj.i) {
        2:(n + 1)
    })))
    j <- array(1:(n^d), rep(n, d))
    j[i] <- rev(j[i])
    return(j)
}

"magicplot" <-
function (m, number = TRUE, do.circuit = FALSE, ...) 
{
    par(pch = 16)
    n <- nrow(m)
    jj <- sort(t(m[n:1, ]), index.return = TRUE)$ix
    x <- process(jj, n)
    y <- (jj - 1)%/%n
    par(pty = "s", xaxt = "n", yaxt = "n")
    plot(x, y, type = "n", asp = 1, xlab = "", ylab = "", frame = FALSE)
    if (number == TRUE) {
        text(x, y, as.character(1:(n * n)))
        if (missing(...)) {
            points(x, y, type = "l")
        }
        else {
            points(x, y, cex = 0, ...)
        }
    }
    else {
        if (missing(...)) {
            points(x, y, type = "o")
        }
        else {
            points(x, y, ...)
        }
    }
    if (do.circuit == TRUE) {
        lines(c(x[1], x[n * n]), c(y[1], y[n * n]), ...)
    }
}

"magic.prime" <-
function (n, i = 2, j = 3) 
{
    a <- matrix(0, n, n)
    return(force.integer(n * (col(a) - i * row(a) + i - 1)%%n + 
        (col(a) - j * row(a) + j - 1)%%n + 1))
}

"magic.product" <-
function (a, b, mat = NULL) 
{
    if (length(a) == 1) {
        a <- magic(a)
    }
    if (length(b) == 1) {
        b <- magic(b)
    }
    if (is.null(mat)) {
        mat <- a * 0
    }
    if (any(dim(mat) != dim(a))) {
        stop("third argument must be same size as a")
    }
    ra <- nrow(a)
    ca <- ncol(a)
    rb <- nrow(b)
    cb <- ncol(b)
    aa <- a
    aa[aa] <- seq_along(a)
    out <- sapply(mat[aa], transf, a = b)
    out <- sweep(out, 2, length(b) * (seq_along(a)-1L), FUN = "+")
    out <- out[, a]
    dim(out) <- c(rb, cb, ra, ca)
    out <- aperm(out, c(1, 3, 2, 4))
    dim(out) <- c(ra * rb, ca * cb)
    return(force.integer(out))
}

"magic.product.fast" <-
function (a, b) 
{
    if ((length(a) == 1) & (length(b) == 1)) {
        return(Recall(magic(a), magic(b)))
    }
    a.l <- nrow(a)
    b.l <- nrow(b)
    return(force.integer(b.l * b.l * (kronecker(a, matrix(1, 
        b.l, b.l)) - 1) + kronecker(matrix(1, a.l, a.l), b)))
}

"minmax" <-
function (x, tol=1e-6) 
{   
    if(is.integer(x)){
      return(identical(max(x), min(x)))
    }
    if(all(x==0)){return(TRUE)}   #special dispensation for all zeros
    if(is.double(x)){
      return(abs(max(x)-min(x))/max(abs(x)) < tol)
    } else {
      return(
             abs(max(Re(x))-min(Re(x)))/max(abs(x)) < tol &
             abs(max(Im(x))-min(Im(x)))/max(abs(x)) < tol)
    }
    
}

"notmagic.2n" <-
function (m) 
{
    options(warn = -1)
    n <- 2 * m
    a <- matrix(NA, n, n)
    s <- seq(from = 2, by = 2, to = m)
    jj.down <- kronecker(rep(1, m), rbind(1:n, n:1))[, 1:m]
    jj.down[, s] <- jj.down[n:1, s]
    jj.down <- cbind(c(1:n, n:1), as.vector(jj.down))
    jj.up <- jj.down
    jj.up[, 2] <- (m + jj.up[, 2])%%n
    jj.up[jj.up == 0] <- n
    jj.both <- rbind(jj.down, jj.up)
    a[jj.both] <- 1:(n^2)
    return(a)
}

"panmagic.4" <-
function (vals = 2^(0:3)) 
{
    a <- rep(1, 2)
    S <- kronecker(a, kronecker(diag(a), t(a)))
    A <- diag(a)[c(1, 2, 1, 2), c(2, 1, 1, 2)]
    N <- t(S)
    C <- t(A)
    jj <- array(c(S, A, N, C), c(4, 4, 4))
    return(force.integer(1 + apply(sweep(jj, 3, vals, "*"), 1:2, 
        sum)))
}

"panmagic.8" <-
function (chosen = 1:6, vals = 2^(0:5)) 
{
    a <- rep(1, 2)
    a.01 <- kronecker(matrix(1, 2, 2), kronecker(diag(a), t(a)))[c(1:4, 
        4:1), ]
    a.03 <- kronecker(t(a), diag(a)[c(1, 2, 1, 2), c(2, 1, 1, 
        2)])[c(1:4, 4:1), ]
    a.05 <- kronecker(a, kronecker(kronecker(a, kronecker(diag(a), 
        t(a))), t(a)))
    a.07 <- kronecker(kronecker(a, diag(a)[c(1, 2, 1, 2), c(2, 
        1, 1, 2)]), t(a))
    a.09 <- kronecker(a, kronecker(kronecker(diag(a), t(c(a, 
        a))), a))
    a.11 <- kronecker(diag(a)[c(1, 2, 1, 2), c(2, 1, 1, 2)], 
        matrix(1, 2, 2))
    a.02 <- t(a.01)
    a.04 <- t(a.03)
    a.06 <- t(a.05)
    a.08 <- t(a.07)
    a.10 <- t(a.09)
    a.12 <- t(a.11)
    jj <- array(c(a.01, a.02, a.03, a.04, a.05, a.06, a.07, a.08, 
        a.09, a.10, a.11, a.12), c(8, 8, 12))
    jj <- jj[, , chosen, drop = FALSE]
    return(force.integer(1 + apply(sweep(jj, 3, vals, "*"), 1:2, 
        sum)))
}

"process" <-
function (x, n) 
{
    x <- x%%n
    x[x == 0] <- as.integer(n)
    return(x)
}

"recurse" <- function (perm,i, start=seq_along(perm)) {
  i <- as.integer(i)
  if (i < 0) {
    invert <- function(perm) {
      perm[perm] <- seq_along(perm)
      return(perm)
    }
    return(Recall(start = start, invert(perm), -i))
  }
  perm.final <- seq_along(perm)
  while (i != 0) {
    perm.final <- perm[perm.final]
    i <- i - 1L
  }
  return(start[perm.final])
}

"shift" <-
function (x, i=1) 
{
    n <- length(x)
    if(n==0){return(x)}
    i <- i%%n
    if (i == 0) {
        return(x)
    }
    return(x[c((n - i + 1):n, 1:(n - i))])
}

"strachey" <-
function (m, square = magic.2np1(m)) 
{
    if(length(m)>1){
      stopifnot(length(m) == length(square))
      funcname <- match.fun(sys.call()[[1]])
      f <- function(i){
        do.call(funcname, list(m=m[i],square=square[[i]]))
      }
      return(lapply(seq_along(m),f))
    }
    
    m <- round(m)
    n <- 4 * m + 2
    r <- 2 * m + 1
    out <- kronecker(matrix(c(0, 3, 2, 1), 2, 2), matrix(1, r, 
        r)) * r^2 + kronecker(matrix(1, 2, 2), square)
    coords.top <- as.matrix(expand.grid(1:r, 1:m))
    coords.top[m + 1, 2] <- m + 1
    if (m > 1) {
        coords.top <- rbind(coords.top, as.matrix(expand.grid(1:r, 
            n:(n - m + 2))))
    }
    coords.low <- sweep(coords.top, 2, c(r, 0), "+")
    jj <- out[coords.top]
    out[coords.top] <- out[coords.low]
    out[coords.low] <- jj
    return(force.integer(out))
}

"subsums" <-
function (a, p, func = "sum", wrap = TRUE, pad = 0) 
{   
    if(length(p)==1){p <- 0*dim(a)+p}
    if (wrap == FALSE) {
        jj <- adiag(array(pad, p - 1), a,pad=pad)
        return(Recall(jj, p, func = func, pad = pad, 
            wrap = TRUE))
    }

    if (is.vector(p)) {
        sub.coords <- 1 - as.matrix(expand.grid(sapply(p, function(i) {
            1:i
        }, simplify = FALSE)))
    }
    else {
        sub.coords <- 1 - p
    }
    out <- apply(sub.coords, 1, function(v) {
        ashift(a, v)
    })
    dim(out) <- c(dim(a), nrow(sub.coords))
    if (nchar(func) == 0) {
        return(out)
    }
    else {
        return(apply(out, seq_along(dim(a)), FUN=func))
    }
    return(out)
}

"transf" <-
function (a, i) 
{
    i <- as.integer(i%%8)
    if (i%%2) {
        a <- t(a)
    }
    if ((i%/%2)%%2) {
        a <- a[nrow(a):1, ]
    }
    if ((i%/%4)%%2) {
        a <- a[, ncol(a):1]
    }
    return(a)
}

"apltake" <- function(a,b, give.indices=FALSE){

  if(is.vector(a)){
    return(as.vector(Recall(as.matrix(a),b=b,give.indices=give.indices)))
  }
  b <- c(b,dim(a)[seq(from=length(b)+1,length=length(dim(a))-length(b),by=1)])

  f <- function(x) {
    if (x[2] <= 0) {
      return(-seq_len(x[1]+x[2]))
    } else {
      return(seq_len(x[2]))
    }
  }
  jj <- apply(cbind(dim(a),b),1,f)
  if(is.matrix(jj)){jj <- as.list(as.data.frame(jj))}
  if(give.indices){
    return(jj)
  } else {
    return(do.call("[",c(list(a),jj ,drop=FALSE)))
  }
}

"apldrop" <- function(a, b, give.indices=FALSE){
    if(is.vector(a)){
    return(as.vector(Recall(as.matrix(a),b=b,give.indices=give.indices)))
  }
  b <- c(b,rep(0,length(dim(a))-length(b),by=1))
  f <- function(x){
    if(x[2] <= 0){
      return(seq(length=x[1]+x[2]))
    } else {
      return(-seq(length=x[2]))
    }
  }
  jj <- apply(cbind(dim(a),b),1,f)
  if(is.matrix(jj)){jj <- as.list(as.data.frame(jj))}
  if(give.indices){
    return(jj)
  } else {
    return(do.call("[",c(list(a),jj ,drop=FALSE)))
  }
}

"apltake<-" <- function(a,b,value){
  do.call("[<-",c(list(a),apltake(a,b,give.indices=TRUE),value))
}

"apldrop<-" <- function(a,b,value){
  do.call("[<-",c(list(a),apldrop(a,b,give.indices=TRUE),value))
}

"fnsd" <- function(a,n=1){
    return(which(dim(a)>1)[seq_len(n)])
}

"apad" <- function(a, l, e=NULL, method="ext", post=TRUE){
  if(is.vector(a)){
    return(drop(Recall(as.matrix(a), l=c(l,0), e=e, method=method,post=post)))
  }
  if(length(l)==1){
    jj <- rep(0,length(dim(a)))
    jj[l] <- e
    l <- jj
  }
  if(post){
    f <-
      switch(method,
             ext    = function(x){c(1:x[1],rep(x[1],x[2]))},
             mirror = function(x){ rep(c(1:x[1],x[1]:1),length=x[1]+x[2])},
             rep    = function(x){ rep(1:x[1],length=x[1]+x[2])}
             )
  } else {
    f <-
      switch(method,
             ext    = function(x){c(rep(1,x[2]), 1:x[1])},
             mirror = function(x){ rev(rep(c(x[1]:1,1:x[1]),length=x[1]+x[2]))},
             rep    = function(x){ rev(rep(x[1]:1,length=x[1]+x[2]))}
             )
  }
  jj <- apply(cbind(dim(a),l),1,f)
  if(is.matrix(jj)){jj <- as.list(as.data.frame(jj))}
  return(do.call("[",c(list(a), jj ,drop=FALSE)))
}

"do.index" <- function(a,f, ...){
       jj <- function(i) {seq_len(dim(a)[i])}
       index <- as.matrix(expand.grid(lapply(seq_len(length(dim(a))), jj),
                                      KEEP.OUT.ATTRS = FALSE) )
       a[index] <- apply(index, 1, f, ...)
       return(a)
}

"sam" <- function(m, u, A=NULL, B=A){
  if(is.null(A)){
    A <- latin(m)
  }
  if(is.null(B)){
    B <- is.latin(m)
  }
  
  if(u%%2){  # u odd
    if(u < 3){
      jj <- NULL
    } else {
      jj <- 8 * seq(from=0 , by=1 , to=round((u-3)/2) )
    }
    JC <- c(0,   6+jj, 13+jj)
    JD <- c(1,   7+jj, 12+jj)
    JS <- c(2,4, 8+jj, 11+jj)
    JT <- c(3,5, 9+jj, 10+jj)
  } else { # u even
    if(u < 4){
      jj <- NULL
    } else {
      jj <- 8 * seq(from=0 , by=1 , to=round((u-4)/2) )
    }
    JC <- c(2,3,   10+jj, 17+jj)
    JD <- c(0,4,   11+jj, 16+jj)
    JS <- c(1,7,9, 12+jj, 15+jj)
    JT <- c(5,6,8, 13+jj, 14+jj)
  }
  
  S <- C <- T <- D <- A*0
  
  i <- row(A)

  for(r in seq_len(u)){
    Ar <- A==r
    Br <- B==r
    S[Br] <-         i[Br] + m*JS[r]
    C[Ar] <- (m+1) - i[Ar] + m*JC[r]
    T[Ar] <-         i[Ar] + m*JT[r]
    D[Br] <- (m+1) - i[Br] + m*JD[r]
  }
  S[B==u+1] <- i[B==u+1] + m*JS[u+1]  # 2
  T[A==u+1] <- i[A==u+1] + m*JT[u+1]  # 3
  
  force.integer(rbind(
                      cbind(C,S),
                      cbind(T,D)
                      )
                )
}

"is.antimagic" <- function(m, give.answers=FALSE, func=sum){
  if (is.list(m)) {
    return(sapply(m, match.fun(sys.call()[[1]]),
                  give.answers=give.answers,
                  func=func
                  ))
  }
  jj <- allsums(m, func=func)
  answer <- all(diff(sort(c(jj$rowsums , jj$colsums)))==1)
  if(give.answers){
    return(c(answer=answer , jj))
  } else {
    return(answer)
  }
}

"is.totally.antimagic" <- function(m, give.answers=FALSE, func=sum){
  if (is.list(m)) {
    return(sapply(m, match.fun(sys.call()[[1]]),
                  give.answers=give.answers,
                  func=func
                  ))
  }
  jj <- allsums(m, func=func)
  answer <-
    all(diff(sort(c(
                    jj$rowsums  ,
                    jj$colsums  ,
                    jj$majors[1],
                    jj$minors[1]
                    )))==1)
  if(give.answers){
    return(c(answer=answer , jj))
  } else {
    return(answer)
  }
}

"is.heterosquare" <- function(m, func = sum){
  if (is.list(m)) {
    return(sapply(m, match.fun(sys.call()[[1]]), func = func))
  }
  sums <- allsums(m, func = func)
  jj <- c(sums$rowsums, sums$colsums)
  if(all(diff(sort(jj)))>0){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

"is.totally.heterosquare" <- function(m, func = sum){
  if (is.list(m)) {
    return(sapply(m, match.fun(sys.call()[[1]]), func = func))
  }
  sums <- allsums(m, func = func)
  jj <- c(sums$rowsums, sums$colsums, sums$majors[1], sums$minors[1])
  if(all(diff(sort(jj)))>0){
    return(TRUE)
  } else {
    return(FALSE)
  }
}
  
"is.sparse" <- function(m){
  if (is.list(m)) {
    return(sapply(m, match.fun(sys.call()[[1]])))
  }
  m <- m[m != 0]
  minmax(c(1,diff(sort(m)))) &  (min(m)==1)
}
  
"is.sam" <- function(m){
  if (is.list(m)) {
    return(sapply(m, match.fun(sys.call()[[1]])))
  }
  is.antimagic(m) & is.sparse(m)
}

"is.stam" <- function(m){
  if (is.list(m)) {
    return(sapply(m, match.fun(sys.call()[[1]])))
  }
  is.totally.antimagic(m) & is.sparse(m)
}

"incidence" <- function(a){
  M <- max(a)
  d <- dim(a)
  sd <- seq_along(d)
  out <- array(0L,dim=c(d,M))
  f <- function(i){out <- rep(0L,M)
                   out[i] <- 1L
                   out
                 }
  aperm(apply(a,sd,f),c(sd+1,1))
}

"is.incidence" <- function(a, include.improper){ 
  f <- function(x){ all(x==0 | x==1) & sum(x)==1 }
  out <- is.semimagichypercube(a, func=f, boolean=TRUE)
  if(include.improper){
    return(out|is.incidence.improper(a))
  } else {
    return(out)
  }
}

"is.incidence.improper" <- function(a){
  f <- function(x){
    (all(x==0 | x==1 | x==(-1)) & sum(x)==1) |
    (all(x==0 | x==1)           & sum(x)==1)
  }
  is.semimagichypercube(a, func=f, boolean=TRUE) & (sum(a == -1) == 1)
}
 
"unincidence" <- function(a){
  stopifnot(is.incidence(a,include.improper=FALSE))
  a <- a>0
  apply(a, seq_len(length(dim(a))-1) , which)
}

"inc_to_inc" <- function(a){ # takes a proper or improper incidence
                                 # array (0/1) and returns an
                                 # incidence array, randomly chosen if
                                 # a is improper.  If 'a' is proper,
                                 # returns an improper array; if
                                 # improper, returns either a proper
                                 # or improper array.

  storage.mode(a) <- "numeric"
  randint <- function(r,n=1){ceiling(runif(n)*r)}
  stopifnot(is.incidence(a, include.improper=TRUE))
  if(is.incidence(a,include.improper=FALSE)){
    proper <- TRUE
  } else {
    proper <- FALSE
  }
  
  if(proper){ # choose a zero
    jj <- which(a==0 , arr.ind=TRUE)
    pivot <- jj[randint(nrow(jj)),,drop=TRUE]
  } else { # choose the (single) -1
    pivot <- which(a == -1, arr.ind=TRUE)
  }

  jj1 <- which(a[        ,pivot[2],pivot[3],drop=TRUE] == 1)
  jj2 <- which(a[pivot[1],        ,pivot[3],drop=TRUE] == 1)
  jj3 <- which(a[pivot[1],pivot[2],        ,drop=TRUE] == 1)
  
  if(!proper){
    jj1 <- jj1[randint(2)]
    jj2 <- jj2[randint(2)]
    jj3 <- jj3[randint(2)]
  }
   
  kk1 <- c(jj1     , pivot[2], pivot[3])
  kk2 <- c(pivot[1], jj2     , pivot[3])
  kk3 <- c(pivot[1], pivot[2], jj3     )    # a[kk[123]] == TRUE
  
  ll1 <- c(pivot[1],     jj2,     jj3)
  ll2 <- c(jj1     ,pivot[2],     jj3)
  ll3 <- c(jj1     ,     jj2,pivot[3])
  
  mm1 <- c(jj1,jj2,jj3)
  
  increment <- rbind(pivot, ll1,ll2,ll3)
  decrement <- rbind(kk1,kk2,kk3,mm1)

  a[increment] <- a[increment] + 1L
  a[decrement] <- a[decrement] - 1L
  return(a)
}

"another_latin" <- function(a){ #given a latin square, returns a _different_ one

  i <- incidence(a)
  anew <- unincidence(i)   #inefficient but clear; anew==a
  while(all(a == anew)){  # iterate until a different one is found
    i <- inc_to_inc(i)
    if(is.incidence(i,FALSE)){
      anew <- unincidence(i)
    }
  }
  return(anew)
}

"another_incidence" <- function(i){ # given a _proper_ incidence
                                    # array, returns a different
                                    # _proper_ incidence array


  out <- i
  while(all(out==i) | !is.incidence(out,FALSE)){
    out <- inc_to_inc(out)
  }
  return(out)
}

"rlatin" <- function(n,size=NULL,start=NULL,burnin=NULL){

  if(is.null(size) & is.null(start)){
    size <- n
    n <- 1
  }
  
  if(is.null(start)){
    start <- latin(size)
  } else {
    stopifnot(is.latin(start))
  }

  if(is.null(burnin)){
    burnin <- prod(dim(start))
  }
  
  out <- array(0L,c(dim(start),n))
  inc <- incidence(start)
  for(i in seq_len(burnin)){inc <- another_incidence(inc)}
  for(i in seq_len(n)){
    out[,,i] <- unincidence(inc)
    inc <- another_incidence(inc)
  }
  return(drop(out))
}

"sylvester" <- function(k){
  stopifnot(k==round(k))
  if(k==0){
    return(matrix(1L,1,1))
  } else {
    return(kronecker(Recall(k-1),matrix(c(1L,1L,1L,-1L),2,2)))
  }
}

"is.hadamard" <- function(m){
  is.matrix(m)              &
  nrow(m)==ncol(m)          &
  all( (m==1)|(m== -1))     &
  all(crossprod(m)==diag(nrow(m),nrow=nrow(m)))
}

"cilleruelo" <- function(n,m){
  matrix(c(
           (n+2)*(m+0), (n+3)*(m+3), (n+1)*(m+2), (n+0)*(m+1),
           (n+1)*(m+1), (n+0)*(m+2), (n+2)*(m+3), (n+3)*(m+0),
           (n+0)*(m+3), (n+1)*(m+0), (n+3)*(m+1), (n+2)*(m+2),
           (n+3)*(m+2), (n+2)*(m+1), (n+0)*(m+0), (n+1)*(m+3)
           ),nrow=4,ncol=4,byrow=TRUE)
}

"bernhardssonA" <- function(n){
  if(n%%2==1){return(adiag(1,Recall(n-1)))}
  out <- matrix(0L,n,n)
  m <- n/2
  j <- seq_len(m)
  out[cbind(j,2*j)] <- 1L
  out[cbind(m + j, 2*j-1 )] <- 1L
  return(out)
}
  
"bernhardssonB" <- function(n){
  if(n%%2==1){return(adiag(1,Recall(n-1)))}
  out <- matrix(0L,n,n)
  m <- n/2
  j <- seq_len(m)
  out[cbind(j,(1+(2*(j-1)+m-1))%%n)] <- 1L
  out[cbind(n+1-j,n - (2*(j-1)+m-1)%%n)] <- 1L
  return(out)
}

"bernhardsson" <- function(n){
  if( (n%%6) %in% 0:1){
    return(bernhardssonA(n))
  } else {
    return(bernhardssonB(n))
  }
}

"is.alicehypercube" <- function(a, ndim, give.answers=FALSE, func=sum, boolean=FALSE){
  stopifnot(minmax(dim(a)))
  n <- dim(a)[1]
  d <- length(dim(a))

  jj <- d-ndim
  
  out <- apply(combn(d,jj),2,function(i){apply(a,i,func)})

  if(boolean){
    answer <- all(out)
  } else {
    answer <- minmax(out)
  }

  if(give.answers){
    dim(out) <- c(rep(n,jj),ncol(out))
    return(list(answer=answer, alice.sums=out))
  } else {
    return(answer)
  }
}




"eq" <- function (m1, m2) { all(m1 == m2) }
"ne" <- function (m1, m2) { any(m1 != m2) }

"gt" <- function (m1, m2) {
  jj <- m1 - m2
  return(ne(m1, m2) && jj[min(which(jj != 0))] > 0)
}

"lt" <- function (m1, m2) {
    jj <- m1 - m2
    return(ne(m1, m2) && jj[min(which(jj != 0))] < 0)
}

"ge" <- function (m1, m2) { eq(m1, m2) || gt(m1, m2) }
"le" <- function (m1, m2) { eq(m1, m2) || lt(m1, m2) }

"%eq%" <- function (m1, m2) { return(eq(m1, m2)) }
"%ne%" <- function (m1, m2) { return(ne(m1, m2)) }
"%gt%" <- function (m1, m2) { return(gt(m1, m2)) }
"%lt%" <- function (m1, m2) { return(lt(m1, m2)) }
"%ge%" <- function (m1, m2) { return(ge(m1, m2)) }
"%le%" <- function (m1, m2) { return(le(m1, m2)) }


panmagic.6npm1 <- function(n){
   if (length(n) > 1) {
        return(sapply(n, match.fun(sys.call()[[1]])))
    }
  apx <- kronecker(t(seq(from=0,by=n-2,len=n)),rep(1,n)) + kronecker(1:n,t(rep(1,n)))
  jj <- process(apx%%n, n)
  return(force.integer(jj+n*t(jj)-n))
}

panmagic.6np1 <- function(m){ panmagic.6npm1(n=6*m+1)}
panmagic.6nm1 <- function(m){ panmagic.6npm1(n=6*m-1)}

panmagic.4n <- function(m){ # returns a square of size [4n x 4n]
  if (length(m) > 1) {
    return(sapply(m, match.fun(sys.call()[[1]])))
  }
  jj <- kronecker(rep(1,m*2),rbind(1:(2*m), (4*m):(2*m+1)))
  jj <- cbind(jj,ashift(jj,v=c(1,0)))
  return(force.integer(jj + 4*m*(arot(jj)-1)))
}


  
