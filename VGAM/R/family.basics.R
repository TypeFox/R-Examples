# These functions are
# Copyright (C) 1998-2015 T.W. Yee, University of Auckland.
# All rights reserved.











subsetcol <-
Select <-
  function(
           data = list(),
           prefix = "y",

           lhs  = NULL,
           rhs  = NULL,  # Can be "0" to suppress an intercept, else "".
           rhs2 = NULL,  # Can be "0" to suppress an intercept, else "".
           rhs3 = NULL,  # Can be "0" to suppress an intercept, else "".

           as.character = FALSE,
           as.formula.arg = FALSE,
           tilde = TRUE,
           exclude = NULL,
           sort.arg = TRUE) {





  if (is.character(exclude))
    if (any(nchar(prefix) == 0))
      stop("bad input for argument 'exclude'")
  if (!is.logical(sort.arg) || 
      length(sort.arg) != 1)
    stop("bad input for argument 'sort.arg'")


  col.names <- colnames(data)
  if (is.logical(prefix)) {
    index <- if (prefix) 1:length(col.names) else
             stop("cannot have 'prefix = FALSE'")
  } else {
    index <- NULL
    for (ii in 1:length(prefix)) {
      small.col.names <- substr(col.names, 1, nchar(prefix[ii]))

      index <- c(index, grep(prefix[ii], small.col.names))
    }
  }





  temp.col.names <- col.names[index]


  if (length(exclude)) {
    exclude.index <- NULL
    for (ii in 1:length(exclude)) {
      exclude.index <- c(exclude.index,
                         (1:length(col.names))[exclude[ii] == col.names])
    }
    exclude.index <- unique(sort(exclude.index))
    index <- setdiff(index, exclude.index)
    temp.col.names <- col.names[index]
  }




  if (sort.arg) {
    ooo <- order(temp.col.names)
    index <- index[ooo]
    temp.col.names <- temp.col.names[ooo]
  }



 ltcn.positive <- (length(temp.col.names) > 0)



  if (as.formula.arg) {
    form.string <- paste0(ifelse(length(lhs), lhs, ""),
                          ifelse(tilde, " ~ ", ""),
                          if (ltcn.positive)
                            paste(temp.col.names, collapse = " + ") else
                            "",
                          ifelse(ltcn.positive && length(rhs ), " + ", ""),
                          ifelse(length(rhs ), rhs, ""),
                          ifelse(length(rhs2), paste(" +", rhs2), ""),
                          ifelse(length(rhs3), paste(" +", rhs3), ""))

    if (as.character) {
      form.string
    } else {
      as.formula(form.string)
    }
  } else {
    if (as.character) {
      paste0("cbind(",
             paste(temp.col.names, collapse = ", "),
             ")")
    } else {
      ans <- if (is.matrix(data)) data[, index] else
             if (is.list(data)) data[index] else
             stop("argument 'data' is neither a list or a matrix")
      if (length(ans)) {
        as.matrix(ans)
      } else {
        NULL
      }
    }
  }
}








if (FALSE)
subsetc <-
  function(x, select,
           prefix = NULL,
           subset = TRUE, drop = FALSE,
           exclude = NULL,
           sort.arg = !is.null(prefix),
           as.character = FALSE) {

  if (!is.null(prefix)) {
    if (!missing(select))
      warning("overwriting argument 'select' by something ",
              "using 'prefix'")
    select <- grepl(paste0("^", prefix), colnames(x))
  }

  if (missing(select)) {
    vars <- TRUE
  } else {
    nl <- as.list(seq_along(x))  # as.list(1L:ncol(x))
    names(nl) <- names(x)  # colnames(x)
    vars <- eval(substitute(select), nl, parent.frame())
  }

  ans <- x[subset & !is.na(subset), vars, drop = drop]
  if (sort.arg) {
    cna <- colnames(ans)
    ooo <- order(cna)
    ans <- ans[, ooo, drop = drop]
  }

  if (!is.null(exclude)) {
    cna <- colnames(ans)
    ooo <- match(exclude, cna)
    ans <- ans[, -ooo, drop = drop]
  }

  if (as.character) {
    cna <- colnames(ans)
    paste0("cbind(", paste0(cna, collapse = ", "), ")")
  } else {
    ans
  }
}







 grid.search <- function(vov, objfun, y, x, w, extraargs = NULL,
                         maximize = TRUE, abs.arg = FALSE,
                         ret.objfun = FALSE, ...) {
  if (!is.vector(vov))
    stop("argument 'vov' must be a vector")
  objvals <- vov
  for (ii in 1:length(vov))
    objvals[ii] <- objfun(vov[ii], y = y, x = x, w = w,
                          extraargs = extraargs, ...)
  try.this <- if (abs.arg) {
               if (maximize) vov[abs(objvals) == max(abs(objvals))] else
               vov[abs(objvals) == min(abs(objvals))]
             } else {
               if (maximize) vov[objvals == max(objvals)] else
               vov[objvals == min(objvals)]
             }
  if (!length(try.this))
    stop("something has gone wrong!")
  ans <- if (length(try.this) == 1)
    try.this else sample(try.this, size = 1)


  myvec <- objvals[ans == vov]  # Could be a vector
  if (ret.objfun) c(Value = ans, ObjFun = myvec[1]) else ans
}








 getind <- function(constraints, M, ncolx) {



  if (!length(constraints)) {

    constraints <- vector("list", ncolx)
    for (ii in 1:ncolx)
      constraints[[ii]] <- diag(M)
  }

  ans <- vector("list", M+1)
  names(ans) <- c(paste("eta", 1:M, sep = ""), "ncolX.vlm")

  temp2 <- matrix(unlist(constraints), nrow = M)
  for (kk in 1:M) {
    ansx <- NULL
    for (ii in 1:length(constraints)) {
      temp <- constraints[[ii]]
      isfox <- any(temp[kk, ] != 0)
      if (isfox) {
        ansx <- c(ansx, ii)
      }
    }
    ans[[kk]] <- list(xindex = ansx,
                  X.vlmindex = (1:ncol(temp2))[temp2[kk,] != 0])
  }
  ans[[M+1]] <- ncol(temp2)

  ans
}




 cm.VGAM <-
  function(cm, x, bool, constraints,
           apply.int = FALSE, 
           cm.default = diag(nrow(cm)),  # 20121226
           cm.intercept.default = diag(nrow(cm))  # 20121226
          ) {




  if (is.null(bool))
    return(NULL)

  if (!is.matrix(cm))
    stop("argument 'cm' is not a matrix")
  M <- nrow(cm)
  asgn <- attr(x, "assign")
  if (is.null(asgn))
    stop("the 'assign' attribute is missing from 'x'; this ",
         "may be due to some missing values")  # 20100306
  nasgn <- names(asgn)
  ninasgn <- nasgn[nasgn != "(Intercept)"]

  if (!length(constraints)) {
    constraints <- vector("list", length(nasgn))
    for (ii in 1:length(nasgn)) {
      constraints[[ii]] <- cm.default  # diag(M)
    }
    names(constraints) <- nasgn


    if (any(nasgn == "(Intercept)"))
      constraints[["(Intercept)"]] <- cm.intercept.default
  } 

  if (!is.list(constraints))
    stop("argument 'constraints' must be a list")

  if (length(constraints) != length(nasgn) ||
      any(sort(names(constraints)) != sort(nasgn))) {
    cat("\nnames(constraints)\n")
    print(names(constraints) )
    cat("\nnames(attr(x, 'assign'))\n")
    print( nasgn )
    stop("The above do not match; 'constraints' is half-pie")
  }




  if (is.logical(bool)) {
    if (bool) {
      if (any(nasgn == "(Intercept)") && apply.int)
        constraints[["(Intercept)"]] <- cm


      if (length(ninasgn))
        for (ii in ninasgn)
          constraints[[ii]] <- cm
    } else {
      return(constraints)
    }
  } else {
    tbool <- terms(bool)
    if (attr(tbool, "response")) {
      ii <- attr(tbool, "factors")
      default <- dimnames(ii)[[1]]
      default <- default[1]
      default <- if (is.null(default[1])) {
        t.or.f <- attr(tbool, "variables")

        t.or.f <- as.character( t.or.f )
        if (t.or.f[1] == "list" && length(t.or.f) == 2 &&
           (t.or.f[2] == "TRUE" || t.or.f[2] == "FALSE")) {
          t.or.f <- as.character( t.or.f[2] )
          parse(text = t.or.f)[[1]]
        } else {
          stop("something gone awry")
        }
      } else {
        parse(text = default[1])[[1]]  # Original
      }
      default <- as.logical(eval(default))
    } else {
      default <- TRUE
    }
    tl <- attr(tbool, "term.labels")
    if (attr(tbool, "intercept"))
      tl <- c("(Intercept)", tl)

    for (ii in nasgn) {
      if ( default &&  any(tl == ii))
        constraints[[ii]] <- cm
      if (!default && !any(tl == ii))
        constraints[[ii]] <- cm
    }
  }

  constraints
}




cm.nointercept.VGAM <- function(constraints, x, nointercept, M) {

  asgn <- attr(x, "assign")
  nasgn <- names(asgn)
  if (is.null(constraints)) {
    constraints <- vector("list", length(nasgn))  # list()
    names(constraints) <- nasgn
  }
  if (!is.list(constraints))
    stop("'constraints' must be a list")

  for (ii in 1:length(asgn))
    constraints[[nasgn[ii]]] <- if (is.null(constraints[[nasgn[ii]]]))
      diag(M) else eval(constraints[[nasgn[ii]]])

  if (is.null(nointercept))
    return(constraints)
  if (!is.numeric(nointercept))
    stop("'nointercept' must be numeric")

  nointercept <- unique(sort(nointercept))
  if (length(nointercept) == 0 || length(nointercept) >= M)
    stop("too few or too many values")

  if (any(nointercept < 1 | nointercept > M))
    stop("'nointercept' out of range")
  if (nasgn[1] != "(Intercept)" || M == 1)
    stop("Need an (Intercept) constraint matrix with M>1")
  if (!all.equal(constraints[["(Intercept)"]], diag(M)))
    warning("Constraint matrix of (Intercept) not diagonal")

  temp <- constraints[["(Intercept)"]]
  temp <- temp[, -nointercept, drop = FALSE]
  constraints[["(Intercept)"]] <- temp 
  constraints
}




 cm.zero.VGAM <- function(constraints, x, zero = NULL, M = 1,
                          predictors.names, M1 = 1) {


  dotzero <- zero  # Transition

  if (is.character(dotzero)) {




  which.numeric.all <- NULL
  for (ii in 1:length(dotzero)) {
    which.ones <-
        grep(dotzero[ii], predictors.names, fixed = TRUE)
    if (length(which.ones)) {
      which.numeric.all <- c(which.numeric.all, which.ones)
    } else {
      warning("some values of argument 'zero' are unmatched. Ignoring them")
    }
  }
  which.numeric <- unique(sort(which.numeric.all))

  if (!length(which.numeric)) {
    warning("No values of argument 'zero' were matched.")
    which.numeric <- NULL
  } else if (length(which.numeric.all) > length(which.numeric)) {
    warning("There were redundant values of argument 'zero'.")
  }

    dotzero <- which.numeric
  }



  posdotzero <-  dotzero[dotzero > 0]
  negdotzero <-  dotzero[dotzero < 0]


  zneg.index <- if (length(negdotzero)) {

    if (!is.Numeric(-negdotzero, positive = TRUE,
                    integer.valued = TRUE) ||
        max(-negdotzero) > M1)
        stop("bad input for argument 'zero'")

    bigUniqInt <- 1080
    zneg.index <- rep(0:bigUniqInt, rep(length(negdotzero),
                      1 + bigUniqInt)) * M1 + abs(negdotzero)
    sort(intersect(zneg.index, 1:M))
  } else {
    NULL
  }

  zpos.index <- if (length(posdotzero)) posdotzero else NULL
  z.Index <- if (!length(dotzero))
               NULL else
               unique(sort(c(zneg.index, zpos.index)))


  zero <- z.Index  # Transition



  asgn <- attr(x, "assign")
  nasgn <- names(asgn)
  if (is.null(constraints)) {
    constraints <- vector("list", length(nasgn))  # list()
    names(constraints) <- nasgn
  }
  if (!is.list(constraints))
    stop("'constraints' must be a list")

  for (ii in 1:length(asgn))
    constraints[[nasgn[ii]]] <- if (is.null(constraints[[nasgn[ii]]]))
      diag(M) else eval(constraints[[nasgn[ii]]])

  if (is.null(zero))
    return(constraints)

  if (any(zero < 1 | zero > M))
    stop("argument 'zero' out of range")
  if (nasgn[1] != "(Intercept)")
    stop("cannot fit an intercept to a no-intercept model")

  if (2 <= length(constraints))
    for (ii in 2:length(constraints)) {
      Hmatk <- constraints[[nasgn[ii]]]
      Hmatk[zero, ] <- 0
      index <- NULL
      for (kk in 1:ncol(Hmatk))
        if (all(Hmatk[, kk] == 0)) index <- c(index, kk)
      if (length(index) == ncol(Hmatk)) 
        stop("constraint matrix has no columns!")
      if (!is.null(index))
        Hmatk <- Hmatk[, -index, drop = FALSE]
      constraints[[nasgn[ii]]] <- Hmatk 
    }
  constraints
}



 process.constraints <-
  function(constraints, x, M,
           by.col = TRUE, specialCM = NULL,
           Check.cm.rank = TRUE  # 20140626
          ) {







    asgn <- attr(x, "assign")
    nasgn <- names(asgn)

  if (is.null(constraints)) {
    constraints <- vector("list", length(nasgn))
    for (ii in 1:length(nasgn))
      constraints[[ii]] <- diag(M)
    names(constraints) <- nasgn
  }

  if (is.matrix(constraints))
    constraints <- list(constraints)

  if (!is.list(constraints))
    stop("'constraints' must be a list")

  lenconstraints <- length(constraints)
  if (lenconstraints > 0)
    for (ii in 1:lenconstraints) {
      list.elt <- constraints[[ii]]

      if (is.function(list.elt)) {
        list.elt <- list.elt()
      }

      constraints[[ii]] <- eval(list.elt)
      if (!is.null  (constraints[[ii]]) &&
          !is.matrix(constraints[[ii]]))
        stop("'constraints[[", ii, "]]' is not a matrix")
    }

  if (is.null(names(constraints))) 
    names(constraints) <- rep(nasgn, length.out = lenconstraints) 

  temp <- vector("list", length(nasgn))
  names(temp) <- nasgn
  for (ii in 1:length(nasgn))
    temp[[nasgn[ii]]] <-
      if (is.null(constraints[[nasgn[ii]]])) diag(M) else
             eval(constraints[[nasgn[ii]]])

  for (ii in 1:length(asgn)) {
    if (!is.matrix(temp[[ii]])) {
      stop("not a constraint matrix")
    }
    if (ncol(temp[[ii]]) > M)
      stop("constraint matrix has too many columns")
  }

  if (!by.col)
    return(temp)

  constraints <- temp
  Hlist <- vector("list", ncol(x))
  for (ii in 1:length(asgn)) {
    cols <- asgn[[ii]]
    ictr <- 0
    for (jay in cols) {
      ictr <- ictr + 1
      cm <- if (is.list(specialCM) &&
                any(nasgn[ii] == names(specialCM))) {
              slist <- specialCM[[(nasgn[ii])]]
              slist[[ictr]]
            } else {
              constraints[[ii]]
            }
      Hlist[[jay]] <- cm 
    }
  }
  names(Hlist) <- dimnames(x)[[2]]



  if (Check.cm.rank) {
    all.svd.d <- function(x) svd(x)$d
    mylist <- lapply(Hlist, all.svd.d)

    if (max(unlist(lapply(mylist, length))) > M)
      stop("some constraint matrices have more than ", M,
           "columns")

    MyVector <- unlist(mylist)
    if (min(MyVector) < 1.0e-10)
      stop("some constraint matrices are not of ",
           "full column-rank: ",
           paste(names(MyVector)[MyVector < 1.0e-10], collapse = ", "))
  }

  Hlist
}





 trivial.constraints <- function(Hlist, target = diag(M)) {


  if (is.null(Hlist))
    return(1)

  if (is.matrix(Hlist))
    Hlist <- list(Hlist)
  M <- dim(Hlist[[1]])[1]

  if (!is.matrix(target)) 
    stop("target is not a matrix")
  dimtar <- dim(target) 

  trivc <- rep(1, length(Hlist))
  names(trivc) <- names(Hlist)
  for (ii in 1:length(Hlist)) {
    d <- dim(Hlist[[ii]])
    if (d[1] != dimtar[1]) trivc[ii] <- 0
    if (d[2] != dimtar[2]) trivc[ii] <- 0
    if (d[1] != M)         trivc[ii] <- 0
    if (length(Hlist[[ii]]) != length(target))
      trivc[ii] <- 0
    if (trivc[ii] == 0) next
    if (!all(c(Hlist[[ii]]) == c(target)))
      trivc[ii] <- 0
    if (trivc[ii] == 0) next
  }
  trivc
}



 add.constraints <- function(constraints, new.constraints,
                             overwrite = FALSE, check = FALSE) {

  empty.list <- function(l)
    (is.null(l) || (is.list(l) && length(l) == 0))

  if (empty.list(constraints))
    if (is.list(new.constraints))
      return(new.constraints) else 
      return(list())  # Both NULL probably

  constraints <- as.list(constraints)
  new.constraints <- as.list(new.constraints)
  nc <- names(constraints)         # May be NULL
  nn <- names(new.constraints)     # May be NULL

  if (is.null(nc) || is.null(nn))
    stop("lists must have names")
  if (any(nc == "") || any(nn == ""))
    stop("lists must have names")

  if (!empty.list(constraints) && !empty.list(new.constraints)) {
    for (ii in nn) {
      if (any(ii == nc)) {
        if (check &&
        (!(all(dim(constraints[[ii]]) == dim(new.constraints[[ii]])) &&
           all(    constraints[[ii]]  ==     new.constraints[[ii]]))))
          stop("apparent contradiction in the specification ",
               "of the constraints")
        if (overwrite)
          constraints[[ii]] <- new.constraints[[ii]]
      } else 
        constraints[[ii]] <- new.constraints[[ii]]
    }
  } else {
    if (!empty.list(constraints))
      return(as.list(constraints)) else
      return(as.list(new.constraints))
  }

  constraints
}















 iam <- function(j, k, M,  # hbw = M,
                 both = FALSE, diag = TRUE) {


  jay <- j 
  kay <- k

  if (M == 1)
    if (!diag) stop("cannot handle this") 

  if (M == 1)
    if (both) return(list(row.index = 1, col.index = 1)) else return(1)

  upper <- if (diag) M else M - 1
  i2 <- as.list(upper:1)
  i2 <- lapply(i2, seq)
  i2 <- unlist(i2)

  i1 <- matrix(1:M, M, M) 
  i1 <- if (diag) c(i1[row(i1) >= col(i1)]) else
                  c(i1[row(i1) >  col(i1)])

  if (both) {
    list(row.index = i2, col.index = i1)
  } else {
    if (jay > M || kay > M || jay < 1 || kay < 1)
      stop("range error in j or k")
    both <- (i1 == jay & i2 == kay) |
            (i1 == kay & i2 == jay)
    (1:length(i2))[both]
  }
}





 dimm <- function(M, hbw = M) {

  if (!is.numeric(hbw))
    hbw <- M

  if (hbw > M || hbw < 1)
    stop("range error in argument 'hbw'")
  hbw * (2 * M - hbw +1) / 2 
}









 m2a <- function(m, M, upper = FALSE, allow.vector = FALSE) {


  if (!is.numeric(m))
      stop("argument 'm' is not numeric")

  if (!is.matrix(m))
    m <- cbind(m)
  n <- nrow(m)
  dimm <- ncol(m)
  index <- iam(NA, NA, M = M, both = TRUE, diag = TRUE)
  if (dimm > length(index$row.index))
    stop("bad value for 'M'; it is too small") 
  if (dimm < M) {
    stop("bad value for 'M'; it is too big") 
  }

  fred <- .C("m2a", as.double(t(m)), ans=double(M*M*n),
      as.integer(dimm),
      as.integer(index$row-1),  
      as.integer(index$col-1),  
      as.integer(n),  as.integer(M),  
      as.integer(as.numeric(upper)), NAOK = TRUE)
  dim(fred$ans) <- c(M, M, n)
  alpn <- NULL
  dimnames(fred$ans) <- list(alpn, alpn, dimnames(m)[[1]])
  fred$a
}



 a2m <- function(a, hbw = M) {



  if (is.matrix(a) && ncol(a) == nrow(a))
    a <- array(a, c(nrow(a), ncol(a), 1))
  if (!is.array(a))
    dim(a) <- c(1,1,length(a))

  M <- dim(a)[1]
  n <- dim(a)[3]
  dimm.value <- dimm(M, hbw)
  index <- iam(NA, NA, M, both = TRUE, diag = TRUE)


  fred <- .C("a2m",
             as.double(a), m = double(dimm.value * n),
      as.integer(dimm.value),
      as.integer(index$row-1),  as.integer(index$col-1),  
      as.integer(n),  as.integer(M), NAOK = TRUE)
  dim(fred$m) <- c(dimm.value,n)
  fred$m <- t(fred$m)

  if (hbw != M) 
    attr(fred$m, "hbw") <- hbw
  if (length(lpn <- dimnames(a)[[1]]) != 0)
    attr(fred$m, "predictors.names") <- lpn
  fred$m
}



 vindex <- function(M, row.arg = FALSE, col.arg = FALSE,
                    length.arg = M * (M + 1) / 2) {



  if ((row.arg + col.arg) != 1)
    stop("only one of row and col must be TRUE") 
  if (M == 1) {
    ans <- 1
  } else {
    if (row.arg) {
      i1 <- matrix(1:M, M, M)
      ans <- c(i1[row(i1) + col(i1) <= (M + 1)])
    } else {
      i1 <- matrix(1:M, M, M) 
      ans <- c(i1[row(i1) >= col(i1)])
    }
  }
  if (length.arg > length(ans))
    stop("argument 'length.arg' too big")
  rep(ans, length.out = length.arg) 
}





 wweights <- function(object, matrix.arg = TRUE, deriv.arg = FALSE,
                      ignore.slot = FALSE, checkwz = TRUE) {







  if (length(wz <- object@weights) && !ignore.slot && !deriv.arg) { 
    return(wz) 
  }

  M <- object@misc$M  # Done below
  n <- object@misc$n  # Done below

  if (any(slotNames(object) == "extra")) {
    extra <- object@extra
    if (length(extra) == 1 && !length(names(extra))) {
      extra <- extra[[1]]
    }
  }
  mu <- object@fitted.values
  if (any(slotNames(object) == "predictors"))
    eta <- object@predictors
  mt <- terms(object)  # object@terms$terms; 20030811
  Hlist <- constraints <- object@constraints 
  new.coeffs <- object@coefficients
  if (any(slotNames(object) == "iter"))
    iter <- object@iter

  w <- rep(1, n)
  if (any(slotNames(object) == "prior.weights"))
    w <- object@prior.weights
  if (!length(w))
    w <- rep(1, n)

  x <- object@x
  if (!length(x))
    x <- model.matrixvlm(object, type = "lm")

  y <- object@y
  if (!length(y))
    y <- depvar(object)



  if (length(object@misc$form2)) {
    Xm2 <- object@Xm2
    if (!length(Xm2))
      Xm2 <- model.matrix(object, type = "lm2")
    Ym2 <- object@Ym2
  }




  if (any(slotNames(object) == "family")) {
    infos.list <- object@family@infos()
    if (length(infos.list))
      for (ii in names(infos.list)) {
        assign(ii, infos.list[[ii]])
      }
  }



  if (any(slotNames(object) == "control"))
    for (ii in names(object@control)) {
      assign(ii, object@control[[ii]])
    } 

  if (length(object@misc))
    for (ii in names(object@misc)) {
      assign(ii, object@misc[[ii]]) 
    } 

  if (any(slotNames(object) == "family")) {
    expr <- object@family@deriv
    deriv.mu <- eval(expr)
    if (!length(wz)) {
      expr <- object@family@weight
      wz <- eval(expr)


      if (M > 1) 
        dimnames(wz) <- list(dimnames(wz)[[1]], NULL)  # Remove colnames
      wz <- if (matrix.arg) as.matrix(wz) else c(wz) 
    }
    if (deriv.arg) list(deriv = deriv.mu, weights = wz) else wz
  } else {
    NULL 
  }
}




 pweights <- function(object, ...) {
  ans <- object@prior.weights
  if (length(ans)) {
    ans 
  } else {
    temp <- object@y
    ans <- rep(1, nrow(temp))  # Assumed all equal and unity.
    names(ans) <- dimnames(temp)[[1]]
    ans 
  }
}





procVec <- function(vec, yn, Default) {








  if (any(is.na(vec)))
    stop("vec cannot contain any NAs")
  L <- length(vec)
  nvec <- names(vec)  # vec[""] undefined
  named <- length(nvec)  # FALSE for c(1,3)
  if (named) {
    index <- (1:L)[nvec == ""]
    default <- if (length(index)) vec[index] else Default
  } else {
    default <- vec
  }

  answer <- rep(default, length.out = length(yn))
  names(answer) <- yn
  if (named) {
    nvec2 <- nvec[nvec != ""]
    if (length(nvec2)) {
      if (any(!is.element(nvec2, yn)))
          stop("some names given which are superfluous")
      answer <- rep(NA_real_, length.out = length(yn))
      names(answer) <- yn
      answer[nvec2] <- vec[nvec2]
      answer[is.na(answer)] <-
        rep(default, length.out <- sum(is.na(answer)))
    }
  }

  answer
}



if (FALSE) {


}



 weightsvglm <- function(object, type = c("prior", "working"),
                        matrix.arg = TRUE, ignore.slot = FALSE,
                        deriv.arg = FALSE, ...) {
  weightsvlm(object, type = type, matrix.arg = matrix.arg,
              ignore.slot = ignore.slot,
              deriv.arg = deriv.arg, ...)
}



 weightsvlm <- function(object, type = c("prior", "working"),
                       matrix.arg = TRUE, ignore.slot = FALSE,
                       deriv.arg = FALSE, ...) {
  if (mode(type) != "character" && mode(type) != "name")
    type <- as.character(substitute(type))
  type <- match.arg(type, c("prior", "working"))[1]

  if (type == "working") {
    wweights(object = object,
             matrix.arg = matrix.arg, deriv.arg = deriv.arg,
             ignore.slot = ignore.slot, ...)
  } else {
    if (deriv.arg)
      stop("cannot set 'deriv = TRUE' when 'type=\"prior\"'")
    ans <- pweights(object)
    if (matrix.arg) as.matrix(ans) else c(ans)
  }
}


if (!isGeneric("weights"))
    setGeneric("weights", function(object, ...)
  standardGeneric("weights"))


setMethod("weights", "vlm",
         function(object, ...)
         weightsvlm(object, ...))


setMethod("weights", "vglm",
         function(object, ...)
         weightsvglm(object, ...))













qnupdate <- function(w, wzold, dderiv, deta, M, keeppd = TRUE, 
                    trace = FALSE, reset = FALSE,
                    effpos=.Machine$double.eps^0.75) {


  if (M == 1) {
    dderiv <- cbind(dderiv)
    deta <- cbind(deta)
  }
  Bs <- mux22(t(wzold), deta, M = M,
              upper = FALSE, as.matrix = TRUE)  # n x M
  sBs <- c( (deta * Bs) %*% rep(1, M) )  # should have positive values
  sy <- c( (dderiv * deta) %*% rep(1, M) )
  wznew <- wzold
  index <- iam(NA, NA, M = M, both = TRUE)
  index$row.index <- rep(index$row.index, len=ncol(wzold))
  index$col.index <- rep(index$col.index, len=ncol(wzold))
  updateThese <- if (keeppd) (sy > effpos) else rep(TRUE, len=length(sy))
  if (!keeppd || any(updateThese)) {
    wznew[updateThese,] <- wznew[updateThese,] -
        Bs[updateThese,index$row] *
        Bs[updateThese,index$col] / sBs[updateThese] +
        dderiv[updateThese,index$row] *
        dderiv[updateThese,index$col] / sy[updateThese]
    notupdated <- sum(!updateThese)
    if (notupdated && trace)
      cat(notupdated,
          "weight matrices not updated out of", length(sy), "\n")
  } else {
    warning("no BFGS quasi-Newton update made at all")
    cat("no BFGS quasi-Newton update made at all\n")
    flush.console()
  }
  wznew
}







mbesselI0 <- function(x, deriv.arg = 0) {
  if (!is.Numeric(deriv.arg, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) &&
      deriv.arg != 0)
    stop("argument 'deriv.arg' must be a single non-negative integer")
  if (!(deriv.arg == 0 || deriv.arg == 1 || deriv.arg == 2))
    stop("argument 'deriv' must be 0, 1, or 2")
  if (!is.Numeric(x))
    stop("bad input for argument 'x'")
  nn <- length(x)

  if (FALSE) {
    }

    ans <- matrix(NA_real_, nrow = nn, ncol = deriv.arg+1)
    ans[, 1] <- besselI(x, nu = 0)
    if (deriv.arg>=1) ans[,2] <- besselI(x, nu = 1) 
    if (deriv.arg>=2) ans[,3] <- ans[,1] - ans[,2] / x
    ans
}



VGAM.matrix.norm <- function(A, power = 2, suppressWarning = FALSE) {
  if ((nrow(A) != ncol(A)) && !suppressWarning)
    warning("norms should be calculated for square matrices; ",
            "'A' is not square")
  if (power == "F") {
    sqrt(sum(A^2)) 
  } else if (power == 1) {
    max(colSums(abs(A)))
  } else if (power == 2) {
    sqrt(max(eigen(t(A) %*% A)$value))
  } else if (!is.finite(power)) {
    max(colSums(abs(A)))
  } else {
    stop("argument 'power' not recognized")
  }
}






rmfromVGAMenv <- function(varnames, prefix = "") {
  evarnames <- paste(prefix, varnames, sep = "")
  for (ii in evarnames) {
    mytext1 <- "exists(x = ii, envir = VGAMenv)"
    myexp1 <- parse(text = mytext1)
    is.there <- eval(myexp1)
    if (is.there) {
      rm(list = ii, envir = VGAMenv)
    }
  }
}




existsinVGAMenv <- function(varnames, prefix = "") {
  evarnames <- paste(prefix, varnames, sep = "")
  ans <- NULL
  for (ii in evarnames) {
    mytext1 <- "exists(x = ii, envir = VGAMenv)"
    myexp1 <- parse(text = mytext1)
    is.there <- eval(myexp1)
    ans <- c(ans, is.there)
  }
  ans
}


assign2VGAMenv <- function(varnames, mylist, prefix = "") {
  evarnames <- paste(prefix, varnames, sep = "")
  for (ii in 1:length(varnames)) {
    assign(evarnames[ii], mylist[[(varnames[ii])]],
           envir = VGAMenv)
  }
}






getfromVGAMenv <- function(varname, prefix = "") {
  varname <- paste(prefix, varname, sep = "")
  if (length(varname) > 1)
    stop("'varname' must be of length 1")
  get(varname, envir = VGAMenv)
}

 

lerch <- function(x, s, v, tolerance = 1.0e-10, iter = 100) {
  if (!is.Numeric(x) || !is.Numeric(s) || !is.Numeric(v))
    stop("bad input in 'x', 's', and/or 'v'")
  if (is.complex(c(x,s,v)))
    stop("complex arguments not allowed in 'x', 's' and 'v'")
  if (!is.Numeric(tolerance, length.arg = 1, positive = TRUE) ||
      tolerance > 0.01)
    stop("bad input for argument 'tolerance'")
  if (!is.Numeric(iter, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'iter'")

  L <- max(length(x), length(s), length(v))
  x <- rep(x, length.out = L);
  s <- rep(s, length.out = L);
  v <- rep(v, length.out = L);
  xok <- abs(x) < 1 & !(v <= 0 & v == round(v))
  x[!xok] <- 0  # Fix this later

  ans <- .C("lerchphi123",
           err = integer(L), as.integer(L),
           as.double(x), as.double(s), as.double(v),
           acc=as.double(tolerance), result=double(L),
           as.integer(iter))

  ifelse(ans$err == 0 & xok , ans$result, NA)
}




negzero.expression.VGAM <- expression({











  if (is.character(dotzero)) {




  which.numeric.all <- NULL
  for (ii in 1:length(dotzero)) {
    which.ones <-
        grep(dotzero[ii], predictors.names, fixed = TRUE)
    if (length(which.ones)) {
      which.numeric.all <- c(which.numeric.all, which.ones)
    } else {
      warning("some values of argument 'zero' are unmatched. Ignoring them")
    }
  }
  which.numeric <- unique(sort(which.numeric.all))

  if (!length(which.numeric)) {
    warning("No values of argument 'zero' were matched.")
    which.numeric <- NULL
  } else if (length(which.numeric.all) > length(which.numeric)) {
    warning("There were redundant values of argument 'zero'.")
  }

    dotzero <- which.numeric
  }



  posdotzero <-  dotzero[dotzero > 0]
  negdotzero <-  dotzero[dotzero < 0]

  zneg.index <- if (length(negdotzero)) {

    if (!is.Numeric(-negdotzero, positive = TRUE,
                    integer.valued = TRUE) ||
        max(-negdotzero) > M1)
        stop("bad input for argument 'zero'")

    bigUniqInt <- 1080
    zneg.index <- rep(0:bigUniqInt, rep(length(negdotzero),
                      1 + bigUniqInt)) * M1 + abs(negdotzero)
    sort(intersect(zneg.index, 1:M))
  } else {
    NULL
  }

  zpos.index <- if (length(posdotzero)) posdotzero else NULL
  z.Index <- if (!length(dotzero))
               NULL else
               unique(sort(c(zneg.index, zpos.index)))

  constraints <- cm.zero.VGAM(constraints, x = x, z.Index, M = M)
})






is.empty.list <- function(mylist) {
  is.list(mylist) &&
  length(unlist(mylist)) == 0
}









  interleave.VGAM  <- function(.M, M1, inverse = FALSE) {
  if (inverse) {
    NRs <- (.M)/M1
    if (round(NRs) != NRs)
      stop("Incompatible number of parameters")
    c(matrix(1:(.M), nrow = NRs, byrow = TRUE))
  } else {
    c(matrix(1:(.M), nrow = M1, byrow = TRUE))
  }
}






interleave.cmat <- function(cmat1, cmat2) {
  ncol1 <- ncol(cmat1)
  ncol2 <- ncol(cmat2)
  if (ncol1 == 1) {
    return(cbind(cmat1, cmat2))
  } else {  # ncol1 > 1
    if (ncol2 == 1) {
      return(cbind(cmat1[, 1], cmat2, cmat1[, -1]))
    } else
    if (ncol1 != ncol2) {
      warning("this function is confused. Returning cbind(cmat1, cmat2)")
      return(cbind(cmat1[, 1], cmat2, cmat1[, -1]))
    } else {  # ncol1 == ncol2 and both are > 1.
      kronecker(cmat1, cbind(1, 0)) +
      kronecker(cmat2, cbind(0, 1))        
    }
  }
}






w.wz.merge <- function(w, wz, n, M, ndepy,
                       intercept.only = FALSE) {





  wz <- as.matrix(wz)

  if (ndepy == 1)
    return( c(w) * wz)


  if (intercept.only)
    warning("yettodo: support intercept.only == TRUE")

  if (ncol(as.matrix(w)) > ndepy)
    stop("number of columns of 'w' exceeds number of responses")

  w  <- matrix(w, n, ndepy)
  w.rep <- matrix(0, n, ncol(wz))
  M1 <- M / ndepy
  all.indices <- iam(NA, NA, M = M, both = TRUE)



  if (FALSE)
  for (ii in 1:ncol(wz)) {

    if ((ind1 <- ceiling(all.indices$row[ii] / M1)) ==
                 ceiling(all.indices$col[ii] / M1)) {
      w.rep[, ii] <- w[, ind1]
    }


  }  # ii


  res.Ind1 <- ceiling(all.indices$row.index / M1)
  Ind1 <- res.Ind1 == ceiling(all.indices$col.index / M1)

  LLLL <- min(ncol(wz), length(Ind1))
  Ind1 <- Ind1[1:LLLL]
  res.Ind1 <- res.Ind1[1:LLLL]

  for (ii in 1:ndepy) {
    sub.ind1 <- (1:LLLL)[Ind1 & (res.Ind1 == ii)]
    w.rep[, sub.ind1] <- w[, ii]
  }  # ii

  w.rep * wz
}






w.y.check <- function(w, y,
                      ncol.w.max = 1, ncol.y.max = 1,
                      ncol.w.min = 1, ncol.y.min = 1,
                      out.wy = FALSE,
                      colsyperw = 1,
                      maximize = FALSE,
                      Is.integer.y = FALSE,
                      Is.positive.y = FALSE,
                      Is.nonnegative.y = FALSE,
                      prefix.w = "PriorWeight",
                      prefix.y = "Response") {



  if (!is.matrix(w))
    w <- as.matrix(w)
  if (!is.matrix(y))
    y <- as.matrix(y)
  n.lm <- nrow(y)
  rn.w <- rownames(w)
  rn.y <- rownames(y)
  cn.w <- colnames(w)
  cn.y <- colnames(y)


  if (Is.integer.y && any(y != round(y)))
    stop("response variable 'y' must be integer-valued")
  if (Is.positive.y && any(y <= 0))
    stop("response variable 'y' must be positive-valued")
  if (Is.nonnegative.y && any(y < 0))
    stop("response variable 'y' must be 0 or positive-valued")

  if (nrow(w) != n.lm)
    stop("nrow(w) should be equal to nrow(y)")

  if (ncol(w) > ncol.w.max)
    stop("prior-weight variable 'w' has too many columns")
  if (ncol(y) > ncol.y.max)
    stop("response variable 'y' has too many columns; ",
         "only ", ncol.y.max, " allowed")

  if (ncol(w) < ncol.w.min)
    stop("prior-weight variable 'w' has too few columns")
  if (ncol(y) < ncol.y.min)
    stop("response variable 'y' has too few columns; ",
         "at least ", ncol.y.max, " needed")

  if (min(w) <= 0)
    stop("prior-weight variable 'w' must contain positive values only")

  if (is.numeric(colsyperw) && ncol(y) %% colsyperw != 0)
    stop("number of columns of the response variable 'y' is not ",
         "a multiple of ", colsyperw)


  if (maximize) {
    Ncol.max.w <- max(ncol(w), ncol(y) / colsyperw)
    Ncol.max.y <- max(ncol(y), ncol(w) * colsyperw)
  } else {
    Ncol.max.w <- ncol(w)
    Ncol.max.y <- ncol(y)
  }

  if (out.wy && ncol(w) < Ncol.max.w) {
    nblanks <- sum(cn.w == "")
    if (nblanks > 0)
      cn.w[cn.w == ""] <- paste(prefix.w, 1:nblanks, sep = "")
    if (length(cn.w) < Ncol.max.w)
      cn.w <- c(cn.w, paste(prefix.w, (length(cn.w)+1):Ncol.max.w,
                            sep = ""))
    w <- matrix(w, n.lm, Ncol.max.w, dimnames = list(rn.w, cn.w))
  }
  if (out.wy && ncol(y) < Ncol.max.y) {
    nblanks <- sum(cn.y == "")
    if (nblanks > 0)
      cn.y[cn.y == ""] <- paste(prefix.y, 1:nblanks, sep = "")
    if (length(cn.y) < Ncol.max.y)
      cn.y <- c(cn.y, paste(prefix.y, (length(cn.y)+1):Ncol.max.y,
                            sep = ""))
    y <- matrix(y, n.lm, Ncol.max.y, dimnames = list(rn.y, cn.y))
  }
       
  list(w = if (out.wy) w else NULL,
       y = if (out.wy) y else NULL)
}





arwz2wz <- function(arwz, M = 1, M1 = 1) {



  if (length(dim.arwz <- dim(arwz)) != 3)
    stop("dimension of 'arwz' should be of length 3")
  n       <- dim.arwz[1]
  ndepy   <- dim.arwz[2]
  dim.val <- dim.arwz[3]

  if (ndepy == 1) {
    dim(arwz) <- c(n, dim.val)
    return(arwz)
  }

  wz <- matrix(0.0, nrow = n, ncol = sum(M:(M-M1+1)))
  ind1 <- iam(NA, NA, M = M1, both = TRUE, diag = TRUE)
  len.ind1 <- dim.val # length(ind1$col.index)

  for (ii in 1:ndepy) {
    for (jlocal in 1:len.ind1) {
      wz[, iam(M1 * (ii - 1) + ind1$row[jlocal],
               M1 * (ii - 1) + ind1$col[jlocal],
               M = M)] <- arwz[, ii, jlocal]
    }
  }

  colind <- ncol(wz)
  while (all(wz[, colind] == 0))
    colind <- colind - 1

  if (colind < ncol(wz))
    wz <- wz[, 1:colind, drop = FALSE]

  wz
}








param.names <- function(string, S) {
  if (S == 1) string else paste(string, 1:S, sep = "")
}





vweighted.mean.default <- function (x, w, ..., na.rm = FALSE) {
  temp5 <- w.y.check(w = w, y = x, ncol.w.max = Inf, ncol.y.max = Inf,
                     out.wy = TRUE,
                     colsyperw = 1,
                     maximize = TRUE,
                     Is.integer.y = FALSE,
                     Is.positive.y = FALSE,
                     Is.nonnegative.y = FALSE,
                     prefix.w = "PriorWeight",
                     prefix.y = "Response")

  x <- temp5$y
  w <- temp5$w

  ans <- numeric(ncol(w))
  for (ii in 1:ncol(w))
    ans[ii] <- weighted.mean(x[, ii], w = w[, ii], ..., na.rm = na.rm)
  ans
}






familyname.vlm <- function(object, all = FALSE, ...) {
  ans <- object@family@vfamily
  if (all) ans else ans[1]
}


familyname.vglmff <- function(object, all = FALSE, ...) {
  ans <- object@vfamily
  if (all) ans else ans[1]
}



if (!isGeneric("familyname"))
    setGeneric("familyname",
  function(object, ...) standardGeneric("familyname"))


setMethod("familyname", "vglmff",
         function(object, ...)
         familyname.vglmff(object, ...))



setMethod("familyname", "vlm",
         function(object, ...)
         familyname.vlm(object, ...))









