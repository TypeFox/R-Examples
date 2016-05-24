# Ichar -------------------------------------------------------------------

Ichar <- function (i, model, npar, ymean, yvar) {
  Ichar.res <- matrix(NA, npar, npar)
  ymean <- gsub("i", i, ymean)
  yvar <- gsub("i", i, yvar)
  inverse.yvar <- paste("1/(", yvar, ")", sep="")
  for (s in 1:npar) {
    for (t in 1:npar) {
      respect1 <- paste("b", s, sep="")
      respect2 <- paste("b", t, sep="")
      charac1 <- paste("D(expression(", ymean, "),", "'", respect1, "'", ")", sep="")
      deriv1 <- eval(parse(text=(charac1)))
      deriv1 <- as.character(as.expression(deriv1))
      charac2 <- paste("D(expression(", ymean, "),", "'", respect2, "'", ")", sep="")
      deriv2 <- eval(parse(text=(charac2)))
      deriv2 <- as.character(as.expression(deriv2))
      deriv3 <- paste(inverse.yvar, deriv1, deriv2, sep="*")
      Ichar.res[s,t] <- deriv3
    }
  }
  return(Ichar.res)
}

# Mchar --------------------------------------------------------------------

Mchar <- function (ndpoints, model , npar ,ymean, yvar, equally = F) {
  Mchar.res = ""
  for (i in 1:ndpoints) {    
    weight = ifelse(equally == F, paste("w", i, sep=""), 1/npar)
    Ichar.res <- paste(weight, Ichar(i = i,model = model, npar= npar,ymean = ymean,
                                     yvar = yvar), sep ="*")
    Mchar.res <- paste(Mchar.res, Ichar.res, sep="+")
  }
  return(Mchar.res)
}

# Cons.M.Des --------------------------------------------------------------

Cons.M.Des <- function (ndpoints, model, vpar, niv, npar, ymean, yvar, r = 5000, lb, ub,
                        equally, prec) {
  detM.Des <- M.Des <- M.Des.mpfr  <- NULL
  Mchar.res <- Mchar(ndpoints = ndpoints, model = model, npar = npar, ymean = ymean, 
                     yvar = yvar, equally = equally)
  for (j in npar:1) {
    b.par <- paste("b", j, sep="")
    Mchar.res <- gsub(b.par, vpar[j] ,Mchar.res)
  }
  Mchar.res.mpfr <- Mchar.res
  if (model == "All") {
    k <- (niv) * ndpoints
    for (i in ndpoints:1) {
      for (j in (niv:1)) {
        x <- paste("x", i, j, sep = "")
        q.des <- paste("q[", k, "]", sep = "")
        q.des.mpfr <- paste ("mpfr(q[", k,"]", ",precBits = ", prec, ")", sep = "")
        Mchar.res <- gsub(x, q.des, Mchar.res)
        Mchar.res.mpfr <- gsub(x, q.des.mpfr,   Mchar.res.mpfr)
        k <- k-1
      }
    }
    if (equally == T) { 
      charac1 <- paste("detM.Des<-function(q){Mat<-matrix(c(", 
                       paste(Mchar.res[1:((npar)*(npar))], sep = ',' ,collapse = ",") , "),", npar, 
                       ",", npar,");d = determinant(Mat, logarithm = TRUE);
                       if(d$sign <= 0) d = -Inf else d = -d$modulus;return(d)}",
                       sep = "") 
      eval(parse(text = charac1))
  }    
    if (equally == F) {
      k <- (niv) * ndpoints
      counter <- ndpoints
      for (i in (k + ndpoints):(k + 1) ) {
        w <- paste("w", counter, sep = "")
        q2 <- paste("q[", i, "]", sep = "")
        q2.mpfr <- paste ("mpfr(q[", i,"]", ",precBits = ", prec, ")", sep = "")
        Mchar.res <- gsub(w, q2, Mchar.res)
        Mchar.res.mpfr <- gsub(w, q2.mpfr, Mchar.res.mpfr)
        counter <- counter - 1
      } 
      pen <- ""
      for (ii in (k + ndpoints):(k + 1)) {
        weigh <- paste("q[", ii, "]", sep = "")
        pen <- paste(pen, weigh, sep = "+")
      }
      pen <- paste(r, "*(", pen, "-1)^2", sep = "")
      charac1 <- paste("detM.Des<-function(q){Mat=matrix(c(", paste(Mchar.res[1:((npar) * (npar))], sep = ',', collapse = ",")
                       , "),", npar, ",", npar, ");d = determinant(Mat, logarithm = TRUE);
                       if(d$sign <= 0) d = -Inf else d = -d$modulus", "; d = d +", pen, ";return(d)}", sep = "")
      eval(parse(text = charac1))
  }
    charac2 <- paste("M.Des<-function(q){Mat=matrix(c(", paste(Mchar.res[1:((npar) * (npar))], sep = ',', collapse = ",")
                     , "),", npar, ",", npar, ");d=Mat;return(d)}", sep="")
    eval(parse(text = charac2))
    charac3 <- paste("M.Des.mpfr<-function(q){Mat=new('mpfrMatrix', c(", paste(Mchar.res.mpfr[1:((npar) * (npar))], sep = ',', collapse = ",")
                     , "),Dim=c(", npar, "L,", npar, "L));d=Mat;return(d)}", sep="")
    eval(parse(text = charac3))
  } else
  { 
    if (model == "Weibull") {
      Mchar.res <- gsub("x11", lb, Mchar.res)
      Mchar.res <- gsub("x21", "q[1]", Mchar.res)
      Mchar.res <- gsub("x31", "q[2]", Mchar.res)
      Mchar.res <- gsub("x41", ub, Mchar.res)
    }
    if (model == "Richards") {
      Mchar.res <- gsub("x11", "q[1]", Mchar.res)
      Mchar.res <- gsub("x21", "q[2]", Mchar.res)
      Mchar.res <- gsub("x31", "q[3]", Mchar.res)
      Mchar.res <- gsub("x41", ub, Mchar.res)
    }
    if (model == "MM") { 
      Mchar.res <- gsub("x11", "q[1]", Mchar.res)
      Mchar.res <- gsub("x21", ub, Mchar.res)
    }
    if (model == "Poisson.dose") { 
      Mchar.res <- gsub("x11", lb, Mchar.res)
      Mchar.res <- gsub("x21", "q[1]", Mchar.res)
    }
    if (model == "Log.lin" || model == "Emax.dosResp" || model == "IQ3.theorem" || model == "Exp.dose")
    {
      Mchar.res <- gsub("x11", lb, Mchar.res)
      Mchar.res <- gsub("x21", "q[1]", Mchar.res)
      Mchar.res <- gsub("x31", ub, Mchar.res)
    }
    if (model == "IQ1.theorem") {
      Mchar.res <- gsub("x11", lb, Mchar.res)
      Mchar.res <- gsub("x21", "q[1]", Mchar.res)
      Mchar.res <- gsub("x31", "q[2]", Mchar.res)
    }
    if (model == "IQ2.theorem") {
      Mchar.res <- gsub("x11", "q[1]", Mchar.res)
      Mchar.res <- gsub("x21", "q[2]", Mchar.res)
      Mchar.res <- gsub("x31", ub, Mchar.res)
    }
    charac1 <- paste("detM.Des<-function(q){Mat<-matrix(c(", 
                     paste(Mchar.res[1:((npar)*(npar))], sep=',' ,collapse = ",") , "),", npar, 
                     ",", npar,");d = determinant(Mat, logarithm = TRUE);
                     if(d$sign <= 0) d = -Inf else d = -d$modulus;return(d)}",
                     sep="") 
    eval(parse(text = charac1))
    charac2 <- paste("M.Des<-function(q){Mat=matrix(c(", paste(Mchar.res[1:((npar) * (npar))], sep = ',', collapse = ",")
                     , "),", npar, ",", npar, ");d=Mat;return(d)}", sep =  "")
    eval(parse(text = charac2)) 
    }
  if (model == "All")
    return(list (M.Des, detM.Des, M.Des.mpfr)) else
      return(list (M.Des, detM.Des))
}

# Mchar.x -----------------------------------------------------------------

Mchar.x <- function (model, niv, vpar, npar, ymean, yvar) {
  Ichar2 <- Ichar(1, model = model, npar = npar, ymean = ymean, yvar = yvar)
  for (i in npar:1) {
    parcharac <- paste("b", i, sep = "")  
    Ichar2 <- gsub(parcharac, vpar[i], Ichar2)    
  }
  k <- niv
  x <- c()
  for (j in (niv:1)) {
    x <- paste("x", 1, j, sep = "")
    q.des <- paste("x", k, sep = "")
    Ichar2 <- gsub(x, q.des, Ichar2)
    k <- k-1
  }
  return(Ichar2)
}

# Cons.M.Par --------------------------------------------------------------

Cons.M.Par <- function (ndpoints, model, MDes, niv, npar, ymean, yvar, equally) {
  Mchar.res <- Mchar(ndpoints = ndpoints, model = model, npar = npar, ymean = ymean,
                     yvar = yvar, equally = equally)
  for (i in ndpoints:1) { 
    weight <- paste("w", i, sep="")
    Mchar.res <- gsub(weight, MDes[i,(niv+1)], Mchar.res)
    for (j in (niv):1) {
      x <- paste("x", i, j, sep="")
      Mchar.res <- gsub(x, MDes[i,j], Mchar.res)
    }
  }  
  for (k in npar:1) {
    par <- paste("b", k, sep="")
    b <- paste("b[", k, "]", sep="")
    Mchar.res <- gsub(par, b, Mchar.res)
  }
  charac2 <- paste("M.Par<-function(b){Mat=matrix(c(", paste(Mchar.res[1:((npar)*(npar))], sep=',', collapse=",")
                   , "),", npar, ",", npar, ");d=Mat;return(d)}", sep="")
  M.Par <- NULL
  eval(parse(text=charac2))
  return(M.Par)
}

# Cons.FreD.x -------------------------------------------------------------


Cons.FreD.x <- function (ndpoints, model, niv, MDes, vpar, npar, ymean, yvar, equally) {
  M.Par <- Cons.M.Par(ndpoints = ndpoints, model = model, MDes = MDes, niv = niv, npar = npar, ymean = ymean, yvar = yvar, equally = equally) 
  try(Minverse <- qr.solve(M.Par(b = vpar),tol = .Machine$double.eps), silent = T)
  singularity <- (!exists("Minverse", inherits = F))
  if (singularity == F) {
    mmatrix <- Mchar.x(model, niv, vpar = vpar, npar, ymean, yvar)
    prodMm <- matrix(NA, dim(Minverse)[1], dim(mmatrix)[2])
    for (k in 1:dim(Minverse)[1]) {
      for (j in 1:dim(mmatrix)[2]) {
        prodMm[k,j] <- paste(Minverse[k,], mmatrix[,j], sep = "*", collapse = "+")
      }
    }  
    charac1 <- paste(diag(prodMm), collapse = "+")
    X.vec <- sapply(c(1:niv), FUN = function (i) paste("x", i, sep = ""), simplify = T)
    X.vec <- paste (X.vec, collapse = ",")
    charac2 <- paste("FreD.x <-function(", X.vec, "){d=", charac1, "-", dim(Minverse)[1], ";return(d)}", sep = "")
    FreD.x <- NULL
    eval(parse(text = charac2))
  }
  if (singularity == T)
    ## this massage will not appear for the user
    FreD.x = "inverse of Fisher information matrix is singular"
  return(list(FreD.x, singularity))
}

# plot.FreD.x -------------------------------------------------------------

plot.FreD.x <-function (ndpoints, model, niv, MDes , vpar, npar, ymean, yvar, lb, ub,
                        equally, arg.curve.user) { 
  Cons.FreD.x.res <- Cons.FreD.x (ndpoints = ndpoints, model = model, niv = niv,
                                  MDes = MDes, vpar = vpar, npar = npar, ymean = ymean,
                                  yvar = yvar, equally = equally) 
  FreD.x <- Cons.FreD.x.res[[1]]
  
  singularity <- Cons.FreD.x.res[[2]]    
  if(singularity == F) { 
    if (niv == 1) {
      arg.curve.default <- list(col = "blue", n = (ub-lb)*1000, main = "",
                                xlab = "design interval", ylab = "derivative function"
                                , from = lb, to = ub)
      arg.curve.user$expr <-  substitute(FreD.x,  env = parent.frame())
      arg.curve <- modifyList (arg.curve.default, arg.curve.user)
      do.call("curve", arg.curve)
      abline(h = 0, v = c(MDes[,1]) ,col = "grey", lty = 3)
      points(x = MDes[,1],  y = rep(0, dim(MDes)[1]), col = "red" ,pch = 16)
    }
  }
  return(singularity)
}

# LocallyD ----------------------------------------------------------------

LocallyD <- function (ndpoints, niv, npar, model, ymean, yvar, lb, ub, vpar, 
                      method="gosolnp", equally, n.restarts, n.sim, tol, prec,
                      rseed, efficiency, user.points, user.weights, arg.curve.user) {
  detM.Des <- Cons.M.Des(ndpoints = ndpoints, model = model, vpar = vpar, niv = niv,
                         npar = npar, ymean = ymean, yvar = yvar,
                         r=5000, lb = lb, ub = ub, equally = equally, prec = prec)[[2]]
  lb.design <- c()
  ub.design <- c()
  if (model == "All") {
    if (equally == F) {
      lb.design <- c(rep(lb, ndpoints), rep(0, ndpoints))
      ub.design <- c(rep(ub, ndpoints), rep(1, ndpoints))      
    }
    if (equally == T) {
      lb.design <- rep(lb, ndpoints) 
      ub.design <- rep(ub, ndpoints)
    }
  }else{ 
    if (model == "Weibull" || model == "Log.lin" || model == "Emax.dosResp" ||
          model == "IQ3.theorem" || model == "Exp.dose") {
      lb.design <- rep(lb, ndpoints-2) 
      ub.design <- rep(ub, ndpoints-2)
    }
    if (model == "Richards" || model == "MM" || model == "IQ1.theorem" ||
          model == "IQ2.theorem" || model == "Poisson.dose") {
      lb.design <- rep(lb, ndpoints-1) 
      ub.design <- rep(ub, ndpoints-1)
    }
  }
  if(method == "gosolnp"){
    gosolnp.res <- suppressWarnings(gosolnp(fun=detM.Des, LB = lb.design, UB = ub.design, n.restarts = n.restarts,
                                            n.sim = n.sim, control = list(tol = tol),
                                            rseed = rseed))
    param1 <- gosolnp.res$par
  }
  if(method == "multistart") {
    multistart.method <- multistart(func=detM.Des, lb=lb.design, ub=ub.design, niv, ndpoints,sum=T)
    min.indice <- which.min(multistart.method$value)
    param1 <-  multistart.method$wset[min.indice,]  
  }
  NA.par <- F
  if (sum(!is.na(param1)) != length(param1)){
    NA.par <- T
    warning("NA value for design points is detected. The calculation is out of R default precision capacity. changing the input arguments can be beneficial.")
  }
  if (model == "All") {
    if (equally == F)
      MDes <- matrix(param1, ndpoints, niv+1)
    if (equally == T) {
      MDes <- matrix(c(param1, rep(1/ndpoints,ndpoints)), ndpoints, niv+1)
      MDes <- MDes[order(MDes[,1]),]
    }
  }else 
  {
    if (model == "Weibull" || model == "Log.lin" || model == "Emax.dosResp" || 
          model == "IQ3.theorem" || model == "Exp.dose") {
      param1 <- sort(c(lb, param1, ub))
    }
    if (model == "Richards" || model == "MM" || model == "IQ2.theorem")
      param1 <- sort(c(param1, ub))
    if (model == "IQ1.theorem" || model == "Poisson.dose")
      param1 <- sort(c(lb, param1))
    MDes <- matrix(c(param1, rep(1/ndpoints,ndpoints)), ndpoints, niv+1)
  }
  if(niv == 1 &&  NA.par == F){
    singular<- plot.FreD.x (ndpoints = ndpoints, model = model, niv = niv, MDes = MDes, vpar = vpar
                            , npar = npar, ymean = ymean, yvar = yvar, lb = lb, ub = ub, equally = equally
                            , arg.curve.user = arg.curve.user)
    if(singular == T)
      warning("inverse of information matrix can not be calculated because of computational singularity at obtained design\n",
              call. = F)
  }
  result.LocallyD <- list(points = MDes[,-(niv+1)], weights =  MDes[, niv+1],
                          det.value = exp(-1 * gosolnp.res$values[length(gosolnp.res$values)]))
  if(efficiency == T){
    D.Efficiency.user <- D.Efficiency (opt.design = c(MDes[,1], MDes[,2])
                                       , user.weights = user.weights, user.points = user.points,
                                       ndpoints = ndpoints, model = model, 
                                       vpar = vpar, niv = niv, npar = npar, ymean = ymean,
                                       yvar = yvar, r = 5000, lb = lb, ub = ub, prec = prec)
    result.LocallyD  <- c(result.LocallyD[1:3], user.eff = round(as.numeric(D.Efficiency.user[1]), 5))
  } 
  return(result.LocallyD)
}

# multistart --------------------------------------------------------------

multistart <- function (func, lb, ub, niv, ndpoints, sum = F) {
  if(sum == T) {
    w.set <- matrix(,, niv * ndpoints + ndpoints)
    lb.design <- c()
    ub.design <- c()
    counter1 <- 0
    counter2 <- 0
    for (k1 in 1:ndpoints) {
      for (k2 in 1:niv) {
        counter1 <- counter1 + 1
        lb.design[counter1] <- lb[k2] 
        ub.design[counter1] <- ub[k2]   
      }
      counter2 <- counter2 + 1
      lb.design[ndpoints * niv + counter2] <- 0
      ub.design[ndpoints * niv + counter2] <- 1
    } 
  }
  if(sum == F) {
    w.set <- matrix(,,niv * ndpoints)
  }
  
  value.w.set <- c()
  step <- 0
  w <- 0
  stop <- F
  while(stop == F) {
    check <- 1
    step <- step + 1
    if(sum == T) {
      p0 <- sapply(X=c(1:length(lb.design)), FUN=function (i) runif(1,lb.design[i],ub.design[i]), simplify = T)
      nweigh <- (niv) * ndpoints
      weight <- c()
      for(i in (nweigh+ndpoints):(nweigh+1))
        weight[i-nweigh] <-p0[i]
      p0 <- sapply(X=c((nweigh+ndpoints):(nweigh+1)),FUN=function (i) p0[i]*(1/sum(weight)), simplify = T)
      LS <- solnp(p0, func, LB=lb.design, UB=ub.design)
    }
    if(sum==F) {
      p0 <- sapply(X=c(1:length(lb)), FUN=function (i) runif(1,lb[i],ub[i]), simplify = T)
      LS <- solnp(p0,func,LB=lb,UB=ub)
    }
    LSvalue <- (LS$value)[3]
    if(step == 1) {
      w.set[1,] <- LS$par
      value.w.set[1] <- LSvalue
      w <- w + 1
    }
    if(step > 1) {
      for (i in 1:dim(w.set)[1]) {
        if(sum(abs(LS$par - w.set[i,]) <= .01) == dim(w.set)[2])
          check <- 0
      }
      if(check == 1) {
        w <- w+1
        value.w.set[w] <- LSvalue
        w.set <- rbind(w.set,LS$par)
      }
      print(w.set)
    }  
    stop <- ((2 * w^2 + 3 * w + 2 < step) + (w * (w + 1) < .01 * (step-1) * step) == 2)    
  }
  return(list( wset= w.set, value = value.w.set)) 
}

# D.Efficiency ------------------------------------------------------------

D.Efficiency <- function (opt.design, user.weights, user.points, ndpoints, model, 
                          vpar, niv, npar, ymean, yvar, r = 5000, lb, ub, prec) {
  M.Des.mpfr <- Cons.M.Des(ndpoints = ndpoints, model = "All", vpar = vpar, niv = niv
                           , npar = npar, ymean = ymean, yvar = yvar, r = r, 
                           lb = lb, ub = ub, equally = F, prec = prec)[[3]]
  det.recur.opt <- det.recur(M.Des.mpfr(opt.design))
  M.Des.mpfr <- Cons.M.Des(ndpoints = length(user.points), model = "All",
                           vpar = vpar, niv = niv, npar, ymean = ymean, yvar = yvar, r=5000,
                           lb = lb, ub = ub, equally = F, prec = prec)[[3]]
  if ((sum(user.weights)-1) >= .Machine$double.eps) 
    user.weights <- user.weights * 1/(sum(user.weights))
  user.design <- c(user.points, user.weights) 
  
  det.recur.user <- det.recur(M.Des.mpfr(user.design))
  D.efficiency1 <- (det.recur.user/ det.recur.opt)^(1/npar)
  return(list <- c(D.efficiency1, det.recur.user))
}

# Determinant : recursive method --------------------------------------------------------

det.recur <- function (Mat) {
  if (dim (Mat)[1] == 2)
    det.recur1 <- Mat[1, 1] * Mat[2, 2] - Mat[2, 1] * Mat[1, 2]
  if (dim (Mat)[1] == 3)
    det.recur1 <-
    Mat[1, 1] * Mat[2, 2] * Mat[3, 3] +
    Mat[2, 1] * Mat[3, 2] * Mat[1, 3] +
    Mat[3, 1] * Mat[1, 2] * Mat[2, 3] -
    Mat[1, 1] * Mat[3, 2] * Mat[2, 3] -
    Mat[3, 1] * Mat[2, 2] * Mat[1, 3] -
    Mat[2, 1] * Mat[1, 2] * Mat[3, 3]
  if (dim (Mat)[1] == 4)
    det.recur1 <-
    Mat[1, 1] * Mat[2, 2] * Mat[3, 3] * Mat[4, 4] +
    Mat[1, 1] * Mat[2, 3] * Mat[3, 4] * Mat[4, 2] +
    Mat[1, 1] * Mat[2, 4] * Mat[3, 2] * Mat[4, 3] +
    Mat[1, 2] * Mat[2, 1] * Mat[3, 4] * Mat[4, 3] +
    Mat[1, 2] * Mat[2, 3] * Mat[3, 1] * Mat[4, 4] +
    Mat[1, 2] * Mat[2, 4] * Mat[3, 3] * Mat[4, 1] +
    Mat[1, 3] * Mat[2, 1] * Mat[3, 2] * Mat[4, 4] +
    Mat[1, 3] * Mat[2, 2] * Mat[3, 4] * Mat[4, 1] +
    Mat[1, 3] * Mat[2, 4] * Mat[3, 1] * Mat[4, 2] +
    Mat[1, 4] * Mat[2, 1] * Mat[3, 3] * Mat[4, 2] +
    Mat[1, 4] * Mat[2, 2] * Mat[3, 1] * Mat[4, 3] +
    Mat[1, 4] * Mat[2, 3] * Mat[3, 2] * Mat[4, 1] -
    Mat[1, 1] * Mat[2, 2] * Mat[3, 4] * Mat[4, 3] -
    Mat[1, 1] * Mat[2, 3] * Mat[3, 2] * Mat[4, 4] -
    Mat[1, 1] * Mat[2, 4] * Mat[3, 3] * Mat[4, 2] -
    Mat[1, 2] * Mat[2, 1] * Mat[3, 3] * Mat[4, 4] -
    Mat[1, 2] * Mat[2, 3] * Mat[3, 4] * Mat[4, 1] -
    Mat[1, 2] * Mat[2, 4] * Mat[3, 1] * Mat[4, 3] -
    Mat[1, 3] * Mat[2, 1] * Mat[3, 4] * Mat[4, 2] -
    Mat[1, 3] * Mat[2, 2] * Mat[3, 1] * Mat[4, 4] -
    Mat[1, 3] * Mat[2, 4] * Mat[3, 2] * Mat[4, 1] -
    Mat[1, 4] * Mat[2, 1] * Mat[3, 2] * Mat[4, 3] -
    Mat[1, 4] * Mat[2, 2] * Mat[3, 3] * Mat[4, 1] -
    Mat[1, 4] * Mat[2, 3] * Mat[3, 1] * Mat[4, 2] 
  return(det.recur1)
}

# Whole number ------------------------------------------------------------

is.naturalnumber <-
  function (x, tol = .Machine$double.eps^0.5)  (abs(x - round(x)) < tol & x >0)

# ldiq.theorem ------------------------------------------------------------

ldiq.theorem <- function (a, b, c, form, lb, ub) {
  s <- lb
  t <- ub
  if (form == 1) {
    gamma.var <- (sqrt (a * c))^{-1}
    delta.var <- (1/2) * (gamma.var + 1 + sqrt (gamma.var^2 + 6 * gamma.var +33 ))
    rho.var <- (delta.var + sqrt(delta.var^2 - 4) ) / (2)
  } 
  if (form == 2) {
    gamma.var <- b * (sqrt (a * c))^{-1}
    delta.var <- (1/2) * (gamma.var + 1 + sqrt (gamma.var^2 + 6 * gamma.var +33 ))
    rho.var <- (delta.var + sqrt(delta.var^2 - 4) ) / (2)
    
  }
  if (  s >= rho.var^{-1} * sqrt (b / c) & t > rho.var * sqrt (b / c) )
    model = "IQ1.theorem" ##lower
  if (  s < rho.var^{-1} * sqrt (b / c) & t <= rho.var * sqrt (b / c) )
    model = "IQ2.theorem" ##upper
  if (  s >= rho.var^{-1} * sqrt (b / c) & t <= rho.var * sqrt (b / c) )
    model = "IQ3.theorem" else model = "All" 
  return (model)  
}  



################################# Export functions

# ldiq --------------------------------------------------------------------

ldiq <- function (a, b, c, form, lb, ub, user.points = NULL, user.weights = NULL,  ...,
                  n.restarts = 1, n.sim = 1, tol = 1e-8, prec = 53, rseed = NULL) {
  arg.curve.user <- as.list(substitute(list(...)))
  if (!is.numeric(a) || length(a) != 1)
    stop("argument 'a' must be numeric of the length one\n")
  if (!is.numeric(b) || length(b) != 1)
    stop("argument 'b' must be numeric of the length one\n")
  if (!is.numeric(c) || length(c) != 1)
    stop("argument 'c' must be numeric of the length one\n")
  if (form != 1 & form != 2)
    stop("argument 'form' must be '1' or '2'\n")
  if (form == 1) {
    if (a <= 0 ||  b <= 0 || c <= 0)
      stop("arguments 'a', 'b' and 'c' must be greater than 0\n")
    if (2 * sqrt(b * c) <= 1)
      stop("value of '2*sqrt(b*c)' must be greater than 1\n")
  }    
  if (form == 2) {
    if (a <= 0 || c <= 0)
      stop("arguments 'a' and 'c' must be greater than 0\n")
    if (abs(b) >= 2 * sqrt(a * c))
      stop("value of '2*sqrt(a*c)' must be greater than '|b|'\n")
  }  
  if (!is.numeric(lb) || length(lb) != 1 )
    stop("argument 'lb' must be numeric of the length one\n")  
  if (lb < 0)
    stop("argument 'lb' must be greater than or equal to  0\n")
  if (!is.numeric(ub) || length(ub) != 1 )
    stop("argument 'ub' must be numeric of the length one\n") 
  if (ub <= lb)
    stop("'lower' > 'upper'\n")
  if (!is.null(user.points) & is.null(user.weights) || is.null(user.points) & !is.null(user.weights))
    stop("argument 'user.points' or 'user.weights' is missing\n" )
  if (!is.null(user.points) & !is.null(user.weights)){
    if (!is.numeric(user.points) || !is.numeric(user.weights) || length(user.weights) != length(user.points))
      stop("arguments 'user.points' and 'user.weights' must be numeric as the same length\n")
    if (sum(user.points >= 0) != length(user.points)  || sum(user.weights >= 0) != length(user.weights))
      stop("arguments 'user.points' and 'user.weights' must be positive\n")
    if (sum(user.points >= lb) != length(user.points)  || sum(user.points <= ub) != length(user.points))
      stop("elements of 'user.points' are out of design interval\n")
    efficiency = T
  }else
    efficiency = F
  if (form == 1)
    output <- LocallyD(ndpoints = 3, niv = 1, npar = 3,
                       model = ldiq.theorem (a = a, b = b, c = c, form = form, lb = lb, ub = ub),
                       ymean = "(b1*xi1)/(b2+xi1+b3*(xi1)^2)", yvar = 1, lb = lb, ub = ub,
                       vpar = c(a, b, c), method="gosolnp", equally = T, n.restarts = n.restarts,
                       n.sim = n.sim, tol = tol, prec = prec, rseed = rseed,
                       efficiency = efficiency,
                       user.points = user.points, user.weights = user.weights, arg.curve.user = arg.curve.user)
  if (form == 2)
    output <- LocallyD(ndpoints = 3, niv = 1, npar = 3,
                       model = ldiq.theorem (a = a, b = b, c = c, form = form, lb = lb, ub = ub),
                       ymean = "xi1/(b1+b2*xi1+b3*(xi1)^2)", yvar = 1, lb = lb, ub = ub,
                       vpar = c(a, b, c), method="gosolnp", equally = T, n.restarts = n.restarts,
                       n.sim = n.sim, tol = tol, prec = prec, rseed = rseed,
                       efficiency = efficiency,
                       user.points = user.points, user.weights = user.weights, arg.curve.user = arg.curve.user)
  return(output)  
}

# ldlogistic  -------------------------------------------------------------


ldlogistic <- function (a, b, form = 1 , lb, ub, user.points = NULL, user.weights = NULL,  ...,
                        n.restarts = 1, n.sim = 1, tol = 1e-8, prec = 53, rseed = NULL) {
  arg.curve.user <- as.list(substitute(list(...)))
  if (!is.numeric(a) || length(a) != 1)
    stop("argument 'a' must be numeric of the length one\n")
  if (!is.numeric(b) || length(b) != 1)
    stop("argument 'b' must be numeric of the length one\n")
  if (form != 1 & form != 2)
    stop("argument 'form' must be 1 or 2\n")
  if (!is.numeric(lb) || length(lb) != 1 )
    stop("argument 'lb' must be numeric of the length one\n") 
  if (!is.numeric(ub) || length(ub) != 1 )
    stop("argument 'ub' must be numeric of the length one\n") 
  if (ub <= lb)
    stop("'lower' > 'upper'\n")
  if (!is.null(user.points) & is.null(user.weights) || is.null(user.points) & !is.null(user.weights))
    stop("argument 'user.points' or 'user.weights' is missing\n" )
  if (!is.null(user.points) & !is.null(user.weights)){
    if (!is.numeric(user.points) || !is.numeric(user.weights) || length(user.weights) != length(user.points))
      stop("arguments 'user.points' and 'user.weights' must be numeric as the same length\n")
    if (sum(user.weights >= 0) != length(user.weights))
      stop("'user.weights' must be positive\n")
    if (sum(user.points >= lb) != length(user.points)  || sum(user.points <= ub) != length(user.points))
      stop("elements of 'user.points' are out of design interval\n")
    efficiency = T
  }else
    efficiency = F
  if (form == 1) 
    output <- LocallyD (ndpoints = 2, niv = 1, npar = 2, model = "All", ymean = "(1/(exp(-b1-b2*xi1)+1))"
                        , yvar = "(1/(exp(-b1-b2*xi1)+1))*(1-(1/(exp(-b1-b2*xi1)+1)))", lb = lb, ub = ub, 
                        vpar = c(a, b), method="gosolnp", equally = T,  n.restarts = n.restarts, n.sim = n.sim,
                        tol = tol, prec = prec, rseed = rseed,
                        efficiency = efficiency, user.points = user.points, user.weights, arg.curve.user = arg.curve.user)
  if (form == 2) 
    output <- LocallyD (ndpoints = 2, niv = 1, npar = 2, model = "All", ymean = "(1/(exp(-b2*(xi1-b1))+1))"
                        , yvar = "(1/(exp(-b2*(xi1-b1))+1))*(1-(1/(exp(-b2*(xi1-b1))+1)))", lb = lb, ub = ub, 
                        vpar = c(a, b), method="gosolnp", equally = T, n.restarts = n.restarts, n.sim = n.sim,
                        tol = tol, prec = prec, rseed = rseed,
                        efficiency = efficiency, user.points = user.points, user.weights = user.weights, arg.curve.user = arg.curve.user)
  return(output)
}

# ldpoisson ---------------------------------------------------------------


ldpoisson <-  function (a, b, form = 1, lb, ub, user.points = NULL, user.weights = NULL, ...,
                        n.restarts = 1, n.sim = 1, tol = 1e-8, prec = 53, rseed = NULL) {
  arg.curve.user <- as.list(substitute(list(...)))
  if (!is.numeric(a) || length(a) != 1)
    stop("argument 'a' must be numeric of the length one\n")
  if (!is.numeric(b) || length(b) != 1)
    stop("argument 'b' must be numeric of the length one\n")
  if (form != 1 & form != 2)
    stop("argument 'form' must be '1' or '2'\n")
  if (!is.numeric(lb) || length(lb) != 1 )
    stop("argument 'lb' must be numeric of the length one\n")  
  if (!is.numeric(ub) || length(ub) != 1 )
    stop("argument 'ub' must be numeric of the length one\n") 
  if (ub <= lb)
    stop("'lower' > 'upper'\n")
  if (!is.null(user.points) & is.null(user.weights) || is.null(user.points) & !is.null(user.weights))
    stop("argument 'user.points' or 'user.weights' is missing\n" )
  if (!is.null(user.points) & !is.null(user.weights)){
    if (!is.numeric(user.points) || !is.numeric(user.weights) || length(user.weights) != length(user.points))
      stop("arguments 'user.points' and 'user.weights' must be numeric as the same length\n")
    if (sum(user.weights >= 0) != length(user.weights))
      stop("'user.weights' must be positive\n")
    if (sum(user.points >= lb) != length(user.points)  || sum(user.points <= ub) != length(user.points))
      stop("elements of 'user.points' are out of design interval\n")
    efficiency = T
  }else
    efficiency = F
  if (form == 1)
    output <- LocallyD (ndpoints = 2, niv = 1, npar = 2, model = "All", ymean = "exp(b1+b2*xi1)"
                        , yvar = "exp(b1+b2*xi1)", lb = lb, ub = ub, 
                        vpar = c(a, b), method="gosolnp", equally = T, n.restarts = n.restarts, n.sim = n.sim,
                        tol = tol, prec = prec, rseed = rseed,
                        efficiency = efficiency, user.points = user.points, user.weights = user.weights, arg.curve.user = arg.curve.user)
  if (form == 2)
    output <- LocallyD (ndpoints = 2, niv = 1, npar = 2, model = "Poisson.dose", ymean = "b1*exp(-b2*xi1)"
                        , yvar = "b1*exp(-b2*xi1)", lb = lb, ub = ub, 
                        vpar = c(a, b), method="gosolnp", equally = T, n.restarts = n.restarts, n.sim = n.sim,
                        tol = tol, prec = prec, rseed = rseed,
                        efficiency = efficiency, user.points = user.points, user.weights = user.weights, arg.curve.user = arg.curve.user)
  return(output)
}

# ldweibull ---------------------------------------------------------------


ldweibull <- function (a, b, lambda, h, lb, ub, user.points = NULL,
                       user.weights = NULL,  ...,  n.restarts = 1, n.sim = 1,
                       tol = 1e-8, prec = 53, rseed = NULL) {
  arg.curve.user <- as.list(substitute(list(...)))
  if (!is.numeric(lambda) || length(lambda) != 1 )
    stop("'lambda' must be numeric of the length one\n")
  if (!is.numeric(lb) || length(lb) != 1 )
    stop("argument 'lb' must be numeric of the length one\n")  
  if (lb <= 0)
    stop("argument 'lb' must be greater than 0\n")
  if (!is.numeric(ub) || length(ub) != 1 )
    stop("argument 'ub' must be numeric of the length one\n") 
  if (ub <= lb)
    stop("'lower' > 'upper'\n")
  if (!is.numeric(a) || length(a) != 1)
    stop("argument 'a' must be numeric of the length one\n")
  if (!is.numeric(b) || length(b) != 1)
    stop("argument 'b' must be numeric of the length one\n")
  if (!is.numeric(h) || length(h) != 1)
    stop("argument 'h' must be numeric of the length one\n")
  if (!is.null(user.points) & is.null(user.weights) || is.null(user.points) & !is.null(user.weights))
    stop("argument 'user.points' or 'user.weights' is missing\n" )
  if (!is.null(user.points) & !is.null(user.weights)){
    if (!is.numeric(user.points) || !is.numeric(user.weights) || length(user.weights) != length(user.points))
      stop("arguments 'user.points' and 'user.weights' must be numeric as the same length\n")
    if (sum(user.points >= 0) != length(user.points)  || sum(user.weights >= 0) != length(user.weights))
      stop("arguments 'user.points' and 'user.weights' must be positive\n")
    if (sum(user.points >= lb) != length(user.points)  || sum(user.points <= ub) != length(user.points))
      stop("elements of 'user.points' are out of design interval\n")
    efficiency = T
  }else
    efficiency = F  
  output <- LocallyD(ndpoints = 4, niv = 1, npar = 4, model = "Weibull", ymean = "b1-b2*exp(-b3*xi1^b4)"
                     , yvar = 1, lb = lb, ub = ub, vpar = c(a, b, lambda, h),
                     method="gosolnp", equally = T, n.restarts = n.restarts, n.sim = n.sim,
                     tol = tol, prec = prec, rseed = rseed,
                     efficiency = efficiency, user.points = user.points, user.weights = user.weights, arg.curve.user = arg.curve.user)
  return(output)
}

# ldrichards --------------------------------------------------------------

ldrichards <- function (a, b, lambda, h, lb, ub, user.points = NULL, user.weights = NULL, ...,
                        n.restarts = 1, n.sim = 1, tol = 1e-8, prec = 53, rseed = NULL) {
  arg.curve.user <- as.list(substitute(list(...)))
  if (!is.numeric(lambda) || length(lambda) != 1 )
    stop("'lambda' must be numeric of the length one\n")
  if (!is.numeric(b) || length(b) != 1)
    stop("argument 'b' must be numeric of the length one\n")
  if (!is.numeric(h) || length(h) != 1)
    stop("argument 'h' must be numeric of the length one\n")
  if (!is.numeric(lb) || length(lb) != 1 )
    stop("argument 'lb' must be numeric of the length one\n")  
  if (lb < 0)
    stop("argument 'lb' must be greater than or equal to  0\n")
  if (!is.numeric(ub) || length(ub) != 1 )
    stop("argument 'ub' must be numeric of the length one\n") 
  if (ub <= lb)
    stop("'lower' > 'upper'\n")
  if (!is.numeric(a) || length(a) != 1)
    stop("argument 'a' must be numeric of the length one\n")
  if (!is.null(user.points) & is.null(user.weights) || is.null(user.points) & !is.null(user.weights))
    stop("argument 'user.points' or 'user.weights' is missing\n" )
  if (!is.null(user.points) & !is.null(user.weights)){
    if (!is.numeric(user.points) || !is.numeric(user.weights) || length(user.weights) != length(user.points))
      stop("arguments 'user.points' and 'user.weights' must be numeric as the same length\n")
    if (sum(user.points >= 0) != length(user.points)  || sum(user.weights >= 0) != length(user.weights))
      stop("arguments 'user.points' and 'user.weights' must be positive\n")
    if (sum(user.points >= lb) != length(user.points)  || sum(user.points <= ub) != length(user.points))
      stop("elements of 'user.points' are out of design interval\n")
    efficiency = T
  }else
    efficiency = F 
  output <- LocallyD(ndpoints = 4, niv = 1, npar = 4, model = "Richards",
                     ymean = "b1/(1+b2*exp(-b3*xi1))^b4", yvar = 1, lb = lb, ub = ub, vpar = c(a, b, lambda, h),
                     method="gosolnp", equally = T, n.restarts = n.restarts,
                     n.sim = n.sim, tol = tol, prec = prec, rseed = rseed,
                     efficiency = efficiency, user.points = user.points, user.weights = user.weights, arg.curve.user = arg.curve.user)
  return(output)
}

# ldmm --------------------------------------------------------

ldmm <- function (a, b, form = 1, lb, ub, user.points = NULL, user.weights = NULL,  ...,
                  n.restarts = 1, n.sim = 1, tol = 1e-8, prec = 53, rseed = NULL) {
  arg.curve.user <- as.list(substitute(list(...)))
  if (!is.numeric(b) || length(b) != 1)
    stop("argument 'b' must be numeric of the length one\n")
  if (form != 1 & form != 2 & form != 3)
    stop("argument 'form' must be '1' or '2' or '3'\n")
  if (!is.numeric(lb) || length(lb) != 1 )
    stop("argument 'lb' must be numeric of the length one\n")  
  if (lb < 0)
    stop("argument 'lb' must be greater than or equal to  0\n")
  if (!is.numeric(ub) || length(ub) != 1 )
    stop("argument 'ub' must be numeric of the length one\n") 
  if (ub <= lb)
    stop("'lower' > 'upper'\n")
  if (!is.numeric(a) || length(a) != 1)
    stop("argument 'a' must be numeric of the length one\n")
  if (!is.null(user.points) & is.null(user.weights) || is.null(user.points) & !is.null(user.weights))
    stop("argument 'user.points' or 'user.weights' is missing\n" )
  if (!is.null(user.points) & !is.null(user.weights)){
    if (!is.numeric(user.points) || !is.numeric(user.weights) || length(user.weights) != length(user.points))
      stop("arguments 'user.points' and 'user.weights' must be numeric as the same length\n")
    if (sum(user.points >= 0) != length(user.points)  || sum(user.weights >= 0) != length(user.weights))
      stop("arguments 'user.points' and 'user.weights' must be positive\n")
    if (sum(user.points >= lb) != length(user.points)  || sum(user.points <= ub) != length(user.points))
      stop("elements of 'user.points' are out of design interval\n")
    efficiency = T
  }else
    efficiency = F 
  if (form == 1) 
    output <- LocallyD(ndpoints = 2, niv = 1, npar = 2, model = "MM",
                       ymean = "(b1*xi1)/(1+b2*xi1)", yvar = 1, lb = lb, ub = ub, vpar = c(a, b),
                       method="gosolnp", equally = T, n.restarts = n.restarts,
                       n.sim = n.sim, tol = tol, prec = prec, rseed = rseed,
                       efficiency = efficiency, user.points = user.points, user.weights = user.weights, arg.curve.user = arg.curve.user)
  if (form == 2) 
    output <- LocallyD(ndpoints = 2, niv = 1, npar = 2, model = "MM",
                       ymean = "(b1*xi1)/(b2+xi1)", yvar = 1, lb = lb, ub = ub, vpar = c(a, b),
                       method="gosolnp", equally = T, n.restarts = n.restarts, n.sim = n.sim, tol = tol,
                       prec = prec, rseed = rseed, efficiency = efficiency,
                       user.points = user.points, user.weights = user.weights,
                       arg.curve.user = arg.curve.user)
  if (form == 3) 
    output <- LocallyD(ndpoints = 2, niv = 1, npar = 2, model = "MM",
                       ymean = "xi1/(b1+b2*xi1)", yvar = 1, lb = lb, ub = ub, vpar = c(a, b),
                       method="gosolnp", equally = T, n.restarts = n.restarts,
                       n.sim = n.sim, tol = tol, prec = prec, rseed = rseed,
                       efficiency = efficiency, user.points = user.points, user.weights = user.weights, arg.curve.user = arg.curve.user)
  return(output)
}

# ldloglin ---------------------------------------------------------------

ldloglin <- function (a, b, c, lb, ub, user.points = NULL, user.weights = NULL, ...,
                       n.restarts = 1, n.sim = 1, tol = 1e-8, prec = 53, rseed = NULL) {
  arg.curve.user <- as.list(substitute(list(...)))
  if (!is.numeric(a) || length(a) != 1)
    stop("argument 'a' must be numeric of the length one\n")
  if (a <= 0)
    stop("'a' must be greater than 0\n")
  if (!is.numeric(b) || length(b) != 1)
    stop("argument 'b' must be numeric of the length one\n")
  if (b <= 0)
    stop("'b' must be greater than 0\n")
  if (!is.numeric(c) || length(c) != 1)
    stop("argument 'c' must be numeric of the length one\n")
  if (c <= 0)
    stop("'c' must be greater than 0\n")
  if (!is.numeric(lb) || length(lb) != 1 )
    stop("argument 'lb' must be numeric of the length one\n")  
  if (lb < 0)
    stop("argument 'lb' must be greater than or equal to  0\n")
  if (!is.numeric(ub) || length(ub) != 1 )
    stop("argument 'ub' must be numeric of the length one\n") 
  if (ub <= lb)
    stop("'lower' > 'upper'\n")
  if (!is.null(user.points) & is.null(user.weights) || is.null(user.points) & !is.null(user.weights))
    stop("'user.points' or 'user.weights' is missing\n" )
  if (!is.null(user.points) & !is.null(user.weights)){
    if (!is.numeric(user.points) || !is.numeric(user.weights) || length(user.weights) != length(user.points))
      stop("arguments 'user.points' and 'user.weights' must be numeric as the same length\n")
    if (sum(user.points >= 0) != length(user.points)  || sum(user.weights >= 0) != length(user.weights))
      stop("arguments 'user.points' and 'user.weights' must be positive\n")
    if (sum(user.points >= lb) != length(user.points)  || sum(user.points <= ub) != length(user.points))
      stop("elements of 'user.points' are out of design interval\n")
    efficiency = T
  }else
    efficiency = F 
  output <- LocallyD(ndpoints = 3, niv = 1, npar = 3, model = "Log.lin",
                     ymean = "b1 + b2 * log(xi1 + b3)", yvar = 1, lb = lb, ub = ub, vpar = c(a, b, c),
                     method="gosolnp", equally = T, n.restarts = n.restarts,
                     n.sim = n.sim, tol = tol, prec = prec, rseed = rseed,
                     efficiency = efficiency, user.points = user.points, user.weights = user.weights, 
                     arg.curve.user = arg.curve.user)
  return(output)
}

# ldnbinom -------------------------------------------------------------

ldnbinom <- function (a, b, theta, lb, ub, user.points = NULL, user.weights = NULL,  ...,
                         n.restarts = 1, n.sim = 1, tol = 1e-8, prec = 53, rseed = NULL) {
  arg.curve.user <- as.list(substitute(list(...)))
  if (!is.numeric(a) || length(a) != 1)
    stop("argument 'a' must be numeric of the length one\n")
  if (!is.numeric(b) || length(b) != 1)
    stop("argument 'b' must be numeric of the length one\n")
  if (!is.numeric(theta) || !is.naturalnumber(theta))
    stop ("'theta' must be a Natural number\n")
  if (!is.numeric(lb) || length(lb) != 1 )
    stop("argument 'lb' must be numeric of the length one\n")  
  if (!is.numeric(ub) || length(ub) != 1 )
    stop("argument 'ub' must be numeric of the length one\n") 
  if (ub <= lb)
    stop("'lower' > 'upper'\n")
  if (!is.null(user.points) & is.null(user.weights) || is.null(user.points) & !is.null(user.weights))
    stop("'user.points' or 'user.weights' is missing\n" )
  if (!is.null(user.points) & !is.null(user.weights)){
    if (!is.numeric(user.points) || !is.numeric(user.weights) || length(user.weights) != length(user.points))
      stop("arguments 'user.points' and 'user.weights' must be numeric as the same length\n")
    if (sum(user.weights >= 0) != length(user.weights))
      stop("'user.weights' must be positive\n")
    if (sum(user.points >= lb) != length(user.points)  || sum(user.points <= ub) != length(user.points))
      stop("elements of 'user.points' are out of design interval\n")
    efficiency = T
  }else
    efficiency = F 
  output <- LocallyD (ndpoints = 2, niv = 1, npar = 2, model = "All", ymean = "b1 * exp(-b2 * xi1)"
                      , yvar = paste ("b1 * exp(-b2 * xi1) * (1 + (1/", theta, ") * b1 * exp(-b2 * xi1))" , sep = "")
                      , lb = lb, ub = ub, vpar = c(a, b), method="gosolnp", equally = T
                      , n.restarts = n.restarts, n.sim = n.sim,
                      tol = tol, prec = prec, rseed = rseed, efficiency = efficiency
                      , user.points = user.points, user.weights = user.weights, arg.curve.user = arg.curve.user)
  return(output)
}

# ldexpdose -----------------------------------------------------------

ldexpdose <- function (a, b, c, lb, ub, user.points = NULL, user.weights = NULL,  ...,
                        n.restarts = 1, n.sim = 1, tol = 1e-8, prec = 53, rseed = NULL) {
  arg.curve.user <- as.list(substitute(list(...)))
  if (!is.numeric(a) || length(a) != 1)
    stop("argument 'a' must be numeric of the length one\n")
  if (!is.numeric(b) || length(b) != 1)
    stop("argument 'b' must be numeric of the length one\n")
  if (!is.numeric(c) || length(c) != 1)
    stop("argument 'c' must be numeric of the length one\n")
  if (lb < 0)
    stop("argument 'lb' must be greater than or equal to  0\n")
  if (!is.numeric(ub) || length(ub) != 1 )
    stop("argument 'ub' must be numeric of the length one\n") 
  if (ub <= lb)
    stop("'lower' > 'upper'\n")
  if (!is.null(user.points) & is.null(user.weights) || is.null(user.points) & !is.null(user.weights))
    stop("'user.points' or 'user.weights' is missing\n" )
  if (!is.null(user.points) & !is.null(user.weights)){
    if (!is.numeric(user.points) || !is.numeric(user.weights) || length(user.weights) != length(user.points))
      stop("arguments 'user.points' and 'user.weights' must be numeric as the same length\n")
    if (sum(user.points >= 0) != length(user.points)  || sum(user.weights >= 0) != length(user.weights))
      stop("arguments 'user.points' and 'user.weights' must be positive\n")
    if (sum(user.points >= lb) != length(user.points)  || sum(user.points <= ub) != length(user.points))
      stop("elements of 'user.points' are out of design interval\n")
    efficiency = T
  }else
    efficiency = F 
  output <- LocallyD(ndpoints = 3, niv = 1, npar = 3, model = "Exp.dose",
                     ymean = "b1 + b2 * exp(xi1/b3)", yvar = 1, lb = lb, ub = ub, vpar = c(a, b, c),
                     method="gosolnp", equally = T, n.restarts = n.restarts, n.sim = n.sim, tol = tol,
                     prec = prec, rseed = rseed,
                     efficiency = efficiency, user.points = user.points,
                     user.weights = user.weights, arg.curve.user = arg.curve.user)
  return(output)
}

# ldemax ------------------------------------------------------------------

ldemax <- function (a, b, c, lb, ub, user.points = NULL, user.weights = NULL,  ...,
                    n.restarts = 1, n.sim = 1, tol = 1e-8, prec = 53, rseed = NULL) {
  arg.curve.user <- as.list(substitute(list(...)))
  if (!is.numeric(b) || length(b) != 1)
    stop("argument 'b' must be numeric of the length one\n")
  if (!is.numeric(c) || length(c) != 1)
    stop("argument 'h' must be numeric of the length one\n")
  if (!is.numeric(lb) || length(lb) != 1 )
    stop("argument 'lb' must be numeric of the length one\n")  
  if (!is.numeric(ub) || length(ub) != 1 )
    stop("argument 'ub' must be numeric of the length one\n") 
  if (ub <= lb)
    stop("'lower' > 'upper'\n")
  if (!is.numeric(a) || length(a) != 1)
    stop("argument 'a' must be numeric of the length one\n")
  if (!is.null(user.points) & is.null(user.weights) || is.null(user.points) & !is.null(user.weights))
    stop("'user.points' or 'user.weights' is missing\n" )
  if (!is.null(user.points) & !is.null(user.weights)){
    if (!is.numeric(user.points) || !is.numeric(user.weights) || length(user.weights) != length(user.points))
      stop("arguments 'user.points' and 'user.weights' must be numeric as the same length\n")
    if (sum(user.weights >= 0) != length(user.weights))
      stop("'user.weights' must be positive\n")
    if (sum(user.points >= lb) != length(user.points)  || sum(user.points <= ub) != length(user.points))
      stop("elements of 'user.points' are out of design interval\n")
    efficiency = T
  }else
    efficiency = F 
  output <- LocallyD (ndpoints = 3 , niv = 1, npar = 3, model = "Emax.dosResp",
                      ymean  = "b1 + (b2 * xi1)/(xi1 + b3)", yvar =  1, lb = lb, ub = ub, 
                      vpar = c(a, b, c), method = "gosolnp", equally = T, n.restarts = n.restarts,
                      n.sim = n.sim, tol = tol, prec = prec, rseed = rseed,
                      efficiency = efficiency, user.points = user.points,
                      user.weights = user.weights, arg.curve.user = arg.curve.user)
  return(output)
}

# cfisher -----------------------------------------------------------------

cfisher <- function (ymean, yvar, ndpoints, prec = 53){
  M.Des <- M.Des.mpfr <- NULL
  if (!is.character(ymean) || !is.character(yvar))
    stop("'ymean' and 'yvar' must be character strings\n")
  stop.while <-  F
  j <- 0
  niv <- 0
  while (stop.while == F){
    j <- j + 1
    x <- paste("x", j, sep = "")
    x.match <- gregexpr(x, ymean)[[1]]
    stop.while  <- (length(x.match) == 1 && x.match == -1)
    niv <- ifelse(stop.while == T, break, niv + 1)
  }
  if (niv == 0)
    stop("there are no characters as independent variables in 'ymean'\n")
  for (k in niv:1){
    xi <- paste("xi", k, sep = "")
    x <- paste("x", k, sep = "")
    ymean <- gsub(x, xi, ymean)
    yvar <- gsub(x, xi, yvar)
  }
  npar <- length(strsplit(ymean, "b")[[1]]) - 1
  if (!is.naturalnumber(ndpoints))
    stop("number of design points must be a natural number\n")
  Mchar.res <- Mchar(ndpoints = ndpoints, model = "All", npar = npar, ymean = ymean, 
                     yvar = yvar, equally = F)
  Mchar.res.mpfr <- Mchar.res
  bnumber <- 0 
  for (k in npar:1) {
    par <- paste("b", k, sep="")
    if (sum(grepl(par, Mchar.res )) != 0)
      bnumber <- bnumber + 1
    b <- paste("b[", k, "]", sep="")
    b.mpfr <- paste ("mpfr(b[", k,"]", ", precBits = ", prec, ")", sep = "")
    Mchar.res <- gsub(par, b, Mchar.res)
    Mchar.res.mpfr <- gsub(par, b.mpfr, Mchar.res.mpfr)
  }
  k <- (niv) * ndpoints
  for (i in ndpoints:1) {
    for (j in (niv:1)) {
      x <- paste("x", i, j, sep = "")
      x.vec <- paste("x[", k, "]", sep = "")
      x.vec.mpfr <- paste ("mpfr(x[", k,"]", ", precBits = ", prec, ")", sep = "")
      Mchar.res <- gsub(x, x.vec, Mchar.res)
      Mchar.res.mpfr <- gsub(x, x.vec.mpfr,   Mchar.res.mpfr)
      k <- k-1
    }
  }
  for (i in ndpoints:1) {
    w <- paste("w", i, sep = "")
    w.vec <- paste("w[", i, "]", sep = "")
    w.vec.mpfr <- paste ("mpfr(w[", i,"]", ", precBits = ", prec, ")", sep = "")
    Mchar.res <- gsub(w, w.vec, Mchar.res)
    Mchar.res.mpfr <- gsub(w, w.vec.mpfr, Mchar.res.mpfr)
  }
  xij <- c()
  counter <- 1
  
  for(i in 1:ndpoints){
    for(j in 1:niv){
      xij[counter]<- paste("x", i, j, sep = "")
      counter <- counter + 1
    }
  }
  x.argument <- paste(xij, collapse = ", ")
  x.argument <- paste( "x = c(",x.argument ,")", sep = "")
  wi <- c()
  for(i in 1:ndpoints)
    wi[i] <- paste("w", i, sep ="")
  w.argument <- paste(wi, collapse = ", ")
  w.argument <- paste("w = c(", w.argument, ")", sep= "")
  if( bnumber == 0){ 
    charac1 <- paste("M.Des <- function(", x.argument, ", ",  w.argument,") {Mat = matrix(c(", paste(Mchar.res[1:((npar) * (npar))], sep = ',', collapse = ",")
                     , "),", npar, ",", npar, "); d=Mat; return(d)}", sep =  "")
    eval(parse(text = charac1)) 
    charac2 <- paste("M.Des.mpfr <- function(", x.argument, ", ",  w.argument, ") {Mat = new('mpfrMatrix', c(", paste(Mchar.res.mpfr[1:((npar) * (npar))], sep = ',', collapse = ",")
                     , "),Dim=c(", npar, "L,", npar, "L));d = Mat; return(d)}", sep="")
    eval(parse(text = charac2))
  }else{
    bi <- c()
    for(i in 1:bnumber)
      bi[i] <- paste("b", i, sep ="")
    b.argument <- paste(bi, collapse = ", ")
    b.argument <- paste("b = c(", b.argument, ")", sep= "")
    charac1 <- paste("M.Des <- function(",x.argument, ", ", w.argument, ", ", b.argument, ") {Mat = matrix(c(", paste(Mchar.res[1:((npar) * (npar))], sep = ',', collapse = ",")
                     , "),", npar, ",", npar, "); d = Mat; return(d)}", sep =  "")
    eval(parse(text = charac1)) 
    charac2 <- paste("M.Des.mpfr <- function(",x.argument, ", ", w.argument, ", ", b.argument, ") {Mat = new('mpfrMatrix', c(", paste(Mchar.res.mpfr[1:((npar) * (npar))], sep = ',', collapse = ",")
                     , "),Dim=c(", npar, "L,", npar, "L));d = Mat; return(d)}", sep="")
    eval(parse(text = charac2))
  }
  return(list(fim = M.Des, fim.mpfr = M.Des.mpfr))
}

# eff ---------------------------------------------------------------------


eff <- function(ymean, yvar,  param, points1, points2, weights1, weights2, prec = 53){
  if (!is.character(ymean) || !is.character(yvar))
    stop("'ymean' and 'yvar' must be character strings\n")
  stop.while <-  F
  j <- 0
  niv <- 0
  while (stop.while == F){
    j <- j + 1
    x <- paste("x", j, sep = "")
    x.match <- gregexpr(x, ymean)[[1]]
    stop.while  <- (length(x.match) == 1 && x.match == -1)
    niv <- ifelse(stop.while == T, break, niv + 1)
  }
  if (niv == 0)
    stop("there are no characters as independent variables in 'ymean'\n")
  for (k in 1:niv){
    xi <- paste("xi", k, sep = "")
    x <- paste("x", k, sep = "")
    ymean <- gsub(x, xi, ymean)
    yvar <- gsub(x, xi, yvar)
  }
  npar <- length(strsplit(ymean, "b")[[1]]) - 1
  if (!is.numeric(param))
    stop("argument 'param' must be a numeric vector\n")
  if (length(param) > 4)
    stop("number of parameters can not be more than 4\n")  
  if (npar > length(param)){
    ## param is missing
    param.missing <- c()
    for(i in (length(param) + 1):npar)
      param.missing[i - length(param)] <- paste("b", i, sep = "")
    l.param.missing <- length(param.missing)
    param.missing <- paste(param.missing, collapse = ", ")
    if (l.param.missing == 1)
      stop("value of ", param.missing, " is missing in 'param'") else
        stop("values of ", param.missing, " are missing in 'param'")
  }
  if (npar < length(param)){
    b.missing <- c()
    for(i in (npar + 1):length(param))
      b.missing[i - npar]<- paste("b", i, sep = "")
    l.b.missing <- length(b.missing)
    b.missing <- paste(b.missing, collapse = ", ")
    if (l.b.missing == 1)
      stop("character ", b.missing, " is missing in 'ymean'") else
        stop("characters ", b.missing, " are missing in 'ymean'")
  }
  if (length(points1) %/% niv < npar || length(points2) %/% niv < npar)
    stop("to avoid exact singularity, length of 'points1' and 'points2' must be at least ", npar * niv)
  if (length(points1) %% niv != 0){
    points1.missing <- c()
    for (i in (length(points1) %% niv + 1):niv)
      points1.missing[i - length(points1) %% niv] <- paste ("x", length(points1) %/% niv + 1, i, sep = "")
    l.points1.missing <- length(points1.missing)
    points1.missing <- paste(points1.missing, collapse = ", ")
    if (l.points1.missing == 1)
      stop("value of ", points1.missing, " is missing in 'points1'") else
        stop("values of ", points1.missing, " are missing in 'points1'")
  }
  if (length(points2) %% niv != 0){
    points2.missing <- c()
    for (i in (length(points2) %% niv + 1):niv)
      points2.missing[i - length(points2) %% niv] <- paste ("x", length(points2) %/% niv + 1, i, sep = "")
    l.points2.missing <- length(points2.missing)
    points2.missing <- paste(points2.missing, collapse = ", ")
    if (l.points2.missing == 1)
      stop("value of ", points2.missing, " is missing in 'points2'") else
        stop("values of ", points2.missing, " are missing in 'points2'")
  }
  ndpoints1 <- length(points1) %/% niv
  if (length(weights1) < ndpoints1){
    weights1.missing <- c()
    for (i in (length(weights1) + 1):ndpoints1)
      weights1.missing[i - length(weights1)] <- paste("w", i, sep = "")
    l.weights1.missing <- length(weights1.missing)
    weights1.missing <- paste (weights1.missing, collapse = ", ")
    if (l.weights1.missing == 1)
      stop(weights1.missing, " is missing in 'weights1'") else
        stop(weights1.missing, " are missing in 'weights1'")    
  }
  if (length(weights1) > ndpoints1){
    weights1.unused <- c()
    for (i in (ndpoints1 + 1):length(weights1))
      weights1.unused[i - ndpoints1] <- paste("w", i, sep = "")
    l.weights1.unused <- length(weights1.unused)
    weights1.unused <- paste(weights1.unused, collapse = ", ")
    if (l.weights1.unused == 1)
      warning(weights1.unused,  " in weights1 is unused") else
        warning(weights1.unused,  " in weights1 are unused")
  }
  if ((sum(weights1)-1) >= .Machine$double.eps) 
    weights1 <- weights1 * 1/(sum(weights1))
  ndpoints2 <- length(points2) %/% niv
  if (length(weights2) < ndpoints2){
    weights2.missing <- c()
    for (i in (length(weights2) + 1):ndpoints2)
      weights2.missing[i - length(weights2)] <- paste("w", i, sep = "")
    l.weights2.missing <- length(weights2.missing)
    weights2.missing <- paste (weights2.missing, collapse = ", ")
    if (l.weights2.missing == 1)
      stop(weights2.missing, " is missing in 'weights2'") else
        stop(weights2.missing, " are missing in 'weights2'")    
  }
  if (length(weights2) > ndpoints2){
    weights2.unused <- c()
    for (i in (ndpoints2 + 1):length(weights2))
      weights2.unused[i - ndpoints2] <- paste("w", i, sep = "")
    l.weights2.unused <- length(weights2.unused)
    weights2.unused <- paste(weights2.unused, collapse = ", ")
    if (l.weights2.unused == 1)
      warning(weights2.unused,  " in weights2 is unused") else
        warning(weights2.unused,  " in weights2 are unused")
  }
  if ((sum(weights2)-1) >= .Machine$double.eps) 
    weights2 <- weights2 * 1/(sum(weights2))
  Mdesign1 <- Cons.M.Des(ndpoints = length(points1), model = "All", vpar = param, niv = niv
                         , npar = npar, ymean = ymean, yvar = yvar, r = 5000, 
                         lb = 0, ub = 1, equally = F, prec = prec)[[3]]
  Mdesign2 <- Cons.M.Des(ndpoints = length(points2), model = "All", vpar = param, niv = niv
                         , npar = npar, ymean = ymean, yvar = yvar, r = 5000, 
                         lb = 0, ub = 1, equally = F, prec = prec)[[3]]
  D.eff <- (det.recur(Mdesign1(c(points1, weights1))) / det.recur(Mdesign2(c(points2, weights2))))^(1/npar)
  return(D.eff)
}

# cfderiv ------------------------------------------------------------------

cfderiv <- function(ymean, yvar, param, points, weights){
  if (!is.character(ymean) || !is.character(yvar))
    stop("'ymean' and 'yvar' must be character strings\n")
  stop.while <-  F
  j <- 0
  niv <- 0
  while (stop.while == F){
    j <- j + 1
    x <- paste("x", j, sep = "")
    x.match <- gregexpr(x, ymean)[[1]]
    stop.while  <- (length(x.match) == 1 && x.match == -1)
    niv <- ifelse(stop.while == T, break, niv + 1)
  }
  if (niv == 0)
    stop("there are no characters as independent variables in 'ymean'\n")
  for (k in 1:niv){
    xi <- paste("xi", k, sep = "")
    x <- paste("x", k, sep = "")
    ymean <- gsub(x, xi, ymean)
    yvar <- gsub(x, xi, yvar)
  }
  npar <- length(strsplit(ymean, "b")[[1]]) - 1
  if (!is.numeric(param))
    stop("argument 'param' must be a numeric vector\n")
  if (npar > length(param)){
    param.missing <- c()
    for(i in (length(param) + 1):npar)
      param.missing[i - length(param)] <- paste("b", i, sep = "")
    l.param.missing <- length(param.missing)
    param.missing <- paste(param.missing, collapse = ", ")
    if (l.param.missing == 1)
      stop("value of ", param.missing, " is missing in 'param'") else
        stop("values of ", param.missing, " are missing in 'param'")
  }
  if (npar < length(param)){
    b.missing <- c()
    for(i in (npar + 1):length(param))
      b.missing[i - npar]<- paste("b", i, sep = "")
    l.b.missing <- length(b.missing)
    b.missing <- paste(b.missing, collapse = ", ")
    if (l.b.missing == 1)
      stop("character ", b.missing, " is missing in 'ymean'") else
        stop("characters ", b.missing, " are missing in 'ymean'")
  }
  if (length(points) %/% niv < npar)
    stop("to avoid exact singularity, length of 'points' must be at least ", npar * niv)
  if (length(points) %% niv != 0){
    points.missing <- c()
    for (i in (length(points) %% niv + 1):niv)
      points.missing[i - length(points) %% niv] <- paste ("x", length(points) %/% niv + 1, i, sep = "")
    l.points.missing <- length(points.missing)
    points.missing <- paste(points.missing, collapse = ", ")
    if (l.points.missing == 1)
      stop("value of ", points.missing, " is missing in 'points'") else
        stop("values of ", points.missing, " are missing in 'points'")
  }
  ndpoints <- length(points) %/% niv
  if (length(weights) < ndpoints){
    weights.missing <- c()
    for (i in (length(weights) + 1):ndpoints)
      weights.missing[i - length(weights)] <- paste("w", i, sep = "")
    l.weights.missing <- length(weights.missing)
    weights.missing <- paste (weights.missing, collapse = ", ")
    if (l.weights.missing == 1)
      stop(weights.missing, " is missing in 'weights'") else
        stop(weights.missing, " are missing in 'weights'")    
  }
  if (length(weights) > ndpoints){
    weights.unused <- c()
    for (i in (ndpoints + 1):length(weights))
      weights.unused[i - ndpoints] <- paste("w", i, sep = "")
    l.weights.unused <- length(weights.unused)
    weights.unused <- paste(weights.unused, collapse = ", ")
    if (l.weights.unused == 1)
      warning(weights.unused,  " in weights is unused") else
        warning(weights.unused,  " in weights are unused")
  }
  if ((sum(weights)-1) >= .Machine$double.eps) 
    weights <- weights * 1/(sum(weights))
  MDes <- matrix (c(points, weights), length(points), niv + 1)
  
  M.Par <- Cons.M.Par(ndpoints = length(points), model = "All", MDes = MDes,
                      niv = niv, npar = npar, ymean = ymean, yvar = yvar,
                      equally = F) 
  try(Minverse <- qr.solve(M.Par(b = param),tol = .Machine$double.eps), silent = T)
  if ((!exists("Minverse", inherits = F)))
    stop("Fisher information matrix at the given design is computationally singular\n")
  mmatrix <- Mchar.x(model = "All", niv = niv, vpar = param, npar = npar
                     , ymean = ymean, yvar = yvar)
  prodMm <- matrix(NA, dim(Minverse)[1], dim(mmatrix)[2])
  for (k in 1:dim(Minverse)[1]) {
    for (j in 1:dim(mmatrix)[2]) {
      prodMm[k,j] <- paste(Minverse[k,], mmatrix[,j], sep = "*", collapse = "+")
    }
  }  
  charac1 <- paste(diag(prodMm), collapse = "+")
  if (niv == 1)
    charac1 <- gsub("x1", "x", charac1) else{
      for (j in niv:1){
        x <- paste("x", j, sep = "")
        xvec <- paste("x[", j, "]", sep = "")
        charac1 <- gsub(x, xvec, charac1)
      }
    }
  charac2 <- paste("FreD <- function(x){d=", charac1, "-", dim(Minverse)[1], ";return(-d)}", sep = "")
  FreD <- NULL
  eval(parse(text = charac2))
  return(fderiv = FreD)
}

