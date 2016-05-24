BinAddHaz <- function (formula, data, subset, weights, na.action,
                       model = TRUE, contrasts = NULL,
                       start = FALSE, start.val, attrib = FALSE, attrib.var,
                       type.attrib = "abs", set.seed = FALSE, seed, bootstrap = FALSE,
                       nbootstrap, parallel = FALSE, type.parallel = "snow", ncpus = 4,...){
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"),
             names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- quote(stats::model.frame)
  mf <- eval.parent(mf)
  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  #wt <- model.weights(mf)

  if(is.null(model.weights(mf))){
    w <- 1} else {
      w <- model.weights(mf)
    }

  if(mean(w) != 1){
    w <- w/mean(w)}

  offset <- model.offset(mf)
  x <- model.matrix(mt, mf, contrasts)

  ### LL function
  BinAddHazLogLik <- function(start.val, data, y, x, wgt){

    beta <- start.val[1: length(start.val)]

    eta_i <- as.vector(x %*% beta)
    pi_i <- 1 - exp(-eta_i)

    dev.resid <- function(y, pi_i, wgt){
      2 * wgt * (y * log(ifelse(y == 0, 1, y/pi_i)) + (1 - y) *
                   log(ifelse(y == 1, 1, (1 - y)/(1 - pi_i))))}

    LL <- sum(dev.resid(y, pi_i, wgt))
    return(LL)}

  ### Initial values
  if(start == FALSE){
  init.val <- glm(formula, family = poisson, data = data)
  start.val <- abs(init.val$coefficients)} else{
    start.val <- start.val
  }

  ### Optimization
  BinLL <- constrOptim(theta = start.val, data = data, f = BinAddHazLogLik, ui = x,
                       ci = rep(0, nrow(x)), method = "Nelder-Mead",  y = y, x = x,
                       wgt = w)
  Coeff <- BinLL$par

  ### Attribution
  if(attrib){
    if(missing(attrib.var)){
      att.var <- (x[, "(Intercept)"])
      x.mat <- x[,-att.var]} else {
        att.var <- data[,deparse(substitute(attrib.var))]
        x.mat <- x[,-c(1:length(unique(att.var)))]}

    att.var.levels <- sort(unique(att.var))
    attrib.length <- length(unique(att.var))

    if(!is.factor(att.var)){
      att.var. <- factor(att.var, levels = att.var.levels,
                         labels = att.var.levels)}
    if(is.factor(att.var)){
      att.var. <- att.var
    }

    if(attrib.length == 1){
      attrib.coef = Coeff[1]
      bin.coef = Coeff[-1]} else {
        attrib.coef = Coeff[1: attrib.length]
        bin.coef = Coeff[-c(1: attrib.length)]
      }

    AttribBin <- function(bin.coef, attrib.coef, data, y, x.mat, att.var.,
                           wgt, type.attrib){
      haz.dis <- t(bin.coef * t(x.mat))
      haz.back <- attrib.coef[att.var.]
      eta_i <- haz.back + apply(haz.dis, 1, sum)
      pi_i <- 1 - exp(-eta_i)
      att.x.mat <- (haz.dis/eta_i) * pi_i
      att.back <- (haz.back/eta_i) * pi_i

      # Attribution matrix
      att.mat <- matrix(NA, nrow = (ncol(x.mat) + 3), ncol = nlevels(att.var.))

      att.mat[1,] <-  tapply(wgt, att.var., sum)
      att.mat[2,] <-  tapply(pi_i * wgt, att.var., sum)
      att.mat[3,] <-  tapply(att.back * wgt, att.var., sum)

      for (i in 1: ncol(x.mat)){
        att.mat[3 + i,] <- tapply(att.x.mat[, i] * wgt, att.var., sum)
      }

      dimnames(att.mat) <- list(c("nn","disab","backgrnd", colnames(x.mat)),
                                levels(att.var.))

      attribution <- list()
      attribution2 <- list()
      att.rel <- list()
      att.abs <- list ()
      att.final <- list()

      for(i in 1: ncol(att.mat)){
        attribution[[i]] <- att.mat[, i]
        row_sub <- attribution[[i]] != 0
        attribution2[[i]] <- attribution[[i]][row_sub]

        if(type.attrib =="rel"){
          att.final[[i]] <- attribution2[[i]][3: length(attribution2[[i]])]/
            attribution2[[i]][2]
        }
        if(type.attrib =="abs"){
          att.final[[i]] <- attribution2[[i]][2: length(attribution2[[i]])]/
            attribution2[[i]][1]}

        if(type.attrib =="both"){
          att.rel[[i]] <- attribution2[[i]][3: length(attribution2[[i]])]/
            attribution2[[i]][2]
          att.abs[[i]] <- attribution2[[i]][2: length(attribution2[[i]])]/
            attribution2[[i]][1]
          att.final <- list(att.rel = att.rel, att.abs = att.abs)
        }
      }
      if(type.attrib == "rel" | type.attrib == "abs"){
        Output <- matrix(unlist(att.final), ncol = 1,
                         dimnames = list(names(unlist(att.final)), "Contribution"))}

      if(type.attrib == "both"){
        Out <- list()
        Out1 <- list()
        for(i in 1:length(att.final)){
          Out[[i]] <- unlist(att.final[[i]])
          Out1[[i]] <- matrix(Out[[i]], ncol = 1)
          colnames(Out1[[i]]) <-  "Contribution"
        }
        names(Out1) <- names(att.final)

        if(attrib.length == 1){
          rownames(Out1[[1]]) <- names(att.final[[1]][[1]])
          rownames(Out1[[2]]) <- names(att.final[[2]][[1]])

        } else {
          rownames(Out1[[1]]) <- paste0(names(Out[[1]]), rep(levels(att.var.),
                                        each = length(Out[[1]])/attrib.length))
          rownames(Out1[[2]]) <- paste0(names(Out[[2]]), rep(levels(att.var.),
                                        each = length(Out[[2]])/attrib.length))}
        Output <- Out1
      }
      return(Output)
    }

    Attribution <- AttribBin(bin.coef = bin.coef, attrib.coef = attrib.coef,
                              data = data, y = y, x.mat = x.mat,
                              att.var. = att.var., wgt = w,
                              type.attrib = type.attrib)
  }

  ### CI
  if(bootstrap == FALSE){
  CovAddHaz <- function(start.val, y, x){
    y.vec <- y
    x <- x
    beta <- start.val
    eta_i <- as.vector(x %*% beta)
    pi_i <- 1 - exp(-eta_i)
    W <- diag(drop((1 - pi_i)/pi_i))
    V <- ginv(t(x) %*% W %*% x)
    return(V)}

    vcov <- CovAddHaz(start.val = BinLL$par, x = x, y = y)
    stdError <- sqrt(diag(vcov))

    CILow <- BinLL$par - (1.96 * stdError)
    CIHigh <- BinLL$par + (1.96 * stdError)
    pvalue <- round(2 * pnorm(-abs(BinLL$par/stdError)), 4)
    CI <- cbind(CILow, CIHigh)
    colnames(CI) <- c("CI2.5", "CI97.5")
    colnames(vcov) <- names(BinLL$par)
    rownames(vcov) <- names(BinLL$par)

    if(attrib){
      Results <- list(coefficients = BinLL$par, ci = CI, resDeviance = BinLL$value,
                      df = nrow(data) - length(BinLL$par), pvalue = pvalue,
                      stdError = stdError, vcov = vcov, contribution = Attribution)
    } else {
      Results <- list(coefficients = BinLL$par, ci = CI, resDeviance = BinLL$value,
                      df = nrow(data) - length(BinLL$par), pvalue = pvalue,
                      stdError = stdError, vcov = vcov, contribution = attrib)}
    }

  ### bootstrap CI
  if(bootstrap){
  if(attrib){
  mystat <- function(data, indices) {
  m1 <-  constrOptim(theta = BinLL$par, data = data[indices,], f = BinAddHazLogLik,
                    ui = x[indices,], ci = 0, method = "Nelder-Mead", y = y[indices],
                    x = x[indices,], wgt = w[indices])

  m2 <- AttribBin(bin.coef = bin.coef, attrib.coef = attrib.coef,
                   data = data[indices,], y = y[indices], x.mat = x.mat[indices,],
                   att.var. = att.var.[indices], wgt = w[indices],
                   type.attrib = type.attrib)
  if(is.list(m2)){
    m3 <- unlist(m2)} else {
      m3 <- m2}
  res <- c(m1$par, m3)
  return(res)
  }

} else {
  mystat <- function(data, indices) {
    m <-  constrOptim(theta = BinLL$par, data = data[indices,], f = BinAddHazLogLik,
                      ui = x[indices,], ci = 0, method = "Nelder-Mead", y = y[indices],
                      x = x[indices,], wgt = w[indices])
    return(m$par)
  }
}

  if(set.seed){
    set.seed(seed)}

  if (parallel == TRUE){
  BootResult <- boot(data = data, statistic = mystat, R = nbootstrap, parallel = type.parallel,
                     ncpus = ncpus)} else {
                       BootResult <- boot(data = data, statistic = mystat, R = nbootstrap)
                     }

  BootCI <- matrix(NA, ncol = 2, nrow = length(BootResult$t0))

  for (i in 1:length(BootResult$t0)){
  BootCI[i,] <- boot.ci(BootResult,  conf = 0.95, type = "perc", index = i)[[4]][, 4:5]}
  colnames(BootCI) <- c("CILow", "CIHigh")

  if(attrib){
    CILow.att <- matrix(BootCI[(length(start.val) + 1): length(BootResult$t0), 1], ncol = 1,
                      dimnames = list(rownames(Attribution), "CILow"))
    CIHigh.att <- matrix(BootCI[(length(start.val) + 1): length(BootResult$t0), 2], ncol = 1,
                       dimnames = list(rownames(Attribution), "CIHigh"))

  if(type.attrib == "both"){
    att1 <- list(att.rel = cbind(Attribution[[1]], CILow.att[1: nrow(Attribution[[1]]),],
                                 CIHigh.att[1: nrow(Attribution[[1]]),]),
                 att.abs = cbind(Attribution[[2]],
                                 CILow.att[(nrow(Attribution[[1]]) + 1):
                                             (nrow(Attribution[[1]]) + nrow(Attribution[[2]])),],
                                 CIHigh.att[(nrow(Attribution[[1]]) + 1):
                                              (nrow(Attribution[[1]]) + nrow(Attribution[[2]])),]))
    colnames(att1[[1]]) <- c("Contribution", "CI2.5", "CI97.5")
    colnames(att1[[2]]) <- c("Contribution", "CI2.5", "CI97.5")
  } else {
    att1 <- cbind(Attribution, CILow.att, CIHigh.att)
  }
  Results <- list(coefficients = Coeff, ci = cbind(matrix(BootCI[1: length(start.val),
                                                                 "CILow"], ncol = 1),
                                                   matrix(BootCI[1:length(start.val),
                                                                 "CIHigh"], ncol = 1)),
                  resDeviance = BinLL$value, df = nrow(data) - length(BinLL$par),
                  contribution = att1, bootsRep = BootResult$t)
  colnames(Results$bootsRep) <- names(BootResult$t0)
  colnames(Results$ci) <- c("CI2.5", "CI97.5")
  rownames(Results$ci) <- names(Coeff)
  } else {

  Results <- list(coefficients = Coeff, ci = cbind(matrix(BootCI[1: length(start.val),
                                                                  "CILow"], ncol = 1),
                                                    matrix(BootCI[1:length(start.val),
                                                                  "CIHigh"], ncol = 1)),
                  resDeviance = BinLL$value, df = nrow(data) - length(BinLL$par),
                  contribution = attrib, bootsRep = BootResult$t)
  colnames(Results$bootsRep) <- names(BootResult$t0)
  colnames(Results$coefficients) <- c("CI2.5", "CI97.5")
  rownames(Results$ci) <- names(Coeff)
  }
  }
  Results$bootstrap <- bootstrap
  Results$fitted.values <- as.vector(x %*% Results$coefficients)
  Results$residuals <- y - Results$fitted.values
  Results$call <- match.call()
  class(Results) <- "binaddhazmod"
  return(Results)
  }

# GENERIC FUNCTIONS IN S3 LANGUAGE (PRINT, SUMMARY, PREDICT)
binaddhazmod <- function(x, ...) UseMethod("binaddhazmod")

print.binaddhazmod <- function(x,...){
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
}

summary.binaddhazmod <- function(object, ...){
  if(object$bootstrap == FALSE){
  se <- sqrt(diag(object$vcov))
  tval <- coef(object)/se
  TAB <- cbind(Estimate = coef(object),
               StdErr = se,
               t.value = tval,
               p.value = 2*pt(-abs(tval), df=object$df))} else {
                 TAB <- cbind(Estimate = coef(object), CI = object$ci)
               }
  if(object$bootstrap == FALSE){
    res <- list(call = object$call, bootstrap = object$bootstrap,
                coefficients = capture.output(printCoefmat(TAB, P.values = TRUE,
                                                              has.Pvalue = TRUE)))}
  else {res <- list(coefficients = TAB, bootstrap = object$bootstrap)}
  class(res) <- "summary.binaddhazmod"
  res
}

predict.binaddhazmod <- function(object, newdata = NULL, ...){
  if(is.null(newdata))
    y <- fitted(object)
  else{
    if(!is.null(object$formula)){
      ## model has been fitted using formula interface
      x <- model.matrix(object$formula, newdata)
    }
    else{
      x <- newdata
    }
    y <- as.vector(x %*% coef(object))
  }
  y
}

