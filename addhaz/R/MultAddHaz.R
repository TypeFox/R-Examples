MultAddHaz <- function(formula, data, subset, weights, na.action,
                       model = TRUE, contrasts = NULL, start.val,
                       attrib = FALSE, attrib.var,
                       type.attrib = "abs", set.seed = FALSE, seed,
                       bootstrap = FALSE,
                       nbootstrap, parallel = FALSE, type.parallel = "snow",
                       ncpus = 4,...){

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

  if(is.null(model.weights(mf))){
    w <- 1} else {
      w <- model.weights(mf)
    }

  if(mean(w) != 1){
    w <- w/mean(w)}

  offset <- model.offset(mf)
  x <- model.matrix(mt, mf, contrasts)

  ### Factor y
  y.levels <- sort(unique(y))

  if(is.factor(y)){
    y. <- model.matrix(~y - 1)
  }

  if(!is.factor(y)){
    y.. <- factor(y, levels = y.levels, labels = y.levels)
    y. <- model.matrix(~y.. -1)
  }
  colnames(y.) <- y.levels
  y.resp <- y.[,-1]

  ### LL function
  MultAddHazLogLik <- function(start.val, data, y, x, wgt){
    beta_ij <- matrix(unlist(split(start.val, cut(seq_along(start.val), ncol(y.resp),
                                              labels = FALSE))), ncol = ncol(y.resp))
    eta_ij <- matrix(NA, ncol = ncol(beta_ij), nrow = nrow(x))
    pi_ij <- matrix(NA, ncol = ncol(beta_ij), nrow = nrow(x))

    for (i in 1:ncol(beta_ij)){
      eta_ij[,i] <- as.vector(x %*% beta_ij[,i])}

    pi_i = 1-exp(-(apply(eta_ij, 1, sum)))
    for (i in 1:ncol(beta_ij)){
    pi_ij[,i] <- pi_i *(eta_ij[,i]/apply(eta_ij, 1, sum))}

    sum.y <- apply(y.resp, 1, sum)
    sum.pi_ij <- apply(pi_ij, 1, sum)

     LL <- -sum(2 * wgt * (((1 - sum.y) * log(ifelse(sum.y == 1, 1, (1 - sum.pi_ij)))) +
                          apply((y.resp * log(ifelse(y.resp == 0, 1, (pi_ij)))),
                                1, sum)))
     return(LL)}

  ### Optimization
  sparse.mat <- paste0("bdiag(", paste0(rep("x,", ncol(y.resp)-1), collapse=""), "x)")
  ui.const <- as.matrix(eval(parse(text = sparse.mat)))
  MultLL <- constrOptim(theta = start.val, data = data, f = MultAddHazLogLik, ui = ui.const,
                        ci = rep(0, nrow(ui.const)), control = list(maxit = 1000),
                        method = "Nelder-Mead", y = y, x = x, wgt = w)

  Coeff <- matrix(unlist(split(MultLL$par, cut(seq_along(MultLL$par),
                                               ncol(y.resp)))), ncol = ncol(y.resp))
  colnames(Coeff) <- colnames(y.resp)
  rownames(Coeff) <- colnames(x)
  Coeff2 <- matrix(Coeff, ncol = 1)
  colnames(Coeff2) <- "Coefficients"
  rownames(Coeff2) <- paste0(colnames(x), rep(1:ncol(y.resp), each = ncol(x)))

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
      attrib.coef = matrix(Coeff[1, ], ncol = ncol(y.resp))
      mult.coef = Coeff[-1,]} else {
        attrib.coef = matrix(Coeff[1: attrib.length, ], ncol = ncol(y.resp))
        mult.coef = matrix(Coeff[-c(1: attrib.length), ], ncol = ncol(y.resp))
      }


    AttribMult <- function(mult.coef, attrib.coef, data, y, x.mat, att.var.,
                           wgt, type.attrib){
      id.disab <- list()
      hdis <- list()
      haz.back <- list()
      eta_ij <- list()
      pi_ij <- list()
      att.x <- list()
      att.back <- list()
      att.mat <- list()

      for (i in 1:ncol(y.resp)){
        id.disab[[i]] <- which(y == 0 | y == i)
        hdis[[i]] <- t(mult.coef[, i] * t(x.mat))
        haz.back[[i]] <- attrib.coef[,i][att.var.]
        eta_ij[[i]] <- haz.back[[i]] + apply(hdis[[i]], 1, sum)

        pi_i = 1-exp(-(Reduce("+",eta_ij)))
        sum.eta_ij <- Reduce("+",eta_ij)
        pi_ij[[i]] <- eta_ij[[i]]/sum.eta_ij * pi_i

        att.x[[i]] <- (hdis[[i]]/eta_ij[[i]]) * pi_ij[[i]]
        att.back[[i]] <- (haz.back[[i]]/eta_ij[[i]]) * pi_ij[[i]]

        # Attribution matrix
        att.mat[[i]] <- matrix(NA, nrow = (ncol(x.mat) + 3), ncol = attrib.length)
        att.mat[[i]][1,] <-  tapply(wgt[id.disab[[i]]], att.var.[id.disab[[i]]], sum)
        att.mat[[i]][2,] <-  tapply(pi_ij[[i]][id.disab[[i]]] * wgt[id.disab[[i]]],
                                    att.var.[id.disab[[i]]], sum)
        att.mat[[i]][3,] <-  tapply(att.back[[i]][id.disab[[i]]] * wgt[id.disab[[i]]],
                                    att.var.[id.disab[[i]]], sum)

          for (j in 1: ncol(x.mat)){
          att.mat[[i]][3 + j,] <- tapply(att.x[[i]][, j][id.disab[[i]]] * wgt[id.disab[[i]]],
                                         att.var.[id.disab[[i]]], sum)
        }

        dimnames(att.mat[[i]]) <- list(c("nn","disab","backgrnd", colnames(x.mat)),
                                       levels(att.var.))
      }

      attribution <- rep(list(list()), ncol(y.resp))
      attribution2 <- rep(list(list()), ncol(y.resp))
      att.final <- rep(list(list()), ncol(y.resp))
      att.abs <- rep(list(list()), ncol(y.resp))
      att.rel <- rep(list(list()), ncol(y.resp))
      names(att.final) <- colnames(y.resp)
      final.list <- list()

      for(i in 1: ncol(y.resp)){
        for(j in 1: ncol(att.mat[[1]])){
          attribution[[i]][[j]] <- att.mat[[i]][, j]
          row_sub <- attribution[[i]][[j]] != 0
          attribution2[[i]][[j]] <- attribution[[i]][[j]][row_sub]

          if(type.attrib == "rel"){
            att.final[[i]][[j]] <- attribution2[[i]][[j]][3: length(attribution2[[i]][[j]])]/
              attribution2[[i]][[j]][2]
          }
          if(type.attrib == "abs") {
            att.final[[i]][[j]] <- attribution2[[i]][[j]][2: length(attribution2[[i]][[j]])]/
              attribution2[[i]][[j]][1]
          }
          if(type.attrib == "both") {
            att.rel[[i]][[j]] <- attribution2[[i]][[j]][3: length(attribution2[[i]][[j]])]/
              attribution2[[i]][[j]][2]
            att.abs[[i]][[j]] <- attribution2[[i]][[j]][2: length(attribution2[[i]][[j]])]/
              attribution2[[i]][[j]][1]
            att.final <- list(att.rel = att.rel, att.abs = att.abs)
          }
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
        rownames(Out1[[1]]) <- paste0(names(att.final[[1]][[1]][[1]]),
                                      rep(1: ncol(y.resp),
                                          each = length(names(att.final[[1]][[1]][[1]]))))
        rownames(Out1[[2]]) <- paste0(names(att.final[[2]][[2]][[1]]),
                                      rep(1: ncol(y.resp),
                                          each = length(names(att.final[[2]][[2]][[1]]))))

        } else {
          rownames(Out1[[1]]) <- paste0(names(Out[[1]]),  rep(levels(att.var.)[-1],
                                        each = length(Out[[1]])/ncol(y.resp)),
                                        rep(rep(1:attrib.length,
                                        each = length(names(att.final[[1]][[1]][[1]]))),
                                        ncol(y.resp)))

          rownames(Out1[[2]]) <- paste0(names(Out[[2]]),  rep(levels(att.var.)[-1],
                                        each = length(Out[[2]])/ncol(y.resp)),
                                        rep(rep(1:attrib.length,
                                        each = length(names(att.final[[2]][[2]][[1]]))),
                                        ncol(y.resp)))
        }

        Output <- Out1}
        return(Output)

    }

    Attribution <- AttribMult(mult.coef = mult.coef, attrib.coef = attrib.coef,
                              data = data, y = y, x.mat = x.mat,
                              att.var. = att.var., wgt = w,
                              type.attrib = type.attrib)
  }

  ### CI
  if(bootstrap == FALSE){

    CovMultAddHaz <- function(start.val, x, y){
      beta_ij <- matrix(unlist(split(start.val, cut(seq_along(start.val), ncol(y.resp),
                                                labels = FALSE))), ncol = ncol(y.resp))
      eta_ij <- matrix(NA, ncol = ncol(beta_ij), nrow = nrow(x))
      pi_ij <- matrix(NA, ncol = ncol(beta_ij), nrow = nrow(x))
      for (i in 1:ncol(beta_ij)){
        eta_ij[,i] <- as.vector(x %*% beta_ij[,i])
      }
      pi_i = 1-exp(-(apply(eta_ij, 1, sum)))
      for (i in 1:ncol(beta_ij)){
        pi_ij[,i] <- pi_i *(eta_ij[,i]/apply(eta_ij, 1, sum))}
      sum.y <- apply(y.resp, 1, sum)
      sum.pi_ij <- apply(pi_ij, 1, sum)

      # j = j
      jj.1 <- matrix(NA, nrow = nrow(x), ncol = ncol(y.resp))
      jj.2 <- matrix(NA, nrow = nrow(x), ncol = ncol(y.resp))

      for (i in 1: ncol(y.resp)){
        jj.1[, i] <- (1 - sum.y) * ( (((1 - pi_ij[,i]) * (1 - sum.pi_ij)) - (1 - pi_ij[,i])^2)/(1 - sum.pi_ij)^2)
        jj.2[, i] <- sum.y * ((pi_ij[, i] - 1)/(pi_ij[, i])^2)}
      jjd <- jj.1 + jj.2

      jj <- list()
      for (i in 1:ncol(y.resp)){
        jj[[i]] <- t(x) %*% diag(jjd[,i]) %*% x}

      # j != j'
      jj.prime <- list()

      for (i in 1:ncol(y.resp)){
        jj.prime[[i]] <- (1 - sum.y) * ( ((1 - pi_ij) * (1 - pi_ij[, i]))/(1 - sum.pi_ij)^2)}

      jjp.mat <- matrix(unlist(jj.prime), ncol = ncol(y.resp)*ncol(y.resp))

      # Removing j = j'
      seq <- NULL
      for (i in 0:(ncol(y.resp) -1)){
        seq[i] = 1 + (i * ncol(y.))}

      jjp.mat2 <- jjp.mat[,-c(1,seq)]

      # Removing duplicates
      if(ncol(y.resp) > 2) {
        dup.col <- duplicated(t(jjp.mat2))
        jjp.mat3 <- jjp.mat2[, !dup.col]} else {
        jjp.mat3 <- jjp.mat2}

      # selecting the important columns of the matrix
      jjp.d <- list()
      for (i in 1:ncol(y.resp)){
        jjp.d[[i]] <- t(x) %*% diag(jjp.mat3[,i]) %*% x}

      if (ncol(y.resp) == 2){
        cov1 <- cbind(jj[[1]], jjp.d[[1]])
        cov2 <- cbind(jjp.d[[2]], jj[[2]])
        vcov <- -1 * rbind(cov1, cov2)} else {

        cov.mat <- bdiag(jj)

        used <- 0

        a <- ncol(x)
        b <- ncol(y.resp)

        for (irow in 1:(b - 1)){

          tempMat <- jjp.d[[used + 1]]

          if ((used + 1) != length(jjp.d)){
             for (k in (used + 2):(used + b - irow)){
               tempMat <- cbind(tempMat, jjp.d[[k]])}
          }
          rows <- seq((irow - 1) * a + 1, irow * a, length = a)
          columns <- seq(irow*a+1, a*b, length=a*b-irow*a)
          cov.mat[rows, columns] <- tempMat

          used <- used+b-irow
        }

        usedL <- 0

        for (irow in 2:b){

          tempMat <- jjp.d[[usedL+1]]

          if ((usedL+1) > 1){
            for (k in (usedL+2):(usedL+2-b+irow)){
              tempMat <- cbind(tempMat, jjp.d[[k]])}
          }
          rows <- seq((irow-1)*a+1, irow*a, length=a)
          columns <- seq(1, a*(irow-1), length=a*(irow-1))
          cov.mat[rows, columns] <- tempMat

          usedL <- usedL+2-b+irow}
        vcov <- -1 * cov.mat
      }
      return(vcov)}

    vcov <- CovMultAddHaz(start.val = Coeff, x = x, y = y)
    Invcov <- solve(vcov)
    stdError <- as.vector(sqrt(diag(Invcov)))

    Std <- matrix(unlist(split(stdError, cut(seq_along(stdError), ncol(y.resp)))),
                  ncol = ncol(y.resp))
    colnames(Std) <- colnames(y.resp)
    rownames(Std) <- colnames(x)

    CILow <- Coeff2 - (1.96 * stdError)
    CIHigh <- Coeff2 + (1.96 * stdError)
    pvalue <- round(2 * pnorm(-abs(Coeff/stdError)), 4)

    CI <- as.matrix(cbind(CILow, CIHigh), ncol = 2,
                    dimnames = list(paste0(colnames(x), rep(1:ncol(y.resp),
                                                            each = ncol(x)))))
    colnames(CI) <- c("CI2.5", "CI97.5")

    colnames(pvalue) <- colnames(pvalue)
    rownames(pvalue) <- colnames(x)

    if(attrib){
      Results <- list(coefficients = Coeff2, ci = CI,
                      resDeviance = MultLL$value,
                      df = nrow(data) - length(MultLL$par),
                      pvalue = matrix(pvalue, ncol = 1,
                                      dimnames = list(rownames(Coeff2), "p.value")),
                      stdError = matrix(Std, ncol = 1, dimnames = list(rownames(Coeff2), "StdErr")),
                      vcov = vcov, contribution = Attribution)
    } else {
    Results <- list(coefficients = Coeff2, ci = CI,
                    resDeviance = MultLL$value,
                    df = nrow(data) - length(MultLL$par),
                    pvalue = matrix(pvalue, ncol = 1,
                                    dimnames = list(rownames(Coeff2), "p.value")),
                    stdError = matrix(Std, ncol = 1, dimnames = list(rownames(Coeff2), "StdErr")),
                    vcov = vcov, contribution = attrib) }
    }

  ### bootstrap CI
  if(bootstrap == TRUE){
  if(attrib){
  mystat <- function(data, indices) {
    indices2 <- rep(indices, length(sort(unique(y))) - 1)
    m1 <-  constrOptim(theta = MultLL$par, data = data[indices,], f = MultAddHazLogLik,
                      ui = ui.const[indices2, ], ci = rep(0, nrow(ui.const)),
                      control=list(maxit=1000), method = "Nelder-Mead",
                      y = y[indices], x = x[indices,], wgt = w[indices])

    m2 <- AttribMult(mult.coef = mult.coef, attrib.coef = attrib.coef,
               data = data[indices,], y = y[indices], x.mat = x.mat[indices,],
               att.var. = att.var.[indices], wgt = w[indices],
               type.attrib = type.attrib)

    if(is.list(m2)){
      m3 <- unlist(m2)} else {
      m3 <- m2}
    names(m1$par) <- rownames(Coeff2)
    res <- c(m1$par, m3)
    return(res)
  }
} else {
  mystat <- function(data, indices) {
    indices2 <- rep(indices, length(sort(unique(y))) - 1)
    m1 <-  constrOptim(theta = MultLL$par, data = data[indices,], f = MultAddHazLogLik,
                       ui = ui.const[indices2, ], ci = rep(0, nrow(ui.const)),
                       control=list(maxit=1000), method = "Nelder-Mead",
                       y = y[indices], x = x[indices,], wgt = w[indices])

    return(m1$par)}
}

  if(set.seed){
    set.seed(seed)}

  if (parallel == TRUE){
    BootResult <- boot(data = data, statistic = mystat, R = nbootstrap, parallel = type.parallel,
                       strata = y, ncpus = ncpus)} else {
    BootResult <- boot(data = data, statistic = mystat, R = nbootstrap, strata = y)}

  BootCI <- matrix(NA, ncol = 2, nrow = length(BootResult$t0))

  for (i in 1:length(BootResult$t0)){
    BootCI[i,] <- boot.ci(BootResult,  conf = 0.95, type = "perc", index = i)[[4]][, 4:5]}

  colnames(BootCI) <- c("CILow", "CIHigh")

  if(attrib){
    CILow.att <- matrix(BootCI[(length(start.val) + 1): length(BootResult$t0), 1], ncol = 1,
                        dimnames = list(rownames(Attribution)))
    CIHigh.att <- matrix(BootCI[(length(start.val) + 1): length(BootResult$t0), 2], ncol = 1,
                         dimnames = list(rownames(Attribution)))

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
        att1 <- list(cbind(Attribution, CILow.att, CIHigh.att))
    }
    colnames(att1[[1]]) <- c("Contribution", "CI2.5", "CI97.5")


    Results <- list(coefficients = Coeff2, ci = cbind(matrix(BootCI[1: length(start.val), "CILow"],
                                                      ncol = 1,
                                                      dimnames = list(rownames(Coeff2))),
                                                      matrix(BootCI[1:length(start.val), "CIHigh"], ncol = 1,
                                                      dimnames = list(rownames(Coeff2)))),
                    resDeviance = MultLL$value, df = nrow(data) - length(MultLL$par),
                    contribution = att1, bootsRep = BootResult$t)

    colnames(Results$bootsRep) <- names(BootResult$t0)
    colnames(Results$ci) <- c("CI2.5", "CI97.5")

    if(type.attrib != "both"){
      names(Results$contribution) <- type.attrib}

    } else {

    Results <- list(coefficients = Coeff2, ci = cbind(matrix(BootCI[1: length(start.val), "CILow"],
                                                             ncol = 1,
                                                      dimnames = list(rownames(Coeff2))),
                                                      matrix(BootCI[1:length(start.val), "CIHigh"], ncol = 1,
                                                      dimnames = list(rownames(Coeff2)))),
                      resDeviance = MultLL$value, df = nrow(data) - length(MultLL$par),
                      contribution = attrib, bootsRep = BootResult$t)

    colnames(Results$bootsRep) <- names(BootResult$t0)
    colnames(Results$ci) <- c("CI2.5", "CI97.5")}
  }

  Results$bootstrap <- bootstrap
  Results$call <- match.call()
  class(Results) <- "multaddhazmod"
  return(Results)
}

# GENERIC FUNCTIONS IN S3 LANGUAGE (PRINT, SUMMARY, PREDICT)
multaddhazmod <- function(x, ...) UseMethod("multaddhazmod")


print.multaddhazmod <- function(x,...){
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)}

summary.multaddhazmod <- function(object, ...){
  if(object$bootstrap == FALSE){
    se <- sqrt(diag(solve(object$vcov)))
    tval <- coef(object)/se
    TAB <- cbind(Estimate = as.vector(coef(object)),
                 StdErr = se,
                 t.value = as.vector(tval),
                 p.value = as.vector(2*pt(-abs(tval), df=object$df)))} else {
                   TAB <- cbind(Estimate = coef(object), CI = object$ci)}
    rownames(TAB) <- rownames(object$coefficients)

  if(object$bootstrap == FALSE){
    res <- list(call = object$call, bootstrap = object$bootstrap,
                coefficients = capture.output(printCoefmat(TAB, P.values = TRUE,
                                                           has.Pvalue = TRUE)))}
  else {res <- list(coefficients = TAB, bootstrap = object$bootstrap)}
  class(res) <- "summary.multaddhazmod"
  res}

predict.multaddhazmod <- function(object, newdata = NULL, ...){
  if(is.null(newdata))
    y <- fitted(object)
  else{
    if(!is.null(object$formula)){
      x <- model.matrix(object$formula, newdata)
    } else{
      x <- newdata}
    y <- as.vector(x %*% coef(object))}
  y}

