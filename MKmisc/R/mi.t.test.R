mi.t.test <- function(miData, ...){
  UseMethod("mi.t.test")
}
mi.t.test.default <- function(miData, x, y = NULL, alternative = c("two.sided", "less", "greater"),
                      mu = 0, paired = FALSE, var.equal = FALSE,
                      conf.level = 0.95, subset = NULL, ...){
  alternative <- match.arg(alternative)
  if (!missing(mu) && (length(mu) != 1 || is.na(mu)))
    stop("'mu' must be a single number")
  if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) ||
                               conf.level < 0 || conf.level > 1))
    stop("'conf.level' must be a single number between 0 and 1")
  if(!is.list(miData))
    stop("'miData' must be a list of imputed datasets")
  if(is.null(y) & paired)
    stop("'y' is missing for paired test")

  nrImp <- length(miData)
  Qs <- numeric(nrImp)
  Us <- numeric(nrImp)
  dfs <- numeric(nrImp)
  if(!is.null(y) & !paired){
    Qs2 <- matrix(NA, nrow = nrImp, ncol = 2)
    Us2 <- matrix(NA, nrow = nrImp, ncol = 2)
  }
  for(i in 1:nrImp){
    if(is.null(subset))
      xi <- miData[[i]][,x]
    else
      xi <- miData[[i]][subset, x]
    if(is.null(y))
      yi <- NULL
    else
      yi <- miData[[i]][,y]
    if(paired){
      xi <- xi-yi
      yi <- NULL
    }
    if(is.null(yi)){
      Qs[i] <- mean(xi)
      Us[i] <- var(xi)/length(xi)
      dfs[i] <- length(xi)-1
    }else{
      L <- split(xi, yi)
      Ms <- sapply(L, mean)
      Qs[i] <- Ms[1]-Ms[2]
      Vs <- sapply(L, var)
      Ns <- sapply(L, length)
      if(var.equal){
        dfs[i] <- Ns[1] + Ns[2] - 2
        Us[i] <- (Vs[1]*(Ns[1]-1) + Vs[2]*(Ns[2]-1))/dfs[i]*(1/Ns[1]+1/Ns[2])
      }else{
        Us[i] <- Vs[1]/Ns[1] + Vs[2]/Ns[2]
        stderr1 <- sqrt(Vs[1]/Ns[1])
        stderr2 <- sqrt(Vs[2]/Ns[2])
        stderr <- sqrt(stderr1^2 + stderr2^2)
        dfs[i] <- stderr^4/(stderr1^4/(Ns[1]-1) + stderr2^4/(Ns[2]-1))
      }
      Qs2[i,] <- Ms
      Us2[i,] <- Vs/Ns
    }
  }
  Q <- mean(Qs)
  U <- mean(Us)
  B <- var(Qs)
  T <- U + (1+1/nrImp)*B
  vm <- (nrImp-1)*(1 + U/((1+1/nrImp)*B))^2
  v0 <- mean(dfs)
  gamma <- (1 + 1/nrImp)*B/T
  v.obs <- (v0+1)/(v0+3)*v0*(1-gamma)
  df.mod <- 1/(1/vm + 1/v.obs)
  stderr <- sqrt(T)
  tstat <- (Q-mu)/stderr
  if (alternative == "less") {
    pval <- pt(tstat, df = df.mod)
    cint <- c(-Inf, tstat + qt(conf.level, df = df.mod))
  }
  else if (alternative == "greater") {
    pval <- pt(tstat, df = df.mod, lower.tail = FALSE)
    cint <- c(tstat - qt(conf.level, df = df.mod), Inf)
  }
  else {
    alpha <- 1 - conf.level
    pval <- 2*pt(-abs(tstat), df = df.mod)
    cint <- qt(1-alpha/2, df = df.mod)
    cint <- tstat + c(-cint, cint)
  }
  cint <- mu + cint*stderr
  names(tstat) <- "t"
  names(df.mod) <- "df"
  attr(cint, "conf.level") <- conf.level
  if(is.null(y))
    names(mu) <- "mean"
  else
    names(mu) <- "difference in means"

  if(is.null(yi)){
    est <- c(Q, sqrt(T)*sqrt(length(xi)))
    if(paired){
      method <- "Multiple Imputation Paired t-test"
      dname <- paste("Variables ", x, "and", y)
      names(est) <- c("mean of difference", "SD of difference")
    }else{
      method <- "Multiple Imputation One Sample t-test"
      dname <- paste("Variable ", x, sep = "")
      names(est) <- c("mean", "SD")
    }
    rval <- list(statistic = tstat, parameter = df.mod, p.value = pval,
                 conf.int = cint, estimate = est, null.value = mu,
                 alternative = alternative, method = method,
                 data.name = dname)
  }else{
    Q2 <- colMeans(Qs2)
    names(Q2) <- paste("mean (", levels(yi), ")", sep = "")
    U2 <- colMeans(Us2)
    B2 <- apply(Qs2, 2, var)
    T2 <- sqrt(U2 + (1+1/nrImp)*B2)*sqrt(Ns)
    names(T2) <- paste("SD (", levels(yi), ")", sep = "")
    est <- c(Q2[1], T2[1], Q2[2], T2[2])
    if(var.equal)
      method <- "Multiple Imputation Two Sample t-test"
    else
      method <- "Multiple Imputation Welch Two Sample t-test"
    dname <- paste("Variable ", x, ": ", paste(paste("group", levels(yi)),
                                               collapse = " vs "), sep = "")
    rval <- list(statistic = tstat, parameter = df.mod, p.value = pval,
                 conf.int = cint, estimate = est, null.value = mu,
                 alternative = alternative, method = method,
                 data.name = dname)
  }
  class(rval) <- "htest"
  rval
}
