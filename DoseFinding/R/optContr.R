## functions for calculating optimal contrasts and critical value

optC <- function(mu, Sinv = NULL, placAdj = FALSE){
  ## calculate optimal contrast for given mu and Sinv (Sinv = proportional to inv covariance matrix)
  if(!placAdj){
    aux <- rowSums(Sinv)  # Sinv %*% 1
    mn <- sum(mu * aux)/sum(aux) # formula is: S^(-1)(mu-mu*S^(-1)*1/(1*S^(-1)1)1)
    val <- Sinv %*% (mu - mn)
    ## now center so that sum is 0
    ## and standardize to have norm 1
    val <- val - sum(val)
  } else { # placAdj = TRUE
    val <- Sinv %*% mu     
  }
    val/sqrt(sum(val^2))
}

constOptC <- function(mu, Sinv = NULL, placAdj = FALSE, direction){
  ## calculate optimal contrasts under the additional constraint that
  ## the control and the active treatment groups have a different sign
  ## in the contrast
  S <- solve(Sinv) # ugly fix, we should use S as argument
  if(!placAdj){
    k <- length(mu)
    CC <- cbind(-1,diag(k-1))
    SPa <- CC%*%S%*%t(CC)
    muPa <- as.numeric(CC%*%mu)
  } else {
    k <- length(mu)+1
    SPa <- S
    muPa <- mu
  }
  ## determine direction of effect
  unContr <- solve(SPa)%*%muPa # unconstrained optimal contrast
  mult <- ifelse(direction == "increasing", 1, -1) # 1 increasing, -1 decreasing
  ## prepare call of quadprog::solve.QP
  D <- SPa
  d <- rep(0,k-1)
  tA <- rbind(muPa, 
              mult*diag(k-1))
  A <- t(tA)
  bvec <- c(1,rep(0,k-1))
  contr <- quadprog::solve.QP(D, d, A, bvec, meq=1)$solution
  contr[abs(contr) < 1e-10] <- 0
  if(!placAdj)
    contr <- c(-sum(contr), contr)
  contr/sqrt(sum(contr^2))
}


modContr <- function(means, W = NULL, Sinv = NULL, placAdj = FALSE,
                     type, direction){
  ## call optC on matrix
  ## check whether constant shape was specified and remove (can happen for linInt model)
  if(!placAdj){ 
    ind <- apply(means, 2, function(x){
      length(unique(x)) > 1 
    })
  } else { ## placAdj
    ind <- apply(means, 2, function(x){
      any(x != 0) 
    })
  }
  if(all(!ind))
    stop("All models correspond to a constant shapes, no optimal contrasts calculated.")
  if(any(!ind)){
    nam <- colnames(means)[!ind]
    namsC <- paste(nam, collapse = ", ")
    if(length(nam) == 1){
      message("The ", namsC, " model has a constant shape, cannot
calculate optimal contrasts for this shape.")
    } else {
      message("The ", namsC, " models have a constant shape, cannot
calculate optimal contrasts for these shapes.")
    }
    means <- means[,ind, drop=FALSE]
  }

  if(is.null(Sinv))
    Sinv <- solve(W)
  if(type == "unconstrained"){
    out <- apply(means, 2, optC, Sinv = Sinv, placAdj = placAdj)
  } else { # type == "constrained"
    out <- apply(means, 2, constOptC, Sinv = Sinv,
                 placAdj = placAdj, direction = direction)
  }
  if(!is.matrix(out)){ ## can happen for placAdj=T and only 1 act dose
    nam <- names(out)
    out <- matrix(out, nrow = 1)
    colnames(out) <- nam
  }
  out
}

optContr <-  function(models, doses, w, S, placAdj = FALSE,
                      type = c("unconstrained", "constrained")){
  ## calculate optimal contrasts and critical value
  if(!(inherits(models, "Mods")))
    stop("models needs to be of class Mods")
  if(missing(doses))
    doses <- attr(models, "doses")
  scal <- attr(models, "scal")
  off <- attr(models, "off")
  nodes <- attr(models, "doses")
  direction <- unique(attr(models, "direction"))
  if(length(direction) > 1)
    stop("need to provide either \"increasing\" or \"decreasing\" as direction to optContr")
  mu <- getResp(models, doses)
  if(placAdj){ 
    mu0 <- getResp(models, 0)
    mu <- mu-matrix(mu0[1,], byrow = TRUE,
                    nrow=nrow(mu), ncol=ncol(mu))
  }
  type <- match.arg(type)
  if(type == "constrained"){
    avail <- requireNamespace("quadprog", quietly = TRUE)
    if(!avail)
      stop("Need suggested package quadprog to calculate constrained contrasts")
  }
  if(any(doses == 0) & placAdj)
    stop("If placAdj == TRUE there should be no placebo group in \"doses\"")
  ## check for n and vCov arguments 
  if(!xor(missing(w), missing(S)))
    stop("Need to specify exactly one of \"w\" or \"S\"")
  if(!missing(w)){
    if(length(w) == 1){ # assume equal weights
      S <- Sinv <- diag(length(doses))
    } else {
      if(length(w) != length(doses))
        stop("w needs to be of length 1 or of the same length as doses")
      S <- diag(1/w)
      Sinv <- diag(w)
    }
  } else { 
    if(!is.matrix(S))
      stop("S needs to be a matrix")
    Sinv <- solve(S)
  }
  contMat <- modContr(mu, Sinv=Sinv, placAdj = placAdj,
                      type = type, direction = direction)
  rownames(contMat) <- doses
  corMat <- cov2cor(t(contMat) %*% S %*% contMat)
  res <- list(contMat = contMat, muMat = mu, corMat = corMat)
  attr(res, "type") <- type
  attr(res, "placAdj") <- placAdj
  class(res) <- "optContr"
  res
}

print.optContr <- function(x, digits = 3, ...){
  cat("Optimal contrasts\n")
  print(round(x$contMat, digits))
}

summary.optContr <- function(object, digits = 3, ...){
  class(object) <- "summary.optContr"
  print(object, digits = digits)
}

print.summary.optContr <- function(x, digits = 3, ...){
  cat("Optimal contrasts\n")
  cat("\n","Optimal Contrasts:","\n", sep="")
  print(round(x$contMat, digits))
  cat("\n","Contrast Correlation Matrix:","\n", sep="")
  print(round(x$corMat, digits))  
  cat("\n")
}

plot.optContr <- function (x, superpose = TRUE, xlab = "Dose",
                           ylab = NULL, plotType = c("contrasts", "means"), ...){
  plotType <- match.arg(plotType)
  if (is.null(ylab)) {
    if (plotType == "contrasts") {
      ylab <- "Contrast coefficients"
    } else {
      ylab <- "Normalized model means"
    }
  }
  cM <- x$contMat
  if (plotType == "means")
    cM <- t(t(x$muMat)/apply(x$muMat, 2, max))
  nD <- nrow(cM)
  nM <- ncol(cM)
  cMtr <- data.frame(resp = as.vector(cM),
                     dose = rep(as.numeric(dimnames(cM)[[1]]), nM),
                     model = factor(rep(dimnames(cM)[[2]], each = nD),
                     levels = dimnames(cM)[[2]]))
  if(superpose){
    spL <- trellis.par.get("superpose.line")
    spL$lty <- rep(spL$lty, nM%/%length(spL$lty) + 1)[1:nM]
    spL$lwd <- rep(spL$lwd, nM%/%length(spL$lwd) + 1)[1:nM]
    spL$col <- rep(spL$col, nM%/%length(spL$col) + 1)[1:nM]
    ## number of columns in legend
    nCol <- ifelse(nM < 5, nM, min(4,ceiling(nM/min(ceiling(nM/4),3))))
    key <- list(lines = spL, transparent = TRUE, 
                text = list(levels(cMtr$model), cex = 0.9),
                columns = nCol)
    ltplot <- xyplot(resp ~ dose, data = cMtr, subscripts = TRUE,
                     groups = cMtr$model, panel = panel.superpose,
                     type = "o", xlab = xlab, ylab = ylab,
                     key = key, ...)
  } else {
    ltplot <- xyplot(resp ~ dose | model, data = cMtr, type = "o", 
                     xlab = xlab, ylab = ylab,
                     strip = function(...){
                       strip.default(..., style = 1)
                     }, ...)
  }
  print(ltplot)
}

