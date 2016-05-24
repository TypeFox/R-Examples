aco1arm <- function (data, svtime, event, treatment, BaselineMarker, id, 
          subcohort, esttype = 1, augment = 1, extra=NULL) 
{
  
  if (!is.data.frame(data)) {
    stop("Argument data must be a data.frame object.")
  }else {
    colNames <- colnames(data)
    if (!(svtime %in% colNames)) {
      stop("Survival time variable was not found in the data.")
    }
    if (!(event %in% colNames)) {
      stop("Failure time indicator variable was not found in the data.")
    }
    if (!(treatment %in% colNames)) {
      stop("Treatment variable was not found in the data.")
    }
    else {
      if (any(levels(factor(data[, treatment])) != c("0", 
                                                     "1"))) {
        warning("Treatment variable must be either 0 or 1 only.")
      }
    }
    if (!(BaselineMarker %in% colNames)) {
      stop("BaselineMarker variable was not found in the data.")
    }
    if (!(id %in% colNames)) {
      stop("id variable was not found in the data.")
    }
    if (!(subcohort %in% colNames)) {
      stop("Subcohort indicator variable was not found in the data.")
    }
    if (! augment %in% c(0,1)) {
      stop("augment variable must be either 0 or 1 only.")
    }
    if (!is.null(extra)) {
      if (!any(extra %in% colNames)) {
        extraNotFound <- paste(extra[!(extra %in% colNames)], 
                               sep = "", collapse = ", ")
        stop(paste("Extra variable(s) was not found in the data:", 
                   extraNotFound))
      }
      tmp <- remove_rarevariants(data[,extra])
      if (any(tmp))
      {
        idx <- tmp==TRUE
        toremove=NULL
        for (i in 1:length(idx))
        {
          if (idx[i]) toremove <- c(toremove,extra[idx[i]])
        }
        warnings(paste0(paste(toremove,sep=", "), " were removed due to rare vairant"))
        extra <- extra[!idx]
        if (length(extra)==0) extra <- NULL
      }
    }
  }

  #Remove missing data
  dat0 <- remove_missingdata(data)$data
  #Biomarker variable should be transformed to numeric
  if (! is.numeric(dat0[, BaselineMarker]))
    dat0[, BaselineMarker] <- char2num(dat0[, BaselineMarker])
  if (remove_rarevariants(dat0[, BaselineMarker]))
  {
    warnings("BaselineMarker variable is rare variant")
    tmpResult <- data.frame(beta=rep(NA,length(extra)+3), stder=rep(NA,length(extra)+3), pVal=rep(NA,length(extra)+3))
    rownames(tmpResult)[4:nrow(tmpResult)]=extra
  }else
  {
    idx <- which(colnames(dat0) == id)
    colnames(dat0)[idx] <- "id"
    idx <- which(colnames(dat0) == subcohort)
    colnames(dat0)[idx] <- "subcohort"
    subcohort <- "subcohort"
    cases <- dat0[dat0[, event] == 1, ]
    fit4 <- glm(cases[, treatment] ~ cases[, BaselineMarker], 
                family = binomial, x = TRUE, y = TRUE)
    bread1 <- t(fit4$x) %*% (fit4$x * fit4$fitted * (1 - fit4$fitted))
    b23 <- fit4$coef
    var23 <- (summary(fit4)$coef[, 2])^2
    dat1 <- dat0[dat0[, treatment] == augment, ]
    n1 <- nrow(dat1)
    nx <- 1 + length(extra)
    ww <- rep(1, n1)
    n <- sum(data[, treatment] == augment)
    fmla <- as.formula(paste0("Surv(", svtime, ",", event, ") ~ ", 
                              paste(paste(c(BaselineMarker, extra), collapse = "+"))))
    if (esttype == 1) {
      fit3 <- cch(fmla, data = dat1, subcoh = ~subcohort, id = ~id, 
                  cohort.size = n, method = "SelfPrentice")
      dat2 <- dat1[dat1[, subcohort] == 1, ]
      xx <- dat2[, c(BaselineMarker, extra)]
      yy <- dat2[, svtime]
      id.ss2 <- which(dat1[, subcohort] == 1)
    }
    else {
      fit3 <- cch(fmla, data = dat1, subcoh = ~subcohort, id = ~id, 
                  cohort.size = n, method = "LinYing", robust = TRUE)
      ww[dat1[, event] == 0] <- (n - sum(dat1[, event] == 1))/(sum(dat1[, 
                                                                        subcohort] == 1) - sum(dat1[, subcohort] == 1 & dat1[, 
                                                                                                                             event] == 1))
      xx <- dat1[, c(BaselineMarker, extra)]
      yy <- dat1[, svtime]
      id.ss2 <- 1:n1
    }
    n2 <- length(id.ss2)
    dat1 <- as.matrix(dat1)
    beta <- fit3$coef
    xx <- as.matrix(xx)
    a <- exp(xx %*% beta) * ww[id.ss2]
    s0 <- rep(0, n1)
    s1 <- matrix(0, n1, nx)
    dd1 <- matrix(0, n1, nx)
    for (i in which(dat1[, event] == 1)) {
      b <- 1 * (yy >= dat1[i, svtime])
      if (sum(b) > 0) {
        s0[i] <- sum(b * a)
        s1[i, ] <- apply(matrix(b * a, nrow = n2, ncol = nx, 
                                byrow = FALSE) * xx, 2, sum)
        dd1[i, ] <- dat1[i, c(BaselineMarker, extra)] - s1[i, 
                                                           ]/s0[i]
      }
      else {
        dd1[i, ] <- as.matrix(dat1[i, c(BaselineMarker, extra)])
      }
    }
    dd2 <- matrix(0, n1, nx)
    for (i in id.ss2) {
      tmp <- matrix(0, n1, nx)
      for (j in which(dat1[, event] == 1)) {
        tmp[j, ] <- ww[i] * (dat1[i, svtime] >= dat1[j, svtime]) * 
          exp(dat1[i, c(BaselineMarker, extra)] %*% beta) * 
          (dat1[i, c(BaselineMarker, extra)] - s1[j, ]/s0[j])/s0[j]
      }
      dd2[i, ] <- apply(tmp, 2, sum)
    }
    ss2 <- dd1 - dd2
    infmat <- matrix(0, nx, nx)
    for (i in which(dat1[, event] == 1)) {
      b <- yy >= dat1[i, svtime]
      s2 <- matrix(0, nx, nx)
      for (j in which(b)) {
        s2 <- s2 + as.matrix(a[j] * xx[j, ]) %*% t(as.matrix(xx[j, 
                                                                ]))
      }
      temp <- s2/s0[i] - as.matrix(s1[i, ]) %*% t(as.matrix(s1[i, 
                                                               ]))/(s0[i]^2)
      infmat <- infmat + temp
    }
    v2 <- (solve(infmat) %*% t(ss2) %*% ss2 %*% solve(infmat))
    est.coef <- c(beta[1], b23, beta[-1])
    est.var <- c(v2[1, 1], var23, diag(v2)[-1])
    if (augment == 1) {
      est.coef[1] <- est.coef[1] - fit4$coef[2]
      cov1 <- solve(bread1) %*% t(fit4$x[fit4$y == 1, ] * (fit4$y[fit4$y == 
                                                                    1] - fit4$fitted[fit4$y == 1])) %*% ss2[dat1[, event] == 
                                                                                                              1, ] %*% solve(infmat)
      est.var[1] <- est.var[1] + summary(fit4)$coef[2, 2]^2 - 
        cov1[2, 1] * 2
    }
    
    pVal <- 2*(1-pnorm(abs(est.coef/sqrt(est.var))))
    tmpResult <- data.frame(beta=round(est.coef, 4), stder=round(sqrt(est.var), 4), pVal=pVal)
  }
    
  rownames(tmpResult)[1] <- paste(BaselineMarker, "(BaselineMarker)")
  rownames(tmpResult)[2] <- paste(treatment, "(Treatment)")
  rownames(tmpResult)[3] <- "interatcion"
  
  return(tmpResult)
}
