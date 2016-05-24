systemfitAR <- function(formula, method = "OLS", inst = NULL, data = list(),
  restrict.matrix = NULL, restrict.rhs = NULL, restrict.regMat = NULL, 
  pooled = FALSE, control = systemfit.control( ... ),
  AR1 = FALSE, rho.sel = c("all", "mean"), model = c("static", "dynamic"), ...) 
{  
  rho.sel <- match.arg(rho.sel) 
  model <- match.arg(model)
  # A1. Estimate SUR and get rho estimates
  est <- systemfit(formula=formula, method="SUR", data=data)
  resids_all <- rho_ind <- NULL
  for (i in 1:length(est$eq)) {
    resids <- est$eq[[i]]$residuals
    resids_all <- c(resids_all, resids)
    res <- resids[-length(resids)]
    rho_ind <- c(rho_ind, sum(resids[-1] * res) / sum(res^2)) 
    rho_ste_ind <- sqrt((1 - rho_ind^2) / length(resids)) 
  }
  res_all <- resids_all[-length(resids_all)]
  rho_all <- sum(resids_all[-1] * res_all) / sum(res_all^2) 
  rho_ste_all <- sqrt((1 - rho_all^2) / length(resids_all))
  if (rho.sel == "all") {
    rho <- rho_all; rho_ste <- rho_ste_all
  } else {
    rho <- mean(rho_ind); rho_ste <- rho_ste_ind
  }
  
  # A2. Adjust raw data for static model; intercept is adjusted/included
  if (AR1 & model == 'static') {
    if(inherits(formula, "formula")) {formula <- list(formula)}
    name.y <- name.x <- NULL; formu.adj <- list()
    for (i in 1:length(formula)) {
      name.y <- c(name.y, all.vars(formula[[i]])[ 1])
      name.x <- c(name.x, all.vars(formula[[i]])[-1])
      formu.adj[[i]] <- bsFormu(name.y=all.vars(formula[[i]])[1], 
        name.x=c("intercept_adj", all.vars(formula[[i]])[-1]), intercept=FALSE)
    }
    nam.y <- unique(name.y); yy <- data[, nam.y]
    nam.x <- unique(name.x); xx <- data[, nam.x]
    if(length(nam.y) == 1) {
      adj.yy <- yy[-1]  - rho * yy[-length(yy)]
    } else { 
      adj.yy <- yy[-1,] - rho * yy[-nrow(yy),] # lost first obs
    }
    adj.intercept <- 1 - rho
    adj.xx <- xx[-1,] - rho * xx[-nrow(xx),]
    data.adj <- cbind(adj.yy, adj.intercept, adj.xx)
    names(data.adj) <- c(nam.y, "intercept_adj", nam.x)
    formula <- formu.adj
    data <- data.adj
  }
  
  # A2.2 Adjust raw data for dynamic model, no intercept
  if (AR1 & model == 'dynamic') {
    if(inherits(formula, "formula")) {formula <- list(formula)}
    name.y <- name.x <- NULL
    for (i in 1:length(formula)) {
      name.y <- c(name.y, all.vars(formula[[i]])[ 1])
      name.x <- c(name.x, all.vars(formula[[i]])[-1])
    }
    nam.y <- unique(name.y); yy <- data[, nam.y]
    nam.x <- unique(name.x); xx <- data[, nam.x]
    if(length(nam.y) == 1) {
      adj.yy <- yy[-1]  - rho * yy[-length(yy)]
    } else { 
      adj.yy <- yy[-1,] - rho * yy[-nrow(yy),]
    }
    adj.xx <- xx[-1,] - rho * xx[-nrow(xx),]
    data.adj <- cbind(adj.yy, adj.xx)
    names(data.adj) <- c(nam.y, nam.x)
    data <- data.adj
  }
  
  # B. Estimate SUR on new data
  result <- systemfit(formula=formula, method=method, inst=inst, data=data,
    restrict.matrix=restrict.matrix, restrict.rhs=restrict.rhs, 
    restrict.regMat=restrict.regMat, pooled=pooled, control=control)
  result$rho <- rho
  result$rho_ste <- rho_ste
  result$data <- data
  result$formula <- formula
  return(result)
}