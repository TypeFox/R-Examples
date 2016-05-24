sf <- function(formula, uhet = NULL, vhet = NULL,
               tmean = NULL, prod = TRUE, data, subset, 
               distribution = c("h", "t"), start.val = NULL, 
               alpha = 0.05, marg.eff = FALSE, digits = 4, 
               print.level = 2) {
 
 if (alpha < 0 | alpha > 99.99) {
  stop("'alpha' must be between 0 and 1 inclusive", call. = FALSE)
 }  
 
 if(length(distribution) != 1){
  stop("Distribution of inefficiency term should be specified.")
 } else {
  distribution <- tolower(substr(distribution, 1,1 ))
 }
 
 if( !distribution %in% c("t","h") ){
  stop("'distribution' is invalid")
 }
 
 if(is.null(tmean) == F & distribution == "h"){
  stop("Option 'tmean' can be used only when distribution of inefficiency term is truncated normal.")
 }
 
 YXZ <- .prepareYXZ(formula = formula, uhet = uhet, vhet = vhet, tmean = tmean, data, subset)
 
 y <- YXZ$Y; X <- YXZ$X; Zu <- YXZ$Zu; Zv <- YXZ$Zv; Zdel <- YXZ$Zdel
 n <- YXZ$n; k <- YXZ$k; ku <- YXZ$ku; kv <- YXZ$kv; kdel <- YXZ$kdel
 esample <- YXZ$esample
 
 names_x <- abbreviate(colnames(X), 9, strict = T, dot = F) ; names_zu <-  abbreviate(colnames(Zu), 9, strict = T, dot = F);  names_zv <-  abbreviate(colnames(Zv), 9, strict = T, dot = F)
 
 if(distribution == "t") {names_del <- abbreviate(colnames(Zdel), 9, strict = T, dot = F); kdel = ncol(Zdel)} else {names_del <- NULL; kdel = 0}
 
 # starting values
 ols <- lm(y ~ 0 + X)
 beta0 <- as.matrix(coef(ols), nrow = length(coef(ols)), ncol = 1)
 olsResid <- y - X%*%beta0
 olsSkewness <- .skewness(olsResid)
  
 # Moment estimators for sigma_u/v squared
 if(is.null(start.val) == TRUE){
  m3 <- sum(olsResid^3)/n
  m2 <- sum(olsResid^2)/n
  su2init <- ((m3/(sqrt(2/pi)*(1-4/pi)))^2)^(1/3)
  if(is.nan(su2init) == T | su2init < 0) su2init <- 0.1
  sv2init <- m2 - (1 - 2/pi)*su2init
  if(is.nan(sv2init) == T | sv2init < 0) sv2init <- 0.1
  beta0[1] <- beta0[1] + sqrt(2/pi)*sqrt(su2init)
  y1 <- 0.5*log(((olsResid^2 - sv2init)/(1 - 2/pi))^2)
  reg_hetu <- lm(y1 ~ 0 + Zu)
  gu0 <- as.matrix(coef(reg_hetu), nrow = length(coef(reg_hetu)), ncol = 1)
  y2 <- 0.5*log((olsResid^2 - (1 - 2/pi)* su2init)^2)
  reg_hetv <- lm(y2 ~ 0 + Zv)
  gv0 <- as.matrix(coef(reg_hetv), nrow = length(coef(reg_hetv)), ncol = 1)
  if((is.null(tmean) == T) & (distribution == "t")){
   theta0 <- rbind( beta0, gv0, gu0, 0)
  } else if (is.null(tmean) == F) {
   theta0 <- rbind( beta0, gv0, gu0, as.matrix(rep(0, kdel)))
  } else {
   theta0 <- rbind( beta0, gv0, gu0)
  }
 } else {
  theta0 = start.val
 }
 rownames(theta0) = c(names_x, names_zv, names_zu, names_del)
 
 
   time.05 <- proc.time()
   if(print.level >= 2){
    # cat("__________________________________________________")
    cat("\nOptimization using 'mlmaximize': \n\n", sep = "")
   }
   obj <- eval(parse(text = paste(" tryCatch(.mlmaximize(theta0, .ll.",noquote(distribution),"n, gr = .g.",noquote(distribution),"n, hess = .hess.",noquote(distribution),"n, prod = prod, k = k, kv = kv, ku = ku, y = y, Zv = Zv, Zu = Zu, kdel = kdel, Zdel = Zdel, X = X, print.level = print.level), error = function(e) e )", sep = "")))
   if(inherits(obj, "error")){
     cat("\n=================")
     cat(" Trying to optimize using 'optim':\n\n", sep = "")
     cat("                  ('optim' minimizes, so 'll' is\n", sep = "")
     cat("                   negative of what is printed) \n\n", sep = "") 
     obj <- eval(parse(text = paste(" tryCatch( optim(par = theta0, fn = .ll.",noquote(distribution),"n, gr = .g.",noquote(distribution),"n, prod = prod, k = k, kv = kv, ku = ku, y = y, Zv = Zv, Zu = Zu, kdel = kdel, Zdel = Zdel, X = X, method = c('BFGS'), control = list(reltol  = .Machine$double.eps, maxit = 200, trace = ifelse(print.level >=2, 1, -1), REPORT = ifelse(print.level >=2, 1, -1), fnscale = -1), hessian = TRUE), error = function(e) e )", sep = "")))
     cannot.est.model <- inherits(obj, "error")
     if(obj$convergence == 0){
       obj$vcov = solve(-obj$hessian)
     } else if(cannot.est.model){
       # what is the error?
      # print(obj)
       warning("convergence was not reached or 'optim' cannot estimate the model at provided starting values")  
      } else {
       theta0 <- obj$par
       cat("\n")
       if(print.level >= 2){
         cat("\n__________________________________________________")
         cat("\nOptimization using 'mlmaximize': \n\n", sep = "")
       }
       obj <- eval(parse(text = paste(" tryCatch(.mlmaximize(theta0, .ll.",noquote(distribution),"n, gr = .g.",noquote(distribution),"n, hess = .hess.",noquote(distribution),"n, prod = prod, k = k, kv = kv, ku = ku, y = y, Zv = Zv, Zu = Zu,kdel = kdel, Zdel = Zdel, X = X, print.level = print.level), error = function(e) e )", sep = "")))
       cannot.est.model <- inherits(obj, "error")
       if ( cannot.est.model ){
        # what is the error?
        # print(obj)
         stop("Cannot estimate the model at provided starting values")
       }
     }}
   
   time.06 <- proc.time()
   est.time.sec <- (time.06-time.05)[1]
   names(est.time.sec) <- "sec"
 if(print.level >= 2){
  .timing(est.time.sec, "\nLog likelihood maximization completed in\n")
  cat("___________________________________________________\n")
 }
  
 # SE for sigmau and sigmav
 sigv <- sqrt(exp(obj$par[k+1]))
 sigu <- sqrt(exp(obj$par[k+kv+1]))
 se_sigv <- (0.5*sigv) * sqrt(obj$vcov[k+1,k+1])
 se_sigu <- (0.5*sigu) * sqrt(obj$vcov[k+kv+1,k+kv+1])
 
 # SE for lambda
 lmd <- sqrt(exp(obj$par[k+kv+1] - obj$par[k+1]))
 glmd <- c(0.5*lmd, -0.5*lmd)
 se_lmd <- as.vector( sqrt( t(glmd) %*% obj$vcov[c((k+1),(k+kv+1)), c((k+1),(k+kv+1))] %*% glmd))
 
 # SE for gamma
 gam <- exp(obj$par[k+kv+1])/(exp(obj$par[k+kv+1]) + exp(obj$par[k+1]))
 ggam <- c(gam*(1 - gam), -gam^2/lmd^2)
 se_gam <- as.vector( sqrt( t(ggam) %*% obj$vcov[c((k+1),(k+kv+1)), c((k+1),(k+kv+1))] %*% ggam)) 
 
 if((is.null(uhet) == T) & (is.null(vhet) == T)){
  par <- c(obj$par, sigv, sigu, lmd, gam); se <- c(sqrt(diag(obj$vcov)), se_sigv, se_sigu, se_lmd, se_gam)
  names(par) = c(names_x, names_zv, names_zu, names_del, "sigma_v  ", "sigma_u  ", "lambda   ", "gamma   ")  
 } else if((is.null(uhet) == T) & (is.null(vhet) == F)){
  par <- c(obj$par, sigu); se <- c(sqrt(diag(obj$vcov)), se_sigu)
  names(par) = c(names_x, names_zv, names_zu, names_del, "sigma_u  ") 
 } else if((is.null(uhet) == F) & (is.null(vhet) == T)){
  par <- c(obj$par, sigv); se <- c(sqrt(diag(obj$vcov)), se_sigv)
  names(par) = c(names_x, names_zv, names_zu, names_del, "sigma_v  ") 
 } else {
  par <- c(obj$par); se <- c(sqrt(diag(obj$vcov)))
  names(par) = c(names_x, names_zv, names_zu, names_del)
 }
 
 output <- cbind(round(par, digits = digits), round(se, digits = digits), round(par/se,digits = 2), round(pnorm(abs(par/se), lower.tail = FALSE)*2, digits = digits))
 colnames(output) <- c("Coef.", "SE ", "z ",  "P>|z|")
 
 # estimation results

 if(print.level >= 1){
  if(prod == TRUE){
   cat("\nCross-sectional stochastic (production) frontier model\n",  sep = "")} 
  else {
   cat("\nCross-sectional stochastic (cost) frontier model\n",  sep = "")
  }
  cat("\nDistributional assumptions\n\n", sep = "")
  Assumptions <- rep("heteroskedastic",2)
  if(kv==1){
   Assumptions[1] <- "homoskedastic"
  } 
  if(ku==1){
   Assumptions[2] <- "homoskedastic"
  }
  Distribution = c("normal ", "half-normal ")
  if(distribution == "t"){
   Distribution[2] <- "truncated-normal "
  }
  a1 <- data.frame(
   Component = c("Random noise: ","Inefficiency: "),
   Distribution = Distribution,
   Assumption = Assumptions
  )
  print(a1, quote = FALSE, right = F)
  cat("\nNumber of observations = ", n, "\n", sep = "")
  cat("\n--------------- Estimation results: --------------\n\n", sep = "")
  .printoutcs(output, digits = digits, k = k, kv = kv, ku = ku, kdel = kdel, na.print = "NA", dist = distribution)
 }
 
 # Technical efficiencies
 e <- y - X%*%obj$par[1:k,]
 sigmas_v <- sqrt(exp(Zv%*%obj$par[(k+1):(k+kv)])) 
 sigmas_u <- sqrt(exp(Zu%*%obj$par[(k+kv+1):(k+kv+ku)]))
 
 if(distribution == "h"){
  eff <- round(.u2efftnm(e, sigmas_u, sigmas_v, mu = 0, alpha = alpha, prod = prod), digits = digits)
  mu = NULL
 } else {
  mu <- Zdel%*%obj$par[-c(1:(k+kv+ku))]
  eff <- round(.u2efftnm(e, sigmas_u, sigmas_v, mu = mu, alpha = alpha, prod = prod), digits = digits)
 }
 
 # Marginal effects
 if(marg.eff == TRUE){
  meff = tryCatch(.me(obj$par[-c(1:(k+kv)),1,drop = FALSE], Zu = Zu, Zdel = Zdel, ku = ku, kdel = kdel, n = n, dist = distribution), error = function(e) { cat('In error handler\n'); print(e); e })
  if(inherits(meff, "error"))  {meff = NULL} else { meff = as.data.frame(meff)}
 } else {meff = NULL}
 
 myeff <- ifelse(prod, "technical", "cost")
 
 if(print.level >= 3){
  cat("\n=================")
  cat(" Summary of ",myeff," efficiencies, exp(-E[ui|ei]): \n\n", sep = "")
  .su1(eff[,2, drop = F])
  
  cat("\n=================")
  cat(" Summary of ",myeff," efficiencies, exp(-M[ui|ei]): \n\n", sep = "")
  .su1(eff[,3, drop = F])
  
  cat("\n=================")
  cat(" Summary of ",myeff," efficiencies, E[exp(-ui)|ei]: \n\n", sep = "")
  .su1(eff[,4, drop = F])
  cat("\n\n")
 }
 
 if(print.level == 4){
  cat("Point and interval estimates of unit-specific ",myeff," efficiencies: \n\n", sep = "")
  print(eff)
 }
 
 if(length(unique(sigmas_u)) == 1) sigmas_u = NULL
 if(length(unique(sigmas_v)) == 1) sigmas_v = NULL
 if(length(unique(mu)) == 1) mu = NULL
 
 if((is.null(uhet) == T) & (is.null(vhet) == T)){
 if (prod == T & olsSkewness > 0) {
  warning("OLS residuals are positively skewed, leading to zero estimate of inefficiency variance. Remaining ML estimates are equal to OLS. All units are estimated to be fully technically efficient.")
 } else if(prod == F & olsSkewness < 0) {
  warning("OLS residuals are negatively skewed, leading to zero estimate of inefficiency variance. Remaining ML estimates are equal to OLS. All units are estimated to be fully cost efficient.")
 }
 }
 
 colnames(obj$vcov) = rownames(obj$vcov) = c(names_x, names_zv, names_zu, names_del)
 
 temp = list(coef = par, vcov = obj$vcov, loglik = obj$ll, efficiencies = eff, marg.effects = meff, sigmas_u = sigmas_u, sigmas_v = sigmas_v, mu = mu, esample = esample)
 
 class(temp) <- "npsf"
 
 return(temp)
 
}
