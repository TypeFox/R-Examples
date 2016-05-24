pmlr <- function(formula, data, weights = NULL, alpha = 0.05, penalized = TRUE, method = c("likelihood", "wald")[1], joint = FALSE) {

# change this!
if (penalized) useAstar.wald <- TRUE else useAstar.wald <- FALSE
useAstar.LRCI <- FALSE
   
   call <- match.call()
   mf <- match.call(expand.dots = FALSE)
   m <- match(c("formula", "data", "weights"), names(mf), 0)
   mf <- mf[c(1, m)]
   mf$drop.unused.levels <- TRUE
   mf[[1]] <- as.name("model.frame")
   mf <- eval(mf, parent.frame())
   mt <- attr(mf, "terms")

   y <- model.response(mf); if (is.factor(y)) y <- indicator.matrix(y)  
   x <- model.matrix(mt, mf)
   wt <- as.vector(model.weights(mf)); if (is.null(wt)) wt <- rep(1, times = dim(y)[1])
   
   p <- ncol(x) # p covariates
   J <- ncol(y) # J categories

   if (penalized) {
      fit <- getPLEs(x, y, wt); B <- fit$B; B.inf <- NULL; Ainv <- fit$Ainv; Astarinv <- fit$Astarinv; l.max <- fit$lstar.max
   } else {
      fit <- getMLEs(x, y, wt); B <- fit$B; B.inf <- fit$B.inf; Ainv <- fit$Ainv; l.max <- fit$l.max
   }
   
   separation <- array(data = NA, dim = c(1,p,J))
   dimnames(separation)<-list("", colnames(x), colnames(y))
   if (!penalized) {
      for (i in 1:J) {
         separation[1,,i] <- t(is.infinite(B.inf))[i,]
         separation[1,,i] <- ifelse(separation[1,,i], t(B.inf)[i,], NA)
      }
   } 
   
   coef <- array(data = NA, dim = c(1,p,J))
   dimnames(coef)<-list("", colnames(x), colnames(y))
   for (i in 1:J) coef[1,,i] <- t(B)[i,]
   
   if (useAstar.wald) Ainv <- Astarinv
   var <- array(data = NA, dim = c(p,p,J))
   dimnames(var)<-list(colnames(x), colnames(x), colnames(y))
   for (i in 1:J) var[,,i] <- Ainv[seq(from = i, by = J, length = p), seq(from = i, by = J, length = p)]

   stat <- pval <- CI.lower <- CI.upper <- array(data = NA, dim = c(1,p,J))
   dimnames(stat) <- dimnames(pval) <- dimnames(CI.lower) <- dimnames(CI.upper) <- list("", colnames(x), colnames(y)) 
   stat.joint <- pval.joint <- array(data = NA, dim = c(1,p,2))                                         
   dimnames(stat.joint) <- dimnames(pval.joint) <- list("", colnames(x), c("H_0: b_{i,1} = b_{i,2} = ... = b_{i,J} = 0","H_0: b_{i,1} = b_{i,2} = ... = b_{i,J}"))

   if (method == "likelihood") {                                        
# Likelihood ratio test for H_0: b_{i,j} = 0 (ith covariate, jth category)   
      test <- test.LR(x, y, wt, B, h0 = 1, penalized); statistic <- test$statistic; pvalue <- test$pvalue
# Profile confidence intervals
      profileCI.lower <- profileCIs(x, y, wt, B, B.inf, side = -1, alpha, l.max, step = 0.05, useAstar = useAstar.LRCI, penalized)
      profileCI.upper <- profileCIs(x, y, wt, B, B.inf, side = 1, alpha, l.max, step = 0.05, useAstar = useAstar.LRCI, penalized)
      for (i in 1:J) {
         stat[1,,i] <- t(statistic)[i,]
         pval[1,,i] <- t(pvalue)[i,]
         CI.lower[1,,i] <- t(profileCI.lower)[i,]
         CI.upper[1,,i] <- t(profileCI.upper)[i,]
      }
      if (joint) {
# Likelihood ratio test for H_0: b_{i,1} = b_{i,2} = ... = b_{i,J} = 0 (ith covariate)         
         test <- test.LR(x, y, wt, B, h0 = 2, penalized); statistic <- test$statistic; pvalue <- test$pvalue
         stat.joint[1,,1] <- t(statistic)
         pval.joint[1,,1] <- t(pvalue)
# Likelihood ratio test for H_0: b_{i,1} = b_{i,2} = ... = b_{i,J} (ith covariate)      
         test <- test.LR(x, y, wt, B, h0 = 3, penalized); statistic <- test$statistic; pvalue <- test$pvalue
         stat.joint[1,,2] <- t(statistic)
         pval.joint[1,,2] <- t(pvalue)
      }
   }

   if (method == "wald") {
# Wald test for H_0: b_{i,j} = 0 (ith covariate, jth category)   
      for (i in 1:J) {
         stat[,,i] <- (coef[,,i]/sqrt(diag(var[,,i])))^2
         pval[,,i] <- pchisq(stat[,,i], df = 1, lower.tail = FALSE)
         CI.lower[,,i] <- coef[,,i] - qnorm(p = 1 - alpha/2) * sqrt(diag(var[,,i]))
         CI.upper[,,i] <- coef[,,i] + qnorm(p = 1 - alpha/2) * sqrt(diag(var[,,i]))
      }
      if (joint) {
# Wald test for H_0: b_{i,1} = b_{i,2} = ... = b_{i,J} = 0 (ith covariate)
         for (i in 1:p) {
            C <- matrix(data = 0, nrow = J, ncol = p * J); C[1:J,(((i - 1) * J) + 1):(i * J)] <- diag(J)
            stat.joint[1,i,1] <- test.wald(vec(t(B)), Ainv, C) 
            pval.joint[1,i,1] <- pchisq(stat.joint[1,i,1], df = J, lower.tail = FALSE)   
         }
# Wald test for H_0: b_{i,1} = b_{i,2} = ... = b_{i,J} for the ith covariate   
         for (i in 1:p) {
            C <- matrix(data = 0, nrow = J-1, ncol = p *J); C[1:(J - 1), ((i - 1) * J + 1):((i * J) - 1)] <- diag(J - 1); for (j in (1:J - 1)) C[j,((i - 1) * J) + 1 + j] <- (-1)
            stat.joint[1,i,2] <- test.wald(vec(t(B)), Ainv, C)
            pval.joint[1,i,2] <- pchisq(stat.joint[1,i,2], df = J - 1, lower.tail = FALSE)
         }     
      }     
   }
   
   fit <- list(coefficients = coef, var = var, stat = stat, pval = pval, CI.lower = CI.lower, CI.upper = CI.upper, separation = separation, stat.joint = stat.joint, pval.joint = pval.joint, call = call, method = method, joint = joint)
   attr(fit, "class") <- c("pmlr")   
   fit

}
