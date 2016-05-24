gofIRT.ppar <- function(object, groups.hl = 10, cutpoint = 0.5)
{
#S3 method for computing 3 deviances and hosmer-lemeshow test
#object ... object of class ppar
#ngroups.hl ... number of percentile groups for Hosmer-Lemeshow Test

  if (max(object$X, na.rm = TRUE) > 1) stop("Tests for polytomous models not implemented yet!")
  if (any(is.na(object$X))) stop("Test for data with missings not implemented yet!")
   
  pi.hat <- pmat(object)
  groups.cldev <- "rawscore"

  #---------------- compute test statistics ----------------------------
  res.cl <- unlist(cldeviance(object, groups.gr = groups.cldev, pi.hat = pi.hat))
  res.hl <- unlist(hoslem(object, groups.hl = groups.hl, pi.hat = pi.hat))
  res.rost <- unlist(rostdeviance(object))
  res.cw <- unlist(cwdeviance(object, pi.hat))
  
  res.table <- rbind(res.cl, res.hl, res.rost, res.cw)
  colnames(res.table) <- c("value","df","p-value")
  rownames(res.table) <- c("Collapsed Deviance", "Hosmer-Lemeshow", "Rost Deviance", "Casewise Deviance")
  #------------------- end test statistics ----------------------------
  
  #---------------------- R-squared -----------------------------------
  res.r2 <- Rsquared(object, pi.hat = pi.hat)
  #---------------------- end R-squared -------------------------------
  
  #--------------------------- classifier stuff -----------------------
  pred.X <- predict(object, cutpoint = cutpoint)        #predicted data matrix
  observed <- as.vector(object$X.ex)
  predicted <- as.vector(pred.X)
  confmat <- table(predicted, observed)
  accuracy <- sum(diag(confmat))/sum(confmat)
  sens <- as.vector((confmat[2,2])/(colSums(confmat)[2]))
  spez <- as.vector((confmat[1,1])/(colSums(confmat)[1]))
  cl.list <- list(confmat = confmat, accuracy = accuracy, sensitivity = sens, specificity = spez)
  
  probvec <- as.vector(pi.hat)
  rocpr.res <- prediction(probvec[!is.na(probvec)], observed[!is.na(observed)])
  roc.res <- performance(rocpr.res, "tpr","fpr")                   #produce ROC output
  
  spezvec <- 1-(roc.res@x.values[[1]])         #vector of specificities (different cuts)
  sensvec <- roc.res@y.values[[1]]             #vector of sensitivities (different cuts)
  cutvec <- roc.res@alpha.values[[1]]          #vector with thresholds
  sscmat <- cbind(cutvec, sensvec - spezvec)[order(abs(sensvec-spezvec), decreasing = FALSE),]
  thresh.opt <- mean(sscmat[1:2,1])
   
  auc.all <- performance(rocpr.res, "auc")                      #area under ROC
  auc.res <- auc.all@y.values[[1]]
  gini <- (2*auc.res)-1
  
  #----------------------- end classifier ----------------------------------
 
  result <- list(test.table = res.table, R2 = res.r2, classifier = cl.list, AUC = auc.res, 
                 Gini = gini, ROC = roc.res, opt.cut = thresh.opt, predobj = rocpr.res)
  class(result) <- "gof"
  result
}
