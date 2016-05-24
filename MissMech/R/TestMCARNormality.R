TestMCARNormality <- function(data, del.lesscases = 6, imputation.number = 1, method = "Auto", 
                              imputation.method = "Dist.Free", nrep = 10000, n.min = 30, 
                              seed = 110, alpha = 0.05, imputed.data = NA)
{
if (!is.na(seed))
 set.seed(seed) 
 if (is.data.frame(data)) {
 data <- as.matrix(data)
 } 
 if(!is.na(imputed.data[1]) && imputation.number!=1)
 {
   cat("Warning: No multiple imputation allowed when imputed data is provided.\n")
 }
 if(!is.matrix(data))
  {
    cat("Warning: Data is not a matrix or data frame.\n")
    stop("")
  }
 if(length(data)==0)
 {
   cat("Warning: Data is empty.\n")
   stop("")
 }
 if(ncol(data)<2)
 {
   cat("Warning: More than 1 variable is required.\n")
   stop("")
 }
 allempty <- which(apply(!is.na(data),1,sum) == 0)
   if (length(allempty) != 0) {
   data <- data[apply(!is.na(data), 1, sum) != 0, ]
   cat("Warning:", length(allempty), "Cases with all variables missing have been removed \n
          from the data.\n")
   }
 newdata <- OrderMissing(data, del.lesscases)
 if(length(newdata$data)==0)
 {
   cat("Warning: There are no data sets after deleting insufficient cases.\n")
   stop("")
 }
 
 if(newdata$g == 1)
 {
   cat("Warning: More than one missing data pattern should be present.\n")
   stop("")
 }
 if(sum(newdata$patcnt==1) > 0)
 {
   cat("Warning: At least 2 cases needed in each missing data patterns.\n")
   stop("")
 }

 y <- newdata$data
 patused <- newdata$patused
 patcnt <- newdata$patcnt
 spatcnt <- newdata$spatcnt
 caseorder <- newdata$caseorder
 removedcases <- newdata$removedcases
 n <- nrow(y)
 p <- ncol(y)
 g <- newdata$g
 spatcntz <- c(0, spatcnt)
 pvalsn <- matrix(0, imputation.number, g)
 adistar <- matrix(0, imputation.number, g)
 pnormality <- c()
 x <- vector("list", g)
 n4sim <- vector("list",g)
#------------------------------imputation-----------------------
 mu <- matrix(0, p, 1)
 sig <- diag(1, p)

    emest <- Mls(newdata, mu, sig, 1e-6)
    mu <- emest$mu
    sig <- emest$sig
 if(is.na(imputed.data[1]))
 {
    yimp <- y
    if (imputation.method == "Dist.Free") {
        iscomp <- apply(patused, 1, sum, na.rm = TRUE) == p
      
        cind <- which(iscomp)
        ncomp <- patcnt[cind]
        if (length(ncomp) == 0) ncomp <- 0
        use.normal <- FALSE
        if (ncomp >= 10 && ncomp>=2*p){
          compy <- y[seq(spatcntz[cind] + 1, spatcntz[cind + 1]), ]
          ybar <- matrix(apply(compy, 2, mean))
          sbar <- cov(compy)
          resid <- (ncomp / (ncomp - 1)) ^ .5 * 
                   (compy - matrix(ybar, ncomp, p, byrow = TRUE))
        } else {
          cat("Warning: There is not sufficient number of complete cases.\n  Dist.Free imputation requires a least 10 complete cases\n  or 2*number of variables, whichever is bigger.\n  imputation.method = normal will be used instead.\n")
          use.normal <- TRUE
        }
    }
    for(k in 1:imputation.number)
    {
        #-----------------normal imputation--------------------------------
        if (imputation.method == "Normal" || use.normal){
          yimp <- Impute(data = y, mu, sig, imputation.method = "Normal")
          yimp <- yimp$yimpOrdered
        }
        #-----------------distribution free imputation---------------------------------
        if (imputation.method == "Dist.Free" && !use.normal){
          yimp <- Impute(data = y, ybar, sbar, imputation.method = "Dist.Free", resid)
        yimp <- yimp$yimpOrdered
        }
        if (k == 1) yimptemp <- yimp
        #--------------Hawkin's test on the completed data------------------
        templist <- Hawkins(yimp,spatcnt)
        fij <- templist$fij
        tail <- templist$a
        ni <- templist$ni
        if (method == "Auto" || method == "Hawkins") {
          #Neyman test of uniformity for each group
           for(i in 1:g)
           {
             if (ni[i] < n.min){
                 if (k == 1) {
                     n4sim[[i]] <- SimNey(ni[i], nrep)
                 }
             }
             templist <- TestUNey(tail[[i]], nrep, sim = n4sim[[i]], n.min)
             pn <- templist$pn
             n4 <- templist$n4
             pn <- pn + (pn == 0) / nrep
             pvalsn[k,i] <- pn
           }
        }
 #--------------Anderson darling test for equality of distribution
        if (method == "Auto" || method == "Nonparametric") {
           if(length(ni)<2)
           {
             cat("Warning: Not enough groups for AndersonDarling test.")
             stop("")
           }
           templist <- AndersonDarling(fij, ni)
           p.ad <- templist$pn
           adistar[k, ] <- templist$adk.all
           pnormality <- c(pnormality, p.ad)
        }
    }
 } else {
    yimp <- imputed.data[caseorder, ]
    yimptemp <- yimp
    templist <- Hawkins(yimp,spatcnt)
    fij <- templist$fij
    tail <- templist$a
    ni <- templist$ni
    if (method == "Auto" || method == "Hawkins") {
      #Neyman test of uniformity for each group
       for(i in 1:g)
       {
          if (ni[i] < n.min){
              n4sim[[i]] <- SimNey(ni[i], nrep)
          }
          templist <- TestUNey(tail[[i]], nrep, sim = n4sim[[i]], n.min)
          pn <- templist$pn
          n4 <- templist$n4
          pn <- pn + (pn == 0) / nrep
          pvalsn[1, i] <- pn
       }
    }
    if (method == "Auto" || method == "Nonparametric") {
 #--------------Anderson darling test for equality of distribution
        templist <- AndersonDarling(fij, ni)
        p.ad <- templist$pn
        adistar[1, ] <- templist$adk.all
        pnormality <- c(pnormality, p.ad)
    }
 }
 adstar <- apply(adistar,1,sum)
 #combine p-values of test of uniformity
 combp <- -2 * apply(log(pvalsn), 1, sum)
 pvalcomb <- pchisq(combp, 2*g, lower.tail = FALSE)
 if (method == "Hawkins") {
    pnormality <- NULL
    adstar <- NULL
    adistar <- NULL
 }
 if (method == "Nonparametric") {
    pvalcomb = NULL
    combp = NULL
    pvalsn = NULL
 }
 yimptemp <- yimptemp[order(caseorder), ]
 if (length(removedcases) == 0) {
    dataused <- data
    }else {dataused <- data[-1 * removedcases, ]}
 homoscedastic <- list(analyzed.data = dataused, imputed.data = yimptemp,
     ordered.data =  y, caseorder = caseorder,
     pnormality = pnormality, adstar = adstar, adistar = adistar,  
     pvalcomb = pvalcomb, combp = combp, pvalsn = pvalsn, g = g, alpha = alpha,
     patused = patused, patcnt = patcnt, imputation.number = imputation.number, mu = mu, sigma = sig)
 homoscedastic$call <- match.call()
 class(homoscedastic) <- "testhomosc"
 homoscedastic
}
#---------------------------------------------------------------------
#testmcar <- function(x, ...) UseImputationMethod("testmcar")
#testmcar.default <- function(data, ncases = 6, imputation.number = 10,
#                             imputation.method = "Normal", nrep = 10000)
#{
#test <- TestMCARNormality(data, ncases = 6, imputation.number = 10,
#                         imputation.method = "Normal", nrep = 10000)
#test$call <- match.call()
#class(test) <- "testmcar"
#test
#}
#---------------------------------------------------------------------
# printing format for the class "testhomosc"
print.testhomosc <- function(x, ...) {
 cat("Call:\n")
 print(x$call)
 #cat("\nNumber of imputation:\n")
 #print(x$imputation.number)
 ni <- x$patcnt
 cat("\nNumber of Patterns: ", x$g,"\n\nTotal number of cases used in the analysis: ", sum(ni),"\n")
 cat("\n Pattern(s) used:\n")
 alpha <- x$alpha 
 disp.patt <- cbind(x$patused, ni)
 colnames(disp.patt)[ncol(disp.patt)] <- "Number of cases"
 rownames(disp.patt) <- rownames(disp.patt, do.NULL = FALSE, prefix = "group.")
 print(disp.patt, print.gap = 3) 
 method <- "Auto"
 if (is.null(x$pnormality)) method <- "Hawkins"
 if (is.null(x$pvalcomb)) method <- "Nonparametric"
 cat("\n\n    Test of normality and Homoscedasticity:\n  -------------------------------------------\n")
 if (method == "Auto") {
    cat("\nHawkins Test:\n")
    cat("\n    P-value for the Hawkins test of normality and homoscedasticity: ", x$pvalcomb[1],"\n")
    if (x$pvalcomb[1] > alpha){
       cat("\n    There is not sufficient evidence to reject normality
    or MCAR at", alpha,"significance level\n")
       }else {
       cat("\n    Either the test of multivariate normality or homoscedasticity (or both) is rejected.\n    Provided that normality can be assumed, the hypothesis of MCAR is 
    rejected at",alpha,"significance level. \n")
    cat("\nNon-Parametric Test:\n")
       cat("\n    P-value for the non-parametric test of homoscedasticity: ", x$pnormality[1],"\n")
       if (x$pnormality[1] > alpha){
          cat("\n    Reject Normality at",alpha,"significance level.
    There is not sufficient evidence to reject MCAR at",alpha,"significance level.\n")
          }else {
          cat("\n    Hypothesis of MCAR is rejected at ",alpha,"significance level.
    The multivariate normality test is inconclusive. \n")
          }
       }
 }
 if (method == "Hawkins"){
 cat("\nHawkins Test:\n")
 cat("\n    P-value for the Hawkins test of normality and homoscedasticity: ", x$pvalcomb[1],"\n")
 }
 if (method == "Nonparametric"){
  cat("\nNon-Parametric Test:\n")
  cat("\n    P-value for the non-parametric test of homoscedasticity: ", x$pnormality[1],"\n")
 }
}
#----------------------------------------------------------------------------
summary.testhomosc <- function(object, ...) {
 ni <- object$patcnt
 cat("\nNumber of imputation: ", object$imputation.number,"\n")
 cat("\nNumber of Patterns: ", object$g,"\n\nTotal number of cases used in the analysis: ", sum(ni),"\n")
 cat("\n Pattern(s) used:\n")
 alpha <- object$alpha 
 disp.patt <- cbind(object$patused, ni)
 colnames(disp.patt)[ncol(disp.patt)] <- "Number of cases"
 rownames(disp.patt) <- rownames(disp.patt, do.NULL = FALSE, prefix = "group.")
 print(disp.patt, print.gap = 3) 
 method <- "Auto"
 if (is.null(object$pnormality)) method <- "Hawkins"
 if (is.null(object$pvalcomb)) method <- "Nonparametric"
 cat("\n\n    Test of normality and Homoscedasticity:\n  -------------------------------------------\n")
 if (method == "Auto") {
    cat("\nHawkins Test:\n")
    cat("\n    P-value for the Hawkins test of normality and homoscedasticity: ", object$pvalcomb[1],"\n")
    cat("\nNon-Parametric Test:\n")
    cat("\n    P-value for the non-parametric test of homoscedasticity: ", object$pnormality[1],"\n")
 }
 if (method == "Hawkins"){
 cat("\nHawkins Test:\n")
 cat("\n    P-value for the Hawkins test of normality and homoscedasticity: ", object$pvalcomb[1],"\n")
 }
 if (method == "Nonparametric"){
  cat("\nNon-Parametric Test:\n")
  cat("\n    P-value for the non-parametric test of homoscedasticity: ", object$pnormality[1],"\n")
 }
}
#-----------------------------------------------------------------------------
# Plot "testhomosc"
boxplot.testhomosc <- function(x, ...) {
 if (is.null(x$pnormality)) {
    par(bg = "cornsilk")
    boxplot(x$pvalsn, col="lightcyan", border = "blue", medlwd = .5, medcol = "red")
    title(main = "Boxplots of p-values corresponding to each set of the missing data patterns\n for the Neyman test of Uniformity",
         xlab = "Missing data pattern group", ylab = "P-value", font.main = 4, 
         col.main = "blue4", cex.main = 1, font.lab = 4, cex.lab = 0.8, 
         col.lab = "blue4")
    abline(h = x$alpha / x$g, col = "red", lty = 2)
 }
 if (is.null(x$pvalcomb)) {
    par(bg = "cornsilk")
    boxplot(x$adistar, col="lightcyan", border = "blue", medlwd = .5, medcol = "red")
    title(main = "Boxplots of the T-value test statistics corresponding to each set of missing\n data patterns for the non-parametric test",
         xlab = "Missing data pattern group", ylab = expression(T[i]), 
         font.main = 4, col.main = "blue4", cex.main = 1, font.lab = 4, 
         cex.lab = 0.8, col.lab = "blue4")
 }
 if (!is.null(x$pvalcomb) && !is.null(x$pnormality)) {
    par(mfrow=c(2,1), bg = "cornsilk")
    boxplot(x$pvalsn, col="lightcyan", border = "blue", medlwd = .5, medcol = "red")
    title(main = "Boxplots of p-values corresponding to each set of the missing data patterns\n for the Neyman test of Uniformity",
         xlab = "Missing data pattern group", ylab = "P-value", font.main = 4, 
         col.main = "blue4", cex.main = 1, font.lab = 4, cex.lab = 0.8, 
         col.lab = "blue4")
    abline(h = x$alpha / x$g, col = "red", lty = 2)
    boxplot(x$adistar, col="lightcyan", border = "blue", medlwd = .5, medcol = "red")
    title(main = "Boxplots of the T-value test statistics corresponding to each set of missing\n data patterns for the non-parametric test",
         xlab = "Missing data pattern group", ylab = expression(T[i]), 
         font.main = 4, col.main = "blue4", cex.main = 1, font.lab = 4, 
         cex.lab = 0.8, col.lab = "blue4")
 }
}