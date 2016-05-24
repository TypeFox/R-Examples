
is.testFriedman <- function(x) {
  is(x, "testFriedman")
}

print.testFriedman <- function (x, ...) {
  cat( sprintf("Friedman test, objetive %s output variable %s. Obtained p-value: %.4e\n",
               x$tags$objetive, x$tags$target, x$tags$pvalue))
}

summary.testFriedman <- function (x, ...) {
  cat( sprintf("Friedman test, objetive %s output variable %s. Obtained p-value: %.4e\n", 
               x$tags$objetive, x$tags$target, x$tags$pvalue))
  cat( sprintf("Chi squared with %d degrees of freedom statistic: %.4f\n", 
               x$tags$distribution, x$tags$statistic) )
  if(x$tags$pvalue < x$tags$alpha)
    cat( sprintf("Test rejected: p-value: %.4e < %.4f\n", 
                 x$tags$pvalue, x$tags$alpha) )
  else
    cat( sprintf("Test accepted: p-value: %.4e >= %.4f\n", 
                 x$tags$pvalue, x$tags$alpha) )
}

.testFriedman <- function(ranks, pvalue, distribution, statistic, objetive, alpha, target, tags) {

  title <- sprintf("Friedman Test for output \"%s\"", target)
  if (pvalue < alpha)
    outcome <- "Rejected"
  else
    outcome <- "Not Rejected"
  
  newTags <- .metaTags(pvalue       = pvalue,
                       distribution = distribution,
                       statistic    = statistic,
                       objetive     = objetive,
                       alpha        = alpha,
                       target       = target,
                       title        = title,
                       outcome      = outcome,
                       alias        = "FriedmanTest")
  
  tags <- .updateTags(tags, newTags)
  
  f <- list(
    "ranks"               = ranks,
    "tags"                = tags
  )
  class(f) <- c("testFriedman", "reportable")
  f
}

.doFriedmanTest <- function(e, output, rankOrder, alpha) {
  #Extracts the problem with only the interested output variable
  data <- e$data[,c(e$method,e$problem,output)]
  #Now transform the problem into a method-versus-problem one.
  data <- reshape2::dcast(data, paste(e$method,"~",e$problem),value.var = output)
  #First column is now the rownames (and then is removed).
  rownames(data) <- data[,1]
  #The parameter drop = FALSE avoid to auto-convert the data.frame into
  #a single array when only one column data.frame is the result of the subseting.
  data <- data[,-1, drop=FALSE]
  
  # Check if it only exists two methods (in that case, 
  # a wilcoxon test would be performed instead)
  if( nrow(data) < 3 )
    stop(.callErrorMessage("friedman2MethodsError"))
  
  fri                         <- .statisticsFriedman(data, rankOrder)
  pvalue                      <- fri$pvalue
  distribution                <- fri$degreesOfFreedom
  statistic                   <- fri$statistic
  ranks                       <- fri$ranks
  
  f <- .testFriedman(ranks  = ranks,
                     pvalue = pvalue,
                     distribution = distribution,
                     statistic = statistic,
                     objetive = .mapId2Name[[rankOrder]],
                     alpha = alpha,
                     target = output,
                     tags   = e$tags)
  
  f
}