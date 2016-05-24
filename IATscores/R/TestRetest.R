TestRetest <- function(IATdata, ...)
{
  
  sess <- unique(IATdata$session)
  test <- RobustScores(IATdata = filter(IATdata, session == sess[1]), ...)
  retest <- RobustScores(IATdata = filter(IATdata, session == sess[2]), ...)
  
  # compute the test-retest correlations
  testretest <-
    data.frame(matrix(ncol=2, nrow = 0,
                      dimnames = list(c(), c("algorithm", "testretest"))))
  algos <- names(test)[names(test) != "subject"]
  
  for(i in 1:length(algos))
  {
    splitdata <- left_join(test[,c("subject", algos[i])],
                           retest[,c("subject", algos[i])], by = "subject")
    
    # correlation
    splitcor <- cor(select(splitdata, -subject),
                    use = "pairwise.complete.obs")[1,2]
    
    # spearman-brown prophetic formula
    testretest[i, "testretest"] <- splitcor
    testretest[i, "algorithm"] <- algos[i]
  }
 
  testretest
}