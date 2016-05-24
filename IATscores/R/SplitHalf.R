SplitHalf <- function(IATdata, ...)
{
  
  # Split the data in two halves
  evenrows <- rep(c(FALSE, FALSE, TRUE, TRUE), length.out = nrow(IATdata))
  oddrows <- !evenrows
  IATdata1 <- IATdata[evenrows,]
  IATdata2 <- IATdata[oddrows,]
  
  # Compute the IAT scores for each half
  split1 <- RobustScores(IATdata1, ...)
  split2 <- RobustScores(IATdata2, ...)
  
  # compute the split-half correlations
  splithalf <- data.frame(
    matrix(ncol=2, nrow = 0, dimnames = list(c(), c("algorithm", "splithalf"))))
  algos <- names(split1)[names(split1)!= "subject"]
  for(i in 1:length(algos))
  {
    splitdata <- left_join(split1[,c("subject", algos[i])],
                           split2[,c("subject", algos[i])], by = "subject")
    
    # correlation
    splitcor <- cor(select(splitdata, -subject),
                    use = "pairwise.complete.obs")[1,2]
    
    # spearman-brown prophetic formula
    splithalf[i, "algorithm"] <- algos[i]
    splithalf[i, "splithalf"] <- 2 * splitcor / (1 + splitcor)
  }
  
  splithalf
}
