alg2param <- function(x)
{
  # convert the name of an algorithm, a string in the form pxxxx,  or a vector
  # of algorithm names into a set of four parameters.
  
  # the four parameters, the order of the arguments is important and must be
  # the same as in the help file of RobustScore and in documentation
  P1 <- c("none", "fxtrim", "fxwins", "trim10", "wins10", "inve10")
  P2 <- c("ignore", "exclude", "recode", "separate", "recode600")
  P3 <- c("dscore", "gscore", "wpr90", "minid", "minid_t10", "minid_w10",
          "minid_i10")
  P4 <- c("nodist", "dist")
  
  p1 <- P1[as.numeric(str_sub(x, 2, 2))]
  p2 <- P2[as.numeric(str_sub(x, 3, 3))]
  p3 <- P3[as.numeric(str_sub(x, 4, 4))]
  p4 <- P4[as.numeric(str_sub(x, 5, 5))]
  
  data.frame(cbind("algorithm" = x, "P1" = p1, "P2" = p2, "P3" = p3, "P4" = p4),
             stringsAsFactors = FALSE)
}
