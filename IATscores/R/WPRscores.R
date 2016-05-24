WPRscores <- function(IATlong)
{
  qnt <- .9 # the quantile for the WPR score
  
  # Compute the Quantiles on all the variants of treated RTs
  Quant <- group_by(IATlong, subject, variable, blockcode) %>%
    filter(variable != "pxxxx") %>%
    summarize(Quant = quantile(RT, probs = qnt, na.rm = TRUE, names = FALSE))
  Quant <- dcast(Quant, subject*variable ~ blockcode, value.var = "Quant")
  
  # pooled-SD
  SDs <- group_by(IATlong, subject, variable) %>%
    filter(variable != "pxxxx") %>%
    summarize(SD = sd(RT, na.rm = TRUE))
  
  # Compute WPR scores
  WPRdata <- left_join(Quant, SDs, by = c("subject", "variable")) %>%
    mutate(WPR = (pair2-pair1)/SD)
  
  # Put data in wide format
  WPRsc <- dcast(WPRdata, subject ~ variable, value.var = "WPR")
  names(WPRsc) <- str_replace(names(WPRsc), "xx", "3x")
  WPRsc
}
