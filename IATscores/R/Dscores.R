Dscores <- function(IATlong)
{
  # Function to compute the Dscores from data in (a certain) long format
  
  # Compute the Means on all the variants of treated RTs
  Means <- group_by(IATlong, subject, variable, blockcode) %>%
    filter(variable != "pxxxx") %>%
    summarize(Mean = mean(RT, na.rm = TRUE))
  Means <- dcast(Means, subject*variable ~ blockcode, value.var = "Mean")
  
  # pooled-SD
  SDs <- group_by(IATlong, subject, variable) %>%
    filter(variable != "pxxxx") %>%
    summarize(SD = sd(RT, na.rm = TRUE))
  
  # Compute D scores
  Ddata <- left_join(Means, SDs, by = c("subject", "variable")) %>%
    mutate(Dscore = (pair2-pair1)/SD)
  
  # Put data in wide format
  Dsc <- dcast(Ddata, subject ~ variable, value.var = "Dscore")
  names(Dsc) <- str_replace(names(Dsc), "xx", "1x")
  Dsc
}