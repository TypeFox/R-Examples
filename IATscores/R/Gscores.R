Gscores <- function(IATlong)
{
  # Comments in "quotes" from Nosek, Bar-Anan, Sriram, Greenwald (2013),
  # "Understanding and using the brief implicit association test: I.
  # recommended scoring procedures". Table 9. http://ssrn.com/abstract=2196002
  
  Gaussianranks <- function(x)
  {
    # It handles NA values by leaving them in the same place as found
    y <- x[!is.na(x)]    
    
    # "1. Assign fractional ranks to N latencies. The longest latency will be
    # assigned a value of 1.0 and the shortest will be assigned a value of 1/N.
    # In the case of ties, ranks are averaged across tied values"
    N <- length(y)
    Fr <- rank(y)/N
    
    # "2. Subtract 1/2N from each fractional rank. Assuming untied values, the
    # largest latency will now have a value of 1-1/2N or (2N-1)/2N.
    # The 1/2N downward adjustment applies generally, even when tied values
    # exist".
    Fr <- Fr - 1/(2*N)
    
    # "3. For each of the N observations, compute the standard normal deviate
    # (mean = 0 and standard deviation = 1) corresponding to the adjusted
    # fractinal rank latency".
    Gr <- scale(Fr)
    
    # Remove the attributes given by the scale() function
    attr(Gr, "scaled:center") <- NULL
    attr(Gr, "scaled:scale") <- NULL
    x[!is.na(x)] <- Gr
    x
  }
  
  # "4. G1 is the mean of the normal deviates in condition 1. G2 is the mean
  # of the normal deviates in condition 2"
  Mranks <- filter(IATlong, variable != "pxxxx") %>%
    group_by(subject, variable) %>% # here do NOT group by blockcode
    mutate(Gr = Gaussianranks(RT)) %>% # compute gaussian ranks
    group_by(subject, variable, blockcode) %>% # here do group by blockcode
    summarize(Mean = mean(Gr, na.rm = TRUE)) # mean rank in each block
  
  Mranks <- dcast(Mranks, subject*variable ~ blockcode, value.var = "Mean")
  
  # "5. G = G2-G1"
  Mranks <- mutate(Mranks, Gscore = pair2 - pair1)
  
  # Put data in wide format
  Gsc <- dcast(Mranks, subject ~ variable, value.var = "Gscore")
  names(Gsc) <- str_replace(names(Gsc), "xx", "2x")
  Gsc  
}