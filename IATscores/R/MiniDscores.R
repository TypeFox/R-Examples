MiniDscores <- function(IATlong, P3 = c("minid", "minid_t10",
                                        "minid_w10", "minid_i10"),
                        maxMemory = 1000)
{
  # Compute all the mini differences between any latency in pair2 and
  # any latency in pair 1. Several kind of scores are then computed by function
  # computeMinid.
  
  # The following code is mainly data manipulation, to get a long format
  # dataframe which includes all the mini differences.
  # Some effort has been made to avoid RAM overflow.
  
  IATlong <- filter(IATlong, variable != "pxxxx")
  
  # Create two dataframes, IL1 and IL2. IL1 includes the data from block
  # pair1 and IL2 from block pair2. NAs are excluded.
  IL1 <- filter(IATlong, blockcode == "pair1") %>%
    select(subject, variable, "pair1" =  RT) %>%
    filter(!is.na(pair1))
  IL2 <- filter(IATlong, blockcode == "pair2") %>%
    select(subject, variable, "pair2" = RT) %>%
    filter(!is.na(pair2))
  rm(IATlong)
  
  # IL3 includes all the combinations of all the RTs in IL1 and IL2
  # for each values of subject and variable. This can impact on memory
  # and determine an overflow. To prevent this, I slice the
  # dataframes IL1 and IL2, so that the size of IL3 is never bigger than
  # parameter maxMemory.
  
  # 1. predict how large it will be IL3, in MB
  nIL1 <- group_by(IL1, subject, variable) %>% summarize(n1 = n())
  nIL2 <- group_by(IL2, subject, variable) %>% summarize(n2 = n())
  nIL <- left_join(nIL1, nIL2, by = c("subject", "variable"))
  nIL <- mutate(nIL, nIL3 = n1*n2)
  rows <- sum(nIL$nIL3, na.rm = TRUE) # the number of rows of IL3
  # estimated size of a single row of IL3
  arow <- (object.size(IL1) + object.size(IL1[,3]*2)) / nrow(IL1)
  rm(nIL1, nIL2, nIL)
  dimIL3 <- as.numeric(arow * rows * 10^-6) # predicted size of IL3, in MB
  
  nslices <- ceiling(dimIL3 / maxMemory)
  suppressWarnings( # the behavior of split is the one we desired, no warnings.
    slices <- split(unique(IL1$subject), 1:nslices)
  )
  
  # feed the slices to function computeMinid
  # and collect the output(s) in the dataframe Score_m
  for(i in 1:length(slices))
  {
    IL3 <- inner_join(filter(IL1, subject %in% slices[[i]]),
                      filter(IL2, subject %in% slices[[i]]),
                      by = c("subject", "variable"))
    IL3 <- mutate(IL3, minid = pair2 - pair1)
    IL3$pair1 <- NULL
    IL3$pair2 <- NULL
    newminid <- computeMinid(IL3, P3)
    
    if(i == 1)
    {
      Score_m <- newminid
    } else {
      Score_m <- rbind(Score_m, newminid)
    }
  }
  Score_m
}



computeMinid <- function(IL3, P3)
{
  # The output variable that will include all the minid scores
  Score_m <- data.frame("subject" = unique(IL3$subject))
  
  if("minid" %in% P3)
  {
    Msc <- group_by(IL3, subject, variable) %>%
      summarize(
        minid = mean(minid) / sd(minid)
      )
    Msc <- dcast(Msc, subject ~ variable, value.var = "minid")
    names(Msc) <- str_replace(names(Msc), "xx", "4x")
    Score_m <- left_join(Score_m, Msc, by = "subject")
  }
  
  if("minid_t10" %in% P3)
  {
    Msc <- group_by(IL3, subject, variable) %>%
      summarize(
        minid = specialmean(minid, type = "trm", tr = .1) / sd(minid)
      )
    Msc <- dcast(Msc, subject ~ variable, value.var = "minid")
    names(Msc) <- str_replace(names(Msc), "xx", "5x")
    Score_m <- left_join(Score_m, Msc, by = "subject")
  }
  
  if("minid_w10" %in% P3)
  {
    Msc <- group_by(IL3, subject, variable) %>%
      summarize(
        minid = specialmean(minid, type = "wns", tr = .1) / sd(minid)
      )
    Msc <- dcast(Msc, subject ~ variable, value.var = "minid")
    names(Msc) <- str_replace(names(Msc), "xx", "6x")
    Score_m <- left_join(Score_m, Msc, by = "subject")
  }
  
  if("minid_i10" %in% P3)
  {
    Msc <- group_by(IL3, subject, variable) %>%
      summarize(
        minid = specialmean(minid, type = "inv", tr = .1) / sd(minid)
      )
    Msc <- dcast(Msc, subject ~ variable, value.var = "minid")
    names(Msc) <- str_replace(names(Msc), "xx", "7x")
    Score_m <- left_join(Score_m, Msc, by = "subject")
  }
  Score_m  
}