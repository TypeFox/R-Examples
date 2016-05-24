doP1P2 <- function(IATdata,
                   P1 = c("none", "fxtrim", "fxwins", "trim10", "wins10",
                          "inve10"),
                   P2 = c("ignore", "exclude", "recode", "separate",
                          "recode600"),
                   verbose = TRUE)
{
  if(verbose) print(paste0(Sys.time(), ": Applying parameters P1 and P2"))
  
  # Apply parameters 1 and 2 and return new columns in the dataframe
  
  k10 <- 10000 # maximum allowed latency. Careful, if change here,
  # consider changing k10 also in function RobustScores
  lofxtrim <- 400 # lower bound for fixed trimminig
  upfxtrim <- 10000 # upper bound for fixed trimming
  lofxwins <- 300 # lower bound for fixed winsorizing
  upfxwins <- 3000 # upper bound for fixed winsorizing
  tr <- .1 # amount of trimming / winsorizing
  k <- 2 # if P2 has option replace, errors are replaced by M+k*SD of correct
  # latencies.
  
  # Exclude all latencies >10 sec
  rm10s <- function(x)
  {
    x[x > k10]  <- NA
    x
  }
  IATdata <- mutate(IATdata, pxxxx = rm10s(latency))
  
  # P2 = 1. ERROR TREATMENT = IGNORE.
  # all possible values of P1
  if("ignore" %in% P2)
  {
    if("none" %in% P1)
    {
      IATdata <- mutate(IATdata, p11xx = pxxxx)
    }  
    if("fxtrim" %in% P1)
    {
      IATdata <- mutate(IATdata,
                        p21xx = fxtrim(pxxxx, lo = lofxtrim, up = upfxtrim))
    }  
    if("fxwins" %in% P1)
    {
      IATdata <- mutate(IATdata,
                        p31xx = fxwins(pxxxx, lo = lofxwins, up = upfxwins))
    }  
    if("trim10" %in% P1)
    {
      IATdata <- group_by(IATdata, subject) %>%
        mutate(p41xx = trim_or_win(pxxxx, type = "trm", tr = tr))
    }  
    
    if("wins10" %in% P1)
    {
      IATdata <- group_by(IATdata, subject) %>%
        mutate(p51xx = trim_or_win(pxxxx, type = "wns", tr = tr))
    }  
    if("inve10" %in% P1)
    {
      IATdata <- group_by(IATdata, subject) %>%
        mutate(p61xx = trim_or_win(pxxxx, type = "inv", tr = tr))
    }  
  }
  
  # P2 = 2. ERRORTREATMENT = EXCLUDE 
  if("exclude" %in% P2)
  {
    
    if("none" %in% P1)
    {
      tmp1 <- filter(IATdata, correct == 1) %>% mutate(p12xx = pxxxx)
      tmp2 <- filter(IATdata, correct == 0) %>% mutate(p12xx = NA)
      IATdata <- rbind_list(tmp1, tmp2) %>% arrange(index)
    }  
    if("fxtrim" %in% P1)
    {
      tmp1 <- filter(IATdata, correct == 1) %>%
        mutate(p22xx = fxtrim(pxxxx, lo = lofxtrim, up = upfxtrim))
      tmp2 <- filter(IATdata, correct == 0) %>% mutate(p22xx = NA)
      IATdata <- rbind_list(tmp1, tmp2) %>% arrange(index)
    }  
    if("fxwins" %in% P1)
    {
      tmp1 <- filter(IATdata, correct == 1) %>%
        mutate(p32xx = fxwins(pxxxx, lo = lofxwins, up = upfxwins))
      tmp2 <- filter(IATdata, correct == 0) %>% mutate(p32xx = NA)
      IATdata <- rbind_list(tmp1, tmp2) %>% arrange(index)
    }  
    if("trim10" %in% P1)
    {
      tmp1 <- filter(IATdata, correct == 1) %>% group_by(subject) %>%
        mutate(p42xx = trim_or_win(pxxxx, type = "trm", tr = tr))
      tmp2 <- filter(IATdata, correct == 0) %>% mutate(p42xx = NA)
      IATdata <- rbind_list(tmp1, tmp2) %>% arrange(index)
    }  
    
    if("wins10" %in% P1)
    {
      tmp1 <- filter(IATdata, correct == 1) %>% group_by(subject) %>%
        mutate(p52xx = trim_or_win(pxxxx, type = "wns", tr = tr))
      tmp2 <- filter(IATdata, correct == 0) %>% mutate(p52xx = NA)
      IATdata <- rbind_list(tmp1, tmp2) %>% arrange(index)
    }  
    if("inve10" %in% P1)
    {
      tmp1 <- filter(IATdata, correct == 1) %>% group_by(subject) %>%
        mutate(p62xx = trim_or_win(pxxxx, type = "inv", tr = tr))
      tmp2 <- filter(IATdata, correct == 0) %>% mutate(p62xx = NA)
      IATdata <- rbind_list(tmp1, tmp2) %>% arrange(index)
    }  
  }
  
  # P2 = 3. ERROR TREATMENT = RECODE with M+2SD
  if("recode" %in% P2)
  {
    ErrReplace <- filter(IATdata, correct == TRUE) %>% # only on correct resp.
      group_by(subject, blockcode) %>% # for each subject and each block
      summarize(ErrReplace = mean(pxxxx, na.rm = TRUE) +
               + k*sd(pxxxx, na.rm = TRUE)) # compute M+2SD of latencies
    
    IATdata <- left_join(IATdata, ErrReplace, by = c("subject", "blockcode"))
    
    if("none" %in% P1)
    {
      tmp1 <- filter(IATdata, correct == 1) %>% mutate(p13xx = pxxxx)
      tmp2 <- filter(IATdata, correct == 0) %>% mutate(p13xx = ErrReplace)
      IATdata <- rbind_list(tmp1, tmp2) %>% arrange(index)
    }  
    if("fxtrim" %in% P1)
    {
      tmp1 <- filter(IATdata, correct == 1) %>%
        mutate(p23xx = fxtrim(pxxxx, lo = lofxtrim, up = upfxtrim))
      tmp2 <- filter(IATdata, correct == 0) %>% mutate(p23xx = ErrReplace)
      IATdata <- rbind_list(tmp1, tmp2) %>% arrange(index)
    }  
    if("fxwins" %in% P1)
    {
      tmp1 <- filter(IATdata, correct == 1) %>%
        mutate(p33xx = fxwins(pxxxx, lo = lofxwins, up = upfxwins))
      tmp2 <- filter(IATdata, correct == 0) %>% mutate(p33xx = ErrReplace)
      IATdata <- rbind_list(tmp1, tmp2) %>% arrange(index)
    }  
    if("trim10" %in% P1)
    {
      tmp1 <- filter(IATdata, correct == 1) %>% group_by(subject) %>%
        mutate(p43xx = trim_or_win(pxxxx, type = "trm", tr = tr))
      tmp2 <- filter(IATdata, correct == 0) %>% mutate(p43xx = ErrReplace)
      IATdata <- rbind_list(tmp1, tmp2) %>% arrange(index)
    }  
    
    if("wins10" %in% P1)
    {
      tmp1 <- filter(IATdata, correct == 1) %>% group_by(subject) %>%
        mutate(p53xx = trim_or_win(pxxxx, type = "wns", tr = tr))
      tmp2 <- filter(IATdata, correct == 0) %>% mutate(p53xx = ErrReplace)
      IATdata <- rbind_list(tmp1, tmp2) %>% arrange(index)
    }  
    if("inve10" %in% P1)
    {
      tmp1 <- filter(IATdata, correct == 1) %>% group_by(subject) %>%
        mutate(p63xx = trim_or_win(pxxxx, type = "inv", tr = tr))
      tmp2 <- filter(IATdata, correct == 0) %>% mutate(p63xx = ErrReplace)
      IATdata <- rbind_list(tmp1, tmp2) %>% arrange(index)
    }  
    IATdata$ErrReplace <- NULL
  }
  
  # P2 = 4. ERROR TREATMENT = SEPARATE.
  if("separate" %in% P2)
  {
    if("none" %in% P1)
    {
      IATdata <- mutate(IATdata, p14xx = pxxxx)
    }  
    if("fxtrim" %in% P1)
    {
      IATdata <- mutate(IATdata,
                        p24xx = fxtrim(pxxxx, lo = lofxtrim, up = upfxtrim))
    }  
    if("fxwins" %in% P1)
    {
      IATdata <- mutate(IATdata,
                        p34xx = fxwins(pxxxx, lo = lofxwins, up = upfxwins))
    }  
    if("trim10" %in% P1)
    {
      tmp1 <- filter(IATdata, correct == 1) %>% group_by(subject) %>%
        mutate(p44xx = trim_or_win(pxxxx, type = "trm", tr = tr))
      tmp2 <- filter(IATdata, correct == 0) %>% group_by(subject) %>%
        mutate(p44xx = trim_or_win(pxxxx, type = "trm", tr = tr))
      IATdata <- rbind_list(tmp1, tmp2) %>% arrange(index)
    }  
    
    if("wins10" %in% P1)
    {
      tmp1 <- filter(IATdata, correct == 1) %>% group_by(subject) %>%
        mutate(p54xx = trim_or_win(pxxxx, type = "wns", tr = tr))
      tmp2 <- filter(IATdata, correct == 0) %>% group_by(subject) %>%
        mutate(p54xx = trim_or_win(pxxxx, type = "wns", tr = tr))
      IATdata <- rbind_list(tmp1, tmp2) %>% arrange(index)
    }  
    if("inve10" %in% P1)
    {
      tmp1 <- filter(IATdata, correct == 1) %>% group_by(subject) %>%
        mutate(p64xx = trim_or_win(pxxxx, type = "inv", tr = tr))
      tmp2 <- filter(IATdata, correct == 0) %>% group_by(subject) %>%
        mutate(p64xx = trim_or_win(pxxxx, type = "inv", tr = tr))
      IATdata <- rbind_list(tmp1, tmp2) %>% arrange(index)
    }  
  }
  
  # P2 = 5. ERROR TREATMENT = RECODE with M+600
  if("recode600" %in% P2)
  {
    ErrReplace <- filter(IATdata, correct == TRUE) %>% # only on correct resp.
      group_by(subject, blockcode) %>% # for each subject and each block
      summarize(ErrReplace = mean(pxxxx, na.rm = TRUE) + 600) # M+600
    
    IATdata <- left_join(IATdata, ErrReplace, by = c("subject", "blockcode"))
    
    if("none" %in% P1)
    {
      tmp1 <- filter(IATdata, correct == 1) %>% mutate(p15xx = pxxxx)
      tmp2 <- filter(IATdata, correct == 0) %>% mutate(p15xx = ErrReplace)
      IATdata <- rbind_list(tmp1, tmp2) %>% arrange(index)
    }  
    if("fxtrim" %in% P1)
    {
      tmp1 <- filter(IATdata, correct == 1) %>%
        mutate(p25xx = fxtrim(pxxxx, lo = lofxtrim, up = upfxtrim))
      tmp2 <- filter(IATdata, correct == 0) %>% mutate(p25xx = ErrReplace)
      IATdata <- rbind_list(tmp1, tmp2) %>% arrange(index)
    }  
    if("fxwins" %in% P1)
    {
      tmp1 <- filter(IATdata, correct == 1) %>%
        mutate(p35xx = fxwins(pxxxx, lo = lofxwins, up = upfxwins))
      tmp2 <- filter(IATdata, correct == 0) %>% mutate(p35xx = ErrReplace)
      IATdata <- rbind_list(tmp1, tmp2) %>% arrange(index)
    }  
    if("trim10" %in% P1)
    {
      tmp1 <- filter(IATdata, correct == 1) %>% group_by(subject) %>%
        mutate(p45xx = trim_or_win(pxxxx, type = "trm", tr = tr))
      tmp2 <- filter(IATdata, correct == 0) %>% mutate(p45xx = ErrReplace)
      IATdata <- rbind_list(tmp1, tmp2) %>% arrange(index)
    }  
    
    if("wins10" %in% P1)
    {
      tmp1 <- filter(IATdata, correct == 1) %>% group_by(subject) %>%
        mutate(p55xx = trim_or_win(pxxxx, type = "wns", tr = tr))
      tmp2 <- filter(IATdata, correct == 0) %>% mutate(p55xx = ErrReplace)
      IATdata <- rbind_list(tmp1, tmp2) %>% arrange(index)
    }  
    if("inve10" %in% P1)
    {
      tmp1 <- filter(IATdata, correct == 1) %>% group_by(subject) %>%
        mutate(p65xx = trim_or_win(pxxxx, type = "inv", tr = tr))
      tmp2 <- filter(IATdata, correct == 0) %>% mutate(p65xx = ErrReplace)
      IATdata <- rbind_list(tmp1, tmp2) %>% arrange(index)
    }  
    IATdata$ErrReplace <- NULL
  }
  
  IATdata
}