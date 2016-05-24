doP1P2P3P4 <- function(IATdata,
                       P1 = c("none", "fxtrim", "fxwins", "trim10", "wins10",
                              "inve10"),
                       P2 = c("ignore", "exclude", "recode", "separate"),
                       P3 = c("dscore", "gscore", "wpr90", "minid", "minid_t10",
                              "minid_w10", "minid_i10"),
                       P4 = c("nodist", "dist"),
                       maxMemory = 1000, verbose = TRUE)
{
  Scores <- data.frame("subject" = unique(IATdata$subject))
  
  # nodist in P4: do not distinguish practice and critical blocks
  if("nodist" %in% P4)
  {
    if(verbose) print(paste0(Sys.time(), ": Applying parameter P4 = nodist"))
    IATdata2 <- doP1P2(IATdata, P1 = P1, P2 = P2, verbose = verbose)
    Ndist <- doP3(IATdata2, P3 = P3, maxMemory = maxMemory, verbose = verbose)
    names(Ndist) <- str_replace(names(Ndist), "x", "1")
    Scores <- left_join(Scores, Ndist, by = "subject")
  }
  
  # dist in P4: distinguish practice and critical blocks
  if("dist" %in% P4)
  {
    if(verbose) print(paste0(Sys.time(), ": Applying parameter P4 = dist"))
    # practice block
    IATdata_prac <- filter(IATdata, praccrit == "prac")
    IATdata_prac2 <- doP1P2(IATdata_prac, P1 = P1, P2 = P2, verbose = verbose)
    dist_prac <- doP3(IATdata_prac2, P3 = P3, maxMemory = maxMemory,
                      verbose = verbose)
    dist_prac$praccrit <- "prac"

    # critical block
    IATdata_crit <- filter(IATdata, praccrit == "crit")
    IATdata_crit2 <- doP1P2(IATdata_crit, P1 = P1, P2 = P2, verbose = verbose)
    dist_crit <- doP3(IATdata_crit2, P3 = P3, maxMemory = maxMemory,
                      verbose = verbose)
    dist_crit$praccrit <- "crit"
    
    # average practice and critical
    dist_l <- rbind(dist_prac, dist_crit)
    dist <- data.frame("subject" = unique(IATdata$subject))
    
    for(i in names(dist_l)[!names(dist_l) %in% c("subject", "praccrit")])
    {
      dist_w <- dcast(dist_l, subject ~ praccrit, value.var = i)
      dist_w <- mutate(dist_w, dscore = (prac+crit)/ 2) %>%
        select(-crit, -prac)
      names(dist_w) <- c("subject", str_replace(i, "x", 2))
      dist <- left_join(dist, dist_w, by = "subject") 
    }
    
    Scores <- left_join(Scores, dist, by = "subject")
  }
  Scores
}
