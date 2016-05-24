RobustScores <- function(IATdata,
                         P1 = c("none", "fxtrim", "fxwins", "trim10", "wins10",
                                "inve10"),
                         P2 = c("ignore", "exclude", "recode", "separate",
                                "recode600"),
                         P3 = c("dscore", "gscore", "wpr90", "minid",
                                "minid_t10", "minid_w10", "minid_i10"),
                         P4 = c("nodist", "dist"),
                         maxMemory = 1000, verbose = TRUE)
{
   
  mincor <- 3 # minimum number of correct responses with lat < k10
  # required to be included int he analyses
  
  # maximum allowed latency. Careful, if change here,
  # consider changing the followint fixed parameters also in function doP1P2
  k10 <- 10000
  upfxtrim <- 10000
  lofxtrim <- 400
  k10 <- min(k10, upfxtrim)
  
  # CHECK THE INPUT
  # column subject must be present and must be numeric
  if(!"subject" %in% names(IATdata)) 
    stop('Bad input IATdata: Column "subject" is missing') 
  if(!is.numeric(IATdata$subject))
    warning('Bad input IATdata: Column "subject" must be numeric',
            immediate. = TRUE)
  
  # column subject must be present and must be numeric
  if(!"latency" %in% names(IATdata)) 
  {
    stop('Bad input IATdata: Column "latency" is missing') 
  } else if(!is.numeric(IATdata$latency))
  {
    stop('Bad input IATdata: Column "latency" must be numeric')
  }
  
  # column correct must be present and binary numerical (0,1) or logical
  if(!"correct" %in% names(IATdata)) 
  {
    stop('Bad input IATdata: Column "correct" is missing') 
  } else if(!is.logical(IATdata$correct) &! is.numeric(IATdata$correct))
  {
    stop('Bad input IATdata: Column "correct" must be logical or binary')
  } else if (is.numeric(IATdata$correct) & !all(IATdata$correct %in% c(0,1)))
  {
    stop('Bad input IATdata: Column "correct" must be logical or binary')
  }
  
  # column blockcode must be present and include only values pair1 and pair2
  if(!"blockcode" %in% names(IATdata)) 
  {
    stop('Bad input IATdata: Column "blockcode" is missing') 
  } 
  
  # praccrit is optional, however if absent, P4 can only be "nodist"
  if(!"praccrit" %in% names(IATdata) & ("dist" %in% P4))
  {
    P4 <- "nodist"
    warning('PARAMETER P4 HAS BEEN SET TO "nodist".
            Parameter P4 includes option "dist", distinction between practice
            and critical blocks. However column praccrit, which would allow
            to distinguish between practice and critical blocks, is not
            specified in the input IATdata.', immediate. = TRUE)
    IATdata$praccrit <- NA
  }
    
  # SELECT COLUMNS
  # drop any irrelevant (and potentially dangerous) column ...
  IATdata <- select(IATdata, subject, latency, correct, blockcode, praccrit) 
  # ... and row
  IATdata <- filter(IATdata, blockcode == "pair1" | blockcode == "pair2")
  # define a useful univocal index by row
  IATdata$index <- 1:nrow(IATdata)
  
  # Exclude participants with less than 3 correct valid latencies (< 10s and
  # > 400ms), in each block
  ncor <- group_by(IATdata, subject, blockcode) %>%
   summarize(ncor = sum(!is.na(correct) & correct == TRUE &
                        !is.na(latency) & latency < k10 &
                          latency >= lofxtrim)) %>%
    filter(ncor < mincor)
  
  if(nrow(ncor) != 0)
  {
    IATdata <- filter(IATdata, !subject %in% ncor$subject)
    warning(paste("The following subjects have been removed because they
            have too few correct responses to compute IAT scores, i.e.,
            less than", mincor, "correct responses with latency less than",
                  k10, "ms and more than", lofxtrim, "ms in at least one block:
                  Subjects =", str_join(ncor$subject, collapse = ", ")),
            immediate. = TRUE)
  }
  
  # COMPUTE THE ROBUST IAT SCORES
  Scores <- doP1P2P3P4(IATdata, P1 = P1, P2 = P2, P3 = P3, P4 = P4,
                       maxMemory = maxMemory, verbose = verbose)
 
  if(verbose) print(paste0(Sys.time(), ": IAT scores have been computed"))
  Scores
}