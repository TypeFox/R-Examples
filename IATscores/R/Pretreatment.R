Pretreatment <- function(IATdata,
                         label_subject = "subject",
                         label_latency = "latency",
                         label_accuracy = "correct",
                         label_block = "blockcode",
                         block_pair1 = c("pair1_left", "pair1_right"),
                         block_pair2 = c("pair2_left", "pair2_right"),
                         label_trial = NA,
                         trial_left = NA,
                         trial_right = NA,
                         label_praccrit = NA,
                         block_prac = NA,
                         block_crit = NA,
                         label_stimulus = NA)
{
  # determine which optional columns are desired
  trials <- !any(is.na(c(label_trial, trial_left, trial_right)))
  praccrit <- !any(is.na(c(label_praccrit, block_prac, block_crit)))
  stimuli <- !is.na(label_stimulus)
  
  #IAT1: subset of IATdata, only the relevant columns are kept, a subset
  # of these: 1)subject, 2)latency, 3)accuracy, 4)block label, 5)trial label,
  # 6) praccrit, 7) stimulus
  cols <- c(label_subject, label_latency, label_accuracy, label_block)
  # columns 5, 6, and 7 are optional
  if(trials) cols <- c(cols, label_trial)
  if(stimuli) cols <- c(cols, label_stimulus)
  if(praccrit) cols <- c(cols, label_praccrit)
  
  IAT1 <- IATdata[,cols]
  rm(IATdata)
  
  # Keep only critical blocks
  IAT1 <- subset(IAT1, (IAT1[,label_block]%in%c(block_pair1, block_pair2)))
  
  # Convert block names to simpler labels
  IAT2 <- data.frame(matrix(ncol=0, nrow=nrow(IAT1)))
  IAT2$subject <- IAT1[,label_subject]
  IAT2$latency <- IAT1[,label_latency]
  IAT2$correct <- IAT1[, label_accuracy]
  IAT2[IAT1[,label_block] %in% block_pair1, "blockcode"] <- "pair1"
  IAT2[IAT1[,label_block] %in% block_pair2, "blockcode"] <- "pair2"
  IAT2$blockcode <- as.factor(IAT2$blockcode)
  
  # optional: if trial codes are specified, convert them to simpler labels
  if(trials)
  {
    IAT2[IAT1[,label_trial] %in% trial_right, "trialcode"] <- "right"
    IAT2[IAT1[,label_trial] %in% trial_left, "trialcode"] <- "left"
    IAT2$trialcode<-as.factor(IAT2$trialcode)
  }

  # optional: if practice / critical double categorization blocks must be
  # distinguished, convert them to simpler labels
  if(praccrit)
  {
    IAT2[IAT1[,label_praccrit] %in% block_prac,"praccrit"] <- "prac"
    IAT2[IAT1[,label_praccrit] %in% block_crit,"praccrit"] <- "crit"
    IAT2[IAT2$blockcode=="practice","praccrit"] <- "null"
    IAT2$praccrit<-as.factor(IAT2$praccrit)
  }
  
  # optional: do you want to keep the stimuli column?
  if(stimuli)
    IAT2$stimulus <- str_trim(IAT1[,label_stimulus])
 
  # clean the NA rows, if any
  IAT2 <- remove.na.rows(IAT2)
  
  # values that are not allowed are removed
  # e.g., (negative latencies, correct != 0 or 1)
  IAT2 <- IAT2 [IAT2$latency >= 0 & IAT2$correct %in% c(0,1), ]
  IAT2
}
