`modelChooseAdd` <-
function(modelShow, criteria="LL") {
 nSubjects      <- dim(modelShow)[1]/8
 mChoose        <- modelChoose(modelShow, criteria=criteria)
 res1           <- data.frame(matrix(NA, nrow=nSubjects,   ncol=2))
 res2           <- data.frame(matrix(NA, nrow=nSubjects*8, ncol=2))
 colnames(res1) <- c("ID","MODEL")
 crit           <- paste("crit",criteria,sep="")
 colnames(res2) <- c("ID",crit)
 for (i in 1:nSubjects) {
  res1[i,1] <- i
  res1[i,2] <- mChoose[i]   ## L'erreur est ici
  res2[which(modelShow$ID == i & modelShow$MODEL == res1[i,2]),2] = TRUE
  }
 modelShow                            <- data.frame(modelShow,res2[,2])
 colnames(modelShow)[ncol(modelShow)] <- crit
 modelShow[crit]                      <- !is.na(modelShow[crit])
 return(modelShow)
 }
 

