.Data.Processing <- function(Dataset=Dataset, Surr=Surr, True=True, Treat=Treat, Trial.ID=Trial.ID, Pat.ID=Pat.ID, Min.Trial.Size=Min.Trial.Size){
  
  Dataset <- data.frame(cbind(Surr, True, Treat, Trial.ID, Pat.ID))
  Dataset <- na.exclude(Dataset[order(Dataset$Trial.ID),])
  if (length(unique(Dataset$Treat))!=2) stop("Please make sure that the treatment variable has only 2 levels.")
  if (((sort(unique(Dataset$Treat))[1]==c(0) & sort(unique(Dataset$Treat))[2]==c(1))==FALSE) & 
        ((sort(unique(Dataset$Treat))[1]==c(-1) & sort(unique(Dataset$Treat))[2]==c(1))==FALSE))
    stop("Please make sure that the treatment is either coded as control = -1 and experimental = 1, or as control = 0 and experimental = 1.")
  keep <- (table(Dataset$Treat, Dataset$Trial.ID)[1,] * table(Dataset$Treat, Dataset$Trial.ID)[2,])!=0  
  mult.Treat <- unique(Dataset$Trial.ID)[keep]
  Dataset$ok <- rep(x=NA, dim(Dataset)[1])
  for (i in 1: dim(Dataset)[1]){for (j in 1: length(mult.Treat)){if (Dataset$Trial.ID[i] == mult.Treat[j]) {Dataset$ok[i] <- c(1)}}}  
  wide <- (na.exclude(Dataset))[,1:5]
  Data.analyze <- reshape(data=wide, direction="long", varying=c("Surr", "True"), v.names="outcome", timevar="endpoint", times=c("-1", "1"), new.row.names = 1:((dim(wide)[1])*2))[,1:5]
  row.names(Data.analyze) <- NULL
  dataS <- data.frame(Data.analyze[Data.analyze$endpoint==-1,]) 
  dataT <- data.frame(Data.analyze[Data.analyze$endpoint==1,]) 
  trialNames <- unique(dataS$Trial.ID)
  
  for (i in 1: length(unique(wide$Trial.ID))){ # remove trials with constant outcome within trial
    if (var(dataS$outcome[dataS$Trial.ID==trialNames[i]])==0 | var(dataT$outcome[dataS$Trial.ID==trialNames[i]])==0) {
      dataS$outcome[dataS$Trial.ID==trialNames[i]] <- dataT$outcome[dataS$Trial.ID==trialNames[i]] <- Data.analyze[Data.analyze$Trial.ID==trialNames[i], ]$outcome <- NA
      wide[wide$Trial.ID==trialNames[i], ]$Surr <- wide[wide$Trial.ID==trialNames[i], ]$True <- NA} 
  }
  for (i in 1: length(unique(wide$Trial.ID))){ # remove trials with less than spec number of patients
    Trial.size.temp <- length(dataS$outcome[dataS$Trial.ID==trialNames[i]])
    if (Trial.size.temp < Min.Trial.Size) {
      dataS$outcome[dataS$Trial.ID==trialNames[i]] <- dataT$outcome[dataS$Trial.ID==trialNames[i]] <- Data.analyze[Data.analyze$Trial.ID==trialNames[i], ]$outcome <- NA
      wide[wide$Trial.ID==trialNames[i], ]$Surr <- wide[wide$Trial.ID==trialNames[i], ]$True <- NA}
  } 
  
  dataS <- na.exclude(dataS)
  dataT <- na.exclude(dataT)
  Data.analyze <- na.exclude(Data.analyze[order(Data.analyze$Pat.ID),])
  rownames(Data.analyze) <- NULL
  wide <- na.exclude(wide)
  rownames(wide) <- NULL
  N.total <- nrow(wide)
  N.trial <- length(unique(wide$Trial.ID))
  Obs.per.trial <- data.frame(cbind(unique(wide$Trial.ID), table(wide$Trial.ID, wide$Treat), table(wide$Trial.ID, wide$Treat)[,1]+table(wide$Trial.ID, wide$Treat)[,2]))
  colnames(Obs.per.trial) <- c("Trial", "Number.cont.Treat", "Number.exp.Treat", "Obs.per.trial")  
  rownames(Obs.per.trial) <- NULL
  
  fit <- list(wide=wide, dataS=dataS, dataT=dataT, Data.analyze=Data.analyze, N.total=N.total, N.trial=N.trial, Obs.per.trial=Obs.per.trial)
  
}
