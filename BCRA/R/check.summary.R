check.summary <-
function(dat, iloop=1, Raw_Ind=1){
    ### list mean summary of "Error_Ind, AbsRsk, RR_Star1 and RR_Star2"
    mean_table<-array(NA, dim=c(4,6))
    colnames(mean_table) <- c("Variable", "Label", "Mean", "StdDev", "N", "NMiss")
    check_cov <- recode.check(dat, Raw_Ind)
    Error_Ind <- as.numeric(as.character(check_cov$Error_Ind))
    RR_Star <- relative.risk(dat,Raw_Ind)
    RR_Star1 <- RR_Star$RR_Star1
    RR_Star2 <- RR_Star$RR_Star2


    AbsRisk <- absolute.risk(dat, iloop, Raw_Ind)
    if (iloop==1){
    mean_table[1:4,1] <- c("Error_Ind", "AbsRisk", "RR_Star1", "RR_Star2")
    }
    if (iloop==2){
    mean_table[1:4,1] <- c("Error_Ind", "AbsRisk_Avg", "RR_Star1", "RR_Star2")
    }
    mean_table[1:4,2] <- c("If mean not 0, implies ERROR in file", "Abs risk(%) of BrCa in age interval [T1,T2)", "Relative risk age lt 50", "Relative risk age ge 50")
    mean_table[1:4,3] <- c(mean(Error_Ind), mean(AbsRisk[which(!is.na(AbsRisk))]), mean(RR_Star1[which(!is.na(RR_Star1))]), mean(RR_Star2[which(!is.na(RR_Star2))]))
    mean_table[1:4,4] <- c(sd(Error_Ind),sd(AbsRisk[which(!is.na(AbsRisk))]), sd(RR_Star1[which(!is.na(RR_Star1))]), sd(RR_Star2[which(!is.na(RR_Star2))]))
    mean_table[1:4,5] <- c(length(which(!is.na(Error_Ind))),length(which(!is.na(AbsRisk))),length(which(!is.na(RR_Star1))),length(which(!is.na(RR_Star2))))
    mean_table[1:4,6] <- c(length(which(is.na(Error_Ind))),length(which(is.na(AbsRisk))),length(which(is.na(RR_Star1))),length(which(is.na(RR_Star2))))
    mean_table<-data.frame(mean_table, row.names=NULL)
    return(mean_table)
}
