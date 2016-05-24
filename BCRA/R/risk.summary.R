risk.summary <-
function(dat, iloop=1, Raw_Ind=1){
    ### for further analysis, create output in a data frame that included age, duration of the projection time interval, covariates and the projected risk
    check_cov <- recode.check(dat, Raw_Ind)
    CharRace <- as.character(check_cov$CharRace)
    AbsRisk <- absolute.risk(dat, iloop, Raw_Ind)
    time_intvl <- dat$T2 - dat$T1
    time_intvl <- round(time_intvl,3)
    risk_table <- cbind(dat$ID, dat$T1, dat$T2, time_intvl, dat$N_Biop, dat$HypPlas, dat$AgeMen, dat$Age1st, dat$N_Rels, dat$Race, CharRace, AbsRisk)
    if(iloop==1){
       colnames(risk_table) <- c("ID", "T1", "T2", "Proj_Intvl", "N_Biop", "HypPlas", "AgeMen", "Age1st", "N_Rels", "Race", "CharRace", "AbsRisk")
    }
    if(iloop==2){
       colnames(risk_table) <- c("ID", "T1", "T2", "Proj_Intvl", "N_Biop", "HypPlas", "AgeMen", "Age1st", "N_Rels", "Race", "CharRace", "AbsRisk_Avg")
    }
    risk_table<-data.frame(risk_table, row.names=NULL)
    write.csv(risk_table, "risk_summary.csv", row.names=F)
    return(risk_table)     
}
