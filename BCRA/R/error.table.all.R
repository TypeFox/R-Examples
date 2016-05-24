error.table.all <-
function(dat, Raw_Ind=1){
    ### list all possible errors for BrCa absolute risk projections
    check_cov <- recode.check(dat, Raw_Ind)
    Error_Ind <- check_cov$Error_Ind

    R_Hyp <- as.numeric(as.character(check_cov$R_Hyp))
    NB_Cat <- as.character(check_cov$NB_Cat)
    AM_Cat <- as.character(check_cov$AM_Cat)
    AF_Cat <- as.character(check_cov$AF_Cat)
    NR_Cat <- as.character(check_cov$NR_Cat)
    CharRace <- as.character(check_cov$CharRace)

    set_T1_missing     <- as.numeric(as.character(check_cov$set_T1_missing))
    set_T2_missing     <- as.numeric(as.character(check_cov$set_T2_missing)) 
    set_HyperP_missing <- as.character(check_cov$set_HyperP_missing)
    set_R_Hyp_missing  <- as.character(check_cov$set_R_Hyp_missing)
    if (length(which(Error_Ind!=0))==0){
        error_table <- NULL
        print("NO ERROR!") 
    }
    if (length(which(Error_Ind!=0))!=0){
        error_table <- array("*", dim=c(length(dat$ID)*2,10))
        origin <- cbind(dat$ID, dat$T1, dat$T2, dat$N_Biop, dat$HypPlas, R_Hyp, dat$AgeMen, dat$Age1st, dat$N_Rels, dat$Race)
        compare <- cbind(dat$ID, set_T1_missing, set_T2_missing, NB_Cat, set_HyperP_missing, set_R_Hyp_missing, AM_Cat, AF_Cat, NR_Cat, CharRace)
        for (i in 1:length(dat$ID)){
             error_table[2*i-1,] <- origin[i,]
             error_table[2*i,] <- compare[i,]
        }
        colnames(error_table) <- c("ID", "T1", "T2", "N_Biop", "HypPlas", "R_Hyp", "AgeMen", "Age1st", "N_Rels", "Race")
        error_table<-data.frame(error_table, row.names=NULL)
        return(error_table)
    } 
}
