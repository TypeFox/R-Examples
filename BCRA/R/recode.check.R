recode.check <-
function(dat, Raw_Ind=1){
    ### set error indicator to default value of 0 for each subject
    ## if mean not 0, implies ERROR in file
    Error_Ind <- rep(0, dim(dat)[1])
 
    ### test for consistency of T1 (initial age) and T2 (projection age)
    set_T1_missing <- dat$T1
    set_T2_missing <- dat$T2
 
    set_T1_missing[which((dat$T1 < 20 | dat$T1 >= 90) | dat$T1 >= dat$T2)] <- NA
    set_T2_missing[which(dat$T2 > 90 | dat$T1 >= dat$T2)] <- NA
 
    Error_Ind[is.na(set_T1_missing)] <- 1
    Error_Ind[is.na(set_T2_missing)] <- 1
 
    ### RR covariates are in raw/original format  
    if (Raw_Ind == 1) {
        ### test for consistency of NumBiop (#biopsies) and Hyperplasia   
        ## set NB_Cat to default value of -1
        NB_Cat <- rep(-1, dim(dat)[1])
 
        ## REQUIREMENT (A)
        NB_Cat[which((dat$N_Biop == 0 | dat$N_Biop == 99) & dat$HypPlas != 99)] <- "A"
        Error_Ind[which(NB_Cat == "A")] <- 1
 
        ## REQUIREMENT (B)
        NB_Cat[which((dat$N_Biop > 0 & dat$N_Biop < 99) & (dat$HypPlas != 0 & dat$HypPlas != 1 & dat$HypPlas != 99))] <- "B"
        Error_Ind[which(NB_Cat == "B")] <- 1
     
        ### editing and recoding for N_Biop
        NB_Cat[which(NB_Cat == -1 & (dat$N_Biop == 0 | dat$N_Biop == 99))] <- 0
        NB_Cat[which(NB_Cat == -1 & dat$N_Biop == 1)] <- 1
        NB_Cat[which(NB_Cat == -1 & (dat$N_Biop >= 2 | dat$N_Biop != 99))] <- 2
        NB_Cat[which(NB_Cat == -1)] <- NA
 
        ### editing and recoding for AgeMen
        AM_Cat <- rep(NA, dim(dat)[1])
        AM_Cat[which((dat$AgeMen >= 14 & dat$AgeMen <= dat$T1) | dat$AgeMen ==99 )] <- 0
        AM_Cat[which(dat$AgeMen >= 12 & dat$AgeMen < 14)] <- 1
        AM_Cat[which(dat$AgeMen > 0 & dat$AgeMen < 12)] <- 2
        AM_Cat[which(dat$AgeMen > dat$T1 & dat$AgeMen !=99)] <- NA
        ## for African-Americans AgeMen code 2 (age <= 11) grouped with code 1(age == 12 or 13)
        AM_Cat[which(dat$Race == 2 & AM_Cat ==2)] <- 1 
 
        ### editing and recoding for Age1st
        AF_Cat <- rep(NA, dim(dat)[1])
        AF_Cat[which(dat$Age1st < 20 | dat$Age1st == 99)] <- 0
        AF_Cat[which(dat$Age1st >= 20 & dat$Age1st < 25)] <- 1
        AF_Cat[which((dat$Age1st >= 25 & dat$Age1st < 30) | dat$Age1st == 98)] <- 2
        AF_Cat[which(dat$Age1st >= 30 & dat$Age1st < 98)] <- 3
        AF_Cat[which(dat$Age1st < dat$AgeMen & dat$AgeMen != 99)] <- NA
        AF_Cat[which(dat$Age1st > dat$T1 & dat$Age1st <98)] <- NA
        ## for African-Americans Age1st is not a RR covariate and not in RR model, set to 0
        AF_Cat[which(dat$Race == 2)] <- 0 
        
        ### editing and recoding for N_Rels
        NR_Cat <- rep(NA, dim(dat)[1])
        NR_Cat[which(dat$N_Rels == 0 | dat$N_Rels == 99)] <- 0
        NR_Cat[which(dat$N_Rels == 1)] <- 1
        NR_Cat[which(dat$N_Rels >= 2 & dat$N_Rels < 99)] <- 2
        ## for Asian-American NR_Cat=2 is pooled with NR_Cat=1
        NR_Cat[which((dat$Race >= 6 & dat$Race <= 11) & NR_Cat == 2)] <- 1
    }
 
    ### Raw_Ind=0 means RR covariates have already been re-coded to 0, 1, 2 or 3 (when necessary)
    ### edit/consistency checks for all relative four risk covariates not performed when Raw_Ind=0. (use this option at your own risk)
    if (Raw_Ind == 0){
        NB_Cat <- dat$N_Biop
        AM_Cat <- dat$AgeMen
        AF_Cat <- dat$Age1st
        NR_Cat <- dat$N_Rels
    }
 
    ### setting RR multiplicative factor for atypical hyperplasia
    R_Hyp <- rep(NA, dim(dat)[1])
    R_Hyp[which(NB_Cat == 0)] <- 1.00
    R_Hyp[which((NB_Cat != "A" & NB_Cat > 0) & dat$HypPlas == 0)] <- 0.93
    R_Hyp[which((NB_Cat != "A" & NB_Cat > 0) & dat$HypPlas == 1)] <- 1.82
    R_Hyp[which((NB_Cat != "A" & NB_Cat > 0) & dat$HypPlas == 99)] <- 1.00
 
    set_HyperP_missing <- dat$HypPlas
    set_R_Hyp_missing <- R_Hyp
    set_HyperP_missing[which(NB_Cat == "A")] <- "A"
    set_R_Hyp_missing[which(NB_Cat == "A")] <- "A"
    set_HyperP_missing[which(NB_Cat == "B")] <- "B"
    set_R_Hyp_missing[which(NB_Cat == "B")] <- "B"
 
    set_Race_missing <- dat$Race
    Race_range<-seq(1,11)
    set_Race_missing[-which(dat$Race %in% Race_range)]<-"U"
 
    Error_Ind[which(is.na(NB_Cat) | is.na(AM_Cat) | is.na(AF_Cat) | is.na(NR_Cat) | set_Race_missing == "U")] <- 1
 
    ### african-american RR model from CARE study:(1) eliminates Age1st from model;(2) groups AgeMen=2 with AgeMen=1;
    ## setting AF_Cat=0 eliminates Age1st and its interaction from RR model;
    AF_Cat[which(dat$Race == 2)] <- 0 
    ## group AgeMen RR level 2 with 1;
    AM_Cat[which(dat$Race == 2 & AM_Cat ==2)] <- 1 
 
    ### for asian-americans NR_Cat=2 is pooled with NR_Cat=1; 
    NR_Cat[which((dat$Race >= 6 & dat$Race <= 11) & NR_Cat == 2)] <- 1
 
    CharRace <- rep(NA, dim(dat)[1])
    CharRace[which(dat$Race == 1)] <- "Wh"      #white SEER 1983:87 BrCa Rate
    CharRace[which(dat$Race == 2)] <- "AA"      #african-american
    CharRace[which(dat$Race == 3)] <- "Hi"      #hispanic
    CharRace[which(dat$Race == 4)] <- "NA"      #native american
    CharRace[which(dat$Race == 5)] <- "Wo"      #white SEER 1995:2003 Rate
    CharRace[which(dat$Race == 6)] <- "Ch"      #chinese
    CharRace[which(dat$Race == 7)] <- "Ja"      #japanese
    CharRace[which(dat$Race == 8)] <- "Fi"      #filipino
    CharRace[which(dat$Race == 9)] <- "Hw"      #hawaiian
    CharRace[which(dat$Race == 10)] <- "oP"     #other pacific islander
    CharRace[which(dat$Race == 11)] <- "oA"     #other asian
    CharRace[which(is.na(CharRace))] <- "??"    #non-applicable race code
 
    recode_check<- cbind(Error_Ind, set_T1_missing, set_T2_missing, NB_Cat, AM_Cat, AF_Cat, NR_Cat, R_Hyp, set_HyperP_missing, set_R_Hyp_missing, set_Race_missing, CharRace)
    recode_check <- data.frame(recode_check, row.names=NULL)
    return(recode_check)
}
