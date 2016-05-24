relative.risk <-
function(dat, Raw_Ind=1){
    ## LN_RR, beta=lnRR, beta for NB, AM, AF, NR, AC*NB and AF*NR, from Gail/CARE model
    White_Beta <- c(0.5292641686, 0.0940103059, 0.2186262218, 0.9583027845, -0.2880424830, -0.1908113865)
    Black_Beta <- c(0.1822121131, 0.2672530336, 0.0, 0.4757242578, -0.1119411682, 0.0)
    Hspnc_Beta <- c(0.5292641686, 0.0940103059, 0.2186262218, 0.9583027845, -0.2880424830, -0.1908113865)
    Other_Beta <- c(0.5292641686, 0.0940103059, 0.2186262218, 0.9583027845, -0.2880424830, -0.1908113865)
    Asian_Beta <- c(0.55263612260619, 0.07499257592975, 0.27638268294593, 0.79185633720481, 0.0, 0.0)
    Wrk_Beta_all <- rbind(White_Beta, Black_Beta, Hspnc_Beta, Other_Beta, White_Beta, Asian_Beta, Asian_Beta, Asian_Beta, Asian_Beta, Asian_Beta, Asian_Beta)

    ### define LP1 = Linear Predictor for woman of interest at ages < 50; LP2 = Linear Predictor for woman of interest at ages >= 50
    LP1 <- rep(NA, dim(dat)[1])
    LP2 <- rep(NA, dim(dat)[1])
   
    ### obtain covariates
    check_cov<-recode.check(dat, Raw_Ind)

    NB_Cat<-check_cov$NB_Cat 
    NB_Cat[which(NB_Cat=="A" | NB_Cat=="B")]<-NA
    NB_Cat <- as.numeric(as.character(NB_Cat))

    AM_Cat    <- as.numeric(as.character(check_cov$AM_Cat)) 
    AF_Cat    <- as.numeric(as.character(check_cov$AF_Cat))
    NR_Cat    <- as.numeric(as.character(check_cov$NR_Cat))
    R_Hyp     <- as.numeric(as.character(check_cov$R_Hyp))
    CharRace  <- check_cov$CharRace

    ### define pattern number when NB_Cat, AM_Cat, AF_Cat, NR_Cat are meaningful
    ### NB_Cat(3 levels), AM_Cat(3 levels), AF_Cat(4 levels), NR_Cat(3 levels), 3*3*4*3 = 108 patterns in total
    ## let PNID be the ID numbers when all "_Cat" variables are numerical
    PatternNumber <- rep(NA, dim(dat)[1])
    PNID <- which(NB_Cat!="A" & NB_Cat!="B" & !is.na(AM_Cat) & !is.na(AF_Cat) & !is.na(NR_Cat))   
    PatternNumber[PNID] <- NB_Cat[PNID]*36+AM_Cat[PNID]*12+AF_Cat[PNID]*3+NR_Cat[PNID]*1+1

    for (i in PNID){
         if (CharRace[i]!="??"){
             Beta <- Wrk_Beta_all[dat$Race[i],]
             ## for woman at ages < 50
             LP1[i] <- NB_Cat[i]*Beta[1]+AM_Cat[i]*Beta[2]+AF_Cat[i]*Beta[3]+NR_Cat[i]*Beta[4]+AF_Cat[i]*NR_Cat[i]*Beta[6]+log(R_Hyp[i])
             LP2[i] <- LP1[i]+NB_Cat[i]*Beta[5]
         }
    }

    ### define RR_Star1 = relative risk for woman of interest at ages < 50; RR_Star2 = relative risk for woman of interest at ages >= 50
    RR_Star1 <- exp(LP1)
    RR_Star2 <- exp(LP2)  
    RR_Star <- cbind(RR_Star1, RR_Star2, PatternNumber)
    RR_Star <- data.frame(RR_Star, row.names=NULL) 
    return(RR_Star)
}
