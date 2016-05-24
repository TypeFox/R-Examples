absolute.risk <-
function(dat, iloop=1, Raw_Ind=1){
    ### set up lambda1*, lambda2, beta & F(t) with known constants used in the nci brca risk disk
    ## lambda1_Star, BrCa composite incidences
    # SEER BrCa incidence rates (current) non-hispanic white women, SEER white 1983:87
    White_lambda1 <- c(0.00001000, 0.00007600, 0.00026600, 0.00066100, 0.00126500, 0.00186600, 0.00221100, 
                       0.00272100, 0.00334800, 0.00392300, 0.00417800, 0.00443900, 0.00442100, 0.00410900)
    # SEER BrCa incidence rates for "avg" non-hispanic white women and "avg" other (native american) women, SEER white 1992:96
    White_lambda1Avg <- c(0.00001220, 0.00007410, 0.00022970, 0.00056490, 0.00116450, 0.00195250, 0.00261540, 
                          0.00302790, 0.00367570, 0.00420290, 0.00473080, 0.00494250, 0.00479760, 0.00401060)
    # SEER BrCa indicdence rates (under study) for non-hispanic white women, SEER white 1995:2003
    White_nlambda1 <- c(0.0000120469, 0.0000746893, 0.0002437767, 0.0005878291, 0.0012069622, 0.0019762053, 0.0026200977, 
                        0.0033401788, 0.0039743676, 0.0044875763, 0.0048945499, 0.0051610641, 0.0048268456, 0.0040407389)
    # SEER black 1994-98
    Black_lambda1 <- c(0.00002696, 0.00011295, 0.00031094, 0.00067639, 0.00119444, 0.00187394, 0.00241504, 
                       0.00291112, 0.00310127, 0.00366560, 0.00393132, 0.00408951, 0.00396793, 0.00363712)
    # SEER hspan 1990:96
    Hspnc_lambda1 <- c(0.00002000, 0.00007100, 0.00019700, 0.00043800, 0.00081100, 0.00130700, 0.00157400, 
                       0.00185700, 0.00215100, 0.00251200, 0.00284600, 0.00275700, 0.00252300, 0.00203900)
    # SEER white 1983:87
    Other_lambda1 <- c(0.00001000, 0.00007600, 0.00026600, 0.00066100, 0.00126500, 0.00186600, 0.00221100, 
                       0.00272100, 0.00334800, 0.00392300, 0.00417800, 0.00443900, 0.00442100, 0.00410900)
    # seer18 chinese  1998:02
    Chnes_lambda1 <- c(0.000004059636, 0.000045944465, 0.000188279352, 0.000492930493, 0.000913603501,
                       0.001471537353, 0.001421275482, 0.001970946494, 0.001674745804, 0.001821581075,
                       0.001834477198, 0.001919911972, 0.002233371071, 0.002247315779)
    # seer18 japanese 1998:02
    Japns_lambda1 <- c(0.000000000001, 0.000099483924, 0.000287041681, 0.000545285759, 0.001152211095,
                       0.001859245108, 0.002606291272, 0.003221751682, 0.004006961859, 0.003521715275,
                       0.003593038294, 0.003589303081, 0.003538507159, 0.002051572909)
    # seer18 filipino 1998:02
    Filip_lambda1 <- c(0.000007500161, 0.000081073945, 0.000227492565, 0.000549786433, 0.001129400541,
                       0.001813873795, 0.002223665639, 0.002680309266, 0.002891219230, 0.002534421279,
                       0.002457159409, 0.002286616920, 0.001814802825, 0.001750879130)
    # seer18 hawaiian 1998:02
    Hawai_lambda1 <- c(0.000045080582, 0.000098570724, 0.000339970860, 0.000852591429, 0.001668562761,
                       0.002552703284, 0.003321774046, 0.005373001776, 0.005237808549, 0.005581732512,
                       0.005677419355, 0.006513409962, 0.003889457523, 0.002949061662)
    # seer18 otr pac isl 1998:02
    OtrPI_lambda1 <- c(0.000000000001, 0.000071525212, 0.000288799028, 0.000602250698, 0.000755579402,
                       0.000766406354, 0.001893124938, 0.002365580107, 0.002843933070, 0.002920921732,
                       0.002330395655, 0.002036291235, 0.001482683983, 0.001012248203)
    # seer18 otr asian 1998:02
    OtrAs_lambda1 <- c(0.000012355409, 0.000059526456, 0.000184320831, 0.000454677273, 0.000791265338,
                       0.001048462801, 0.001372467817, 0.001495473711, 0.001646746198, 0.001478363563,
                       0.001216010125, 0.001067663700, 0.001376104012, 0.000661576644)
    ## lambda2, Competing hazards
    #nchs competing mortality (current) for non-hispanic white women, NCHS white 1985:87
    White_lambda2 <- c(0.00049300, 0.00053100, 0.00062500, 0.00082500, 0.00130700, 0.00218100, 0.00365500, 
                       0.00585200, 0.00943900, 0.01502800, 0.02383900, 0.03883200, 0.06682800, 0.14490800)
    # nchs competing mortality for "avg" non-hispanic white women and "avg" other (native american) women, NCHS white 1992:96
    White_lambda2Avg <- c(0.00044120, 0.00052540, 0.00067460, 0.00090920, 0.00125340, 0.00195700, 0.00329840, 
                          0.00546220, 0.00910350, 0.01418540, 0.02259350, 0.03611460, 0.06136260, 0.14206630)
    # nchs competing mortality (under study) for non-hispanic white women, NCHS white 1995:2003
    White_nlambda2 <- c(0.0004000377, 0.0004280396, 0.0005656742, 0.0008474486, 0.0012752947, 0.0018601059, 0.0028780622, 
                        0.0046903348, 0.0078835252, 0.0127434461, 0.0208586233, 0.0335901145, 0.0575791439, 0.1377327125)
    # NCHS black 1996-00
    Black_lambda2 <- c(0.00074354, 0.00101698, 0.00145937, 0.00215933, 0.00315077, 0.00448779, 0.00632281, 
                       0.00963037, 0.01471818, 0.02116304, 0.03266035, 0.04564087, 0.06835185, 0.13271262)
    # NCHS hspan 1990:96
    Hspnc_lambda2 <- c(0.00043700, 0.00053300, 0.00070000, 0.00089700, 0.00116300, 0.00170200, 0.00264600,
                       0.00421600, 0.00696000, 0.01086700, 0.01685800, 0.02515600, 0.04186600, 0.08947600)
    # NCHS white 1985:87
    Other_lambda2 <- c(0.00049300, 0.00053100, 0.00062500, 0.00082500, 0.00130700, 0.00218100, 0.00365500, 
                       0.00585200, 0.00943900, 0.01502800, 0.02383900, 0.03883200, 0.06682800, 0.14490800)
    # NCHS mortality chinese  1998:02
    Chnes_lambda2 <- c(0.000210649076, 0.000192644865, 0.000244435215, 0.000317895949, 0.000473261994,
                       0.000800271380, 0.001217480226, 0.002099836508, 0.003436889186, 0.006097405623,
                       0.010664526765, 0.020148678452, 0.037990796590, 0.098333900733) 
    # NCHS mortality japanese 1998:02
    Japns_lambda2 <- c(0.000173593803, 0.000295805882, 0.000228322534, 0.000363242389, 0.000590633044,
                       0.001086079485, 0.001859999966, 0.003216600974, 0.004719402141, 0.008535331402,
                       0.012433511681, 0.020230197885, 0.037725498348, 0.106149118663)
    # NCHS mortality filipino 1998:02
    Filip_lambda2 <- c(0.000229120979, 0.000262988494, 0.000314844090, 0.000394471908, 0.000647622610,
                       0.001170202327, 0.001809380379, 0.002614170568, 0.004483330681, 0.007393665092,
                       0.012233059675, 0.021127058106, 0.037936954809, 0.085138518334)
    # NCHS mortality hawaiian 1998:02
    Hawai_lambda2 <- c(0.000563507269, 0.000369640217, 0.001019912579, 0.001234013911, 0.002098344078,
                       0.002982934175, 0.005402445702, 0.009591474245, 0.016315472607, 0.020152229069,
                       0.027354838710, 0.050446998723, 0.072262026612, 0.145844504021)
    # NCHS mortality otr pac isl 1998:02
    OtrPI_lambda2 <- c(0.000465500812, 0.000600466920, 0.000851057138, 0.001478265376, 0.001931486788,
                       0.003866623959, 0.004924932309, 0.008177071806, 0.008638202890, 0.018974658371,
                       0.029257567105, 0.038408980974, 0.052869579345, 0.074745721133)
    # NCHS mortality otr asian 1998:02
    OtrAs_lambda2 <- c(0.000212632332, 0.000242170741, 0.000301552711, 0.000369053354, 0.000543002943,
                       0.000893862331, 0.001515172239, 0.002574669551, 0.004324370426, 0.007419621918,
                       0.013251765130, 0.022291427490, 0.041746550635, 0.087485802065)
    ## F(t), 1-Attributable Risk=F(t) 
    White_1_AR <- c(0.5788413, 0.5788413)
    Black_1_AR <- c(0.72949880, 0.74397137)
    Hspnc_1_AR <- c(0.5788413, 0.5788413) 
    Other_1_AR <- c(0.5788413, 0.5788413)
    Asian_1_AR <- c(0.47519806426735, 0.50316401683903)

    ### intialize "avg white women" and "avg" other (native american women) rate for each year in the 5yr age cat
    Avg_lambda1 <- array(0, dim=c(14,5))
    Avg_lambda1[,1:(dim(Avg_lambda1)[2])] <- White_lambda1Avg
    Avg_lambda2 <- array(0, dim=c(14,5))
    Avg_lambda2[,1:(dim(Avg_lambda2)[2])] <- White_lambda2Avg

    ### initialize rate vectors with the correct rates for each woman under study based on her race
    ## for i=1 to 11, when Race=i, Wrk_lambda1=Wrk_lambda1_all[i], Wrk_lambda2=Wrk_lambda2_all[i], Wrk_Beta<-Wrk_Beta_all[i], Wrk_1_AR=Wrk_1_AR_all[i]
    Wrk_lambda1_all <- rbind(White_lambda1, Black_lambda1, Hspnc_lambda1, Other_lambda1, White_nlambda1, Chnes_lambda1, Japns_lambda1, Filip_lambda1, Hawai_lambda1, OtrPI_lambda1, OtrAs_lambda1)
    Wrk_lambda2_all <- rbind(White_lambda2, Black_lambda2, Hspnc_lambda2, Other_lambda2, White_nlambda2, Chnes_lambda2, Japns_lambda2, Filip_lambda2, Hawai_lambda2, OtrPI_lambda2, OtrAs_lambda2) 
    Wrk_1_AR_all <- rbind(White_1_AR, Black_1_AR, Hspnc_1_AR, Other_1_AR, White_1_AR, Asian_1_AR, Asian_1_AR, Asian_1_AR, Asian_1_AR, Asian_1_AR, Asian_1_AR)

    AbsRisk <- rep(NA,dim(dat)[1])
    ## obtain IDs without any error
    check_cov <- recode.check(dat, Raw_Ind)
    Error_Ind <- check_cov$Error_Ind
    IDwoERR <- which(Error_Ind==0)
    for (i in IDwoERR){
         obs <- dat[i,]
         RR_Star <- relative.risk(dat,Raw_Ind)
         rrstar1 <- RR_Star$RR_Star1[i]
         rrstar2 <- RR_Star$RR_Star2[i]

         One_AR_RR <- rep(NA, 70)    
         Strt_Intvl <- floor(obs$T1)-20+1
         End_Intvl <- ceiling(obs$T2)-20+0
         NumbrIntvl <- ceiling(obs$T2)-floor(obs$T1)
         RskWrk <- 0
         Cum_lambda <- 0
         lambda1.temp <- array(0, dim=c(14,5))
         lambda2.temp <- array(0, dim=c(14,5))

         ## calculate abs risk  
         if (iloop == 1){
             One_AR1 <- Wrk_1_AR_all[obs$Race,1]
             One_AR2 <- Wrk_1_AR_all[obs$Race,2]
             # (1-AR)*RR at ages < 50
             One_AR_RR1 <- One_AR1*rrstar1
             # (1-AR)*RR at ages >= 50
             One_AR_RR2 <- One_AR2*rrstar2
             # define One_AR_RR
             One_AR_RR[1:30] <- One_AR_RR1
             One_AR_RR[31:70] <- One_AR_RR2
             lambda1.temp[,1:(dim(lambda1.temp)[2])] <- Wrk_lambda1_all[obs$Race,]
             lambda2.temp[,1:(dim(lambda2.temp)[2])] <- Wrk_lambda2_all[obs$Race,]
             lambda1 <- c(t(lambda1.temp))
             lambda2 <- c(t(lambda2.temp))
         }
         ## calculate avg abs risk
         if (iloop == 2){
             # define One_AR_RR
             One_AR_RR <- rep(1, 70)
             lambda1.temp[,1:(dim(lambda1.temp)[2])] <- Wrk_lambda1_all[obs$Race,]
             lambda2.temp[,1:(dim(lambda2.temp)[2])] <- Wrk_lambda2_all[obs$Race,]
             # if Race=1, 4 or 5, lambda1.temp[race,1:5]<-Avg_lambda1[race,],lambda2.temp[race,1:5]<-Avg_lambda2[race,]
             if (obs$Race==1 | obs$Race==4 | obs$Race==5){
                 lambda1.temp <- Avg_lambda1
                 lambda2.temp <- Avg_lambda2
             }
             lambda1 <- c(t(lambda1.temp))
             lambda2 <- c(t(lambda2.temp))
         }
         for (j in 1:NumbrIntvl){
              j_intvl <- Strt_Intvl+j-1 
              if (NumbrIntvl>1 & j>1 & j<NumbrIntvl){
                  IntgrlLngth <- 1
              }
              if (NumbrIntvl>1 & j==1){
                  IntgrlLngth <- 1-(obs$T1-floor(obs$T1))        
              }
              if (NumbrIntvl>1 & j==NumbrIntvl){
                  z1 <- ifelse((obs$T2>floor(obs$T2)), 1, 0)
                  z2 <- ifelse((obs$T2==floor(obs$T2)), 1, 0)
                  IntgrlLngth <- (obs$T2-floor(obs$T2))*z1+z2   
              }
              if (NumbrIntvl==1){
                  IntgrlLngth <- obs$T2-obs$T1
              }
              lambdaj <- lambda1[j_intvl]*One_AR_RR[j_intvl]+lambda2[j_intvl]
              PI_j <- ((One_AR_RR[j_intvl]*lambda1[j_intvl]/lambdaj)*exp(-Cum_lambda))*(1-exp(-lambdaj*IntgrlLngth))
              RskWrk <- RskWrk+PI_j 
              Cum_lambda <- Cum_lambda+lambdaj*IntgrlLngth     
          } 
          AbsRisk[i] <- 100*RskWrk    
     }
     return(AbsRisk)
}
