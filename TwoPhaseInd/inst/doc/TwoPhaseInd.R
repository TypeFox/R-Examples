### R code from vignette source 'TwoPhaseInd.Rnw'

###################################################
### code chunk number 1: loadLibrary
###################################################

library(TwoPhaseInd)


###################################################
### code chunk number 2: data
###################################################
data(acodata)
dim(acodata)
str(acodata)


###################################################
### code chunk number 3: cfit
###################################################
cfit=caseonly(data=acodata[acodata[,2]==1,], ##dataset
              treatment="f_treat",  ##treatment variable
              BaselineMarker="fcgr2a.3")  ##biomarker
cfit


###################################################
### code chunk number 4: whiBioMarker
###################################################
data(whiBioMarker)

dim(whiBioMarker)
str(whiBioMarker)


###################################################
### code chunk number 5: spmleNonIndExtra
###################################################
spmleNonIndExtra <- spmle(data=whiBioMarker,  ## dataset
               response="stroke",  ## response variable
               treatment="hrtdisp",	## treatment variable
               BaselineMarker="papbl",	## biomarker
               extra=c(
                       "age"
                        , "dias" 	
                        , "hyp" 	
                        , "syst" 	
                        , "diabtrt" 
                        , "lmsepi" 
                            ),	## extra variable(s)
               phase="phase",	## phase indicator
               ind=FALSE	## independent or non-indepentent
)

spmleNonIndExtra


###################################################
### code chunk number 6: spmleIndExtra
###################################################
spmleIndExtra <- spmle(data=whiBioMarker,	## dataset
            response="stroke",	## response variable
            treatment="hrtdisp",	## treatment variable
            BaselineMarker="papbl",	## biomarker
            extra=c(
               "age"  
		   , "dias"	
              , "hyp" 
              , "syst" 
              , "diabtrt"	
              , "lmsepi" 
                 ),	## extra variable(s)
            phase="phase", ## phase indicator
            ind=TRUE ## independent or non-indepentent
)

spmleIndExtra


###################################################
### code chunk number 7: melIndExtra
###################################################
melIndExtra <- mele(data=whiBioMarker,	## dataset
          response="stroke",	## response variable
          treatment="hrtdisp",		## treatment variable
          BaselineMarker="papbl",		## biomarker
          extra=c(
             "age" 	
              , "dias"  
              , "hyp" ## 
              , "syst" 	
              , "diabtrt"	
              , "lmsepi" 
              ),	## extra variable(s)
          phase="phase",	## phase indicator
          ind=TRUE	## independent or non-indepentent
)
melIndExtra


###################################################
### code chunk number 8: melNoIndExtra
###################################################
melNoIndExtra <- mele(data=whiBioMarker,	## dataset
            response="stroke",	## response variable
            treatment="hrtdisp",	## treatment variable
            BaselineMarker="papbl",	## biomarker
            extra=c(
                "age"
                , "dias" 	
                , "hyp" 	
                , "syst" 	
                , "diabtrt"	
                , "lmsepi" 
                ),	## extra variable(s)
            phase="phase",	## phase indicator
            ind=FALSE	## independent or non-indepentent
)
melNoIndExtra


###################################################
### code chunk number 9: data
###################################################
data(acodata)

dim(acodata)
str(acodata)


###################################################
### code chunk number 10: rfit0
###################################################
rfit0 <- acoarm(data=acodata,  ## dataset
                 svtime="vacc1_evinf", ## survival time
                 event="f_evinf",  ## event
                 treatment="f_treat", ## treatment
                 BaselineMarker="fcgr2a.3",  #biomarker
                 id="ptid",  #participant id
                 subcohort="subcoh", #subcohort
                 esttype=1, ## use Self-Prentice method
                 augment=0, ## augment from placebo arm
                 extra=c("f_agele30"
                         ,"f_hsv_2"
                         ,"f_ad5gt18"
                         ,"f_crcm"
                         ,"any_drug"
                         ,"num_male_part_cat"
                         ,"uias"
                         ,"uras")) ## extra varibles
rfit0


###################################################
### code chunk number 11: rfit1
###################################################
rfit1 <- acoarm(data=acodata,  ## dataset
                 svtime="vacc1_evinf",  ## survival time
                 event="f_evinf",  ## event
                 treatment="f_treat", ## treatment
                 BaselineMarker="fcgr2a.3",  #biomarker
                 id="ptid",  #participant id
                 subcohort="subcoh", #subcohort
                 esttype=1, ## use Self-Prentice method
                 augment=1,## augment from active arm
                 extra=c("f_agele30"
                         ,"f_hsv_2"
                         ,"f_ad5gt18"
                         ,"f_crcm"
                         ,"any_drug"
                         ,"num_male_part_cat"
                         ,"uias"
                         ,"uras")) ## extra varibles
rfit1


###################################################
### code chunk number 12: rfit2
###################################################
rfit2 <- acoarm(data=acodata,  ## dataset
                 svtime="vacc1_evinf",  ## survival time
                 event="f_evinf",  ## event
                 treatment="f_treat", ## treatment
                 BaselineMarker="fcgr2a.3",  #biomarker
                 id="ptid",  #participant id
                 subcohort="subcoh", #subcohort
                 esttype=1, ## use Self-Prentice method
                 augment=2,## augment from both arms
                 extra=c("f_agele30"
                         ,"f_hsv_2"
                         ,"f_ad5gt18"
                         ,"f_crcm"
                         ,"any_drug"
                         ,"num_male_part_cat"
                         ,"uias"
                         ,"uras")) ## extra varibles
rfit2


###################################################
### code chunk number 13: TwoPhaseInd.Rnw:295-300
###################################################

cat('\\begin{figure}[h]\n')
file = "./figure1new.png"
cat('\\includegraphics{', file, '}\n', sep = '')
cat('\\end{figure}\n')


###################################################
### code chunk number 14: sessionInfo
###################################################
sessionInfo()


