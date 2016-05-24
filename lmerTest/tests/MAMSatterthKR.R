require(lmerTest)
testType1 <- TRUE


lm.pred <- lm(as.formula(paste("Lightlevel", "~", 
                               paste(c("TVset","Picture"), collapse="*"), sep="")),
              data=TVbo)
TVbo$x <- scale(predict(lm.pred), scale=FALSE)
lmerTVpic <- lmer(Lightlevel ~ TVset*Picture +   Assessor:x  + (1|Assessor) +
                    (1|TVset:Assessor) + (1|Picture:Assessor) + 
                    (1|TVset:Picture:Assessor), data=TVbo)

## TODO: check with SAS dfs for Satterthwaite and KR to agree
if(testType1){
tools::assertWarning(anova(lmerTVpic, type=1)) ## warning: ddf=0 for TVset
anova(lmerTVpic, type=1, ddf="Kenward-Roger")
}

## check MAM for BO data
# 
# m.bo <- lmer(att1 ~ Track*SPL*Car + (1|Assessor) + (1|SPL:Assessor) + 
#                      (1|Track:SPL:Assessor) + (1|Car:SPL:Assessor), 
#                    data = sound_data_balanced)
# 
# lm.pred <- lm(as.formula(paste("att1", "~", 
#                                paste(c("Track", "SPL", "Car"), collapse="*"), sep="")),
#               data=sound_data_balanced)
# sound_data_balanced$x <- scale(predict(lm.pred), scale=FALSE)
# 
# m.bo.mam <- lmer(att1 ~ Track*SPL*Car + Assessor:x + (1|Assessor) + (1|SPL:Assessor) + 
#                (1|Track:SPL:Assessor) + (1|Car:SPL:Assessor), 
#              data = sound_data_balanced)
# 
# sound_data_balanced$y <- rep(1, nrow(sound_data_balanced))
# 
# m.bo.mam <- lmer(att1 ~ (Track + SPL+ Car)^2 + Assessor:x:y + (1|Assessor) + (1|SPL:Assessor) + 
#                    (1|Track:SPL:Assessor) + (1|Car:SPL:Assessor), 
#                  data = sound_data_balanced)
# 
# anova(m.bo.mam, type=1)
# 
# 
# ## non standard order
# m.bo <- lmer(att1 ~ Track:SPL:Car + Track + Car + (1|Assessor) + (1|SPL:Assessor) + 
#                (1|Track:SPL:Assessor) + (1|Car:SPL:Assessor), 
#              data = sound_data_balanced)
