require(SensMixed)
load(system.file("testdata","sensBObalanc.RData",package="SensMixed"))
load(system.file("testdata","sensBO.RData",package="SensMixed"))
testBO <- FALSE
###########################################
## check for TVbo without replication
###########################################



TVbo <- convertToFactors(TVbo, c("Assessor", "Repeat", "Picture"))
result <- sensmixed(c("Noise", "Elasticeffect"),
                    Prod_effects = c("TVset", "Picture"),
                    replication = "Repeat", 
                    individual = "Assessor", data = TVbo,
                    calc_post_hoc = TRUE)

m.test.noise <- lm(Noise ~ TVset + Picture + Picture:Assessor + 
                  TVset:Picture + Assessor, data = TVbo)

stopifnot(all.equal( sqrt(anova(m.test.noise)[1,4]*2/64), 
                     result$fixed$dprimeav[1, "Noise"], tol = 1e-4))

stopifnot(all.equal( sqrt(anova(m.test.noise)[2,4]*2/48), 
                     result$fixed$dprimeav[2, "Noise"], tol = 1e-4))

deltas <- tapply(TVbo$Noise,factor(TVbo$TVset:TVbo$Picture),mean)-
  rep(tapply(TVbo$Noise,factor(TVbo$TVset),mean),rep(4,3))-
  rep(tapply(TVbo$Noise,factor(TVbo$Picture),mean),3)+
  rep(mean(TVbo$Noise),12)

deltadifs=matrix(rep(0,144),ncol=12)
for (i in 1:12) for (j in 1:i) deltadifs[i,j]=deltas[i]-deltas[j]


stopifnot(all.equal( sqrt(sum(deltadifs^2/anova(m.test.noise)[6,3])/(66)), 
                     result$fixed$dprimeav[3, "Noise"], tol = 1e-4))

stopifnot(all.equal( sqrt(anova(m.test.noise)[5,4]*(1/8)*(6/11)), 
                     result$fixed$dprimeav[3, "Noise"], tol = 1e-4))
  
## without replication
res <- sensmixed(c("Coloursaturation", "Cutting"),
                Prod_effects = c("TVset"), replication="Repeat", 
                individual="Assessor", data=TVbo, parallel = FALSE,
                calc_post_hoc = TRUE)

m.test.fixed <- lm(Coloursaturation ~ TVset + TVset:Assessor + Assessor, 
                     data = TVbo)
stopifnot(all.equal( sqrt(anova(m.test.fixed)[1,4]*2/64), 
           res$fixed$dprimeav[, "Coloursaturation"], tol = 1e-4))

m.test.fixed.cut <- lm(Cutting ~ TVset + TVset:Assessor + Assessor, 
                   data = TVbo)
stopifnot(all.equal( sqrt(anova(m.test.fixed.cut)[1,4]*2/64), 
           res$fixed$dprimeav[, "Cutting"], tol = 1e-4))


## with replication
res.rep <- sensmixed(c("Coloursaturation", "Cutting"),
                 Prod_effects = c("TVset"), replication="Repeat", 
                 individual="Assessor", data=TVbo, error_structure = "3-WAY", 
                 parallel = FALSE, calc_post_hoc = TRUE)



m.testrep.cut <- lm(Cutting ~ TVset + TVset:Assessor + TVset:Repeat + Assessor, 
                   data = TVbo)
stopifnot(all.equal( sqrt(anova(m.testrep.cut)[1,4]*2/64), 
           res.rep$fixed$dprimeav[, "Cutting"], tol = 1e-4))

## with interaction
res.inter <- sensmixed(c("Colourbalance", "Cutting"),
                 Prod_effects = c("TVset", "Picture"), replication="Repeat", 
                 individual="Assessor", data=TVbo, parallel = FALSE,
                 calc_post_hoc = TRUE)

m.testinter.cut <- lm(Cutting ~ TVset*Picture + TVset:Assessor + Assessor, 
                    data = TVbo)

stopifnot(all.equal( sqrt(anova(m.testinter.cut)[1,4]*2/64),
           res.inter$fixed$dprimeav[1, "Cutting"], tol = 1e-4))
stopifnot(all.equal( sqrt(anova(m.testinter.cut)[2,4]*2/48),
           res.inter$fixed$dprimeav[2, "Cutting"], tol = 1e-4))

m.testinter.col <- lm(Colourbalance ~ TVset*Picture + TVset:Assessor + Assessor+
                        Picture:Assessor, 
                      data = TVbo)
stopifnot(all.equal( sqrt(anova(m.testinter.col)[1,4]*2/64),
           res.inter$fixed$dprimeav[1, "Colourbalance"], tol = 1e-3))
stopifnot(all.equal( sqrt(anova(m.testinter.col)[2,4]*2/48),
           res.inter$fixed$dprimeav[2, "Colourbalance"], tol = 1e-3))

## check interaction term for Cutting
deltas <- tapply(TVbo$Cutting,factor(TVbo$TVset:TVbo$Picture),mean)-
  rep(tapply(TVbo$Cutting,factor(TVbo$TVset),mean),rep(4,3))-
  rep(tapply(TVbo$Cutting,factor(TVbo$Picture),mean),3)+
  rep(mean(TVbo$Cutting),12)

deltadifs=matrix(rep(0,144),ncol=12)
for (i in 1:12) for (j in 1:i) deltadifs[i,j]=deltas[i]-deltas[j]


stopifnot(all.equal( sqrt(sum(deltadifs^2/anova(m.testinter.cut)[6,3])/(66)), 
           res.inter$fixed$dprimeav[3, "Cutting"], tol = 1e-4))

stopifnot(all.equal( sqrt(anova(m.testinter.cut)[4,4]*(1/8)*(6/11)), 
           res.inter$fixed$dprimeav[3, "Cutting"], tol = 1e-4))

## check interaction term for Colourbalance
deltas <- tapply(TVbo$Colourbalance,factor(TVbo$TVset:TVbo$Picture),mean)-
  rep(tapply(TVbo$Colourbalance,factor(TVbo$TVset),mean),rep(4,3))-
  rep(tapply(TVbo$Colourbalance,factor(TVbo$Picture),mean),3)+
  rep(mean(TVbo$Colourbalance),12)

deltadifs=matrix(rep(0,144),ncol=12)
for (i in 1:12) for (j in 1:i) deltadifs[i,j]=deltas[i]-deltas[j]


stopifnot(all.equal( sqrt(sum(deltadifs^2/anova(m.testinter.col)[7,3])/(66)), 
           res.inter$fixed$dprimeav[3, "Colourbalance"], tol = 1e-3))

stopifnot(all.equal( sqrt(anova(m.testinter.col)[4,4]*(1/8)*(6/11)), 
           res.inter$fixed$dprimeav[3, "Colourbalance"], tol = 1e-3))


## compare examples from the d -prime presentation
res2 <- sensmixed(c("Colourbalance", "Cutting"),
                       Prod_effects = c("TVset", "Picture"), replication="Repeat", 
                       individual="Assessor", data=TVbo, parallel = FALSE,
                       calc_post_hoc = TRUE, reduce.random = FALSE)

## check for cutting
m.testinter.cut <- lm(Cutting ~ TVset*Picture + TVset:Assessor + Assessor +
                        Picture:Assessor + TVset:Picture:Assessor, 
                      data = TVbo)

stopifnot(all.equal( sqrt(anova(m.testinter.cut)[1,4]*2/64),
                     res2$fixed$dprimeav[1, "Cutting"], tol = 1e-2))
stopifnot(all.equal( sqrt(anova(m.testinter.cut)[2,4]*2/48),
                     res2$fixed$dprimeav[2, "Cutting"], tol = 1e-2))

## check interaction term for Cutting
deltas <- tapply(TVbo$Cutting,factor(TVbo$TVset:TVbo$Picture),mean)-
  rep(tapply(TVbo$Cutting,factor(TVbo$TVset),mean),rep(4,3))-
  rep(tapply(TVbo$Cutting,factor(TVbo$Picture),mean),3)+
  rep(mean(TVbo$Cutting),12)

deltadifs=matrix(rep(0,144),ncol=12)
for (i in 1:12) for (j in 1:i) deltadifs[i,j]=deltas[i]-deltas[j]


stopifnot(all.equal( sqrt(sum(deltadifs^2/anova(m.testinter.cut)[8,3])/(66)), 
                     res2$fixed$dprimeav[3, "Cutting"], tol = 1e-2))

stopifnot(all.equal( sqrt(anova(m.testinter.cut)[4,4]*(1/8)*(6/11)), 
                     res.inter$fixed$dprimeav[3, "Cutting"], tol = 1e-1))

## check for colourbalance
m.testinter.col <- lm(Colourbalance ~ TVset*Picture + TVset:Assessor + Assessor +
                        Picture:Assessor + TVset:Picture:Assessor, 
                      data = TVbo)

stopifnot(all.equal( sqrt(anova(m.testinter.col)[1,4]*2/64),
                     res2$fixed$dprimeav[1, "Colourbalance"], tol = 1e-4))
stopifnot(all.equal( sqrt(anova(m.testinter.col)[2,4]*2/48),
                     res2$fixed$dprimeav[2, "Colourbalance"], tol = 1e-4))



## check for 3-way product structure for BO data

# sound_data_balanced <- read.csv("C:/Users/alku/SkyDrive/Documents/Work/work PhD/SensMixed/BO data/sound_data_balanced.csv", sep=";")
# str(sound_data_balanced)
# sound_data_balanced <- convertToFactors(sound_data_balanced, c("Assessor", "Rep", 
#                                                               "Clip", "Track", 
#                                                               "Car", "SPL"))
# 

if(testBO){
  result_bo <- sensmixed(c("att1", "att2"),
                         Prod_effects=c("Track", "Car", "SPL"),
                        replication="Rep", individual="Assessor", 
                         product_structure=3, data=sound_data_balanced,
                        calc_post_hoc = TRUE)
  
  #m.bo <- lmer(att1 ~ Track*SPL*Car + (1|Assessor) + (1|SPL:Assessor) + 
  #                     (1|Track:SPL:Assessor) + (1|Car:SPL:Assessor), 
  #                   data = sound_data_balanced)
  
  m.bo.fixed <- lm(att1 ~ Track*SPL*Car + Assessor + SPL:Assessor + 
                       Track:SPL:Assessor + Car:SPL:Assessor, 
                     data = sound_data_balanced)
  
  
  ##sigma <- summary(m.bo, "lme4")$sigma
  
  deltas <- tapply(sound_data_balanced$att1,
                   factor(sound_data_balanced$Track), mean)
  
  deltadifs=matrix(rep(0,25),ncol=5)
  for (i in 1:5) for (j in 1:i) deltadifs[i,j]=deltas[i]-deltas[j]
  
  stopifnot(all.equal( sqrt(sum(deltadifs^2/anova(m.bo.fixed)["Residuals",3])/(10)), 
             result_bo$fixed$dprimeav[1, "att1"], tol = 1e-3))
  
  ##result_bo$fixed
  ## for Track
  stopifnot(all.equal(sqrt(anova(m.bo.fixed)[1,4]*(2/288)), 
            result_bo$fixed$dprimeav["Track", "att1"], tol = 1e-3))
  ## for SPL
  stopifnot(all.equal(sqrt(anova(m.bo.fixed)[2,4]*(2/480)), 
            result_bo$fixed$dprimeav["SPL", "att1"], tol = 1e-3))
  ## for Car
  stopifnot(all.equal(sqrt(anova(m.bo.fixed)[3,4]*(2/240)), 
            result_bo$fixed$dprimeav["Car", "att1"], tol = 1e-3))
  ## for Track:SPL
  stopifnot(all.equal(sqrt(anova(m.bo.fixed)[5,4]*(2/96)*(8/14)), 
            result_bo$fixed$dprimeav["Track:SPL", "att1"], tol = 1e-3))
  ## for Track:Car
  stopifnot(all.equal(sqrt(anova(m.bo.fixed)[6,4]*(2/48)*(20/29)), 
            result_bo$fixed$dprimeav["Track:Car", "att1"], tol = 1e-3))
  ## for SPL:Car
  stopifnot(all.equal(sqrt(anova(m.bo.fixed)[7,4]*(2/80)*(10/17)), 
            result_bo$fixed$dprimeav["Car:SPL", "att1"], tol = 1e-3))
  ## for Track:SPL:Car
  stopifnot(all.equal(sqrt(anova(m.bo.fixed)[9,4]*(2/16)*(40/89)), 
            result_bo$fixed$dprimeav["Track:Car:SPL", "att1"], tol = 1e-3))
}
