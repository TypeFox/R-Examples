require(lmerTest)

# Read in data set
load(system.file("testdata","potdata.RData", package="lmerTest"))

# Mixed model
lmerout <- lmer(biomass ~ CO2*nutrients + (1|chamber),data=potdata)

an.sat <- anova(lmerout)
if(require(pbkrtest))
  an.kr <- anova(lmerout, ddf="Kenward-Roger")

TOL <- 1e-7
stopifnot(all.equal(an.kr[,"Pr(>F)"], c(0.0224955602, 1e-11, 0.020905569) , 
                    tol=TOL),
          all.equal(an.kr[,"DenDF"], 
                    c(2, 10, 10) , tol=TOL),
          TRUE)
## error in Satterthwaite!denom  should be 2 for CO2 like in SAS
TOL <- 1e-5
stopifnot(all.equal(an.sat[,"DenDF"], 
          c(2, 10, 10) , tol=TOL), TRUE)
