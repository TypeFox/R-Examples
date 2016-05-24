library(MCPMod)
set.seed(1)
# Approximately reproduces analysis from Biometrics paper
# (slightly different linlog model used in paper)
data(biom)
models <- list(emax = 0.2, linlog = NULL, linear = NULL, exponential = c(0.5/log(6),0.15), quadratic = c(-1.7485/2.0485,-1))
MCPMod(biom, models, clinRel = 0.4, dePar = .025, pVal = T, doseEst = "MED2", off = 1, alpha = 0.05)
# some variations
MCPMod(biom, models, clinRel = 0.4, dePar = .05, doseEst = "MED3", off = 1, alpha = 0.05)
MCPMod(biom, models, clinRel = 0.4, dePar = .05, doseEst = "MED1", off = 1, 
       selModel = "AIC", alpha = 0.05)
models <- list(sigEmax = c(0.2,2), logistic = c(0.2,0.1), betaMod = c(1,1))
MCPMod(biom, models, clinRel = 0.4, dePar = .5, doseEst = "ED", off = 1, 
       selModel = "aveBIC", alpha = 0.05, scal = 1.5)

# Examples from JBS paper
doses <- c(0,10,25,50,100,150)
models <- list(linear = NULL, emax = 25,                               
               logistic = c(50, 10.88111), exponential= 85,            
               betaMod=matrix(c(0.33,2.31,1.39,1.39), byrow=TRUE, nrow=2))
planMM(models, doses, n = rep(50,6), alpha = 0.05, scal=200)
#powerMM(models, doses, base = 0, maxEff = 0.4, sigma = 1,       
#        lower = 10, upper = 100, step = 20, scal = 200, alpha = 0.05)
#sampSize(models, doses, base = 0, maxEff = 0.4, sigma = 1,             
#         upperN = 80, scal = 200, alpha = 0.05)
#LP(models, model = "emax", type = "both", paramRange = c(10,70),
#    doses = doses, base = 0, maxEff = 0.4, sigma = 1, n = 60,
#    alpha = 0.05, len = 15, scal = 200)
    
# Example from R News 1(2) p. 28, 29
CM <- c(1, 1, 1, 0, 0, -1, 0, 0, 1, 0, 0, -1, 0, 0,
        1, 0, 0, 0, -1, -1, 0, 0, -1, 0, 0)
CM <- t(matrix(CM, ncol = 5))
critVal(CM, n=c(26, 24, 20, 33, 32), alpha = 0.05, 
        twoSide = TRUE)
