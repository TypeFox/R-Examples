data("Examples", package = "mlogit")
library("lmtest")
library("mlogit")
data("Train", package = "mlogit")
Tr <- mlogit.data(Train, shape = "wide", varying = 4:11, 
                  choice = "choice", sep = "", 
                  opposite = c("price", "time", "change", "comfort"),
                  alt.levels=c("choice1", "choice2"), id="id")

# mixed mogit model :
#  - H0 : uncorrelated model
#waldtest(Train.mxlc, hyp = "uncorrelated")
waldtest(Train.mxlc, correlation = FALSE)
scoretest(Train.mxlu, correlation = TRUE)
lrtest(Train.mxlc, Train.mxlu)

#  - H0 : no random effects vs H1 : correlated random effects
#waldtest(Train.mxlc, hyp = "fixed")
waldtest(Train.mxlc, rpar = NULL)
# ou plus simplement
waldtest(Train.mxlc)
scoretest(Train.ml, rpar = c(time = "n", change = "n", comfort = "n"),
          R = 100, correlation = TRUE, halton = NA, panel = TRUE)
lrtest(Train.mxlc, Train.ml)

#  - H0 : no random effects vs H1 : uncorrelated random effects
waldtest(Train.mxlu, rpar = NULL)
# ou plus simplement
waldtest(Train.mxlu)
scoretest(Train.ml, rpar = c(time = "n", change = "n", comfort = "n"),
          R = 100, correlation = FALSE, halton = NA, panel = TRUE)
lrtest(Train.mxlu, Train.ml)

# Nested logit model

data("HC", package = "mlogit")
HC <- mlogit.data(HC, varying = c(2:8, 10:16), choice = "depvar", shape = "wide")
cooling.modes <- HC$alt %in% c("gcc", "ecc", "erc","hpc")
HC$icca[!cooling.modes] <- HC$occa[!cooling.modes] <- 0
ml <- mlogit(depvar~occa+icca+och+ich, HC)
nl <- mlogit(depvar~occa+icca+och+ich, HC, 
             nests = list(cooling = c('ecc', 'erc', 'gcc', 'hpc'), 
               noncool = c('ec', 'gc', 'er')))
nlu <- update(nl, un.nest.el = TRUE)

# - H0 : unique nest elasticity vs different nest elasticities
scoretest(nlu, un.nest.el = FALSE)
#waldtest(nl, hyp = "un.nest.el")
waldtest(nl, un.nest.el = TRUE)
lrtest(nl, nlu)

# - H0 : unique nest elasticity nested logit vs no nests
scoretest(ml, nests = list(cooling = c('ecc', 'erc', 'gcc', 'hpc'), 
                noncool = c('ec', 'gc', 'er')), un.nest.el = TRUE)
waldtest(nlu, nests = NULL)
# ou plus simplement
waldtest(nlu)
lrtest(nlu, ml)

# - H0 : nested logit vs no nests
scoretest(ml, nests = list(cooling = c('ecc', 'erc', 'gcc', 'hpc'), 
                noncool = c('ec', 'gc', 'er')))
waldtest(nl, nests = NULL)
# ou plus simplement
waldtest(nl)
lrtest(nl, ml)
