## ---- echo=FALSE---------------------------------------------------------
set.seed(1)

## ------------------------------------------------------------------------
library(saeSim)
setup <- sim_base() %>% 
  sim_gen_x() %>% 
  sim_gen_e() %>% 
  sim_gen_v() %>%
  sim_resp_eq(y = 100 + 2 * x + v + e) %>% 
  sim_simName("Doku")
setup

## ----eval=FALSE----------------------------------------------------------
#  dataList <- sim(setup, R = 10)

## ----eval=FALSE----------------------------------------------------------
#  simData <- sim_base() %>%
#    sim_gen_x() %>%
#    sim_gen_e() %>%
#    as.data.frame
#  simData

## ---- eval=FALSE---------------------------------------------------------
#  setup <- sim_base() %>%
#    sim_gen_x() %>%
#    sim_gen_e() %>%
#    sim_resp_eq(y = 100 + 2 * x + e)
#  
#  setup1 <- setup %>% sim_sample(sample_fraction(0.05))
#  setup2 <- setup %>% sim_sample(sample_number(5))

## ----echo=FALSE, dpi=100-------------------------------------------------
# This is too large in terms of MB for a vignette inside a package.
# library(DiagrammeR)
# grViz("
# digraph boxes_and_circles {
# node [shape = box]
# sim_base
# sim_gen
# sim_sample
# sim_comp_sample
# sim_agg
# 
# sim_base -> sim_gen [label = '       add variables to data in sim_base']; 
# sim_gen -> sim_sample [label = '    draw a sample'];
# sim_sample -> sim_comp_sample [label = '      make predictions']; 
# sim_comp_sample -> sim_agg [label = '     aggregate your data'];
# 
# }
# ")

## ------------------------------------------------------------------------
setup <- sim_base_lmm()
plot(setup)
autoplot(setup)
autoplot(setup, "e")
autoplot(setup %>% sim_gen_vc())

## ------------------------------------------------------------------------
base_id(2, 3) %>% 
  sim_gen(gen_generic(rnorm, mean = 5, sd = 10, name = "x", groupVars = "idD"))

## ------------------------------------------------------------------------
library(saeSim)
setup <- sim_base() %>% 
  sim_gen_x() %>% # Variable 'x'
  sim_gen_e() %>% # Variable 'e'
  sim_gen_v() %>% # Variable 'v' as a random-effect
  sim_gen(gen_v_sar(name = "vSp")) %>% # random-effect following a SAR(1)
  sim_resp_eq(y = 100 + x + v + vSp + e) # Computing 'y'
setup

## ------------------------------------------------------------------------
contSetup <- setup %>% 
  sim_gen_cont(
    gen_v_sar(sd = 40, name = "vSp"), # defining the model
    nCont = 0.05, # 5 per cent outliers
    type = "area", # whole areas are outliers, i.e. all obs within
    areaVar = "idD", # var name to identify domain
    fixed = TRUE # if in each iteration the same area is an outlier
  )

## ------------------------------------------------------------------------
base_id(3, 4) %>% 
  sim_gen_x() %>% 
  sim_gen_e() %>% 
  sim_gen_ec(mean = 0, sd = 150, name = "eCont", nCont = c(1, 2, 3)) %>%
  as.data.frame

## ------------------------------------------------------------------------
base_id(2, 3) %>% 
  sim_gen_x() %>% 
  sim_gen_e() %>% 
  sim_gen_ec() %>% 
  sim_resp_eq(y = 100 + x + e) %>%
   # the mean in each domain:
  sim_comp_pop(comp_var(popMean = mean(y)), by = "idD")

## ------------------------------------------------------------------------
comp_linearPredictor <- function(dat) {
  dat$linearPredictor <- lm(y ~ x, dat) %>% predict
  dat
}

sim_base_lm() %>% 
  sim_comp_pop(comp_linearPredictor)

## ------------------------------------------------------------------------
sim_base_lm() %>% 
  sim_comp_pop(function(dat) lm(y ~ x, dat)) %>%
  sim(R = 1)

comp_linearModelAsAttr <- function(dat) {
  attr(dat, "linearModel") <- lm(y ~ x, dat)
  dat
}

dat <- sim_base_lm() %>% 
  sim_comp_pop(comp_linearModelAsAttr) %>%
  as.data.frame

attr(dat, "linearModel")

## ------------------------------------------------------------------------
sim_base_lm() %>% 
  sim_sample() %>%
  sim_comp_sample(comp_linearPredictor)

## ------------------------------------------------------------------------
sim_base_lm() %>% 
  sim_sample() %>%
  sim_agg() %>% 
  sim_comp_agg(comp_linearPredictor)

## ------------------------------------------------------------------------
base_id(3, 4) %>% 
  sim_gen_x() %>% 
  sim_sample(sample_number(1L))
base_id(3, 4) %>% 
  sim_gen_x() %>% 
  sim_sample(sample_number(1L, groupVars = "idD"))

## ------------------------------------------------------------------------
# simple random sampling:
sim_base_lm() %>% sim_sample(sample_number(size = 10L))
sim_base_lm() %>% sim_sample(sample_fraction(size = 0.05))
# srs in each domain/cluster
sim_base_lm() %>% sim_sample(sample_number(size = 10L, groupVars = "idD"))
sim_base_lm() %>% sim_sample(sample_fraction(size = 0.05, groupVars = "idD"))

