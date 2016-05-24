## ---- echo=FALSE, include=FALSE------------------------------------------
library(heemod)

## ------------------------------------------------------------------------
# a function to return age-related mortality rate
# given age and sex
param <- define_parameters(
    age_init = 60,
    sex = 0,
    # age increases with cycles
    age = age_init + markov_cycle,
    
    # operative mortality rates
    omrPTHR = .02,
    omrRTHR = .02,
    
    # re-revision mortality rate
    rrr = .04,
    
    # parameters for calculating primary revision rate
    cons = -5.49094,
    ageC = -.0367,
    maleC = .768536,
    lambda = exp(cons + ageC * age_init + maleC * sex),
    gamma = 1.45367786,
    
    rrNP1 = .260677,
    
    # revision probability of primary procedure
    standardRR = 1 - exp(lambda * ((markov_cycle - 1) ^ gamma -
                                     markov_cycle ^ gamma)),
    np1RR = 1 - exp(lambda * rrNP1 * ((markov_cycle - 1) ^ gamma - 
                                        markov_cycle ^ gamma)),
    
    # age-related mortality rate
    sex_cat = ifelse(sex == 0, "FMLE", "MLE"),
    mr = get_who_mr(age, sex_cat, country = "GBR")
)
param

## ------------------------------------------------------------------------
mat_trans <- define_matrix(
    state_names = c(
      "PrimaryTHR",
      "SuccessP",
      "RevisionTHR",
      "SuccessR",
      "Death"
    ),
    0, C, 0,          0, omrPTHR,
    0, C, standardRR, 0, mr,
    0, 0, 0,          C, omrRTHR+mr,
    0, 0, rrr,        C, mr,
    0, 0, 0,          0, 1
)
mat_trans

## ---- fig.width = 6, fig.height=6, fig.align='center'--------------------
plot(mat_trans)

## ------------------------------------------------------------------------
mod_standard <- define_model(
  transition_matrix = mat_trans,
  PrimaryTHR = define_state(
    utility = 0,
    cost = 394
  ),
  SuccessP = define_state(
    utility = discount(.85, .015),
    cost = 0
  ),
  RevisionTHR = define_state(
    utility = discount(.30, .015),
    cost = discount(5294, .06)
  ),
  SuccessR = define_state(
    utility = discount(.75, .015),
    cost = 0
  ),
  Death = define_state(
    utility = 0,
    cost = 0
  )
)
mod_standard

mod_np1 <- define_model(
  transition_matrix = mat_trans,
  PrimaryTHR = define_state(
    utility = 0,
    cost = 579
  ),
  SuccessP = define_state(
    utility = discount(.85, .015),
    cost = 0
  ),
  RevisionTHR = define_state(
    utility = discount(.30, .015),
    cost = discount(5294, .06)
  ),
  SuccessR = define_state(
    utility = discount(.75, .015),
    cost = 0
  ),
  Death = define_state(
    utility = 0,
    cost = 0
  )
)

## ------------------------------------------------------------------------
res_mod <- run_models(
  standard = mod_standard,
  np1 = mod_np1,
  parameters = param,
  cycles = 60,
  cost = cost,
  effect = utility,
  method = "end"
)

## ------------------------------------------------------------------------
summary(res_mod)

## ---- fig.width = 6, fig.align='center'----------------------------------
plot(res_mod, model = "standard", type = "counts")

## ---- fig.width = 6, fig.align='center'----------------------------------
plot(res_mod, model = "np1", type = "counts")

