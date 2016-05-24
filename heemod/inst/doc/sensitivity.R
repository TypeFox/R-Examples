## ---- echo=FALSE, include=FALSE------------------------------------------
library(heemod)

param <- define_parameters(
  rr = .509,
  
  p_AB_base = .202,
  p_AC_base = .067,
  p_AD_base = .010,
  
  p_BC_base = .407,
  p_BD_base = .012,
  
  p_CD_base = .250,
  
  
  p_AB = p_AB_base,
  p_AC = p_AC_base,
  p_AD = p_AD_base,
  
  p_BC = p_BC_base,
  p_BD = p_BD_base,
  
  p_CD = p_CD_base,
  
  
  cost_zido = 2278,
  cost_lami = 2086,
  
  cost_A = 2756,
  cost_B = 3052,
  cost_C = 9007,
  
  dr = .06
)

mat_trans_mono <- define_matrix(
  C,    p_AB, p_AC, p_AD,
  .000, C,    p_BC, p_BD,
  .000, .000, C,    p_CD,
  .000, .000, .000, 1.00
)

mat_trans_comb <- define_matrix(
  C,    p_AB * rr, p_AC * rr, p_AD * rr,
  .000, C,    p_BC * rr, p_BD * rr,
  .000, .000, C,    p_CD * rr,
  .000, .000, .000, 1.00
)

A_mono <- define_state(
  cost_health = cost_A,
  cost_drugs = cost_zido,
  cost_total = discount(cost_health + cost_drugs, dr),
  life_year = 1
)
B_mono <- define_state(
  cost_health = cost_B,
  cost_drugs = cost_zido,
  cost_total = discount(cost_health + cost_drugs, dr),
  life_year = 1
)
C_mono <- define_state(
  cost_health = cost_C,
  cost_drugs = cost_zido,
  cost_total = discount(cost_health + cost_drugs, dr),
  life_year = 1
)
D_mono <- define_state(
  cost_health = 0,
  cost_drugs = 0,
  cost_total = discount(cost_health + cost_drugs, dr),
  life_year = 0
)

A_comb <- define_state(
  cost_health = cost_A,
  cost_drugs = cost_zido + cost_lami,
  cost_total = discount(cost_health + cost_drugs, dr),
  life_year = 1
)
B_comb <- define_state(
  cost_health = cost_B,
  cost_drugs = cost_zido + cost_lami,
  cost_total = discount(cost_health + cost_drugs, dr),
  life_year = 1
)
C_comb <- define_state(
  cost_health = cost_C,
  cost_drugs = cost_zido + cost_lami,
  cost_total = discount(cost_health + cost_drugs, dr),
  life_year = 1
)
D_comb <- define_state(
  cost_health = 0,
  cost_drugs = 0,
  cost_total = discount(cost_health + cost_drugs, dr),
  life_year = 0
)

mod_mono <- define_model(
  transition_matrix = mat_trans_mono,
  A_mono,
  B_mono,
  C_mono,
  D_mono
)

mod_comb <- define_model(
  transition_matrix = mat_trans_comb,
  A_comb,
  B_comb,
  C_comb,
  D_comb
)

res_mod <- run_models(
  mono = mod_mono,
  comb = mod_comb,
  parameters = param,
  cycles = 20,
  cost = cost_total,
  effect = life_year
)

## ------------------------------------------------------------------------
se <- define_sensitivity(
  rr = c(.4, .6),
  
  cost_zido = c(1500, 3000),
  cost_lami = c(1500, 3000),
  
  dr = c(.04, .08)
)

## ------------------------------------------------------------------------
res <- run_sensitivity(
  model = res_mod,
  sensitivity = se
)

## ---- fig.width = 6, fig.align='center'----------------------------------
plot(res, model = "mono")

## ---- fig.width = 6, fig.align='center'----------------------------------
plot(res, model = "comb")

## ---- fig.width = 6, fig.align='center'----------------------------------
plot(res, model = "comb", type = "diff")

