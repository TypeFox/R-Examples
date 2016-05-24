## ---- echo=FALSE, include=FALSE------------------------------------------
library(heemod)

## ------------------------------------------------------------------------
param <- define_parameters(
  rr = .509,
  
  p_AA_base = .721,
  p_AB_base = .202,
  p_AC_base = .067,
  p_AD_base = .010,
  
  p_BC_base = .407,
  p_BD_base = .012,
  
  p_CD_base = .250,
  
  
  p_AB_comb = p_AB_base * rr,
  p_AC_comb = p_AC_base * rr,
  p_AD_comb = p_AD_base * rr,
  
  p_BC_comb = p_BC_base * rr,
  p_BD_comb = p_BD_base * rr,
  
  p_CD_comb = p_CD_base * rr,
  
  p_AA_comb = 1 - (p_AB_comb + p_AC_comb + p_AD_comb),
  
  
  cost_zido = 2278,
  cost_lami = 2086,
  
  cost_A = 2756,
  cost_B = 3052,
  cost_C = 9007
)

## ------------------------------------------------------------------------
mat_trans_mono <- define_matrix(
  p_AA_base, p_AB_base, p_AC_base, p_AD_base,
  .000, C,    p_BC_base, p_BD_base,
  .000, .000, C,    p_CD_base,
  .000, .000, .000, 1.00
)
mat_trans_comb <- define_matrix(
  p_AA_comb, p_AB_comb, p_AC_comb, p_AD_comb,
  .000, C,    p_BC_comb, p_BD_comb,
  .000, .000, C,    p_CD_comb,
  .000, .000, .000, 1.00
)

## ------------------------------------------------------------------------
A_mono <- define_state(
  cost_health = cost_A,
  cost_drugs = cost_zido,
  cost_total = discount(cost_health + cost_drugs, .06),
  life_year = 1
)
B_mono <- define_state(
  cost_health = cost_B,
  cost_drugs = cost_zido,
  cost_total = discount(cost_health + cost_drugs, .06),
  life_year = 1
)
C_mono <- define_state(
  cost_health = cost_C,
  cost_drugs = cost_zido,
  cost_total = discount(cost_health + cost_drugs, .06),
  life_year = 1
)
D_mono <- define_state(
  cost_health = 0,
  cost_drugs = 0,
  cost_total = discount(cost_health + cost_drugs, .06),
  life_year = 0
)

A_comb <- define_state(
  cost_health = cost_A,
  cost_drugs = cost_zido + cost_lami,
  cost_total = discount(cost_health + cost_drugs, .06),
  life_year = 1
)
B_comb <- define_state(
  cost_health = cost_B,
  cost_drugs = cost_zido + cost_lami,
  cost_total = discount(cost_health + cost_drugs, .06),
  life_year = 1
)
C_comb <- define_state(
  cost_health = cost_C,
  cost_drugs = cost_zido + cost_lami,
  cost_total = discount(cost_health + cost_drugs, .06),
  life_year = 1
)
D_comb <- define_state(
  cost_health = 0,
  cost_drugs = 0,
  cost_total = discount(cost_health + cost_drugs, .06),
  life_year = 0
)

## ------------------------------------------------------------------------
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
rsp <- define_distrib(
  rr ~ lognormal(mean = .509, sdlog = .173),
  
  cost_A ~ make_gamma(mean = 2756, sd = sqrt(2756)),
  cost_B ~ make_gamma(mean = 3052, sd = sqrt(3052)),
  cost_C ~ make_gamma(mean = 9007, sd = sqrt(9007)),
  
  p_CD_base ~ prop(prob = .25, size = 40),
  
  p_AA_base +
    p_AB_base +
    p_AC_base + 
    p_AD_base ~ multinom(721, 202, 67, 10)
)

## ------------------------------------------------------------------------
pm <- run_probabilistic(
  model = res_mod,
  resample = rsp,
  N = 100
)

## ---- fig.width = 6, fig.height=4, fig.align='center'--------------------
plot(pm, type = "ce")

## ---- fig.width = 6, fig.align='center'----------------------------------
plot(pm, type = "ac", values = seq(0, 25e3, 1e3))

## ---- fig.align='center', fig.height=4, fig.width=6, message=FALSE-------
library(ggplot2)

plot(pm, type = "ce") +
  xlab("Life-years gained") +
  ylab("Additional cost") +
  scale_color_brewer(
    name = "Strategy",
    palette = "Set1"
  ) +
  theme_minimal()

