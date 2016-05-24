## ---- echo=FALSE, include=FALSE------------------------------------------
library(heemod)
library(ggplot2)

## ------------------------------------------------------------------------
mat_mono <-
  define_matrix(
    .721, .202, .067, .010,
    .000, .581, .407, .012,
    .000, .000, .750, .250,
    .000, .000, .000, 1.00
  )
mat_mono

## ------------------------------------------------------------------------
rr <- .509

mat_comb <-
  define_matrix(
    C,    .202*rr, .067*rr, .010*rr,
    .000, C,       .407*rr, .012*rr,
    .000, .000,    C,       .250*rr,
    .000, .000,    .000,    1.00
  )
mat_comb

## ---- fig.width = 6, fig.height=6, fig.align='center'--------------------
plot(mat_mono)

## ---- fig.width = 6, fig.height=6, fig.align='center'--------------------
plot(mat_comb)

## ------------------------------------------------------------------------
cost_zido <- 2278
cost_lami <- 2086

## ------------------------------------------------------------------------
A_mono <-
  define_state(
    cost_health = 2756,
    cost_drugs = cost_zido,
    cost_total = discount(cost_health + cost_drugs, .06),
    life_year = 1
  )
A_mono

## ------------------------------------------------------------------------
B_mono <-
  define_state(
    cost_health = 3052,
    cost_drugs = cost_zido,
    cost_total = discount(cost_health + cost_drugs, .06),
    life_year = 1
  )
C_mono <-
  define_state(
    cost_health = 9007,
    cost_drugs = cost_zido,
    cost_total = discount(cost_health + cost_drugs, .06),
    life_year = 1
  )
D_mono <-
  define_state(
    cost_health = 0,
    cost_drugs = 0,
    cost_total = discount(cost_health + cost_drugs, .06),
    life_year = 0
  )

## ------------------------------------------------------------------------
A_comb <-
  define_state(
    cost_health = 3052,
    cost_drugs = cost_zido + cost_lami,
    cost_total = discount(cost_health + cost_drugs, .06),
    life_year = 1
  )
B_comb <-
  define_state(
    cost_health = 3052 + cost_lami,
    cost_drugs = cost_zido,
    cost_total = discount(cost_health + cost_drugs, .06),
    life_year = 1
  )
C_comb <-
  define_state(
    cost_health = 9007 + cost_lami,
    cost_drugs = cost_zido,
    cost_total = discount(cost_health + cost_drugs, .06),
    life_year = 1
  )
D_comb <-
  define_state(
    cost_health = 0,
    cost_drugs = 0,
    cost_total = discount(cost_health + cost_drugs, .06),
    life_year = 0
  )

## ------------------------------------------------------------------------
mod_mono <- define_model(
  transition_matrix = mat_mono,
  A_mono,
  B_mono,
  C_mono,
  D_mono
)
mod_mono

## ------------------------------------------------------------------------
mod_comb <- define_model(
  transition_matrix = mat_comb,
  A_comb,
  B_comb,
  C_comb,
  D_comb
)

## ------------------------------------------------------------------------
res_mod <- run_models(
  mono = mod_mono,
  comb = mod_comb,
  cycles = 20,
  cost = cost_total,
  effect = life_year,
  method = "end"
)

## ------------------------------------------------------------------------
summary(res_mod)

## ---- fig.align='center', fig.width=6, message=FALSE---------------------
plot(res_mod, model = "mono", type = "counts") +
  xlab("Time") +
  theme_minimal() +
  scale_color_brewer(
    name = "State",
    palette = "Set1"
  )

## ---- fig.align='center', fig.width=6, message=FALSE---------------------
plot(res_mod, model = "comb", type = "counts") +
  xlab("Time") +
  theme_minimal() +
  scale_color_brewer(
    name = "State",
    palette = "Set1"
  )

