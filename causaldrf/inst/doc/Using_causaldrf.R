## ----include=FALSE-------------------------------------------------------
library(knitr)
opts_chunk$set(
concordance=TRUE
)

## ----echo=FALSE-----------------------------------------------------
options(width=70)

## -------------------------------------------------------------------
library (causaldrf)

## -------------------------------------------------------------------
set.seed(301)
hi_sample <- function(N){
  X1 <- rexp(N)
  X2 <- rexp(N)
  T <- rexp(N, X1 + X2)
  gps <- (X1 + X2) * exp(-(X1 + X2) * T)
  Y <- T + gps + rnorm(N)
  hi_data <- data.frame(cbind(X1, X2, T, gps, Y))
  return(hi_data)
}


hi_sim_data <- hi_sample(1000)
head(hi_sim_data)

## ----echo= FALSE----------------------------------------------------

overlap_list <- overlap_fun(Y = Y,
                            treat = T,
                            treat_formula = T ~ X1 + X2,
                            data_set = hi_sim_data,
                            n_class = 3,
                            treat_mod = "Gamma",
                            link_function = "inverse")

overlap_data <- overlap_list$overlap_dataset


t_mod_list <- t_mod(treat = T,
                    treat_formula = T ~ X1 + X2 + I(X1^2) + I(X2^2),
                    data = overlap_data,
                    treat_mod = "Gamma",
                    link_function = "inverse")
cond_exp_data <- t_mod_list$T_data
full_data <- cbind(overlap_data, cond_exp_data)

## -------------------------------------------------------------------
add_spl_estimate <- add_spl_est(Y = Y,
                                treat = T,
                                treat_formula = T ~ X1 + X2,
                                data = hi_sim_data,
                                grid_val = quantile(hi_sim_data$T,
                                            probs = seq(0, .95, by = 0.01)),
                                knot_num = 3,
                                treat_mod = "Gamma",
                                link_function = "inverse")

## -------------------------------------------------------------------
gam_estimate <- gam_est(Y = Y,
                        treat = T,
                        treat_formula = T ~ X1 + X2,
                        data = hi_sim_data,
                        grid_val = quantile(hi_sim_data$T,
                                    probs = seq(0, .95, by = 0.01)),
                        treat_mod = "Gamma",
                        link_function = "inverse")

## -------------------------------------------------------------------
hi_estimate <- hi_est(Y = Y,
                      treat = T,
                      treat_formula = T ~ X1 + X2,
                      outcome_formula = Y ~ T + I(T^2) +
                        gps + I(gps^2) + T * gps,
                      data = hi_sim_data,
                      grid_val = quantile(hi_sim_data$T,
                                  probs = seq(0, .95, by = 0.01)),
                      treat_mod = "Gamma",
                      link_function = "inverse")

## -------------------------------------------------------------------
iptw_estimate <- iptw_est(Y = Y,
                          treat = T,
                          treat_formula = T ~ X1 + X2,
                          numerator_formula = T ~ 1,
                          data = hi_sim_data,
                          degree = 2,
                          treat_mod = "Gamma",
                          link_function = "inverse")

## ----echo=FALSE-----------------------------------------------------
load("nmes_10192015.RData")

## -------------------------------------------------------------------
data("nmes_data")
dim (nmes_data)
summary(nmes_data)

## -------------------------------------------------------------------
t(p_val_bal_cond)
t(p_val_bal_no_cond)

## -------------------------------------------------------------------
pf_estimate <- reg_est(Y = TOTALEXP,
                       treat = packyears,
                       covar_formula = ~ 1,
                       data = full_data_orig,
                       degree = 2,
                       wt = full_data_orig$HSQACCWT,
                       method = "same")
pf_estimate

## -------------------------------------------------------------------
reg_estimate <- reg_est(Y = TOTALEXP,
                        treat = packyears,
                        covar_formula = ~ LASTAGE + LASTAGE2 +
                          AGESMOKE + AGESMOKE2 + MALE + beltuse +
                          educate + marital + POVSTALB + RACE3,
                        covar_lin_formula = ~ 1,
                        covar_sq_formula = ~ 1,
                        data = full_data_orig,
                        degree = 2,
                        wt = full_data_orig$HSQACCWT,
                        method = "different")
reg_estimate

## -------------------------------------------------------------------
spline_estimate <- prop_spline_est(Y = TOTALEXP,
                                   treat = packyears,
                                   covar_formula = ~ LASTAGE + LASTAGE2 +
                                     AGESMOKE + AGESMOKE2 + MALE + beltuse +
                                     educate + marital + POVSTALB + RACE3,
                                   covar_lin_formula = ~ 1,
                                   covar_sq_formula = ~ 1,
                                   data = full_data_orig,
                                   e_treat_1 = full_data_orig$est_treat,
                                   degree = 2,
                                   wt = full_data_orig$HSQACCWT,
                                   method = "different",
                                   spline_df = 5,
                                   spline_const = 4,
                                   spline_linear = 4,
                                   spline_quad = 4)
spline_estimate

## -------------------------------------------------------------------
ivd_estimate <- prop_spline_est(Y = TOTALEXP,
                                treat = packyears,
                                covar_formula = ~ 1,
                                covar_lin_formula = ~ 1,
                                covar_sq_formula = ~ 1,
                                data = full_data_orig,
                                e_treat_1 = full_data_orig$est_treat,
                                degree = 2,
                                wt = full_data_orig$HSQACCWT,
                                method = "different",
                                spline_df = 5,
                                spline_const = 4,
                                spline_linear = 4,
                                spline_quad = 4)
ivd_estimate

## ----eval = FALSE, echo = TRUE--------------------------------------
#  library(Hmisc)
#  mydata <- sasxport.get("09795-0141-Data-card_image.xpt")
#  data_58 <- mydata[[58]]
#  ihdp_raw <- data_58
#  # restricts data to treated cases
#  treated_raw <- ihdp_raw[which(ihdp_raw$tg == "I"),]
#  # continuous treatment variable
#  treat_value <- treated$cdays.t

## ----eval = FALSE, echo = TRUE--------------------------------------
#  overlap_temp <- overlap_fun(Y = iqsb.36,
#                              treat = ncdctt,
#                              treat_formula = t_formula,
#                              data = data_set,
#                              n_class = 3,
#                              treat_mod = "Normal")
#  
#  median_list <- overlap_temp[[2]]
#  overlap_orig <- overlap_temp[[1]]
#  overlap_3 <- overlap_temp[[3]]
#  fitted_values_overlap <- overlap_3$fitted.values

## ----eval = FALSE, echo = TRUE--------------------------------------
#  bart_estimate <-  bart_est(Y = iqsb.36,
#                               treat = ncdctt,
#                               outcome_formula = iqsb.36 ~ ncdctt + bw +
#                               female + mom.lths +
#                               site1 + site7 + momblack +
#                               workdur.imp,
#                               data = full_data_orig,
#                               grid_val = grid_treat)

## ----eval = FALSE, echo = TRUE--------------------------------------
#  iw_estimate <- iw_est(Y = iqsb.36,
#                         treat = ncdctt,
#                         treat_formula = ncdctt ~ bw + female + mom.lths +
#                          site1 + site7 + momblack +
#                          workdur.imp,
#                         data = full_data_orig,
#                         grid_val = grid_treat,
#                         bandw = 2 * bw.SJ(full_data_orig$ncdctt),
#                         treat_mod = "Normal")

## ----eval = FALSE, echo = TRUE--------------------------------------
#  nw_estimate <- nw_est(Y = iqsb.36,
#                        treat = ncdctt,
#                        treat_formula = ncdctt ~ bw + female + mom.lths +
#                          site1 + site7 + momblack +
#                          workdur.imp,
#                        data = full_data_orig,
#                        grid_val = grid_treat,
#                        bandw = 2 * bw.SJ(full_data_orig$ncdctt),
#                        treat_mod = "Normal")

## ----eval = FALSE, echo = TRUE--------------------------------------
#  spline_estimate <- prop_spline_est(Y = iqsb.36,
#                                     treat = ncdctt,
#                                     covar_formula = ~ bw + female +
#                                       mom.lths + site1 + site7 +
#                                       momblack + workdur.imp,
#                                     covar_lin_formula = ~ 1,
#                                     covar_sq_formula = ~ 1,
#                                     data = full_data_orig,
#                                     e_treat_1 = full_data_orig$est_treat,
#                                     degree = 2,
#                                     wt = NULL,
#                                     method = "different",
#                                     spline_df = 5,
#                                     spline_const = 2,
#                                     spline_linear = 2,
#                                     spline_quad = 2)

