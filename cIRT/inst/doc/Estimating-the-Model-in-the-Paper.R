## ----load_pkg------------------------------------------------------------
library(cIRT)

## ----load_data-----------------------------------------------------------
data(trial_matrix)
data(choice_matrix)

## ----thurstone_design----------------------------------------------------
# Create the Thurstone Design Matrices
hard_items = choice_matrix$hard_q_id
easy_items = choice_matrix$easy_q_id

D_easy = model.matrix(~-1+factor(easy_items))
D_hard = -1*model.matrix(~-1+factor(hard_items))[,-c(5,10,15)]

## ----effect_coding-------------------------------------------------------
# Defining effect-coded contrasts
high_contrasts <- rbind(-1,diag(4))
rownames(high_contrasts) = 12:16
low_contrasts <- rbind(-1,diag(2))
rownames(low_contrasts) = 4:6

# Creating high & low factors
high = factor(choice_matrix[,'high_value'])
low = factor(choice_matrix[,'low_value'])
contrasts(high) = high_contrasts
contrasts(low) = low_contrasts

fixed_effects = model.matrix(~high+low)
fixed_effects_base = fixed_effects[,1]
fixed_effects_int = model.matrix(~high*low)

## ----model_data----------------------------------------------------------
# Model with Thurstone D matrix
system.time({ 
  out_model_thurstone = cIRT(choice_matrix[,'subject_id'], 
                             cbind(fixed_effects[,-1],D_easy,D_hard),
                             c(1:ncol(fixed_effects)), 
                             as.matrix(fixed_effects),
                             as.matrix(trial_matrix), 
                             choice_matrix[,'choose_hard_q'],
                             20000,
                             25000)
})

## ----param_ests----------------------------------------------------------
vlabels_thurstone = colnames(cbind(fixed_effects[,-1],D_easy,D_hard))
G_thurstone = t(apply(out_model_thurstone$gs0, 2, FUN = quantile,probs=c(.5,.025,.975)))
rownames(G_thurstone)=vlabels_thurstone
B_thurstone = t(apply(out_model_thurstone$beta, 2, FUN = quantile,probs=c(.5,0.025,.975)))
rownames(B_thurstone)=colnames(fixed_effects)

S_thurstone = solve(apply(out_model_thurstone$Sigma_zeta_inv, c(1,2), FUN = mean))
inv_sd = diag(1/sqrt(diag(solve(apply(out_model_thurstone$Sigma_zeta_inv, c(1,2), FUN = mean)))))
corrmat = inv_sd%*%S_thurstone%*%inv_sd
as = apply(out_model_thurstone$as, 2, FUN = mean)
bs = apply(out_model_thurstone$bs, 2, FUN = mean)

## ----param_results-------------------------------------------------------
# gs0
G_thurstone

# betas
B_thurstone

# Sigma Thurstone
S_thurstone

# item parameters

# a
as

# b
bs

