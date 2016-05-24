data(corona)

# 1. Setup vectors for ei_est_gen; and other functions
cands <- c("pct_husted", "pct_spiegel","pct_ruth","pct_button","pct_montanez","pct_fox" )
race_group3 <- c("~ pct_hisp", "~ pct_asian","~ pct_white") 
table_names <- c("EI: Pct Lat", "EI: Pct Asian", "EI: Pct White")

#2. Run ei_est_gen for EI Results
results <- ei_est_gen(cand_vector=cands, race_group = race_group3,
                      total = "totvote", data = corona, 
                      table_names = table_names, tomog=F)

# EI:RxC

#3. Generate RxC formula, must add to 1 (100%) for both candidates and racial groups
form <- formula(cbind(pct_husted,pct_spiegel,pct_ruth,pct_button,pct_montanez,pct_fox) ~ cbind(pct_hisp, pct_asian, pct_white)) 

#4. Run RxC Bayesian model
suppressWarnings (
ei_bayes <- ei.reg.bayes(form, data=corona, sample=10000, truncate=T)
)
#5. Generate table names vector for RxC model
table_names <- c("RxC: Pct Hisp", "RxC: Pct Asian", "RxC: Pct White")

#6. Create RxC Table
ei_bayes_res <- bayes_table_make(ei_bayes, cands, table_names)

#7. Put EI and EI:RxC together
ei_rc_combine <- ei_rc_good_table(results, ei_bayes_res, groups= c("Latino", "Asian", "White"))

#8. Plot differences
plot(ei_rc_combine) 

