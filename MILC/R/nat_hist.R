nat_hist <-
function(dat, pred_yrs, gender, status, ts, tq, m, cdiagn, creg, cdist)
 {
temp.MILC.Env <- new.env()
data(ci.lung, current.other, former.other, never.other, envir=temp.MILC.Env)

ci.lung      <- get("ci.lung",        envir=temp.MILC.Env)
curren.other <- get("current.other",  envir=temp.MILC.Env)
former.other <- get("former.other",   envir=temp.MILC.Env)
never.other  <- get("never.other",    envir=temp.MILC.Env)

 uv <- dat[1:5] 
 age_curr <- dat[6]
 T_curr <- age_curr
 d <- dat[7]
 T_mal <- t_mal(uv[1], gender, ts, tq, d)
 co <- c(age_curr, gender, d_grp(d))
 T_do <- tdeath_other(uv[2], uv[3], status, covs_other=co)$time
 Tprog <- t_prog(1, m, cdiagn, creg, cdist)

 T_reg<- Tprog$Treg   + T_mal
 T_dist<- Tprog$Tdist  + T_mal
 T_diagn<- Tprog$Tdiagn + T_mal
 D_diagn<- Tprog$Ddiagn
 stage<- Tprog$stage

 T_pred <- age_curr + pred_yrs
 cl <- c(stage, T_diagn, gender)
 T_dl <- tdeath_lung(uv[4], uv[5], cl)$time

 if(is.na(T_do) | is.na(T_dl)) {T_final <- NA ; lung_inc <- NA ; excl <- NA ; death <- NA ; T_death <- NA; cause<- NA} else {
 if(T_do < T_dl) {cause <- "other" ; T_death <- T_do}  else {cause <- "lung" ; T_death <- T_dl} 
 death <- 1*(T_death < T_pred)
 T_final <- T_death*death + T_pred*(1-death)
 excl <- (1)*(T_death<T_curr)# exclude these erroneous cases
 lung_inc <- 1*(T_diagn<T_final)# indicator variable for the development of lung cancer before 2006
                                                                                 }
 natural_history <- list("T_entry"=age_curr, "T_mal"=T_mal, "T_reg"=T_reg, "T_dist"=T_dist, 
"T_diagn"=T_diagn, "D_diagn"=D_diagn, "stage"=stage,
"T_pred"=T_pred, 
"T_do"=T_do, "T_dl"=T_dl,
"T_final"=T_final,
"lung_inc"=lung_inc, "excl"=excl,
"cause"=cause, "T_death"=T_death, 
"gender"=gender, 
"status"=status, "ts"=ts, "tq"=tq, "intensity"=d)
 return(natural_history)
 }
