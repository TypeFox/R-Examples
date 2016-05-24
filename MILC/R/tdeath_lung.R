tdeath_lung <-
function(u1, u2, covs_lung)
 {
temp.MILC.Env <- new.env()
data(ci.lung, current.other, former.other, never.other, envir=temp.MILC.Env)

#ci.lung      <- get("ci.lung",        envir=temp.MILC.Env)
#curren.other <- get("current.other",  envir=temp.MILC.Env)
#former.other <- get("former.other",   envir=temp.MILC.Env)
#never.other  <- get("never.other",    envir=temp.MILC.Env)

 obj       <- get("ci.lung", envir=temp.MILC.Env)
 stage     <- covs_lung[1]
 age_curr  <- as.numeric(covs_lung[2])# current age
 age_group <- age_grp(age_curr)# current age group
 gender    <- covs_lung[3]
 
 ind2 <- grep(paste(stage, age_group, gender), names(obj))

 if(u1 > max(obj[[ind2]][[2]])){t_int_ld <- age_curr + max(obj[[ind2]][[1]])/12 ; event <- 0} else {
low  <- which( unique(abs(u1 - obj[[ind2]][[2]])) == min(unique(abs(u1 - obj[[ind2]][[2]]))) )
up   <- low+1
if (   up > length( unique(obj[[ind2]][[1]])) ) { t_int_ld <- age_curr + unique(obj[[ind2]][[1]])[low]/12 ; event <- 0  } else {
int_ld <- unique(obj[[ind2]][[1]])[c(low,up)] # if the estimated age at diagnosis is larger than the max(age) 
# for which predictions can be made, the predicted time to death is censored
t_int_ld <- (u2 * (int_ld[2] - int_ld[1]) + int_ld[1])/12 +age_curr  ; event <- 1 }}
 return(list(u1, u2, ind2, "time"=t_int_ld, "event"=event, obj[ind2]))
 }
