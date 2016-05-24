tdeath_other <-
function(u1, u2, status, covs_other)
 {
temp.MILC.Env <- new.env()
data(ci.lung, current.other, former.other, never.other, envir=temp.MILC.Env)

#ci.lung      <- get("ci.lung",        envir=temp.MILC.Env)
#curren.other <- get("current.other",  envir=temp.MILC.Env)
#former.other <- get("former.other",   envir=temp.MILC.Env)
#never.other  <- get("never.other",    envir=temp.MILC.Env)


  if (status == "never")  {obj <- get("never.other", envir=temp.MILC.Env)} else 
  if (status == "former") {obj <- get("former.other", envir=temp.MILC.Env)} else obj <- get("current.other", envir=temp.MILC.Env)
  age_curr  <- as.numeric(covs_other[1])# current age
  age_group <- age_grp(age_curr)# current age group
  gender <- covs_other[2]
  intensity <- covs_other[3]

  if(status!="current") { ind1 <- grep(paste(age_group, gender), names(obj)) 
} else ind1 <- grep(paste(age_group, gender, intensity), names(obj)) 
 
  if(u1 > max(obj[[ind1]][[2]])) {low <- age_curr + max(obj[[ind1]][[1]]) ; up <- 110 ; int_od <- c(low, up)
    t_int_od <- ( u2 * (int_od[2] - int_od[1]) ) + age_curr} else {
 low  <- which( unique(abs(u1 - obj[[ind1]][[2]])) == min(unique(abs(u1 - obj[[ind1]][[2]]))) )
 up   <- low+1
 int_od <- unique(obj[[ind1]][[1]])[c(low,up)] 
  t_int_od <- (u2 * (int_od[2] - int_od[1]) + int_od[1] ) + age_curr  }
 return(list(u1, u2, ind1, int_od, "time"=t_int_od, obj[ind1]))
 }
