penf <- function(est, age, sex, mut, base.dist, agemin){
  H <- cumhaz(base.dist, age-agemin, exp(est[1:2]))*exp(est[3]*sex+est[4]*mut) 
  1-exp(-H)
}