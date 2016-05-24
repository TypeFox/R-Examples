epi.pooled <- function(se, sp, P, m, r){
   # Herd specificity:
   PlSp <- sp^m
   HSp <- (PlSp)^r

   # Herd sensitivity:
   HSe <- 1 - ((1 - (1 - P)^m) * (1 - se) + (1 - P)^m * PlSp)^r
   
   # Herd level apparent prevalence:
   HAPneg <- 1 - HSp

   rval <- list(HAPneg = HAPneg, HSe = HSe, HSp = HSp)
   rval
}

 