Test.Mono <- function(pi1_1_, pi0_1_, pi1_0_, pi_1_1, pi_1_0, pi_0_1){

pi_0_0 <-  1 - (pi_1_1+pi_1_0+pi_0_1)
pi0_0_ <- 1-(pi1_1_+pi0_1_+pi1_0_)

if (pi_0_0 < 0) {stop("\nThe specified marginals for pi_1_1, pi_1_0, and pi_0_1 imply a negative pi_0_0\n")}
if (pi0_0_ < 0) {stop("\nThe specified marginals for pi1_1_, pi0_1_, and pi1_0_ imply a negative pi0_0_\n")}

mono.ST.test <-
  (pi1_1_ <= pi_1_1) & (pi_0_0 <= pi0_0_) & (pi0_1_ <= (pi_0_1+pi_1_1)) & (pi1_0_ <= (pi_1_1+pi_1_0)) & 
  (pi_1_0 <= (pi0_0_+pi1_0_)) & (pi_0_1 <= (pi0_0_+pi0_1_))

mono.S.test <-
  (pi1_1_ <= (pi_0_1 + pi_1_1)) & (pi0_1_ <= (pi_0_1 + pi_1_1)) & (pi_1_0 <= (pi0_0_ + pi1_0_)) & (pi_0_0 <= (pi0_0_ + pi1_0_)) 

mono.T.test <- 
  (pi1_1_ <= (pi_1_0 + pi_1_1)) & (pi1_0_ <= (pi_1_1 + pi_1_0)) & (pi_0_1 <= (pi0_0_ + pi0_1_)) & (pi_0_0 <= (pi0_0_+pi0_1_))

cat("\nThe data are compatible with: \n")
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
cat("\nMonotonicity for both S and T:",  mono.ST.test)
cat("\nMonotonicity for S:",  mono.S.test)
cat("\nMonotonicity for T:",  mono.T.test, "\n\n")

fit <- 
  list(Mono.ST.Test=mono.ST.test, Mono.S.Test=mono.S.test, Mono.T.Test=mono.T.test, Call=match.call())   

class(fit) <- "Test.Mono"
fit

}
