HT_mal <-
function(t, gender, d, smok_status)# Function based on the results from the three cohorts (estimates from CPS-I and CPS-II)
 {
 ### FIXED ###
 X <- 10^7
 ### VARIED ###   
# ----------------- #
# Males/Non-Smokers #
# ----------------- #
 if (gender=="male" & smok_status=="N")
 {
  a0 <- 7.7
  g0 <- 0.09
  v0 <- 7e-8
  p4 <- 0
  p5 <- 0
  p6 <- 0
 }
# ------------------- #
# Females/Non-Smokers #
# ------------------- #
 if (gender=="female" & smok_status=="N")
 {
  a0 <- 15.5
  g0 <- 0.071
  v0 <- 1e-7
  p4 <- 0
  p5 <- 0
  p6 <- 0
 }
# ------------- #
# Males/Smokers #
# ------------- #
 if (gender=="male" & smok_status=="Y")
 {
 a0 <- 7.5
 g0 <- 0.09
 v0 <- 7e-8
 p4 <- 0.01
 p5 <- 0.6
 p6 <- 0.22
 }
# --------------- #
# Females/Smokers #
# --------------- #
 if (gender=="female" & smok_status=="Y")
 {
  a0 <- 15.5
  g0 <- 0.071
  v0 <- 1e-7
  p4 <- 0.02
  p5 <- 0.5
  p6 <- 0.32
 }

 params <- c(a0, g0, v0, p4, p5, p6)

 ### CALCULATED ###
 v    <- v0*(1+p4)
 g    <- g0*(1+p5*d^p6)
 a    <- a0*(1+p5*d^p6)
 m0   <- v0
 m    <- m0
 beta<- a - m - g
 B    <- (1/2)*(-g+sqrt((g^2) + 4*a*m))
 K    <- v*m*X/(g+B)
 H    <- K*(-t + (1/B)*log(g+B+B*exp((2*B+g)*t)))
 return (H)
 }
