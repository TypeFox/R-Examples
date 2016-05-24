risk <-
function(x, alpha = c(0.05), beta = 1, p = 2)
{
 StD = VaR = EL = ELD = ES = SDR = EVaR = DEVaR = ENT = DENT = ML = rep(0, length(alpha))
 expect <- function(x, alpha)
 {
  int <- function(e)
  {
   ind = ifelse(x < e, 1, 0)
   sum(abs(alpha - ind) * ((x - e) ^ 2))
  }
  optimize(int, c(-20, 20)) $ minimum
 }
 for (i in 1 : length(alpha))
  {
   StD[i] = sd(x)
   VaR[i] = -quantile(x, alpha[i])
   EL[i] = -mean(x)
   ELD[i] = EL[i] + beta * (mean(ifelse(x < -EL[i], (x - (-EL[i]))^p, 0)))^(1/p)
   ES[i] = -mean(x[x < -VaR[i]])
   SDR[i] = ES[i] + beta * (mean(ifelse(x < -ES[i], (x - (-ES[i]))^p, 0)))^(1/p)
   EVaR[i] = -expect(x, alpha[i])
   DEVaR[i] = EVaR[i] + beta * (mean(ifelse(x < -EVaR[i], (x - (-EVaR[i]))^p, 0)))^(1/p)
   ENT[i] = (1 / beta) * log(mean(exp(-beta * x)))
   DENT[i] = ENT[i] + beta * (mean(ifelse(x < -ENT[i], (x - (-ENT[i]))^p, 0)))^(1/p)
   ML [i] = -min(x)
  }
 ans <- rbind(StD, VaR, EL, ELD, ES, SDR, EVaR, DEVaR, ENT, DENT, ML)
 colnames(ans) <- paste(round(100*alpha, 2), "%", sep="")
 return(ans)
}
