# Probleme mit der Bewertung der Partial-time-end Down-and-out-Puts
#
# Literatur/Anmerkungen:
# 1. Haug (2007), S.168 Barrier Option Symmetries
# 2. Quellcode Visual Basic vs. Quellcode R-Paket fExoticOptions
# 3. Für Start des Beobachtungszeitraums t1 --> 0 müsste sich der Preis eines normalen Down-and-out-Puts ergeben
# 4. Hui (1997) : Time-Dependent Barrier Option Values, Journal of Futures Markets, pp. 667-688
#library(fExoticOptions)

PTSingleAssetBarrierOption2 <- function (TypeFlag = c("cdoA", "cuoA", "pdoA", "puoA", "coB1", 
    "poB1", "cdoB2", "cuoB2"), S, X, H, time1, Time2, r, b, sigma, 
    title = NULL, description = NULL)
{
  t1 = time1
  T2 = Time2
  
  if (TypeFlag == "cdoA") 
        eta = 1
  if (TypeFlag == "cuoA") 
      eta = -1
  d1 = (log(S/X) + (b + sigma^2/2) * T2)/(sigma * sqrt(T2))
  d2 = d1 - sigma * sqrt(T2)
  f1 = (log(S/X) + 2 * log(H/S) + (b + sigma^2/2) * T2)/(sigma * sqrt(T2))
  f2 = f1 - sigma * sqrt(T2)
  
  e1 = (log(S/H) + (b + sigma^2/2) * t1)/(sigma * sqrt(t1))
  e2 = e1 - sigma * sqrt(t1)
  e3 = e1 + 2 * log(H/S)/(sigma * sqrt(t1))
  e4 = e3 - sigma * sqrt(t1)
  mu = (b - sigma^2/2)/sigma^2
  rho = sqrt(t1/T2)
  g1 = (log(S/H) + (b + sigma^2/2) * T2)/(sigma * sqrt(T2))
  g2 = g1 - sigma * sqrt(T2)
  g3 = g1 + 2 * log(H/S)/(sigma * sqrt(T2))
  g4 = g3 - sigma * sqrt(T2)
    
  z1 = CND(e2) - (H/S)^(2 * mu) * CND(e4)
  z2 = CND(-e2) - (H/S)^(2 * mu) * CND(-e4)
  z3 = CBND(g2, e2, rho) - (H/S)^(2 * mu) * CBND(g4, -e4, -rho)
  z4 = CBND(-g2, -e2, rho) - (H/S)^(2 * mu) * CBND(-g4, e4, -rho)
  z5 = CND(e1) - (H/S)^(2 * (mu + 1)) * CND(e3)
  z6 = CND(-e1) - (H/S)^(2 * (mu + 1)) * CND(-e3)
  z7 = CBND(g1, e1, rho) - (H/S)^(2 * (mu + 1)) * CBND(g3, -e3, -rho)
  z8 = CBND(-g1, -e1, rho) - (H/S)^(2 * (mu + 1)) * CBND(-g3, e3, -rho)
  
  # Nutze die Symmetrie von Haug (2007), S.168
  #if (TypeFlag == "poB1") {
  #      PartialTimeBarrier = PTSingleAssetBarrierOption("coB1", 
  #          X, S, S*X/H, t1, T2, r+b, -b, sigma)@price
  #}
  #if (TypeFlag == "pdoB2") {
  #      PartialTimeBarrier = PTSingleAssetBarrierOption("cuoB2", 
  #          X, S, S*X/H, t1, T2, r+b, -b, sigma)@price
  #}
  
  # Haug Code aus dem Buch
  if (TypeFlag == "cdoB2" && X < H) {
        PartialTimeBarrier = S * exp((b - r) * T2) * (CBND(g1, e1, rho) - (H/S)^(2 * (mu + 1)) * CBND(g3, -e3, -rho)) -X * exp(-r * T2) * (CBND(g2, e2, rho) - (H/S)^(2 * mu) * CBND(g4, -e4, -rho))
  }
  else if (TypeFlag == "cdoB2" && X > H) {
    PartialTimeBarrier = PTSingleAssetBarrierOption2("coB1", S, X, H, t1, T2, r, b, sigma)@price
  }
  else if (TypeFlag == "poB1") {  # put out type B1
    PartialTimeBarrier = PTSingleAssetBarrierOption2("coB1", S, X, H, t1, T2, r, b, sigma)@price - S * exp((b - r) * T2) * z8 + X * exp(-r * T2) * z4 - S * exp((b - r) * T2) * z7 + X * exp(-r * T2) * z3
  }
  else if (TypeFlag == "pdoB2") {  # put down-and-out type B2)
    PartialTimeBarrier = 
      PTSingleAssetBarrierOption2("cdoB2", S, X, H, t1, T2, r, b, sigma)@price - S * exp((b - r) * T2) * z7 + X * exp(-r * T2) * z3
  }
  else if (TypeFlag == "coB1" && X > H) {  # call out type B1
    PartialTimeBarrier = S * exp((b - r) * T2) * (CBND(d1, e1, rho) - (H / S) ^ (2 * (mu + 1)) * CBND(f1, -e3, -rho)) - X * exp(-r * T2) * (CBND(d2, e2, rho) - (H / S) ^ (2 * mu) * CBND(f2, -e4, -rho))
  }
  
  param = list()
  param$TypeFlag = TypeFlag
  param$S = S
  param$X = X
  param$H = H
  param$time1 = time1
  param$Time2 = Time2
  param$r = r
  param$b = b
  param$sigma = sigma
  if (is.null(title)) 
      title = "Partial Time Single Asset Barrier Option"
  if (is.null(description)) 
      description = as.character(date())
  new("fOPTION", call = match.call(), parameters = param, price = PartialTimeBarrier, 
      title = title, description = description)
}   



