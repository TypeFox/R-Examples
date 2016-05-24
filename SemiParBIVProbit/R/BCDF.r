BCDF <- function (u1, u2, family, par1){

if(family %in% c(2:5) ) res <- (u1^-par1 + u2^-par1 - 1)^(-1/par1) # Clayton

if(family %in% c(6:9) ){ # Joe 
  bit1 <- (1 - u1)^par1
  bit2 <- (1 - u2)^par1
  res <- 1 - (bit1 + bit2 - bit1*bit2)^(1/par1)
}

if(family %in% c(10:13) ) res <- exp( - ( (-log(u1))^par1 + (-log(u2))^par1 )^(1/par1))  # Gumbel


if(family == 14 ){ # Frank -(1/par1)*log(  1 + ( (exp(-par1*u1) - 1)*(exp(-par1*u2) - 1) )/( exp(-par1) - 1 ) ) 
   bit <- -expm1(-par1) # 1 - exp(-par1)
   res <- -(1/par1)*log( (bit - (1 - exp(-par1*u1))*(1 - exp(-par1*u2)))/bit ) 
}




if(family == 55) res <- u1 * u2 / (1 - par1 * (1 - u1) * (1 - u2)) # AMH

if(family == 56) res <- u1 * u2 * (1 + par1 * (1 - u1) * (1 - u2)) # FGM



res

}
