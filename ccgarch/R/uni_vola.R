# returning univariate volatility

# for "R_uni_vola.c"
   uni.vola <- function(a,u){
      usq <- u^2
#      .Call("uni_vola", a, usq, PACKAGE="ccgarch")
      .Call("uni_vola", a, usq)

   }
