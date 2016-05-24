# stationariry condition. 
# for details, see He and Ter\"{a}svirta (2004) and Nakatani and Ter\"{a}svirta (2007).
   stationarity <- function(A,B){
      G <- A + B
      max(Mod(eigen(G)$values))
   }
