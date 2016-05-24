# fourth order moment condition 
# for details, see He and Ter\"{a}svirta (2004) and Nakatani and Ter\"{a}svirta (2007).
   fourth <- function(A,B,R){
     AB <- A+B
     AA <- A%x%A
     G <- AB%x%AB + 2*AA*(diag(as.vector(R))^2)
     max(Mod(eigen(G)$values))
   }

