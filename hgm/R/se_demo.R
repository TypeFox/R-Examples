# 2014/03/24 sei
# A demo program; von Mises--Fisher on S^{3-1}
# $OpenXM: OpenXM/src/R/r-packages/hgm/R/se_demo.R,v 1.1 2014/03/25 02:25:26 takayama Exp $

# source("se_hgm.R")

hgm.se.demo1<-function() {

G.exact = function(th){  # exact value by built-in function
  c( sinh(th[1])/th[1], cosh(th[1])/th[1] - sinh(th[1])/th[1]^2 )
}

dG.fun = function(th, G, fn.params=NULL){  # Pfaffian
  dG = array(0, c(1, 2))
  sh = G[1] * th[1]
  ch = G[2] * th[1] + G[1]
  dG[1,1] = G[2] # Pfaffian eq's
  dG[1,2] = sh/th[1] - 2*ch/th[1]^2 + 2*sh/th[1]^3
  dG
}

th0 = 0.5
th1 = 15

G0 = G.exact(th0)
G0

G1 = hgm.se.hgm(th0, G0, th1, dG.fun)  # HGM
G1

G1.exact = G.exact(th1)
G1.exact

return(data.frame(exact=G1.exact,byHGM=G1,start=G0));
}
# EOF
