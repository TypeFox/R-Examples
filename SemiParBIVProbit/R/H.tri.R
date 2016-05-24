H.tri <- function(respvec, VC, TIn, LgTRI){

dst.1 <- dnorm( (TIn$eta2  - TIn$theta12 * TIn$eta1 )/sqrt(1 - TIn$theta12^2) )  
pst.1 <- pnorm( ( ((TIn$eta3 - TIn$theta13 * TIn$eta1)/sqrt(1 - TIn$theta13^2)) - ((TIn$theta23 - TIn$theta12 * TIn$theta13)/sqrt((1 - TIn$theta12^2) * (1 - TIn$theta13^2))) * ((TIn$eta2  - TIn$theta12 * TIn$eta1)/sqrt(1 - TIn$theta12^2)) )/sqrt(1 - ((TIn$theta23 - TIn$theta12 * TIn$theta13)/sqrt((1 - TIn$theta12^2) * (1 - TIn$theta13^2)))^2))
pst.1 <- mm(pst.1)
 st.1 <- -TIn$theta12/sqrt(1 - TIn$theta12^2)

dst.2 <- dnorm((TIn$eta3 - TIn$theta13 * TIn$eta1)/sqrt(1 - TIn$theta13^2))
pst.2 <- pnorm( ( ((TIn$eta2  - TIn$theta12 * TIn$eta1)/sqrt(1 - TIn$theta12^2)) - ((TIn$theta23 - TIn$theta12 * TIn$theta13)/sqrt((1 - TIn$theta12^2) * (1 - TIn$theta13^2))) * ((TIn$eta3 - TIn$theta13 * TIn$eta1)/sqrt(1 - TIn$theta13^2)) )/sqrt(1 - ((TIn$theta23 - TIn$theta12 * TIn$theta13)/sqrt((1 - TIn$theta12^2) * (1 - TIn$theta13^2)))^2))
pst.2 <- mm(pst.2)
 st.2 <- -TIn$theta13/sqrt(1 - TIn$theta13^2)

dst.3 <- dnorm((TIn$eta1 - TIn$theta12 * TIn$eta2)/sqrt(1 - TIn$theta12^2))
pst.3 <- pnorm( ( ((TIn$eta3 - TIn$theta23 * TIn$eta2)/sqrt(1 - TIn$theta23^2)) - ((TIn$theta13 - TIn$theta12 * TIn$theta23)/sqrt((1 - TIn$theta12^2) * (1 - TIn$theta23^2))) * ((TIn$eta1 - TIn$theta12 * TIn$eta2)/sqrt(1 - TIn$theta12^2)) )/sqrt(1 - ((TIn$theta13 - TIn$theta12 * TIn$theta23)/sqrt((1 - TIn$theta12^2) * (1 - TIn$theta23^2)))^2))
pst.3 <- mm(pst.3)
 st.3 <- -TIn$theta12/sqrt(1 - TIn$theta12^2)

dst.4 <- dnorm((TIn$eta3 - TIn$theta23 * TIn$eta2)/sqrt(1 - TIn$theta23^2))
pst.4 <- pnorm( ( ((TIn$eta1 - TIn$theta12 * TIn$eta2)/sqrt(1 - TIn$theta12^2)) - ((TIn$theta13 - TIn$theta12 * TIn$theta23)/sqrt((1 - TIn$theta12^2) * (1 - TIn$theta23^2))) * ((TIn$eta3 - TIn$theta23 * TIn$eta2)/sqrt(1 - TIn$theta23^2)) )/sqrt(1 - ((TIn$theta13 - TIn$theta12 * TIn$theta23)/sqrt((1 - TIn$theta12^2) * (1 - TIn$theta23^2)))^2))
pst.4 <- mm(pst.4)
st.4  <- -TIn$theta23/sqrt(1 - TIn$theta23^2)

dst.5 <- dnorm((TIn$eta1 - TIn$theta13 * TIn$eta3)/sqrt(1 - TIn$theta13^2))
pst.5 <- pnorm( ( ((TIn$eta2  - TIn$theta23 * TIn$eta3)/sqrt(1 - TIn$theta23^2)) - ((TIn$theta12 - TIn$theta13 * TIn$theta23)/sqrt((1 - TIn$theta13^2) * (1 - TIn$theta23^2))) * ((TIn$eta1 - TIn$theta13 * TIn$eta3)/sqrt(1 - TIn$theta13^2)) )/sqrt(1 - ((TIn$theta12 - TIn$theta13 * TIn$theta23)/sqrt((1 - TIn$theta13^2) * (1 - TIn$theta23^2)))^2))
pst.5 <- mm(pst.5)
st.5  <- -TIn$theta13/sqrt(1 - TIn$theta13^2)

dst.6 <- dnorm((TIn$eta2 - TIn$theta23 * TIn$eta3)/sqrt(1 - TIn$theta23^2))
pst.6 <- pnorm( ( ((TIn$eta1 - TIn$theta13 * TIn$eta3)/sqrt(1 - TIn$theta13^2)) - ((TIn$theta12 - TIn$theta13 * TIn$theta23)/sqrt((1 - TIn$theta13^2) * (1 - TIn$theta23^2))) * ((TIn$eta2 - TIn$theta23 * TIn$eta3)/sqrt(1 - TIn$theta23^2)) )/sqrt(1 - ((TIn$theta12 - TIn$theta13 * TIn$theta23)/sqrt((1 - TIn$theta13^2) * (1 - TIn$theta23^2)))^2))
pst.6 <- mm(pst.6)
st.6  <- -TIn$theta23/sqrt(1 - TIn$theta23^2)

dp.1.11.de1 <- dst.1 * pst.1        * st.1      + dst.2 * pst.2       * st.2
dp.1.10.de1 <- dst.1 * (1 - pst.1)  * st.1      + dst.2 * pst.2       * ( -st.2 )
dp.1.00.de1 <- dst.1 * (1 - pst.1)  * ( -st.1 ) + dst.2 * (1 - pst.2) * ( -st.2 )
dp.1.01.de1 <- dst.1 * pst.1        * ( -st.1 ) + dst.2 * (1 - pst.2) * st.2 

dp.2.11.de2 <- dst.3 * pst.3        * st.3      + dst.4 * pst.4       * st.4
dp.2.10.de2 <- dst.3 * (1 - pst.3)  * st.3      + dst.4 * pst.4       * ( -st.4 )
dp.2.00.de2 <- dst.3 * (1 - pst.3)  * ( -st.3 ) + dst.4 * (1 - pst.4) * ( -st.4 )
dp.2.01.de2 <- dst.3 * pst.3        * ( -st.3 ) + dst.4 * (1 - pst.4) * st.4 

dp.3.11.de3 <- dst.5 * pst.5        * st.5      + dst.6 * pst.6       * st.6
dp.3.10.de3 <- dst.5 * (1 - pst.5)  * st.5      + dst.6 * pst.6       * ( -st.6 )
dp.3.00.de3 <- dst.5 * (1 - pst.5)  * ( -st.5 ) + dst.6 * (1 - pst.6) * ( -st.6 )
dp.3.01.de3 <- dst.5 * pst.5        * ( -st.5 ) + dst.6 * (1 - pst.6) * st.6 

d.1 <- dnorm(TIn$eta1)
d.2 <- dnorm(TIn$eta2)
d.3 <- dnorm(TIn$eta3)


d2l.de1.e1 <- respvec$y1.y2.y3  * ( -1/TIn$p111^2 * (d.1 * LgTRI$p.1.11)^2 + 1/TIn$p111 * ( -TIn$eta1 * d.1 * LgTRI$p.1.11 + d.1 * dp.1.11.de1) ) + 
  respvec$y1.y2.cy3 * ( -1/TIn$p110^2 * (d.1 * LgTRI$p.1.10)^2 + 1/TIn$p110 * ( -TIn$eta1 * d.1 * LgTRI$p.1.10 + d.1 * dp.1.10.de1) ) -
  respvec$cy1.y2.y3  * (  1/TIn$p011^2 * (d.1 * LgTRI$p.1.11)^2 + 1/TIn$p011 * ( -TIn$eta1 * d.1 * LgTRI$p.1.11 + d.1 * dp.1.11.de1) ) -
  respvec$cy1.y2.cy3 * (  1/TIn$p010^2 * (d.1 * LgTRI$p.1.10)^2 + 1/TIn$p010 * ( -TIn$eta1 * d.1 * LgTRI$p.1.10 + d.1 * dp.1.10.de1) ) -
  respvec$cy1.cy2.cy3 * (  1/TIn$p000^2 * (d.1 * LgTRI$p.1.00)^2 + 1/TIn$p000 * ( -TIn$eta1 * d.1 * LgTRI$p.1.00 + d.1 * dp.1.00.de1) ) -
  respvec$cy1.cy2.y3 * (  1/TIn$p001^2 * (d.1 * LgTRI$p.1.01)^2 + 1/TIn$p001 * ( -TIn$eta1 * d.1 * LgTRI$p.1.01 + d.1 * dp.1.01.de1) ) +
  respvec$y1.cy2.cy3 * ( -1/TIn$p100^2 * (d.1 * LgTRI$p.1.00)^2 + 1/TIn$p100 * ( -TIn$eta1 * d.1 * LgTRI$p.1.00 + d.1 * dp.1.00.de1) ) +
  respvec$y1.cy2.y3  * ( -1/TIn$p101^2 * (d.1 * LgTRI$p.1.01)^2 + 1/TIn$p101 * ( -TIn$eta1 * d.1 * LgTRI$p.1.01 + d.1 * dp.1.01.de1) )



d2l.de2.e2 <-  respvec$y1.y2.y3 * ( -1/TIn$p111^2 * (d.2 * LgTRI$p.2.11)^2 + 1/TIn$p111 * ( -TIn$eta2   * d.2 * LgTRI$p.2.11 + d.2 * dp.2.11.de2) ) + 
  respvec$y1.y2.cy3 * ( -1/TIn$p110^2 * (d.2 * LgTRI$p.2.10)^2 + 1/TIn$p110 * ( -TIn$eta2  * d.2 * LgTRI$p.2.10 + d.2 * dp.2.10.de2) ) +
  respvec$cy1.y2.y3 * ( -1/TIn$p011^2 * (d.2 * LgTRI$p.2.01)^2 + 1/TIn$p011 * ( -TIn$eta2  * d.2 * LgTRI$p.2.01 + d.2 * dp.2.01.de2) ) +
  respvec$cy1.y2.cy3 * ( -1/TIn$p010^2 * (d.2 * LgTRI$p.2.00)^2 + 1/TIn$p010 * ( -TIn$eta2  * d.2 * LgTRI$p.2.00 + d.2 * dp.2.00.de2) ) -
  respvec$cy1.cy2.cy3 * (  1/TIn$p000^2 * (d.2 * LgTRI$p.2.00)^2 + 1/TIn$p000 * ( -TIn$eta2  * d.2 * LgTRI$p.2.00 + d.2 * dp.2.00.de2) ) -
  respvec$cy1.cy2.y3 * (  1/TIn$p001^2 * (d.2 * LgTRI$p.2.01)^2 + 1/TIn$p001 * ( -TIn$eta2  * d.2 * LgTRI$p.2.01 + d.2 * dp.2.01.de2) ) -
  respvec$y1.cy2.cy3 * (  1/TIn$p100^2 * (d.2 * LgTRI$p.2.10)^2 + 1/TIn$p100 * ( -TIn$eta2  * d.2 * LgTRI$p.2.10 + d.2 * dp.2.10.de2) ) -
  respvec$y1.cy2.y3 * (  1/TIn$p101^2 * (d.2 * LgTRI$p.2.11)^2 + 1/TIn$p101 * ( -TIn$eta2  * d.2 * LgTRI$p.2.11 + d.2 * dp.2.11.de2) )



d2l.de3.e3 <-  respvec$y1.y2.y3 * ( -1/TIn$p111^2 * (d.3 * LgTRI$p.3.11)^2 + 1/TIn$p111 * ( -TIn$eta3 * d.3 * LgTRI$p.3.11 + d.3 * dp.3.11.de3) ) - 
  respvec$y1.y2.cy3 * (  1/TIn$p110^2 * (d.3 * LgTRI$p.3.11)^2 + 1/TIn$p110 * ( -TIn$eta3 * d.3 * LgTRI$p.3.11 + d.3 * dp.3.11.de3) ) +
  respvec$cy1.y2.y3 * ( -1/TIn$p011^2 * (d.3 * LgTRI$p.3.01)^2 + 1/TIn$p011 * ( -TIn$eta3 * d.3 * LgTRI$p.3.01 + d.3 * dp.3.01.de3) ) -
  respvec$cy1.y2.cy3 * (  1/TIn$p010^2 * (d.3 * LgTRI$p.3.01)^2 + 1/TIn$p010 * ( -TIn$eta3 * d.3 * LgTRI$p.3.01 + d.3 * dp.3.01.de3) ) -
  respvec$cy1.cy2.cy3 * (  1/TIn$p000^2 * (d.3 * LgTRI$p.3.00)^2 + 1/TIn$p000 * ( -TIn$eta3 * d.3 * LgTRI$p.3.00 + d.3 * dp.3.00.de3) ) +
  respvec$cy1.cy2.y3 * ( -1/TIn$p001^2 * (d.3 * LgTRI$p.3.00)^2 + 1/TIn$p001 * ( -TIn$eta3 * d.3 * LgTRI$p.3.00 + d.3 * dp.3.00.de3) ) -
  respvec$y1.cy2.cy3 * (  1/TIn$p100^2 * (d.3 * LgTRI$p.3.10)^2 + 1/TIn$p100 * ( -TIn$eta3 * d.3 * LgTRI$p.3.10 + d.3 * dp.3.10.de3) ) +
  respvec$y1.cy2.y3 * ( -1/TIn$p101^2 * (d.3 * LgTRI$p.3.10)^2 + 1/TIn$p101 * ( -TIn$eta3 * d.3 * LgTRI$p.3.10 + d.3 * dp.3.10.de3) )


mean11 <- TIn$theta23 * TIn$eta2 + ((TIn$theta13 - TIn$theta12 * TIn$theta23) * (TIn$eta1 - TIn$theta12 * TIn$eta2))/(1 - TIn$theta12^2)
mean22 <- TIn$theta23 * TIn$eta3 + ((TIn$theta12 - TIn$theta13 * TIn$theta23) * (TIn$eta1 - TIn$theta13 * TIn$eta3))/(1 - TIn$theta13^2)
mean33 <- TIn$theta13 * TIn$eta3 + ((TIn$theta12 - TIn$theta13 * TIn$theta23) * (TIn$eta2 - TIn$theta23 * TIn$eta3))/(1 - TIn$theta23^2)

deno <- 1 - TIn$theta12^2 - TIn$theta13^2 - TIn$theta23^2 + 2 * TIn$theta12 * TIn$theta13 * TIn$theta23

sd11 <- sqrt( deno / ( 1 - TIn$theta12^2 ) )
sd22 <- sqrt( deno / ( 1 - TIn$theta13^2 ) )
sd33 <- sqrt( deno / ( 1 - TIn$theta23^2 ) )

p1 <- mm( pnorm((TIn$eta3 - mean11)/sd11) )
p2 <- mm( pnorm((TIn$eta2 - mean22)/sd22) )
p3 <- mm( pnorm((TIn$eta1 - mean33)/sd33) )

p1.c <- mm(1 - p1)
p2.c <- mm(1 - p2)
p3.c <- mm(1 - p3)
 
d.1.1  <- dnorm(TIn$eta1, mean = TIn$theta12 * TIn$eta2, sd = sqrt(1 - TIn$theta12^2))
d.1.2  <- dnorm(TIn$eta1, mean = TIn$theta13 * TIn$eta3, sd = sqrt(1 - TIn$theta13^2))
d.1.3  <- dnorm(TIn$eta2, mean = TIn$theta23 * TIn$eta3, sd = sqrt(1 - TIn$theta23^2))

d2l.de1.e2 <-  respvec$y1.y2.y3 * ( -1/TIn$p111^2 * d.2 * LgTRI$p.2.11 * d.1 * LgTRI$p.1.11 + 1/TIn$p111 * d.2 * d.1.1 * p1   ) +
  respvec$y1.y2.cy3 * ( -1/TIn$p110^2 * d.2 * LgTRI$p.2.10 * d.1 * LgTRI$p.1.10 + 1/TIn$p110 * d.2 * d.1.1 * p1.c ) -
  respvec$cy1.y2.y3 * ( -1/TIn$p011^2 * d.2 * LgTRI$p.2.01 * d.1 * LgTRI$p.1.11 + 1/TIn$p011 * d.2 * d.1.1 * p1   ) -
  respvec$cy1.y2.cy3 * ( -1/TIn$p010^2 * d.2 * LgTRI$p.2.00 * d.1 * LgTRI$p.1.10 + 1/TIn$p010 * d.2 * d.1.1 * p1.c ) +
  respvec$cy1.cy2.cy3 * ( -1/TIn$p000^2 * d.2 * LgTRI$p.2.00 * d.1 * LgTRI$p.1.00 + 1/TIn$p000 * d.2 * d.1.1 * p1.c ) +
  respvec$cy1.cy2.y3 * ( -1/TIn$p001^2 * d.2 * LgTRI$p.2.01 * d.1 * LgTRI$p.1.01 + 1/TIn$p001 * d.2 * d.1.1 * p1   ) -
  respvec$y1.cy2.cy3 * ( -1/TIn$p100^2 * d.2 * LgTRI$p.2.10 * d.1 * LgTRI$p.1.00 + 1/TIn$p100 * d.2 * d.1.1 * p1.c ) -
  respvec$y1.cy2.y3 * ( -1/TIn$p101^2 * d.2 * LgTRI$p.2.11 * d.1 * LgTRI$p.1.01 + 1/TIn$p101 * d.2 * d.1.1 * p1   ) 

d2l.de1.e3 <-  respvec$y1.y2.y3 * ( -1/TIn$p111^2 * d.3 * LgTRI$p.3.11 * d.1 * LgTRI$p.1.11 + 1/TIn$p111 * d.3 * d.1.2 * p2   ) -
  respvec$y1.y2.cy3 * ( -1/TIn$p110^2 * d.3 * LgTRI$p.3.11 * d.1 * LgTRI$p.1.10 + 1/TIn$p110 * d.3 * d.1.2 * p2   ) -
  respvec$cy1.y2.y3 * ( -1/TIn$p011^2 * d.3 * LgTRI$p.3.01 * d.1 * LgTRI$p.1.11 + 1/TIn$p011 * d.3 * d.1.2 * p2   ) +
  respvec$cy1.y2.cy3 * ( -1/TIn$p010^2 * d.3 * LgTRI$p.3.01 * d.1 * LgTRI$p.1.10 + 1/TIn$p010 * d.3 * d.1.2 * p2   ) +
  respvec$cy1.cy2.cy3 * ( -1/TIn$p000^2 * d.3 * LgTRI$p.3.00 * d.1 * LgTRI$p.1.00 + 1/TIn$p000 * d.3 * d.1.2 * p2.c ) -
  respvec$cy1.cy2.y3 * ( -1/TIn$p001^2 * d.3 * LgTRI$p.3.00 * d.1 * LgTRI$p.1.01 + 1/TIn$p001 * d.3 * d.1.2 * p2.c ) -
  respvec$y1.cy2.cy3 * ( -1/TIn$p100^2 * d.3 * LgTRI$p.3.10 * d.1 * LgTRI$p.1.00 + 1/TIn$p100 * d.3 * d.1.2 * p2.c ) +
  respvec$y1.cy2.y3 * ( -1/TIn$p101^2 * d.3 * LgTRI$p.3.10 * d.1 * LgTRI$p.1.01 + 1/TIn$p101 * d.3 * d.1.2 * p2.c ) 

d2l.de2.e3 <-  respvec$y1.y2.y3 * ( -1/TIn$p111^2 * d.3 * LgTRI$p.3.11 * d.2 * LgTRI$p.2.11 + 1/TIn$p111 * d.3 * d.1.3 * p3   ) -
  respvec$y1.y2.cy3 * ( -1/TIn$p110^2 * d.3 * LgTRI$p.3.11 * d.2 * LgTRI$p.2.10 + 1/TIn$p110 * d.3 * d.1.3 * p3   ) +
  respvec$cy1.y2.y3 * ( -1/TIn$p011^2 * d.3 * LgTRI$p.3.01 * d.2 * LgTRI$p.2.01 + 1/TIn$p011 * d.3 * d.1.3 * p3.c ) -
  respvec$cy1.y2.cy3 * ( -1/TIn$p010^2 * d.3 * LgTRI$p.3.01 * d.2 * LgTRI$p.2.00 + 1/TIn$p010 * d.3 * d.1.3 * p3.c ) +
  respvec$cy1.cy2.cy3 * ( -1/TIn$p000^2 * d.3 * LgTRI$p.3.00 * d.2 * LgTRI$p.2.00 + 1/TIn$p000 * d.3 * d.1.3 * p3.c ) -
  respvec$cy1.cy2.y3 * ( -1/TIn$p001^2 * d.3 * LgTRI$p.3.00 * d.2 * LgTRI$p.2.01 + 1/TIn$p001 * d.3 * d.1.3 * p3.c ) +
  respvec$y1.cy2.cy3 * ( -1/TIn$p100^2 * d.3 * LgTRI$p.3.10 * d.2 * LgTRI$p.2.10 + 1/TIn$p100 * d.3 * d.1.3 * p3   ) -
  respvec$y1.cy2.y3 * ( -1/TIn$p101^2 * d.3 * LgTRI$p.3.10 * d.2 * LgTRI$p.2.11 + 1/TIn$p101 * d.3 * d.1.3 * p3   ) 

d12 <- dnorm( (TIn$eta3 - LgTRI$mean.12)/LgTRI$sd.12 )
d13 <- dnorm( (TIn$eta2  - LgTRI$mean.13)/LgTRI$sd.13 )
d23 <- dnorm( (TIn$eta1 - LgTRI$mean.23)/LgTRI$sd.23 )


d2l.de1.theta12 <- respvec$y1.y2.y3 * ( -1/TIn$p111^2 * d.1 * LgTRI$p.1.11 * LgTRI$d11.12 * LgTRI$p12.g   + 1/TIn$p111 * ( LgTRI$d11.12 * (TIn$theta12*TIn$eta2  - TIn$eta1)/(1 - TIn$theta12^2) * LgTRI$p12.g   - LgTRI$d11.12 * d12/LgTRI$sd.12 * ((TIn$theta13 - TIn$theta12 * TIn$theta23)/(1 - TIn$theta12^2)) )) +
  respvec$y1.y2.cy3 * ( -1/TIn$p110^2 * d.1 * LgTRI$p.1.10 * LgTRI$d11.12 * LgTRI$p12.g.c + 1/TIn$p110 * ( LgTRI$d11.12 * (TIn$theta12*TIn$eta2  - TIn$eta1)/(1 - TIn$theta12^2) * LgTRI$p12.g.c + LgTRI$d11.12 * d12/LgTRI$sd.12 * ((TIn$theta13 - TIn$theta12 * TIn$theta23)/(1 - TIn$theta12^2)) )) -
  respvec$cy1.y2.y3 * (  1/TIn$p011^2 * d.1 * LgTRI$p.1.11 * LgTRI$d11.12 * LgTRI$p12.g   + 1/TIn$p011 * ( LgTRI$d11.12 * (TIn$theta12*TIn$eta2  - TIn$eta1)/(1 - TIn$theta12^2) * LgTRI$p12.g   - LgTRI$d11.12 * d12/LgTRI$sd.12 * ((TIn$theta13 - TIn$theta12 * TIn$theta23)/(1 - TIn$theta12^2)) )) -
  respvec$cy1.y2.cy3 * (  1/TIn$p010^2 * d.1 * LgTRI$p.1.10 * LgTRI$d11.12 * LgTRI$p12.g.c + 1/TIn$p010 * ( LgTRI$d11.12 * (TIn$theta12*TIn$eta2  - TIn$eta1)/(1 - TIn$theta12^2) * LgTRI$p12.g.c + LgTRI$d11.12 * d12/LgTRI$sd.12 * ((TIn$theta13 - TIn$theta12 * TIn$theta23)/(1 - TIn$theta12^2)) )) +
  respvec$cy1.cy2.cy3 * (  1/TIn$p000^2 * d.1 * LgTRI$p.1.00 * LgTRI$d11.12 * LgTRI$p12.g.c + 1/TIn$p000 * ( LgTRI$d11.12 * (TIn$theta12*TIn$eta2  - TIn$eta1)/(1 - TIn$theta12^2) * LgTRI$p12.g.c + LgTRI$d11.12 * d12/LgTRI$sd.12 * ((TIn$theta13 - TIn$theta12 * TIn$theta23)/(1 - TIn$theta12^2)) )) +
  respvec$cy1.cy2.y3 * (  1/TIn$p001^2 * d.1 * LgTRI$p.1.01 * LgTRI$d11.12 * LgTRI$p12.g   + 1/TIn$p001 * ( LgTRI$d11.12 * (TIn$theta12*TIn$eta2  - TIn$eta1)/(1 - TIn$theta12^2) * LgTRI$p12.g   - LgTRI$d11.12 * d12/LgTRI$sd.12 * ((TIn$theta13 - TIn$theta12 * TIn$theta23)/(1 - TIn$theta12^2)) )) -
  respvec$y1.cy2.cy3 * ( -1/TIn$p100^2 * d.1 * LgTRI$p.1.00 * LgTRI$d11.12 * LgTRI$p12.g.c + 1/TIn$p100 * ( LgTRI$d11.12 * (TIn$theta12*TIn$eta2  - TIn$eta1)/(1 - TIn$theta12^2) * LgTRI$p12.g.c + LgTRI$d11.12 * d12/LgTRI$sd.12 * ((TIn$theta13 - TIn$theta12 * TIn$theta23)/(1 - TIn$theta12^2)) )) -
  respvec$y1.cy2.y3 * ( -1/TIn$p101^2 * d.1 * LgTRI$p.1.01 * LgTRI$d11.12 * LgTRI$p12.g   + 1/TIn$p101 * ( LgTRI$d11.12 * (TIn$theta12*TIn$eta2  - TIn$eta1)/(1 - TIn$theta12^2) * LgTRI$p12.g   - LgTRI$d11.12 * d12/LgTRI$sd.12 * ((TIn$theta13 - TIn$theta12 * TIn$theta23)/(1 - TIn$theta12^2)) )) 

d2l.de1.theta13 <- respvec$y1.y2.y3 * ( -1/TIn$p111^2 * d.1 * LgTRI$p.1.11 * LgTRI$d11.13 * LgTRI$p13.g   + 1/TIn$p111 * ( LgTRI$d11.13 * (TIn$theta13 * TIn$eta3 - TIn$eta1)/(1 - TIn$theta13^2) * LgTRI$p13.g   - LgTRI$d11.13 * d13/LgTRI$sd.13 * ((TIn$theta12 - TIn$theta13 * TIn$theta23)/(1 - TIn$theta13^2)) )) -
  respvec$y1.y2.cy3   * ( -1/TIn$p110^2 * d.1 * LgTRI$p.1.10 * LgTRI$d11.13 * LgTRI$p13.g   + 1/TIn$p110 * ( LgTRI$d11.13 * (TIn$theta13 * TIn$eta3 - TIn$eta1)/(1 - TIn$theta13^2) * LgTRI$p13.g   - LgTRI$d11.13 * d13/LgTRI$sd.13 * ((TIn$theta12 - TIn$theta13 * TIn$theta23)/(1 - TIn$theta13^2)) )) -
  respvec$cy1.y2.y3   * (  1/TIn$p011^2 * d.1 * LgTRI$p.1.11 * LgTRI$d11.13 * LgTRI$p13.g   + 1/TIn$p011 * ( LgTRI$d11.13 * (TIn$theta13 * TIn$eta3 - TIn$eta1)/(1 - TIn$theta13^2) * LgTRI$p13.g   - LgTRI$d11.13 * d13/LgTRI$sd.13 * ((TIn$theta12 - TIn$theta13 * TIn$theta23)/(1 - TIn$theta13^2)) )) +
  respvec$cy1.y2.cy3  * (  1/TIn$p010^2 * d.1 * LgTRI$p.1.10 * LgTRI$d11.13 * LgTRI$p13.g   + 1/TIn$p010 * ( LgTRI$d11.13 * (TIn$theta13 * TIn$eta3 - TIn$eta1)/(1 - TIn$theta13^2) * LgTRI$p13.g   - LgTRI$d11.13 * d13/LgTRI$sd.13 * ((TIn$theta12 - TIn$theta13 * TIn$theta23)/(1 - TIn$theta13^2)) )) +
  respvec$cy1.cy2.cy3 * (  1/TIn$p000^2 * d.1 * LgTRI$p.1.00 * LgTRI$d11.13 * LgTRI$p13.g.c + 1/TIn$p000 * ( LgTRI$d11.13 * (TIn$theta13 * TIn$eta3 - TIn$eta1)/(1 - TIn$theta13^2) * LgTRI$p13.g.c + LgTRI$d11.13 * d13/LgTRI$sd.13 * ((TIn$theta12 - TIn$theta13 * TIn$theta23)/(1 - TIn$theta13^2)) )) -
  respvec$cy1.cy2.y3  * (  1/TIn$p001^2 * d.1 * LgTRI$p.1.01 * LgTRI$d11.13 * LgTRI$p13.g.c + 1/TIn$p001 * ( LgTRI$d11.13 * (TIn$theta13 * TIn$eta3 - TIn$eta1)/(1 - TIn$theta13^2) * LgTRI$p13.g.c + LgTRI$d11.13 * d13/LgTRI$sd.13 * ((TIn$theta12 - TIn$theta13 * TIn$theta23)/(1 - TIn$theta13^2)) )) -
  respvec$y1.cy2.cy3  * ( -1/TIn$p100^2 * d.1 * LgTRI$p.1.00 * LgTRI$d11.13 * LgTRI$p13.g.c + 1/TIn$p100 * ( LgTRI$d11.13 * (TIn$theta13 * TIn$eta3 - TIn$eta1)/(1 - TIn$theta13^2) * LgTRI$p13.g.c + LgTRI$d11.13 * d13/LgTRI$sd.13 * ((TIn$theta12 - TIn$theta13 * TIn$theta23)/(1 - TIn$theta13^2)) )) +
  respvec$y1.cy2.y3   * ( -1/TIn$p101^2 * d.1 * LgTRI$p.1.01 * LgTRI$d11.13 * LgTRI$p13.g.c + 1/TIn$p101 * ( LgTRI$d11.13 * (TIn$theta13 * TIn$eta3 - TIn$eta1)/(1 - TIn$theta13^2) * LgTRI$p13.g.c + LgTRI$d11.13 * d13/LgTRI$sd.13 * ((TIn$theta12 - TIn$theta13 * TIn$theta23)/(1 - TIn$theta13^2)) )) 

d2l.de1.theta23 <- respvec$y1.y2.y3 * ( -1/TIn$p111^2 * d.1 * LgTRI$p.1.11 * LgTRI$d11.23 * LgTRI$p23.g   + 1/TIn$p111 * LgTRI$d11.23 * d23/LgTRI$sd.23) -
  respvec$y1.y2.cy3  * ( -1/TIn$p110^2 * d.1 * LgTRI$p.1.10 * LgTRI$d11.23 * LgTRI$p23.g   + 1/TIn$p110 * LgTRI$d11.23 * d23/LgTRI$sd.23 ) -
  respvec$cy1.y2.y3 * ( -1/TIn$p011^2 * d.1 * LgTRI$p.1.11 * LgTRI$d11.23 * LgTRI$p23.g.c + 1/TIn$p011 * LgTRI$d11.23 * d23/LgTRI$sd.23 ) +
  respvec$cy1.y2.cy3 * ( -1/TIn$p010^2 * d.1 * LgTRI$p.1.10 * LgTRI$d11.23 * LgTRI$p23.g.c + 1/TIn$p010 * LgTRI$d11.23 * d23/LgTRI$sd.23 ) -
  respvec$cy1.cy2.cy3 * ( -1/TIn$p000^2 * d.1 * LgTRI$p.1.00 * LgTRI$d11.23 * LgTRI$p23.g.c + 1/TIn$p000 * LgTRI$d11.23 * d23/LgTRI$sd.23 ) +
  respvec$cy1.cy2.y3 * ( -1/TIn$p001^2 * d.1 * LgTRI$p.1.01 * LgTRI$d11.23 * LgTRI$p23.g.c + 1/TIn$p001 * LgTRI$d11.23 * d23/LgTRI$sd.23 ) +
  respvec$y1.cy2.cy3 * ( -1/TIn$p100^2 * d.1 * LgTRI$p.1.00 * LgTRI$d11.23 * LgTRI$p23.g   + 1/TIn$p100 * LgTRI$d11.23 * d23/LgTRI$sd.23 ) -
  respvec$y1.cy2.y3 * ( -1/TIn$p101^2 * d.1 * LgTRI$p.1.01 * LgTRI$d11.23 * LgTRI$p23.g   + 1/TIn$p101 * LgTRI$d11.23 * d23/LgTRI$sd.23 ) 

d2l.de2.theta12 <- respvec$y1.y2.y3 * ( -1/TIn$p111^2 * d.2 * LgTRI$p.2.11 * LgTRI$d11.12 * LgTRI$p12.g   + 1/TIn$p111 * (LgTRI$d11.12 * (TIn$theta12*TIn$eta1 - TIn$eta2)/(1 - TIn$theta12^2) * LgTRI$p12.g   - LgTRI$d11.12 * (d12/LgTRI$sd.12) * (TIn$theta23 - TIn$theta12 * TIn$theta13)/(1 - TIn$theta12^2) ) ) +
  respvec$y1.y2.cy3 * ( -1/TIn$p110^2 * d.2 * LgTRI$p.2.10 * LgTRI$d11.12 * LgTRI$p12.g.c + 1/TIn$p110 * (LgTRI$d11.12 * (TIn$theta12 * TIn$eta1 - TIn$eta2)/(1 - TIn$theta12^2) * LgTRI$p12.g.c + LgTRI$d11.12 * (d12/LgTRI$sd.12) * (TIn$theta23 - TIn$theta12 * TIn$theta13)/(1 - TIn$theta12^2) ) ) -
  respvec$cy1.y2.y3 * ( -1/TIn$p011^2 * d.2 * LgTRI$p.2.01 * LgTRI$d11.12 * LgTRI$p12.g   + 1/TIn$p011 * (LgTRI$d11.12 * (TIn$theta12 * TIn$eta1 - TIn$eta2)/(1 - TIn$theta12^2) * LgTRI$p12.g   - LgTRI$d11.12 * (d12/LgTRI$sd.12) * (TIn$theta23 - TIn$theta12 * TIn$theta13)/(1 - TIn$theta12^2) ) ) -
  respvec$cy1.y2.cy3 * ( -1/TIn$p010^2 * d.2 * LgTRI$p.2.00 * LgTRI$d11.12 * LgTRI$p12.g.c + 1/TIn$p010 * (LgTRI$d11.12 * (TIn$theta12 * TIn$eta1 - TIn$eta2)/(1 - TIn$theta12^2) * LgTRI$p12.g.c + LgTRI$d11.12 * (d12/LgTRI$sd.12) * (TIn$theta23 - TIn$theta12 * TIn$theta13)/(1 - TIn$theta12^2) ) ) +
  respvec$cy1.cy2.cy3 * (  1/TIn$p000^2 * d.2 * LgTRI$p.2.00 * LgTRI$d11.12 * LgTRI$p12.g.c + 1/TIn$p000 * (LgTRI$d11.12 * (TIn$theta12 * TIn$eta1 - TIn$eta2)/(1 - TIn$theta12^2) * LgTRI$p12.g.c + LgTRI$d11.12 * (d12/LgTRI$sd.12) * (TIn$theta23 - TIn$theta12 * TIn$theta13)/(1 - TIn$theta12^2) ) ) +
  respvec$cy1.cy2.y3 * (  1/TIn$p001^2 * d.2 * LgTRI$p.2.01 * LgTRI$d11.12 * LgTRI$p12.g   + 1/TIn$p001 * (LgTRI$d11.12 * (TIn$theta12 * TIn$eta1 - TIn$eta2)/(1 - TIn$theta12^2) * LgTRI$p12.g   - LgTRI$d11.12 * (d12/LgTRI$sd.12) * (TIn$theta23 - TIn$theta12 * TIn$theta13)/(1 - TIn$theta12^2) ) ) -
  respvec$y1.cy2.cy3 * (  1/TIn$p100^2 * d.2 * LgTRI$p.2.10 * LgTRI$d11.12 * LgTRI$p12.g.c + 1/TIn$p100 * (LgTRI$d11.12 * (TIn$theta12 * TIn$eta1 - TIn$eta2)/(1 - TIn$theta12^2) * LgTRI$p12.g.c + LgTRI$d11.12 * (d12/LgTRI$sd.12) * (TIn$theta23 - TIn$theta12 * TIn$theta13)/(1 - TIn$theta12^2) ) ) -
  respvec$y1.cy2.y3 * (  1/TIn$p101^2 * d.2 * LgTRI$p.2.11 * LgTRI$d11.12 * LgTRI$p12.g   + 1/TIn$p101 * (LgTRI$d11.12 * (TIn$theta12 * TIn$eta1 - TIn$eta2)/(1 - TIn$theta12^2) * LgTRI$p12.g   - LgTRI$d11.12 * (d12/LgTRI$sd.12) * (TIn$theta23 - TIn$theta12 * TIn$theta13)/(1 - TIn$theta12^2) ) )

d2l.de2.theta13 <- respvec$y1.y2.y3 * ( -1/TIn$p111^2 * d.2 * LgTRI$p.2.11 * LgTRI$d11.13 * LgTRI$p13.g   + 1/TIn$p111 * LgTRI$d11.13 * d13/LgTRI$sd.13 )-
  respvec$y1.y2.cy3 * ( -1/TIn$p110^2 * d.2 * LgTRI$p.2.10 * LgTRI$d11.13 * LgTRI$p13.g   + 1/TIn$p110 * LgTRI$d11.13 * d13/LgTRI$sd.13 )-
  respvec$cy1.y2.y3 * ( -1/TIn$p011^2 * d.2 * LgTRI$p.2.01 * LgTRI$d11.13 * LgTRI$p13.g   + 1/TIn$p011 * LgTRI$d11.13 * d13/LgTRI$sd.13 ) +
  respvec$cy1.y2.cy3 * ( -1/TIn$p010^2 * d.2 * LgTRI$p.2.00 * LgTRI$d11.13 * LgTRI$p13.g   + 1/TIn$p010 * LgTRI$d11.13 * d13/LgTRI$sd.13 ) -
  respvec$cy1.cy2.cy3 * ( -1/TIn$p000^2 * d.2 * LgTRI$p.2.00 * LgTRI$d11.13 * LgTRI$p13.g.c + 1/TIn$p000 * LgTRI$d11.13 * d13/LgTRI$sd.13 ) +
  respvec$cy1.cy2.y3 * ( -1/TIn$p001^2 * d.2 * LgTRI$p.2.01 * LgTRI$d11.13 * LgTRI$p13.g.c + 1/TIn$p001 * LgTRI$d11.13 * d13/LgTRI$sd.13 ) +
  respvec$y1.cy2.cy3 * ( -1/TIn$p100^2 * d.2 * LgTRI$p.2.10 * LgTRI$d11.13 * LgTRI$p13.g.c + 1/TIn$p100 * LgTRI$d11.13 * d13/LgTRI$sd.13 ) -
  respvec$y1.cy2.y3 * ( -1/TIn$p101^2 * d.2 * LgTRI$p.2.11 * LgTRI$d11.13 * LgTRI$p13.g.c + 1/TIn$p101 * LgTRI$d11.13 * d13/LgTRI$sd.13 ) 

d2l.de2.theta23 <- respvec$y1.y2.y3 * ( -1/TIn$p111^2 * d.2 * LgTRI$p.2.11 * LgTRI$d11.23 * LgTRI$p23.g   + 1/TIn$p111 * (LgTRI$d11.23 * (TIn$theta23 * TIn$eta3 - TIn$eta2)/(1 - TIn$theta23^2) * LgTRI$p23.g   - LgTRI$d11.23 * d23/LgTRI$sd.23 * (TIn$theta12 - TIn$theta13 * TIn$theta23)/(1 - TIn$theta23^2) ) ) -
  respvec$y1.y2.cy3 * ( -1/TIn$p110^2 * d.2 * LgTRI$p.2.10 * LgTRI$d11.23 * LgTRI$p23.g   + 1/TIn$p110 * (LgTRI$d11.23 * (TIn$theta23 * TIn$eta3 - TIn$eta2)/(1 - TIn$theta23^2) * LgTRI$p23.g   - LgTRI$d11.23 * d23/LgTRI$sd.23 * (TIn$theta12 - TIn$theta13 * TIn$theta23)/(1 - TIn$theta23^2) ) ) +
  respvec$cy1.y2.y3 * ( -1/TIn$p011^2 * d.2 * LgTRI$p.2.01 * LgTRI$d11.23 * LgTRI$p23.g.c + 1/TIn$p011 * (LgTRI$d11.23 * (TIn$theta23 * TIn$eta3 - TIn$eta2)/(1 - TIn$theta23^2) * LgTRI$p23.g.c + LgTRI$d11.23 * d23/LgTRI$sd.23 * (TIn$theta12 - TIn$theta13 * TIn$theta23)/(1 - TIn$theta23^2) ) ) -
  respvec$cy1.y2.cy3 * ( -1/TIn$p010^2 * d.2 * LgTRI$p.2.00 * LgTRI$d11.23 * LgTRI$p23.g.c + 1/TIn$p010 * (LgTRI$d11.23 * (TIn$theta23 * TIn$eta3 - TIn$eta2)/(1 - TIn$theta23^2) * LgTRI$p23.g.c + LgTRI$d11.23 * d23/LgTRI$sd.23 * (TIn$theta12 - TIn$theta13 * TIn$theta23)/(1 - TIn$theta23^2) ) ) +
  respvec$cy1.cy2.cy3 * (  1/TIn$p000^2 * d.2 * LgTRI$p.2.00 * LgTRI$d11.23 * LgTRI$p23.g.c + 1/TIn$p000 * (LgTRI$d11.23 * (TIn$theta23 * TIn$eta3 - TIn$eta2)/(1 - TIn$theta23^2) * LgTRI$p23.g.c + LgTRI$d11.23 * d23/LgTRI$sd.23 * (TIn$theta12 - TIn$theta13 * TIn$theta23)/(1 - TIn$theta23^2) ) ) -
  respvec$cy1.cy2.y3 * (  1/TIn$p001^2 * d.2 * LgTRI$p.2.01 * LgTRI$d11.23 * LgTRI$p23.g.c + 1/TIn$p001 * (LgTRI$d11.23 * (TIn$theta23 * TIn$eta3 - TIn$eta2)/(1 - TIn$theta23^2) * LgTRI$p23.g.c + LgTRI$d11.23 * d23/LgTRI$sd.23 * (TIn$theta12 - TIn$theta13 * TIn$theta23)/(1 - TIn$theta23^2) ) ) +
  respvec$y1.cy2.cy3 * (  1/TIn$p100^2 * d.2 * LgTRI$p.2.10 * LgTRI$d11.23 * LgTRI$p23.g   + 1/TIn$p100 * (LgTRI$d11.23 * (TIn$theta23 * TIn$eta3 - TIn$eta2)/(1 - TIn$theta23^2) * LgTRI$p23.g   - LgTRI$d11.23 * d23/LgTRI$sd.23 * (TIn$theta12 - TIn$theta13 * TIn$theta23)/(1 - TIn$theta23^2) ) ) -
  respvec$y1.cy2.y3 * (  1/TIn$p101^2 * d.2 * LgTRI$p.2.11 * LgTRI$d11.23 * LgTRI$p23.g   + 1/TIn$p101 * (LgTRI$d11.23 * (TIn$theta23 * TIn$eta3 - TIn$eta2)/(1 - TIn$theta23^2) * LgTRI$p23.g   - LgTRI$d11.23 * d23/LgTRI$sd.23 * (TIn$theta12 - TIn$theta13 * TIn$theta23)/(1 - TIn$theta23^2) ) ) 

d2l.de3.theta12 <- respvec$y1.y2.y3 * ( -1/TIn$p111^2 * d.3 * LgTRI$p.3.11 * LgTRI$d11.12 * LgTRI$p12.g   + 1/TIn$p111 * LgTRI$d11.12 * d12/LgTRI$sd.12 ) -
  respvec$y1.y2.cy3 * ( -1/TIn$p110^2 * d.3 * LgTRI$p.3.11 * LgTRI$d11.12 * LgTRI$p12.g.c + 1/TIn$p110 * LgTRI$d11.12 * d12/LgTRI$sd.12 ) -
  respvec$cy1.y2.y3 * ( -1/TIn$p011^2 * d.3 * LgTRI$p.3.01 * LgTRI$d11.12 * LgTRI$p12.g   + 1/TIn$p011 * LgTRI$d11.12 * d12/LgTRI$sd.12 ) +
  respvec$cy1.y2.cy3 * ( -1/TIn$p010^2 * d.3 * LgTRI$p.3.01 * LgTRI$d11.12 * LgTRI$p12.g.c + 1/TIn$p010 * LgTRI$d11.12 * d12/LgTRI$sd.12 ) -
  respvec$cy1.cy2.cy3 * ( -1/TIn$p000^2 * d.3 * LgTRI$p.3.00 * LgTRI$d11.12 * LgTRI$p12.g.c + 1/TIn$p000 * LgTRI$d11.12 * d12/LgTRI$sd.12 ) +
  respvec$cy1.cy2.y3 * ( -1/TIn$p001^2 * d.3 * LgTRI$p.3.00 * LgTRI$d11.12 * LgTRI$p12.g   + 1/TIn$p001 * LgTRI$d11.12 * d12/LgTRI$sd.12 ) +
  respvec$y1.cy2.cy3 * ( -1/TIn$p100^2 * d.3 * LgTRI$p.3.10 * LgTRI$d11.12 * LgTRI$p12.g.c + 1/TIn$p100 * LgTRI$d11.12 * d12/LgTRI$sd.12 ) -
  respvec$y1.cy2.y3 * ( -1/TIn$p101^2 * d.3 * LgTRI$p.3.10 * LgTRI$d11.12 * LgTRI$p12.g   + 1/TIn$p101 * LgTRI$d11.12 * d12/LgTRI$sd.12 ) 

d2l.de3.theta13 <- respvec$y1.y2.y3 * ( -1/TIn$p111^2 * d.3 * LgTRI$p.3.11 * LgTRI$d11.13 * LgTRI$p13.g   + 1/TIn$p111 * (LgTRI$d11.13 * ((TIn$theta13*TIn$eta1 - TIn$eta3)/(1 - TIn$theta13^2)) * LgTRI$p13.g   - LgTRI$d11.13 * (d13/LgTRI$sd.13) * (TIn$theta23 - TIn$theta12 * TIn$theta13)/(1 - TIn$theta13^2) ) ) -
  respvec$y1.y2.cy3 * (  1/TIn$p110^2 * d.3 * LgTRI$p.3.11 * LgTRI$d11.13 * LgTRI$p13.g   + 1/TIn$p110 * (LgTRI$d11.13 * ((TIn$theta13 * TIn$eta1 - TIn$eta3)/(1 - TIn$theta13^2)) * LgTRI$p13.g   - LgTRI$d11.13 * (d13/LgTRI$sd.13) * (TIn$theta23 - TIn$theta12 * TIn$theta13)/(1 - TIn$theta13^2) ) ) -
  respvec$cy1.y2.y3 * ( -1/TIn$p011^2 * d.3 * LgTRI$p.3.01 * LgTRI$d11.13 * LgTRI$p13.g   + 1/TIn$p011 * (LgTRI$d11.13 * ((TIn$theta13 * TIn$eta1 - TIn$eta3)/(1 - TIn$theta13^2)) * LgTRI$p13.g   - LgTRI$d11.13 * (d13/LgTRI$sd.13) * (TIn$theta23 - TIn$theta12 * TIn$theta13)/(1 - TIn$theta13^2) ) ) +
  respvec$cy1.y2.cy3 * (  1/TIn$p010^2 * d.3 * LgTRI$p.3.01 * LgTRI$d11.13 * LgTRI$p13.g   + 1/TIn$p010 * (LgTRI$d11.13 * ((TIn$theta13 * TIn$eta1 - TIn$eta3)/(1 - TIn$theta13^2)) * LgTRI$p13.g   - LgTRI$d11.13 * (d13/LgTRI$sd.13) * (TIn$theta23 - TIn$theta12 * TIn$theta13)/(1 - TIn$theta13^2) ) ) +
  respvec$cy1.cy2.cy3 * (  1/TIn$p000^2 * d.3 * LgTRI$p.3.00 * LgTRI$d11.13 * LgTRI$p13.g.c + 1/TIn$p000 * (LgTRI$d11.13 * ((TIn$theta13 * TIn$eta1 - TIn$eta3)/(1 - TIn$theta13^2)) * LgTRI$p13.g.c + LgTRI$d11.13 * (d13/LgTRI$sd.13) * (TIn$theta23 - TIn$theta12 * TIn$theta13)/(1 - TIn$theta13^2) ) ) -
  respvec$cy1.cy2.y3 * ( -1/TIn$p001^2 * d.3 * LgTRI$p.3.00 * LgTRI$d11.13 * LgTRI$p13.g.c + 1/TIn$p001 * (LgTRI$d11.13 * ((TIn$theta13 * TIn$eta1 - TIn$eta3)/(1 - TIn$theta13^2)) * LgTRI$p13.g.c + LgTRI$d11.13 * (d13/LgTRI$sd.13) * (TIn$theta23 - TIn$theta12 * TIn$theta13)/(1 - TIn$theta13^2) ) ) -
  respvec$y1.cy2.cy3 * (  1/TIn$p100^2 * d.3 * LgTRI$p.3.10 * LgTRI$d11.13 * LgTRI$p13.g.c + 1/TIn$p100 * (LgTRI$d11.13 * ((TIn$theta13 * TIn$eta1 - TIn$eta3)/(1 - TIn$theta13^2)) * LgTRI$p13.g.c + LgTRI$d11.13 * (d13/LgTRI$sd.13) * (TIn$theta23 - TIn$theta12 * TIn$theta13)/(1 - TIn$theta13^2) ) ) +
  respvec$y1.cy2.y3 * ( -1/TIn$p101^2 * d.3 * LgTRI$p.3.10 * LgTRI$d11.13 * LgTRI$p13.g.c + 1/TIn$p101 * (LgTRI$d11.13 * ((TIn$theta13 * TIn$eta1 - TIn$eta3)/(1 - TIn$theta13^2)) * LgTRI$p13.g.c + LgTRI$d11.13 * (d13/LgTRI$sd.13) * (TIn$theta23 - TIn$theta12 * TIn$theta13)/(1 - TIn$theta13^2) ) ) 

d2l.de3.theta23 <- respvec$y1.y2.y3 * ( -1/TIn$p111^2 * d.3 * LgTRI$p.3.11 * LgTRI$d11.23 * LgTRI$p23.g   + 1/TIn$p111 * (LgTRI$d11.23 * ((TIn$theta23*TIn$eta2  - TIn$eta3)/(1 - TIn$theta23^2)) * LgTRI$p23.g   - LgTRI$d11.23 * d23/LgTRI$sd.23 * (TIn$theta13 - TIn$theta12 * TIn$theta23)/(1 - TIn$theta23^2) ) ) -
  respvec$y1.y2.cy3 * (  1/TIn$p110^2 * d.3 * LgTRI$p.3.11 * LgTRI$d11.23 * LgTRI$p23.g   + 1/TIn$p110 * (LgTRI$d11.23 * ((TIn$theta23*TIn$eta2  - TIn$eta3)/(1 - TIn$theta23^2)) * LgTRI$p23.g   - LgTRI$d11.23 * d23/LgTRI$sd.23 * (TIn$theta13 - TIn$theta12*TIn$theta23)/(1 - TIn$theta23^2) ) ) +
  respvec$cy1.y2.y3 * ( -1/TIn$p011^2 * d.3 * LgTRI$p.3.01 * LgTRI$d11.23 * LgTRI$p23.g.c + 1/TIn$p011 * (LgTRI$d11.23 * ((TIn$theta23*TIn$eta2  - TIn$eta3)/(1 - TIn$theta23^2)) * LgTRI$p23.g.c + LgTRI$d11.23 * d23/LgTRI$sd.23 * (TIn$theta13 - TIn$theta12*TIn$theta23)/(1 - TIn$theta23^2) ) ) -
  respvec$cy1.y2.cy3 * (  1/TIn$p010^2 * d.3 * LgTRI$p.3.01 * LgTRI$d11.23 * LgTRI$p23.g.c + 1/TIn$p010 * (LgTRI$d11.23 * ((TIn$theta23*TIn$eta2  - TIn$eta3)/(1 - TIn$theta23^2)) * LgTRI$p23.g.c + LgTRI$d11.23 * d23/LgTRI$sd.23 * (TIn$theta13 - TIn$theta12 * TIn$theta23)/(1 - TIn$theta23^2) ) ) +
  respvec$cy1.cy2.cy3 * (  1/TIn$p000^2 * d.3 * LgTRI$p.3.00 * LgTRI$d11.23 * LgTRI$p23.g.c + 1/TIn$p000 * (LgTRI$d11.23 * ((TIn$theta23*TIn$eta2  - TIn$eta3)/(1 - TIn$theta23^2)) * LgTRI$p23.g.c + LgTRI$d11.23 * d23/LgTRI$sd.23 * (TIn$theta13 - TIn$theta12 * TIn$theta23)/(1 - TIn$theta23^2) ) ) -
  respvec$cy1.cy2.y3 * ( -1/TIn$p001^2 * d.3 * LgTRI$p.3.00 * LgTRI$d11.23 * LgTRI$p23.g.c + 1/TIn$p001 * (LgTRI$d11.23 * ((TIn$theta23*TIn$eta2  - TIn$eta3)/(1 - TIn$theta23^2)) * LgTRI$p23.g.c + LgTRI$d11.23 * d23/LgTRI$sd.23 * (TIn$theta13 - TIn$theta12 * TIn$theta23)/(1 - TIn$theta23^2) ) ) +
  respvec$y1.cy2.cy3 * (  1/TIn$p100^2 * d.3 * LgTRI$p.3.10 * LgTRI$d11.23 * LgTRI$p23.g   + 1/TIn$p100 * (LgTRI$d11.23 * ((TIn$theta23*TIn$eta2  - TIn$eta3)/(1 - TIn$theta23^2)) * LgTRI$p23.g   - LgTRI$d11.23 * d23/LgTRI$sd.23 * (TIn$theta13 - TIn$theta12 * TIn$theta23)/(1 - TIn$theta23^2) ) ) -
  respvec$y1.cy2.y3 * ( -1/TIn$p101^2 * d.3 * LgTRI$p.3.10 * LgTRI$d11.23 * LgTRI$p23.g   + 1/TIn$p101 * (LgTRI$d11.23 * ((TIn$theta23*TIn$eta2  - TIn$eta3)/(1 - TIn$theta23^2)) * LgTRI$p23.g   - LgTRI$d11.23 * d23/LgTRI$sd.23 * (TIn$theta13 - TIn$theta12 * TIn$theta23)/(1 - TIn$theta23^2) ) ) 

dd11.12.dtheta12 <- (LgTRI$d11.12/(1 - TIn$theta12^2)) * ( ( (TIn$theta12 * TIn$eta2  - TIn$eta1) * (TIn$theta12 * TIn$eta1 - TIn$eta2)/(1 - TIn$theta12^2) ) + TIn$theta12)
dd11.13.dtheta13 <- (LgTRI$d11.13/(1 - TIn$theta13^2)) * ( ( (TIn$theta13 * TIn$eta3 - TIn$eta1) * (TIn$theta13 * TIn$eta1 - TIn$eta3)/(1 - TIn$theta13^2) ) + TIn$theta13)
dd11.23.dtheta23 <- (LgTRI$d11.23/(1 - TIn$theta23^2)) * ( ( (TIn$theta23 * TIn$eta3 - TIn$eta2) * (TIn$theta23 * TIn$eta2  - TIn$eta3)/(1 - TIn$theta23^2) ) + TIn$theta23)


dmean12.dtheta12 <- ( (-TIn$eta1 * TIn$theta23 - TIn$eta2 * TIn$theta13) * (1 - TIn$theta12^2) + 2 * TIn$theta12 * (TIn$eta1 * (TIn$theta13 - TIn$theta12 * TIn$theta23) + TIn$eta2 * (TIn$theta23 - TIn$theta12 * TIn$theta13)))/(1 - TIn$theta12^2)^2
dmean13.dtheta13 <- ( (-TIn$eta1 * TIn$theta23 - TIn$eta3 * TIn$theta12) * (1 - TIn$theta13^2) + 2 * TIn$theta13 * (TIn$eta1 * (TIn$theta12 - TIn$theta13 * TIn$theta23) + TIn$eta3 * (TIn$theta23 - TIn$theta12 * TIn$theta13)))/(1 - TIn$theta13^2)^2
dmean23.dtheta23 <- ( (-TIn$eta2 * TIn$theta13 - TIn$eta3 * TIn$theta12) * (1 - TIn$theta23^2) + 2 * TIn$theta23 * (TIn$eta2 * (TIn$theta12 - TIn$theta13 * TIn$theta23) + TIn$eta3 * (TIn$theta13 - TIn$theta12 * TIn$theta23)))/(1 - TIn$theta23^2)^2

dmean13.dtheta12 <- ( TIn$eta1 - TIn$eta3 * TIn$theta13 )/( 1 - TIn$theta13^2 )
dmean23.dtheta12 <- ( TIn$eta2 - TIn$eta3 * TIn$theta23 )/( 1 - TIn$theta23^2 )
dmean23.dtheta13 <- ( TIn$eta3 - TIn$eta2  * TIn$theta23 )/( 1 - TIn$theta23^2 )

dvar12.dtheta12 <- ( 2 * (TIn$theta13 * TIn$theta23 - TIn$theta12) * (1 - TIn$theta12^2) + 2 * TIn$theta12 * deno )/(1 - TIn$theta12^2)^2
dvar13.dtheta13 <- ( 2 * (TIn$theta12 * TIn$theta23 - TIn$theta13) * (1 - TIn$theta13^2) + 2 * TIn$theta13 * deno )/(1 - TIn$theta13^2)^2
dvar23.dtheta23 <- ( 2 * (TIn$theta12 * TIn$theta13 - TIn$theta23) * (1 - TIn$theta23^2) + 2 * TIn$theta23 * deno )/(1 - TIn$theta23^2)^2

dvar13.dtheta12 <- ( 2 * (TIn$theta13 * TIn$theta23 - TIn$theta12 ) )/(1 - TIn$theta13^2)
dvar23.dtheta12 <- ( 2 * (TIn$theta13 * TIn$theta23 - TIn$theta12 ) )/(1 - TIn$theta23^2)
dvar23.dtheta13 <- ( 2 * (TIn$theta12 * TIn$theta23 - TIn$theta13 ) )/(1 - TIn$theta23^2)


d2l.dtheta12.theta12 <- respvec$y1.y2.y3 * ( -1/TIn$p111^2 * (LgTRI$d11.12 * LgTRI$p12.g   )^2 + 1/TIn$p111 * (( dd11.12.dtheta12 * LgTRI$p12.g  ) - ((d12/LgTRI$sd.12) * (dmean12.dtheta12  + (((TIn$eta3 - LgTRI$mean.12)/(2 * LgTRI$sd.12^2)) * dvar12.dtheta12)) * LgTRI$d11.12) )) +
  respvec$y1.y2.cy3 * ( -1/TIn$p110^2 * (LgTRI$d11.12 * LgTRI$p12.g.c )^2 + 1/TIn$p110 * (( dd11.12.dtheta12 * LgTRI$p12.g.c) + ((d12/LgTRI$sd.12) * (dmean12.dtheta12  + (((TIn$eta3 - LgTRI$mean.12)/(2 * LgTRI$sd.12^2)) * dvar12.dtheta12)) * LgTRI$d11.12) )) -
  respvec$cy1.y2.y3 * (  1/TIn$p011^2 * (LgTRI$d11.12 * LgTRI$p12.g   )^2 + 1/TIn$p011 * (( dd11.12.dtheta12 * LgTRI$p12.g  ) - ((d12/LgTRI$sd.12) * (dmean12.dtheta12  + (((TIn$eta3 - LgTRI$mean.12)/(2 * LgTRI$sd.12^2)) * dvar12.dtheta12)) * LgTRI$d11.12) )) -
  respvec$cy1.y2.cy3 * (  1/TIn$p010^2 * (LgTRI$d11.12 * LgTRI$p12.g.c )^2 + 1/TIn$p010 * (( dd11.12.dtheta12 * LgTRI$p12.g.c) + ((d12/LgTRI$sd.12) * (dmean12.dtheta12  + (((TIn$eta3 - LgTRI$mean.12)/(2 * LgTRI$sd.12^2)) * dvar12.dtheta12)) * LgTRI$d11.12) )) +
  respvec$cy1.cy2.cy3 * ( -1/TIn$p000^2 * (LgTRI$d11.12 * LgTRI$p12.g.c )^2 + 1/TIn$p000 * (( dd11.12.dtheta12 * LgTRI$p12.g.c) + ((d12/LgTRI$sd.12) * (dmean12.dtheta12  + (((TIn$eta3 - LgTRI$mean.12)/(2 * LgTRI$sd.12^2)) * dvar12.dtheta12)) * LgTRI$d11.12) )) +
  respvec$cy1.cy2.y3 * ( -1/TIn$p001^2 * (LgTRI$d11.12 * LgTRI$p12.g   )^2 + 1/TIn$p001 * (( dd11.12.dtheta12 * LgTRI$p12.g  ) - ((d12/LgTRI$sd.12) * (dmean12.dtheta12  + (((TIn$eta3 - LgTRI$mean.12)/(2 * LgTRI$sd.12^2)) * dvar12.dtheta12)) * LgTRI$d11.12) )) -
  respvec$y1.cy2.cy3 * (  1/TIn$p100^2 * (LgTRI$d11.12 * LgTRI$p12.g.c )^2 + 1/TIn$p100 * (( dd11.12.dtheta12 * LgTRI$p12.g.c) + ((d12/LgTRI$sd.12) * (dmean12.dtheta12  + (((TIn$eta3 - LgTRI$mean.12)/(2 * LgTRI$sd.12^2)) * dvar12.dtheta12)) * LgTRI$d11.12) )) -
  respvec$y1.cy2.y3 * (  1/TIn$p101^2 * (LgTRI$d11.12 * LgTRI$p12.g   )^2 + 1/TIn$p101 * (( dd11.12.dtheta12 * LgTRI$p12.g  ) - ((d12/LgTRI$sd.12) * (dmean12.dtheta12  + (((TIn$eta3 - LgTRI$mean.12)/(2 * LgTRI$sd.12^2)) * dvar12.dtheta12)) * LgTRI$d11.12) )) 

d2l.dtheta13.theta13 <- respvec$y1.y2.y3 * ( -1/TIn$p111^2 * (LgTRI$d11.13 * LgTRI$p13.g)^2 + 1/TIn$p111 * (( dd11.13.dtheta13 * LgTRI$p13.g  ) - ((d13/LgTRI$sd.13) * (dmean13.dtheta13 + (((TIn$eta2 - LgTRI$mean.13)/(2 * LgTRI$sd.13^2)) * dvar13.dtheta13)) * LgTRI$d11.13) )) -
  respvec$y1.y2.cy3 * (  1/TIn$p110^2 * (LgTRI$d11.13 * LgTRI$p13.g)^2 + 1/TIn$p110 * ((dd11.13.dtheta13 * LgTRI$p13.g  ) - ((d13/LgTRI$sd.13) * (dmean13.dtheta13 + (((TIn$eta2 - LgTRI$mean.13)/(2 * LgTRI$sd.13^2)) * dvar13.dtheta13)) * LgTRI$d11.13) )) -
  respvec$cy1.y2.y3 * (  1/TIn$p011^2 * (LgTRI$d11.13 * LgTRI$p13.g)^2 + 1/TIn$p011 * ((dd11.13.dtheta13 * LgTRI$p13.g  ) - ((d13/LgTRI$sd.13) * (dmean13.dtheta13 + (((TIn$eta2 - LgTRI$mean.13)/(2 * LgTRI$sd.13^2)) * dvar13.dtheta13)) * LgTRI$d11.13) )) +
  respvec$cy1.y2.cy3 * ( -1/TIn$p010^2 * (LgTRI$d11.13 * LgTRI$p13.g)^2 + 1/TIn$p010 * ((dd11.13.dtheta13 * LgTRI$p13.g  ) - ((d13/LgTRI$sd.13) * (dmean13.dtheta13 + (((TIn$eta2 - LgTRI$mean.13)/(2 * LgTRI$sd.13^2)) * dvar13.dtheta13)) * LgTRI$d11.13) )) +
  respvec$cy1.cy2.cy3 * ( -1/TIn$p000^2 * (LgTRI$d11.13 * LgTRI$p13.g.c )^2 + 1/TIn$p000 * ((dd11.13.dtheta13 * LgTRI$p13.g.c) + ((d13/LgTRI$sd.13) * (dmean13.dtheta13 + (((TIn$eta2 - LgTRI$mean.13)/(2 * LgTRI$sd.13^2)) * dvar13.dtheta13)) * LgTRI$d11.13) )) -
  respvec$cy1.cy2.y3 * (  1/TIn$p001^2 * (LgTRI$d11.13 * LgTRI$p13.g.c )^2 + 1/TIn$p001 * ((dd11.13.dtheta13 * LgTRI$p13.g.c) + ((d13/LgTRI$sd.13) * (dmean13.dtheta13 + (((TIn$eta2 - LgTRI$mean.13)/(2 * LgTRI$sd.13^2)) * dvar13.dtheta13)) * LgTRI$d11.13) )) -
  respvec$y1.cy2.cy3 * (  1/TIn$p100^2 * (LgTRI$d11.13 * LgTRI$p13.g.c )^2 + 1/TIn$p100 * ((dd11.13.dtheta13 * LgTRI$p13.g.c) + ((d13/LgTRI$sd.13) * (dmean13.dtheta13 + (((TIn$eta2 - LgTRI$mean.13)/(2 * LgTRI$sd.13^2)) * dvar13.dtheta13)) * LgTRI$d11.13) )) +
  respvec$y1.cy2.y3 * ( -1/TIn$p101^2 * (LgTRI$d11.13 * LgTRI$p13.g.c )^2 + 1/TIn$p101 * ((dd11.13.dtheta13 * LgTRI$p13.g.c) + ((d13/LgTRI$sd.13) * (dmean13.dtheta13 + (((TIn$eta2 - LgTRI$mean.13)/(2 * LgTRI$sd.13^2)) * dvar13.dtheta13)) * LgTRI$d11.13) ))   

d2l.dtheta23.theta23 <- respvec$y1.y2.y3 * ( -1/TIn$p111^2 * (LgTRI$d11.23 * LgTRI$p23.g   )^2 + 1/TIn$p111 * ((dd11.23.dtheta23 * LgTRI$p23.g  ) - ((d23/LgTRI$sd.23) * (dmean23.dtheta23 + ((TIn$eta1 - LgTRI$mean.23)/(2*LgTRI$sd.23^2) * dvar23.dtheta23)) * LgTRI$d11.23) )) -
  respvec$y1.y2.cy3 * (  1/TIn$p110^2 * (LgTRI$d11.23 * LgTRI$p23.g   )^2 + 1/TIn$p110 * ((dd11.23.dtheta23 * LgTRI$p23.g  ) - ((d23/LgTRI$sd.23) * (dmean23.dtheta23 + ((TIn$eta1 - LgTRI$mean.23)/(2*LgTRI$sd.23^2) * dvar23.dtheta23)) * LgTRI$d11.23) )) +
  respvec$cy1.y2.y3 * ( -1/TIn$p011^2 * (LgTRI$d11.23 * LgTRI$p23.g.c )^2 + 1/TIn$p011 * ((dd11.23.dtheta23 * LgTRI$p23.g.c) + ((d23/LgTRI$sd.23) * (dmean23.dtheta23 + ((TIn$eta1 - LgTRI$mean.23)/(2*LgTRI$sd.23^2) * dvar23.dtheta23)) * LgTRI$d11.23) )) -
  respvec$cy1.y2.cy3 * (  1/TIn$p010^2 * (LgTRI$d11.23 * LgTRI$p23.g.c )^2 + 1/TIn$p010 * ((dd11.23.dtheta23 * LgTRI$p23.g.c) + ((d23/LgTRI$sd.23) * (dmean23.dtheta23 + ((TIn$eta1 - LgTRI$mean.23)/(2*LgTRI$sd.23^2) * dvar23.dtheta23)) * LgTRI$d11.23) )) +
  respvec$cy1.cy2.cy3 * ( -1/TIn$p000^2 * (LgTRI$d11.23 * LgTRI$p23.g.c )^2 + 1/TIn$p000 * ((dd11.23.dtheta23 * LgTRI$p23.g.c) + ((d23/LgTRI$sd.23) * (dmean23.dtheta23 + ((TIn$eta1 - LgTRI$mean.23)/(2*LgTRI$sd.23^2) * dvar23.dtheta23)) * LgTRI$d11.23) )) -
  respvec$cy1.cy2.y3 * (  1/TIn$p001^2 * (LgTRI$d11.23 * LgTRI$p23.g.c )^2 + 1/TIn$p001 * ((dd11.23.dtheta23 * LgTRI$p23.g.c) + ((d23/LgTRI$sd.23) * (dmean23.dtheta23 + ((TIn$eta1 - LgTRI$mean.23)/(2*LgTRI$sd.23^2) * dvar23.dtheta23)) * LgTRI$d11.23) )) +
  respvec$y1.cy2.cy3 * ( -1/TIn$p100^2 * (LgTRI$d11.23 * LgTRI$p23.g   )^2 + 1/TIn$p100 * ((dd11.23.dtheta23 * LgTRI$p23.g  ) - ((d23/LgTRI$sd.23) * (dmean23.dtheta23 + ((TIn$eta1 - LgTRI$mean.23)/(2*LgTRI$sd.23^2) * dvar23.dtheta23)) * LgTRI$d11.23) )) -
  respvec$y1.cy2.y3 * (  1/TIn$p101^2 * (LgTRI$d11.23 * LgTRI$p23.g   )^2 + 1/TIn$p101 * ((dd11.23.dtheta23 * LgTRI$p23.g  ) - ((d23/LgTRI$sd.23) * (dmean23.dtheta23 + ((TIn$eta1 - LgTRI$mean.23)/(2*LgTRI$sd.23^2) * dvar23.dtheta23)) * LgTRI$d11.23) )) 

d2l.dtheta12.theta13 <- respvec$y1.y2.y3 * ( -1/TIn$p111^2 * LgTRI$d11.12 * LgTRI$p12.g   * LgTRI$d11.13 * LgTRI$p13.g   - 1/TIn$p111 * LgTRI$d11.13 * (d13/LgTRI$sd.13) * (dmean13.dtheta12 + ((TIn$eta2-LgTRI$mean.13)/(2*LgTRI$sd.13^2) * dvar13.dtheta12) )) -
  respvec$y1.y2.cy3 * ( -1/TIn$p110^2 * LgTRI$d11.12 * LgTRI$p12.g.c * LgTRI$d11.13 * LgTRI$p13.g   - 1/TIn$p110 * LgTRI$d11.13 * (d13/LgTRI$sd.13) * (dmean13.dtheta12 + ((TIn$eta2-LgTRI$mean.13)/(2*LgTRI$sd.13^2) * dvar13.dtheta12) )) -
  respvec$cy1.y2.y3 * (  1/TIn$p011^2 * LgTRI$d11.12 * LgTRI$p12.g   * LgTRI$d11.13 * LgTRI$p13.g   - 1/TIn$p011 * LgTRI$d11.13 * (d13/LgTRI$sd.13) * (dmean13.dtheta12 + ((TIn$eta2-LgTRI$mean.13)/(2*LgTRI$sd.13^2) * dvar13.dtheta12) )) +
  respvec$cy1.y2.cy3 * (  1/TIn$p010^2 * LgTRI$d11.12 * LgTRI$p12.g.c * LgTRI$d11.13 * LgTRI$p13.g   - 1/TIn$p010 * LgTRI$d11.13 * (d13/LgTRI$sd.13) * (dmean13.dtheta12 + ((TIn$eta2-LgTRI$mean.13)/(2*LgTRI$sd.13^2) * dvar13.dtheta12) )) +
  respvec$cy1.cy2.cy3 * ( -1/TIn$p000^2 * LgTRI$d11.12 * LgTRI$p12.g.c * LgTRI$d11.13 * LgTRI$p13.g.c + 1/TIn$p000 * LgTRI$d11.13 * (d13/LgTRI$sd.13) * (dmean13.dtheta12 + ((TIn$eta2-LgTRI$mean.13)/(2*LgTRI$sd.13^2) * dvar13.dtheta12) )) -
  respvec$cy1.cy2.y3 * ( -1/TIn$p001^2 * LgTRI$d11.12 * LgTRI$p12.g   * LgTRI$d11.13 * LgTRI$p13.g.c + 1/TIn$p001 * LgTRI$d11.13 * (d13/LgTRI$sd.13) * (dmean13.dtheta12 + ((TIn$eta2-LgTRI$mean.13)/(2*LgTRI$sd.13^2) * dvar13.dtheta12) )) -
  respvec$y1.cy2.cy3 * (  1/TIn$p100^2 * LgTRI$d11.12 * LgTRI$p12.g.c * LgTRI$d11.13 * LgTRI$p13.g.c + 1/TIn$p100 * LgTRI$d11.13 * (d13/LgTRI$sd.13) * (dmean13.dtheta12 + ((TIn$eta2-LgTRI$mean.13)/(2*LgTRI$sd.13^2) * dvar13.dtheta12) )) +
  respvec$y1.cy2.y3 * (  1/TIn$p101^2 * LgTRI$d11.12 * LgTRI$p12.g   * LgTRI$d11.13 * LgTRI$p13.g.c + 1/TIn$p101 * LgTRI$d11.13 * (d13/LgTRI$sd.13) * (dmean13.dtheta12 + ((TIn$eta2-LgTRI$mean.13)/(2*LgTRI$sd.13^2) * dvar13.dtheta12) ))

d2l.dtheta12.theta23 <- respvec$y1.y2.y3 * ( -1/TIn$p111^2 * LgTRI$d11.12 * LgTRI$p12.g   * LgTRI$d11.23 * LgTRI$p23.g   - 1/TIn$p111 * LgTRI$d11.23 * (d23/LgTRI$sd.23) * (dmean23.dtheta12 + ((TIn$eta1-LgTRI$mean.23)/(2*LgTRI$sd.23^2) * dvar23.dtheta12) )) -
  respvec$y1.y2.cy3 * ( -1/TIn$p110^2 * LgTRI$d11.12 * LgTRI$p12.g.c * LgTRI$d11.23 * LgTRI$p23.g   - 1/TIn$p110 * LgTRI$d11.23 * (d23/LgTRI$sd.23) * (dmean23.dtheta12 + ((TIn$eta1-LgTRI$mean.23)/(2*LgTRI$sd.23^2) * dvar23.dtheta12) )) +
  respvec$cy1.y2.y3 * (  1/TIn$p011^2 * LgTRI$d11.12 * LgTRI$p12.g   * LgTRI$d11.23 * LgTRI$p23.g.c + 1/TIn$p011 * LgTRI$d11.23 * (d23/LgTRI$sd.23) * (dmean23.dtheta12 + ((TIn$eta1-LgTRI$mean.23)/(2*LgTRI$sd.23^2) * dvar23.dtheta12) )) -
  respvec$cy1.y2.cy3 * (  1/TIn$p010^2 * LgTRI$d11.12 * LgTRI$p12.g.c * LgTRI$d11.23 * LgTRI$p23.g.c + 1/TIn$p010 * LgTRI$d11.23 * (d23/LgTRI$sd.23) * (dmean23.dtheta12 + ((TIn$eta1-LgTRI$mean.23)/(2*LgTRI$sd.23^2) * dvar23.dtheta12) )) +
  respvec$cy1.cy2.cy3 * ( -1/TIn$p000^2 * LgTRI$d11.12 * LgTRI$p12.g.c * LgTRI$d11.23 * LgTRI$p23.g.c + 1/TIn$p000 * LgTRI$d11.23 * (d23/LgTRI$sd.23) * (dmean23.dtheta12 + ((TIn$eta1-LgTRI$mean.23)/(2*LgTRI$sd.23^2) * dvar23.dtheta12) )) -
  respvec$cy1.cy2.y3 * ( -1/TIn$p001^2 * LgTRI$d11.12 * LgTRI$p12.g   * LgTRI$d11.23 * LgTRI$p23.g.c + 1/TIn$p001 * LgTRI$d11.23 * (d23/LgTRI$sd.23) * (dmean23.dtheta12 + ((TIn$eta1-LgTRI$mean.23)/(2*LgTRI$sd.23^2) * dvar23.dtheta12) )) +
  respvec$y1.cy2.cy3 * (  1/TIn$p100^2 * LgTRI$d11.12 * LgTRI$p12.g.c * LgTRI$d11.23 * LgTRI$p23.g   - 1/TIn$p100 * LgTRI$d11.23 * (d23/LgTRI$sd.23) * (dmean23.dtheta12 + ((TIn$eta1-LgTRI$mean.23)/(2*LgTRI$sd.23^2) * dvar23.dtheta12) )) -
  respvec$y1.cy2.y3 * (  1/TIn$p101^2 * LgTRI$d11.12 * LgTRI$p12.g   * LgTRI$d11.23 * LgTRI$p23.g   - 1/TIn$p101 * LgTRI$d11.23 * (d23/LgTRI$sd.23) * (dmean23.dtheta12 + ((TIn$eta1-LgTRI$mean.23)/(2*LgTRI$sd.23^2) * dvar23.dtheta12) ))

d2l.dtheta13.theta23 <- respvec$y1.y2.y3 * ( -1/TIn$p111^2 * LgTRI$d11.13 * LgTRI$p13.g   * LgTRI$d11.23 * LgTRI$p23.g   - 1/TIn$p111 * LgTRI$d11.23 * (d23/LgTRI$sd.23) * (dmean23.dtheta13 + ((TIn$eta1-LgTRI$mean.23)/(2 * LgTRI$sd.23^2) * dvar23.dtheta13) )) -
  respvec$y1.y2.cy3 * (  1/TIn$p110^2 * LgTRI$d11.13 * LgTRI$p13.g   * LgTRI$d11.23 * LgTRI$p23.g   - 1/TIn$p110 * LgTRI$d11.23 * (d23/LgTRI$sd.23) * (dmean23.dtheta13 + ((TIn$eta1-LgTRI$mean.23)/(2 * LgTRI$sd.23^2) * dvar23.dtheta13) )) +
  respvec$cy1.y2.y3 * (  1/TIn$p011^2 * LgTRI$d11.13 * LgTRI$p13.g   * LgTRI$d11.23 * LgTRI$p23.g.c + 1/TIn$p011 * LgTRI$d11.23 * (d23/LgTRI$sd.23) * (dmean23.dtheta13 + ((TIn$eta1-LgTRI$mean.23)/(2 * LgTRI$sd.23^2) * dvar23.dtheta13) )) -
  respvec$cy1.y2.cy3 * ( -1/TIn$p010^2 * LgTRI$d11.13 * LgTRI$p13.g   * LgTRI$d11.23 * LgTRI$p23.g.c + 1/TIn$p010 * LgTRI$d11.23 * (d23/LgTRI$sd.23) * (dmean23.dtheta13 + ((TIn$eta1-LgTRI$mean.23)/(2 * LgTRI$sd.23^2) * dvar23.dtheta13) )) +
  respvec$cy1.cy2.cy3 * ( -1/TIn$p000^2 * LgTRI$d11.13 * LgTRI$p13.g.c * LgTRI$d11.23 * LgTRI$p23.g.c + 1/TIn$p000 * LgTRI$d11.23 * (d23/LgTRI$sd.23) * (dmean23.dtheta13 + ((TIn$eta1-LgTRI$mean.23)/(2 * LgTRI$sd.23^2) * dvar23.dtheta13) )) -
  respvec$cy1.cy2.y3 * (  1/TIn$p001^2 * LgTRI$d11.13 * LgTRI$p13.g.c * LgTRI$d11.23 * LgTRI$p23.g.c + 1/TIn$p001 * LgTRI$d11.23 * (d23/LgTRI$sd.23) * (dmean23.dtheta13 + ((TIn$eta1-LgTRI$mean.23)/(2 * LgTRI$sd.23^2) * dvar23.dtheta13) )) +
  respvec$y1.cy2.cy3 * (  1/TIn$p100^2 * LgTRI$d11.13 * LgTRI$p13.g.c * LgTRI$d11.23 * LgTRI$p23.g   - 1/TIn$p100 * LgTRI$d11.23 * (d23/LgTRI$sd.23) * (dmean23.dtheta13 + ((TIn$eta1-LgTRI$mean.23)/(2 * LgTRI$sd.23^2) * dvar23.dtheta13) )) -
  respvec$y1.cy2.y3 * ( -1/TIn$p101^2 * LgTRI$d11.13 * LgTRI$p13.g.c * LgTRI$d11.23 * LgTRI$p23.g   - 1/TIn$p101 * LgTRI$d11.23 * (d23/LgTRI$sd.23) * (dmean23.dtheta13 + ((TIn$eta1-LgTRI$mean.23)/(2 * LgTRI$sd.23^2) * dvar23.dtheta13) )) 


dtheta12.dtheta12.st <- 4 * exp( 2 * TIn$theta12.st )/( exp(2 * TIn$theta12.st) + 1 )^2
dtheta13.dtheta13.st <- 4 * exp( 2 * TIn$theta13.st )/( exp(2 * TIn$theta13.st) + 1 )^2
dtheta23.dtheta23.st <- 4 * exp( 2 * TIn$theta23.st )/( exp(2 * TIn$theta23.st) + 1 )^2

d2theta12.theta12.st <- (8 * exp( 2 * TIn$theta12.st ) - 8 * exp( 4 * TIn$theta12.st ))/( exp(2 * TIn$theta12.st) + 1 )^3
d2theta13.theta13.st <- (8 * exp( 2 * TIn$theta13.st ) - 8 * exp( 4 * TIn$theta13.st ))/( exp(2 * TIn$theta13.st) + 1 )^3
d2theta23.theta23.st <- (8 * exp( 2 * TIn$theta23.st ) - 8 * exp( 4 * TIn$theta23.st ))/( exp(2 * TIn$theta23.st) + 1 )^3

d2l.dbe1.be1 <- crossprod( VC$X1 * c(VC$weights*d2l.de1.e1), VC$X1 )
d2l.dbe2.be2 <- crossprod( VC$X2 * c(VC$weights*d2l.de2.e2), VC$X2 )
d2l.dbe3.be3 <- crossprod( VC$X3 * c(VC$weights*d2l.de3.e3), VC$X3 )
d2l.dbe1.be2 <- crossprod( VC$X1 * c(VC$weights*d2l.de1.e2), VC$X2 )
d2l.dbe1.be3 <- crossprod( VC$X1 * c(VC$weights*d2l.de1.e3), VC$X3 )
d2l.dbe2.be3 <- crossprod( VC$X2 * c(VC$weights*d2l.de2.e3), VC$X3 )

d2l.dbe1.theta12.st <- colSums(c(VC$weights*d2l.de1.theta12* dtheta12.dtheta12.st) * VC$X1) 
d2l.dbe1.theta13.st <- colSums(c(VC$weights*d2l.de1.theta13* dtheta13.dtheta13.st) * VC$X1) 
d2l.dbe1.theta23.st <- colSums(c(VC$weights*d2l.de1.theta23* dtheta23.dtheta23.st) * VC$X1) 

d2l.dbe2.theta12.st <- colSums(c(VC$weights*d2l.de2.theta12* dtheta12.dtheta12.st) * VC$X2) 
d2l.dbe2.theta13.st <- colSums(c(VC$weights*d2l.de2.theta13* dtheta13.dtheta13.st) * VC$X2) 
d2l.dbe2.theta23.st <- colSums(c(VC$weights*d2l.de2.theta23* dtheta23.dtheta23.st) * VC$X2) 

d2l.dbe3.theta12.st <- colSums(c(VC$weights*d2l.de3.theta12* dtheta12.dtheta12.st) * VC$X3) 
d2l.dbe3.theta13.st <- colSums(c(VC$weights*d2l.de3.theta13* dtheta13.dtheta13.st) * VC$X3) 
d2l.dbe3.theta23.st <- colSums(c(VC$weights*d2l.de3.theta23* dtheta23.dtheta23.st) * VC$X3) 

d2l.dtheta12.st <- sum( VC$weights*(d2l.dtheta12.theta12 * (dtheta12.dtheta12.st)^2  + LgTRI$dl.dtheta12 * d2theta12.theta12.st) )
d2l.dtheta13.st <- sum( VC$weights*(d2l.dtheta13.theta13 * (dtheta13.dtheta13.st)^2  + LgTRI$dl.dtheta13 * d2theta13.theta13.st) )
d2l.dtheta23.st <- sum( VC$weights*(d2l.dtheta23.theta23 * (dtheta23.dtheta23.st)^2  + LgTRI$dl.dtheta23 * d2theta23.theta23.st) )

d2l.dtheta12.st.theta13.st <- sum(VC$weights*d2l.dtheta12.theta13 * dtheta12.dtheta12.st * dtheta13.dtheta13.st )
d2l.dtheta12.st.theta23.st <- sum(VC$weights*d2l.dtheta12.theta23 * dtheta12.dtheta12.st * dtheta23.dtheta23.st )
d2l.dtheta13.st.theta23.st <- sum(VC$weights*d2l.dtheta13.theta23 * dtheta13.dtheta13.st * dtheta23.dtheta23.st )

h1 <- cbind(d2l.dbe1.be1, d2l.dbe1.be2, d2l.dbe1.be3, d2l.dbe1.theta12.st, d2l.dbe1.theta13.st, d2l.dbe1.theta23.st)
h2 <- cbind(t(d2l.dbe1.be2), d2l.dbe2.be2, d2l.dbe2.be3, d2l.dbe2.theta12.st, d2l.dbe2.theta13.st, d2l.dbe2.theta23.st)
h3 <- cbind(t(d2l.dbe1.be3), t(d2l.dbe2.be3), d2l.dbe3.be3, d2l.dbe3.theta12.st, d2l.dbe3.theta13.st, d2l.dbe3.theta23.st)
h4 <- cbind(t(d2l.dbe1.theta12.st), t(d2l.dbe2.theta12.st), t(d2l.dbe3.theta12.st), d2l.dtheta12.st, d2l.dtheta12.st.theta13.st, d2l.dtheta12.st.theta23.st)
h5 <- cbind(t(d2l.dbe1.theta13.st), t(d2l.dbe2.theta13.st), t(d2l.dbe3.theta13.st), t(d2l.dtheta12.st.theta13.st), d2l.dtheta13.st, d2l.dtheta13.st.theta23.st)
h6 <- cbind(t(d2l.dbe1.theta23.st), t(d2l.dbe2.theta23.st), t(d2l.dbe3.theta23.st), t(d2l.dtheta12.st.theta23.st), t(d2l.dtheta13.st.theta23.st), d2l.dtheta23.st)

HTRIVec <- list(h1 = h1, h2 = h2, h3 = h3, h4 = h4, h5 = h5, h6 = h6)

HTRIVec

}