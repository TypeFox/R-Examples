"med.regressCOP" <-
function(u=seq(0.01,0.99, by=0.01), cop=NULL, para=NULL, ...) {
  return(qua.regressCOP(f=0.5, u=u, cop=cop, para=para, ...))
}
