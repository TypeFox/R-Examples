"med.regressCOP2" <-
function(v=seq(0.01,0.99, by=0.01), cop=NULL, para=NULL, ...) {
  return(qua.regressCOP2(f=0.5, v=v, cop=cop, para=para, ...))
}
