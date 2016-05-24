"MLnp" <-
function(Xr,yr,iv,tp){np <- ncol(Xr)
# Maximum likelihood
beta <- tp$beta; tn <- tp$tn
# coefficients
# z  <- riclls(Xr,yr); th <- z$theta[1:np] 
  z  <- lsfit(Xr, yr, intercept=FALSE)
  th <- z$coef
# scale
  if (iv==1) {rr  <- as.vector(yr-Xr%*%as.matrix(th)); v <- Scalen(rr,tn-np,beta)}
list(th1=th,v1=v)}

