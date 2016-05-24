"MLn1" <-
function(yr,iv,tp){
beta  <- tp$beta; tn <- tp$tn
# location
  th <- mean(yr)
# scale
  if (iv==1) v <- Scalen(yr-th,tn-1,beta)
list(th1=th,v1=v)}

