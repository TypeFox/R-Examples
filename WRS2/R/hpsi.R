hpsi<-function(x,bend=1.28){
#
#   Evaluate Huber`s Psi function for each value in the vector x
#   The bending constant defaults to 1.28.
#
hpsi<-ifelse(abs(x)<=bend,x,bend*sign(x))
hpsi
}