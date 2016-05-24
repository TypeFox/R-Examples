gvar<-function(m){
#
# Compute the generalized variance of a matrix m
#
m<-elimna(m)
temp<-var(m)
gvar<-prod(eigen(temp)$values)
gvar
}
