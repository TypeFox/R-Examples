cell_effect_mult_or <-
function(x,y,z,w,q)
{ q^(-1)*((1+x*y/y*z/z)^(-1)+(1+x*y/y*z)^(-1)*w)*((1+x*y/y*z/z)^(-1)+(1+x*y/y*z)^(-1)*w*z)^(-1)*
((1+x*y)^(-1)+(1+x*y*z*q)^(-1)*w*z*q)*((1+x*y*z/z)^(-1)+(1+x*y*z*q)^(-1)*w)^(-1)}
