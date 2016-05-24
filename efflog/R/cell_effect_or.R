cell_effect_or <-
function(x,y,z,w)
{ ((1+x*y/y*z/z)^(-1)+(1+x*y/y*z)^(-1)*w)*((1+x*y/y*z/z)^(-1)+(1+x*y/y*z)^(-1)*w*z)^(-1)*
((1+x*y*z/z)^(-1)+(1+x*y*z)^(-1)*w*z)*((1+x*y*z/z)^(-1)+(1+x*y*z)^(-1)*w)^(-1)}
