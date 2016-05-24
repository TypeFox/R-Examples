kern <-
function(x0,h,freq) {

vec_fin=exp(-diag(t((t(freq) - x0))%*%(t(freq) - x0))/h^2)    
list(v=vec_fin)

}
