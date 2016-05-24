"dst" <-
function(x,nu=3){
gamma( (nu+1)/2)/(gamma(nu/2)*sqrt(pi)*sqrt(nu-2)*( 1+ (x^2)/(nu-2))^((nu+1)/2))
}
