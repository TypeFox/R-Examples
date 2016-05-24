weights.Bartlett = function(x){ 1 - abs(x) }
weights.trunc = function(x){ x[] = 1 ; x }
weights.Tukey = function(x){ 0.54 + 0.46 * cos(pi*x) }
weights.Parzen = function(x){ res = x ; yes = abs(x) <= 0.5 ; res[yes] = 1 - 6*x[yes]^2+6*abs(x[yes])^3 ; res[!yes] = 2*(1-abs(x[!yes]))^3 ; res }
weights.Bohman = function(x){ ((1-abs(x))^3) *cos(pi*x) + sin(pi*abs(x))/pi }
weights.Daniell = function(x) { res = sin(pi*x) / (pi*x) ; res[res > 1] = 1 ; res[is.nan(res)] = 1 ; res }
#weights.BarlettPriestley = function(x) { min(abs(3*(min(sin(pi*x) / pi*x,1) - cos(pi*x)) / (pi*pi*x*x)),1) }
weights.ParzenCogburnDavis = function(x, r = 1) { 1/(1+abs(x)^(2*r)) }
