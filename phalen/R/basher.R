basher <-
function(X, A, K) {
  
  if (A<=K) {
    if (A<1) {
      E = 1 - A
      X = X + E; A = A + E; K = K + E 
    } else {
      E = 0
    }
    r = 1/(K - A)
    M = A-(log(1/(K * r))/ -r)
    y = ifelse(X<A, X, K * (1 - exp( -r * (X - M)))) - E
  } else { 
    if (K<1) {
      E = 1 - K
      X = X + E; A = A + E; K = K + E
    } else {
      E = 0
    }
    r = 1/(A - K)
    M = -(log((A/K) - 1)/r) + A
    y = ifelse(X>A, X, K * (1 + exp(r * (X - M)))) - E
  }
  structure(list(y = y,
                 A = A - E,
                 K = K - E,
                 r = r,
                 M = M,
                 E = E,
                 penalty = ifelse(A<=K,'growth','decay')), class = 'basher')
}
