prolif_q <-
function(t, x0, q0=0, d, q, ad, aq, tol=1e-6) {
  # Handle limit cases that would otherwise result in division by zero
  if(abs(d-q) < tol)
  {
    if(abs(ad - aq) < tol)
    {
      exp(-aq*t) * (q0 + x0*q*t)
    } else {
      exp(-aq*t) / (ad-aq) * (
        x0 * q * (1-exp((aq-ad)*t)) + 
        (ad - aq)*q0
      )
    }
  } else {
    beta <- d - q + aq - ad
    exp(-aq*t) / beta * (
      q * (x0*(exp(beta*t)-1)-q0) + 
      (beta+q)*q0
    ) 
  }
}
