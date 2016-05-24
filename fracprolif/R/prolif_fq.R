prolif_fq <-
function(t, x0, q0=0, d, q, ad, aq, tol=1e-6) 
{
  prolif_q(t, x0, q0, d, q, ad, aq, tol) / prolif_tot(t, x0, q0, d, q, ad, aq, tol)
}
