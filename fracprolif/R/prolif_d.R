prolif_d <-
function(t, x0, q0=0, d, q, ad, aq, tol=1e-6)
{
    x0 * exp((d-q-ad)*t)
}
