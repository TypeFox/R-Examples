nudge.fit <-
function(data, avg = NULL, weights = NULL, weights.cutoff= -1.345, 
  pi = NULL, mu = NULL, sigma = NULL, tol=1e-5, max.iter=2000, z = NULL)
{
  inudge.fit(data, avg = avg, K = 1, weights = weights, 
    weights.cutoff = weights.cutoff, pi = pi, mu = mu,
    sigma = sigma, tol=tol, max.iter=max.iter, z = z);
}