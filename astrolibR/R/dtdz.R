# Integrand for cosmic age vs. redshift in standard cosmology

dtdz = function(z, lambda0, q0) {
  term1 = (1.0 + z)
  term2 = 2.0 * (q0 + lambda0) * z + 1.0 - lambda0
  term3 = (1.0 + z) * (1.0 +z)
  return(1.0 / (term1 * sqrt(term2 * term3 + lambda0)))
}
