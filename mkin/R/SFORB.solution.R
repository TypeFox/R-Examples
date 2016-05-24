SFORB.solution = function(t, parent.0, k_12, k_21, k_1output) {
  sqrt_exp = sqrt(1/4 * (k_12 + k_21 + k_1output)^2 + k_12 * k_21 - (k_12 + k_1output) * k_21)
  b1 = 0.5 * (k_12 + k_21 + k_1output) + sqrt_exp
  b2 = 0.5 * (k_12 + k_21 + k_1output) - sqrt_exp
  
  parent = parent.0 * 
        (((k_12 + k_21 - b1)/(b2 - b1)) * exp(-b1 * t) +
        ((k_12 + k_21 - b2)/(b1 - b2)) * exp(-b2 * t))
}
