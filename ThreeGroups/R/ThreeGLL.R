.ThreeGLL <-
function(par, s_b, f_b, s_t_c, f_t_c, s_t_nc, f_t_nc, s_p_c, f_p_c, s_p_nc, f_p_nc) {

  alpha <- par[1]
  p_r   <- par[2]
  p_nr  <- par[3]
  tau   <- par[4]

  fvalue <- 
  s_b    * log(     alpha * (p_r-tau) + (1-alpha) * p_nr ) +
  f_b    * log(1 -( alpha * (p_r-tau) + (1-alpha) * p_nr )) +
  s_t_c  * log( alpha * (   p_r) ) +
  f_t_c  * log( alpha * (1-(p_r)) ) +
  s_t_nc * log( (1 -  alpha) * (    p_nr) ) +
  f_t_nc * log( (1 -  alpha) * (1 - p_nr) ) +
  s_p_c  * log( alpha * (   p_r - tau) ) +
  f_p_c  * log( alpha * (1-(p_r - tau)) ) +
  s_p_nc * log( (1 -  alpha) * (    p_nr) ) +
  f_p_nc * log( (1 -  alpha) * (1 - p_nr) ) 

  if(is.na(fvalue)) fvalue <- -9999999
  if(fvalue == Inf | fvalue == -Inf) fvalue <- -9999999
  return(fvalue)

  }
