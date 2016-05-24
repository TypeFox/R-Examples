CVM_pvall <-
function (sn, n)
{
   sn_star <- 0
   star_pval <- 0
   sn_star <- (sn - (0.4 / n) + (0.6 / (n^2))) * (1 + (1 / n))
   if (sn_star >= 2.303) 
       star_pval <- 0.0001
   else if (sn_star >= 2.098) 
            star_pval <- approx(c(2.098, 2.303), c(0.005, 0.001), sn_star)$y
   else if (sn_star >= 2.001) 
            star_pval <- approx(c(2.001, 2.098), c(0.01, 0.005), sn_star)$y
   else if (sn_star >= 1.862) 
            star_pval <- approx(c(1.862, 2.001), c(0.025, 0.01), sn_star)$y
   else if (sn_star >= 1.747) 
            star_pval <- approx(c(1.747, 1.862), c(0.05, 0.025), sn_star)$y
   else if (sn_star >= 1.62)  
            star_pval <- approx(c(1.62, 1.747), c(0.1, 0.05), sn_star)$y
   else if (sn_star >= 1.537) 
            star_pval <- approx(c(1.537, 1.62), c(0.15, 0.1), sn_star)$y
   else if (sn_star >= 1.420) 
            star_pval <- approx(c(1.420, 1.537), c(0.25, 0.15), sn_star)$y
   else if (sn_star >= 0.976) 
            star_pval <- approx(c(0.976, 1.42), c(0.85, 0.25), sn_star)$y
   else if (sn_star >= 0.928) 
            star_pval <- approx(c(0.928, 0.976), c(0.9, 0.85), sn_star)$y
   else if (sn_star >= 0.861) 
            star_pval <- approx(c(0.861, 0.928), c(0.95, 0.9), sn_star)$y
   else if (sn_star >= 0.81)  
            star_pval <- approx(c(0.81, 0.861), c(0.975, 0.95), sn_star)$y
   else if (sn_star >= 0.755) 
            star_pval <- approx(c(0.755, 0.81), c(0.99, 0.975), sn_star)$y
   else star_pval <- 0.9999

   return(star_pval)
}
