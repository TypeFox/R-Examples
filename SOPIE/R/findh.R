findh <-
function (data, h = 1, to = 1) 
{
   #validate whether the choice of h is between 1 and 9 according to Table 1
   if ((h < 1) || (h > 9) || is.na(h)) 
           stop("\n", "invalid smoothing parameter index")

   #ensure that the domain of the data is correct. Removes all duplicate observations while scaling the data to a domain of [0,1] if required
   if (round(to, 2) == round(2 * pi, 2))
       data <- unique(data / (2 * pi))

   #Create Objects of class circular (require package "circular")
   datac <- circular(data * 2 * pi, type = "angles", units = "radians", template = "none", modulo = "asis", zero = 0, rotation = "counter")

   if (h == 1)
   {
       linearSD <- sd(data)
       h_hat_plugin <- 1.06 * linearSD * length(data) ^ (-1 / 5)
   }
   else if (h == 2)
   {
            Circular_SD <- sqrt(var.circular(datac))
            h_hat_plugin <- 1.06 * Circular_SD * length(data) ^ (-1 / 5)
   }
   else if (h == 3)
   {
            Circular_mean_dev <- meandeviation(datac) / (2 * pi)
            h_hat_plugin <- 1.06 * Circular_mean_dev * length(data) ^ (-1 / 5)
   }
   else if (h == 4)
   {
            Circular_median <- as.numeric(median.circular(datac)) / (2 * pi)
            Circular_median_abs_dev <- median(abs(data - Circular_median))
            h_hat_plugin <- 1.06 * Circular_median_abs_dev * length(data) ^ (-1 / 5)
   }
   else if (h == 5)
   {
            Circular_Quantiles <- as.vector(quantile.circular(datac, probs = c(0.25, 0.75)))
            Circular_IQR <- abs((Circular_Quantiles[2]-Circular_Quantiles[1]) / (2 * pi))
            h_hat_plugin <- 1.06 * Circular_IQR * length(data) ^ (-1 / 5)
   }
   else if (h == 6)
   {
            Circular_Quantiles <- as.vector(quantile.circular(datac, probs = c(0.25, 0.75)))
            Circular_IQR <- abs((Circular_Quantiles[2]-Circular_Quantiles[1])/(2 * pi))
            h_hat_plugin <- 0.79 * Circular_IQR * length(data) ^ (-1 / 5)
   }
   else if (h == 7)
   {
            Circular_SD <- sqrt(var.circular(datac))
            h_hat_plugin <- 0.9 * Circular_SD * length(data) ^ (-1 / 5)
   }
   else if (h == 8)
   {
            Circular_Quantiles <- as.vector(quantile.circular(datac, probs=c(0.25, 0.75)))
            Circular_IQR <- abs((Circular_Quantiles[2] - Circular_Quantiles[1])/(2 * pi))
            h_hat_plugin <- 0.9 * Circular_IQR / 1.349 * length(data) ^ (-1 / 5)
   }
   else
   {
      linearSD <- sd(data)
      Circular_var <- var.circular(datac)
      CircularSD <- sqrt(Circular_var)
      Circular_mean_dev <- meandeviation(datac) / (2 * pi)
      Circular_median <- as.numeric(median.circular(datac)) / (2 * pi)
      Circular_median_abs_dev <- median(abs(data-Circular_median))
      Circular_Quantiles <- as.vector(quantile.circular(datac, probs = c(0.25, 0.75)))
      Circular_IQR <- abs((Circular_Quantiles[2] - Circular_Quantiles[1]) / (2 * pi))
      
      h_hat_plugin <- numeric(8)
      h_hat_plugin[1] <- 1.06 * linearSD * length(data) ^ (-1 / 5)
      h_hat_plugin[2] <- 1.06 * CircularSD * length(data) ^ (-1 / 5)
      h_hat_plugin[3] <- 1.06 * Circular_mean_dev * length(data) ^ (-1 / 5)
      h_hat_plugin[4] <- 1.06 * Circular_median_abs_dev * length(data) ^ (-1 / 5)
      h_hat_plugin[5] <- 1.06 * Circular_IQR * length(data) ^ (-1 / 5)
      h_hat_plugin[6] <- 0.79 * Circular_IQR * length(data) ^ (-1 / 5)
      h_hat_plugin[7] <- 0.9 * CircularSD * length(data) ^ (-1 / 5)
      h_hat_plugin[8] <- 0.9 * Circular_IQR / 1.349 * length(data) ^ (-1 / 5)
  }

    return(round(mean(h_hat_plugin), 2))
}
