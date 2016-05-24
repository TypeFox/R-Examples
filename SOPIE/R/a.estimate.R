a.estimate <-
function (data, to = 1, min_points, alpha = 0.05, g = 1, r = 1)
{

   cl <- match.call()

   #ensure that the domain of the data is correct. Removes all duplicate observations while scaling the data to a domain of [0,1] if required
   if (round(max(data), 1) == round(2 * pi, 1))
       data<-unique(data/(2*pi))
   data<-sort(unique(data))

   data <- c(data - 1, data, data + 1)
   close_to_min_index <- unique(findInterval(min_points, data) + 1)

   CVM_rejectionpnt <- numeric(length(close_to_min_index))
   AD_rejectionpnt <- numeric(length(close_to_min_index))
   KS_rejectionpnt <- numeric(length(close_to_min_index))
   rayleigh_rejectionpnt <- numeric(length(close_to_min_index))

   #This loop iterates the algorithm for each of the unique minimum points specified in the arguent "min_points"
   for (k in 1 : length(close_to_min_index))
   {
       #end represents the last value in the complete interval over which uniformity is tested 
       end <- close_to_min_index[k]
       #begin represents the first value in the complete interval over which uniformity is tested
       begin <- close_to_min_index[k] - (length(data) / 3)
       start <- begin
       endpoint <- end
       lengte <- endpoint-start-1

       # Create empty vectors for the different GOF test statistics and p values
       ks <- numeric(trunc(lengte / g))
       ks_p <- numeric(trunc(lengte / g))
       c <- numeric(trunc(lengte / g))
       d <- numeric(trunc(lengte / g))
       asq <- numeric(trunc(lengte / g))
       asq_p <- numeric(trunc(lengte / g))
       CVM_pvalue <- numeric(trunc(lengte / g))
       rayleigh_p <- numeric(trunc(lengte / g))
       tol <- 1e-20



       # j represents the set number
       j <- 1

       #Boolean operators to establish when all GOF tests have been rejected
       CVM_rejection <- FALSE
       KS_rejection <- FALSE
       AD_rejection <- FALSE
       rayleigh_rejection <- FALSE
       rejected <- (CVM_rejection & KS_rejection & AD_rejection & rayleigh_rejection)


       while (!rejected & (j <= trunc(lengte / g)))
       {
              sn<-0
              n<-0
              sn_star<-0
              star_pval<-0
              start <- endpoint - (j * g) - 1
              n <- endpoint - start - 1
      d[j] <- n
              
              # The next lines contains the calculation of the test statistics and p-values for the different GOF tests.

              #Anderson-Darling test statistic and p-value calculation
              if (n < 8)
              {
                  asq[j] <- 0
                  asq_p[j] <-0.999
              }
              else 
              {
                  sn <- ad.test(data[(start + 1) : (endpoint - 1)],punif, data[start]-tol, data[endpoint]+tol)
  asq[j] <-sn$statistic
  asq_p[j] <-sn$p.value

              }

      #Cramer-von Mises test statistic and p-value calculation (Stephens 1970, p101,p105)
              sn <- ((data[(start + 1) : (endpoint - 1)] - data[start]) / (data[endpoint] - data[start]) - (1 : n - 0.5)/(n))^2
              c[j] <- sum(sn) + (1 / (12 * n))
              CVM_pvalue[j] <- CVM_pvall(c[j], n)

      #Kolmogorov-Smirnov test statistic and p-value calculation (ks.test function in package "stats")
              ks[j] <- ks.test(data[(start + 1) : (endpoint - 1)], punif, data[start], data[endpoint])$statistic
              ks_p[j] <- ks.test(data[(start + 1) : (endpoint - 1)], punif, data[start], data[endpoint])$p.value

              

      #Rayleigh p-value calculation
              #The test statistic is n * R ^2 (Mardia and Jupp 2000) and the p values are calculated with the approximation proposed by 
              #Greenwood and Durand (1955) in the circular package
            
               rayleigh_p[j] <- rayleigh.test(circular(((data[(start + 1) : (endpoint - 1)]) - data[(start + 1)]) / (data[(endpoint - 1)] - data[(start + 1)]) * 2 * pi))$p.value
            
               if (rayleigh_p[j] == "NaN")
                  rayleigh_p[j] <- 1


              if (j >= r)
              {
                  #Calulate the point of rejection of uniformity for each GOF test and for each selected minimum point by establishing the point where the
                   #p-value(s) is less than alpha

                  if (CVM_rejection == FALSE)
                  {
                      if (max(CVM_pvalue[(j - r + 1) : j]) < alpha)
                      {
                          CVM_rejection <- TRUE
                          if (endpoint - 1 - ((j - r + 1) * g) >= ((length(data) / 3) + 1))
                              CVM_rejectionpnt[k] <- (data[endpoint - 1 - ((j - r + 1) * g)]) / 1
                          else CVM_rejectionpnt[k] <- (data[endpoint - 1- ((j - r + 1) * g) + (length(data) / 3)])/1
                      }
                  }

                  if (KS_rejection == FALSE)
                  {
                      if (max(ks_p[(j - r + 1) : j]) < alpha)
                      {
                          KS_rejection <- TRUE
                          if (endpoint - 1 - ((j - r + 1) * g) >= ((length(data) / 3) + 1))
                              KS_rejectionpnt[k] <- (data[endpoint - 1 - ((j - r + 1) * g)]) / 1
                          else KS_rejectionpnt[k] <- (data[endpoint - 1 - ((j - r + 1) * g)+(length(data) / 3)]) / 1
                      }
                  }

                  if (AD_rejection == FALSE)
                  {
                      if (max(asq_p[(j - r + 1) : j]) < alpha)
                      {
                          AD_rejection <- TRUE
                          if (endpoint - 1 - ((j - r + 1) * g) >= ((length(data) / 3) + 1))
                              AD_rejectionpnt[k] <- (data[endpoint - 1 - ((j - r + 1) * g)]) / 1
                          else AD_rejectionpnt[k] <- (data[endpoint - 1 - ((j - r + 1) * g) + (length(data) / 3)]) / 1
                      }
                  }

                  if (rayleigh_rejection == FALSE)
                  {
                      if (max(rayleigh_p[(j - r + 1) : j]) < alpha)
                      {
                          rayleigh_rejection <- TRUE
                          if (endpoint - 1 - ((j - r + 1) * g) >= ((length(data) / 3) + 1)) 
                              rayleigh_rejectionpnt[k] <- (data[endpoint - 1 - ((j - r + 1) * g)]) / 1
                          else rayleigh_rejectionpnt[k] <- (data[endpoint - 1 - ((j - r + 1) * g) + (length(data) / 3)]) / 1
                      }
                  }

             }
             # if j > rejection loop end

             rejected <- (CVM_rejection & KS_rejection & AD_rejection & rayleigh_rejection) 
             j <- j + 1
       }
       # while rejection loop end

       #Clean up the vector of test statistics and p-values before the next minimum point is analysed
       if (!(is.na(which(c == 0)[1])))
           CVM_pvalue <- CVM_pvalue[1 : ((which(CVM_pvalue == 0)[1]) - 1)]
       
       #if (!(is.na(which(d == 0)[1])))
       #    d <- d[1 : ((which(d == 0)[1] - 1))]
       
       if (!(is.na(which(ks == 0)[1])))
           ks <- ks[1 : ((which(ks == 0)[1] - 1))]
       
       if (!(is.na(which(ks_p == 0)[1]))) 
           ks_p <- ks_p[1 : ((which(ks_p == 0)[1] - 1))]
       
       if (!(is.na(which(asq == 0)[1])))
           asq <- asq[1 : ((which(asq == 0)[1] - 1))]
       
       if (!(is.na(which(asq_p == 0)[1])))
           asq_p <- asq_p[1 : ((which(asq_p == 0)[1] - 1))]
       
       if (!(is.na(which(rayleigh_p == 0)[1])))
           rayleigh_p<-rayleigh_p[1 : ((which(rayleigh_p == 0)[1] - 1))]
       
   }
   #for k loop end (number of minimum points)

   #Combine results in different data structures for output
   vec <- c(mean(CVM_rejectionpnt), mean(KS_rejectionpnt), mean(AD_rejectionpnt), mean(rayleigh_rejectionpnt))
   sum <- matrix(c(mean(CVM_rejectionpnt), mean(KS_rejectionpnt), mean(AD_rejectionpnt), mean(rayleigh_rejectionpnt), median(vec)), ncol=5, dimnames=list(c("a-hat"),
          c("Cramer von Mises", "Kolmogorov-Smirnoff", "Anderson-Darling", "Rayleigh", "MEDIAN")))

   CVM <- list(rejection = CVM_rejectionpnt, mean = mean(CVM_rejectionpnt))
   KS <- list(rejection = KS_rejectionpnt, mean = mean(KS_rejectionpnt))
   AD <- list(rejection = AD_rejectionpnt, mean = mean(AD_rejectionpnt))
   rayleigh <- list(rejection = rayleigh_rejectionpnt, mean = mean(rayleigh_rejectionpnt))
   general<-list(call = cl, Minimums = data[close_to_min_index], alpha = alpha, grow = g, nr_reject = r, kernel_function = "Epanechnikov")

   comb<-list(summary = sum, General = general)

   return(comb)
}
