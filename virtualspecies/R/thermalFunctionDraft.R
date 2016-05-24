########## WORK IN PROGRESS ##########
# 
# The function "thermalResponse" is currently being improved so that its  
# two parameters (rho and sigma) can be automatically adjusted according to
# user-defined critical temperature thresholds (CTmin and CTmax).
# 
# 
# x <- seq(0, 30, length = 1000)
# 
# thermalResponse <- function(x, To, rho, sigma)
# {
#   exp(-exp(rho * (x - To) - 6) - sigma * (x - To)^2)
# }
# 
# rescaled.thermalResponse <- function(x, To, rho, sigma)
# {
#   (thermalResponse(x, To, rho, sigma) - min(thermalResponse(x, To, rho, sigma))) /
#     max(thermalResponse(x, To, rho, sigma) - min(thermalResponse(x, To, rho, sigma)))
# }
# 
# par(mfrow = c(3, 3))
# resp <- thermalResponse(x = x, To = 5, rho = exp(100), sigma = 0.1)
# 
# if(rescale)
# {
#   tf <- rescaled.thermalResponse
# } else
# {
#   tf <- thermalResponse
# }
# 
# CTmax <- 16
# To = 15
# CTmin <- 0
# 
# precision.CTmax <- 0.01
# precision.CTmin <- 0.01
# 
# 
# rho.min <- -10
# rho.max <- 10
# sigma.min <- 0.0001
# sigma.max <- 10000
# 
# 
# ############# SECTION SIGMA ####################
# x0 <- x
# xmin <- min(range(x)) - 2 * diff(range(x))
# xmax <- max(range(x)) + 2 * diff(range(x))
# x <- seq(xmin, xmax, length = 100000)
# resp.min <- thermalResponse(x = x, To = To, rho = exp(rho.min), sigma = sigma.min)
# while(min(resp.min) > precision.CTmin & sigma.min != sigma.max)
# {
#   sigma.min <- sigma.min * 10
#   resp.min <- thermalResponse(x = x, To = To, rho = exp(rho.min), sigma = sigma.min)
# }
# 
# resp.max <- thermalResponse(x = x, To = To, rho = exp(rho.min), sigma = sigma.max)
# 
# 
# 
# 
# if(CTmin > max(x[resp.max < 0.01 & x < To]))
# {
#   stop(paste("CTmin is too high, try to decrease it (maximum possible value: ", 
#              round(max(x[resp.max < 0.01 & x < To]), 3), ")", sep = ""))
# }
# 
# epsilon.up <- CTmin - max(x[resp.min < 0.01 & x < To])
# epsilon.bot <- CTmin - max(x[resp.max< 0.01 & x < To])
# i <- 1
# while(min(abs(epsilon.up), abs(epsilon.bot)) > precision.CTmin)
# {
#   sigma.mid <- (sigma.min + sigma.max) / 2
#   resp.mid <- thermalResponse(x = x, To = To, rho = exp(rho.min), sigma = sigma.mid)
#   epsilon.mid <- CTmin - max(x[resp.mid < 0.01 & x < To])
#   
#   if(epsilon.mid < 0)
#   {
#     sigma.max <- sigma.mid
#     resp.max <- thermalResponse(x = x, To = To, rho = exp(rho.min), sigma = sigma.max)
#     epsilon.bot <- CTmin - max(x[resp.max< 0.01 & x < To])
#   } else if (epsilon.mid > 0)
#   {
#     sigma.min <- sigma.mid
#     resp.min <- thermalResponse(x = x, To = To, rho = exp(rho.min), sigma = sigma.min)
#     epsilon.up <- CTmin - max(x[resp.min < 0.01 & x < To])
#   }
#   cat(paste(i, "\n"))
#   i <- i+1
# }
# 
# 
# if(min(abs(epsilon.up), abs(epsilon.bot)) == abs(epsilon.bot))
# {
#   sigma <- sigma.max
# } else if(min(abs(epsilon.up), abs(epsilon.bot)) == abs(epsilon.up))
# {
#   sigma <- sigma.min
# }
# 
# sigma
# 
# ############# SECTION RHO ####################
# 
# resp.min <- thermalResponse(x = x0, To = To, rho = exp(rho.min), sigma = sigma)
# resp.max <- thermalResponse(x = x0, To = To, rho = exp(rho.max), sigma = sigma)
# plot(resp.min ~ x0, type = "l")
# plot(resp.max ~ x0, type = "l")
# plot(resp.mid ~ x, type = "l")
# 
# 
# if(min(resp.min[x0 > To]) > 0.01)
# {
#   #change cutoff ?
#   # Investigate the issue between the used cutoff for the y axis (0.01)
#   # and the condition testing in the lines below
# }
# 
# if(CTmax < min(x[resp.max < 0.01 & x > To]))
# {
#   stop(paste("CTmax is too low, try to increase it (minimum possible value: ", 
#              round(min(x[resp.max < 0.01 & x > To]), 3), ")", sep = ""))
# } else if(CTmax > min(x[resp.min < 0.01 & x > To]))
# {
#   stop(paste("CTmax is too high, try to increase sigma for a higher CTmax (maximum possible value: ", 
#              round(min(x[resp.min < 0.01 & x > To]), 3), ")", sep = ""))
# }
# 
# epsilon.up <- CTmax - min(x[resp.min < 0.01 & x > To])
# epsilon.bot <- CTmax - min(x[resp.max < 0.01 & x > To])
# i <- 1
# while(min(abs(epsilon.up), abs(epsilon.bot)) > precision.CTmax)
# {
#   rho.mid <- (rho.min + rho.max) / 2
#   resp.mid <- thermalResponse(x = x, To = To, rho = exp(rho.mid), sigma = sigma)
#   epsilon.mid <- CTmax - min(x[resp.mid < 0.01 & x > To])
#   
#   if(epsilon.mid > 0)
#   {
#     rho.max <- rho.mid
#     resp.max <- thermalResponse(x = x, To = To, rho = exp(rho.max), sigma = sigma)
#     epsilon.bot <- CTmax - min(x[resp.max < 0.01 & x > To])
#   } else if (epsilon.mid < 0)
#   {
#     rho.min <- rho.mid
#     resp.min <- thermalResponse(x = x, To = To, rho = exp(rho.min), sigma = sigma)
#     epsilon.up <- CTmax - min(x[resp.min < 0.01 & x > To])
#   }
#   cat(paste(i, "\n"))
#   i <- i+1
# }
# 
# if(min(abs(epsilon.up), abs(epsilon.bot)) == abs(epsilon.bot))
# {
#   rho <- exp(rho.max)
# } else if(min(abs(epsilon.up), abs(epsilon.bot)) == abs(epsilon.up))
# {
#   rho <- exp(rho.min)
# }
# 
# rho
# 
# 
# 
# 
# 
# CTmax <- 7
# 
# resptot <- thermalResponse(x = x0, To = To, rho = rho, sigma = sigma)
# 
# plot(resptot~x0, type = "l")
# abline(v = CTmax)
# abline(v = CTmin)
# 
# min(x[resp < 0.01])
# 
# optimize(f = thermalResponse, interval = c(15, 100), To = 15, rho = 5, sigma = .01, tol = 1)
