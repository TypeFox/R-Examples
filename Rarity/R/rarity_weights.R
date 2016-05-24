######################################
######   Weight calculation    #######
######################################

setGeneric("rWeights", 
           def = function(occData, 
                          Qmax = max(occData), 
                          Qmin = min(occData), 
                          wMethods = "W", # Possibilities: W, invQ, oldW (attention)
                          rCutoff = "Gaston", # A value, "Leroy" or "Gaston"
                          normalised = T, # Should weights be normalised?
                          assemblages, # Assemblages for meanCutoff calculation
                          extended = F, # Should one-scale weights be provided in the weight table?
                          rounding = 3) # F = no rounding, an integer = the number of digits for rounding
                  {
                    standardGeneric( "rWeights" )
                  }
           )

setMethod("rWeights",
          signature(occData = "vector"),
          function(occData, Qmax = max(occData), Qmin = min(occData), wMethods = "W", rCutoff, normalised = T, assemblages, rounding = 3)
          {   
            if(any(is.na(Qmax)))
            {
              Qmax = max(occData, na.rm = T)
            }
            if(any(is.na(Qmin)))
            {
              Qmin = min(occData, na.rm = T)
            }
            
            if (is.numeric(rCutoff))
            {
              if(rCutoff > 1 | round(rCutoff * Qmax, 10) < Qmin)
              {
                stop("Cutoff value must be between (Qmin/Qmax) and 1")
              }
              
              if(round(rCutoff * Qmax, 10) == Qmin)
              {
                rCutoff <- rCutoff + .000000001
                warning("The cutoff cannot be equal to Qmin for mathematical reasons. Setting a cutoff slightly above Qmin (Qmin + .000000001).")
              }
              cat("Rarity cut-off point:", rCutoff, "/", rCutoff * Qmax, "\n")
            } else if (rCutoff == "Leroy")
            {
              if(!is.matrix(assemblages))
              {
                if(!is.data.frame(assemblages))
                {
                  stop("Please supply a correct assemblage matrix or data.frame for the Leroy rarity cut-off selection; rownames = species names, colnames = assemblage labels")
                } 
                assemblages <- as.matrix(assemblages)
              }
              if(any(apply(assemblages, 2, sum) < 5))
              {
                if(any(apply(assemblages, 2, sum) < 1))
                {
                  stop("Some assemblages do not have any species, please remove them")
                }
                cat("Warning: some assemblages have very low species richnesses, which might bias cutoff calculations\n")
              }
              if(nrow(assemblages) != length(occData))
              {
                cat("Remark: Number of species different between assemblages and occData, calculating cutoff using species provided in assemblages\n")
                if(any(!(rownames(assemblages) %in% names(occData))))
                {
                  stop("Row names (sp. names) of assemblages do not correspond to names (sp. names) of occData")
                }
                occData2 <- occData[match(rownames(assemblages), names(occData))]
              } else
              {
                occData2 <- occData
              }
              assemblages <- occData2 * assemblages[, 1:ncol(assemblages)]
              assemblages[which(assemblages == 0)] <- NA
              rCutoff <- mean(apply(assemblages, 2, quantile, .25, na.rm=T)) / Qmax
              if(round(rCutoff * Qmax, 10) == Qmin)
              {
                rCutoff <- rCutoff + .000000001
                warning("The cutoff cannot be equal to Qmin for mathematical reasons. Setting a cutoff slightly above Qmin (Qmin + .000000001).")
              }
              cat("Rarity cut-off point:", rCutoff, "/", rCutoff * Qmax, "\n")
            } else if (rCutoff == "Gaston")
            {
              rCutoff <- quantile(occData, .25, na.rm = T)/Qmax
              if(round(rCutoff * Qmax, 10) == Qmin)
              {
                rCutoff <- rCutoff + .000000001
                warning("The cutoff cannot be equal to Qmin for mathematical reasons. Setting a cutoff slightly above Qmin (Qmin + .000000001).")
              }
              cat("Rarity cut-off point:", rCutoff, "/", rCutoff * Qmax, "\n")
            } else stop("Please specify a correct cutoff (either a value or a correct method)")
              
              

            
            spWeights <- cbind(Q = occData,
                               R = occData <= rCutoff * Qmax)
            
            if ("W" %in% wMethods)
            {
              spWeights <- cbind(spWeights,
                                 W = exp( -(((occData - Qmin)/(rCutoff * Qmax - Qmin)) * 0.97 + 1.05)^2))
              if(normalised == T)
              {
                spWeights[, "W"] <- spWeights[, "W"] / exp( -(((Qmin - Qmin)/(rCutoff * Qmax - Qmin)) * 0.97 + 1.05)^2)
              }
            }
            if ("invQ" %in% wMethods)
            {
              spWeights <- cbind(spWeights,
                                 invQ = 1 / occData)
              if(normalised == T)
              {
                spWeights[, "invQ"] <- spWeights[, "invQ"] / (1 / Qmin)
              }
            }
            if ("oldW" %in% wMethods)
            {
              n <- .nFind(Qmin = Qmin, Qmax = Qmax, r = rCutoff)
              spWeights <- cbind(spWeights,
                                 oldW = exp( -(((occData / Qmax) * n + 1)^2) ))
              if(normalised == T)
              {
                spWeights[, "oldW"] <- spWeights[, "oldW"] / exp( -(((Qmin / Qmax) * n + 1)^2) )
              }
            }
            if(any(c("oldW", "W") %in% wMethods))
            {
              spWeights <- cbind(spWeights,
                                 cut.off = rCutoff)
            }
            if(rounding)
            {
              if(!is.numeric(rounding))
              {
                stop("Either provide an integer value (if you want rounding) or FALSE (if you do not want rounding) for argument 'rounding'")
              }
              rounding <- as.integer(rounding)
              spWeights[, which(!(colnames(spWeights) %in% c("Q", "R", "cut.off")))] <- round(spWeights[, which(!(colnames(spWeights) %in% c("Q", "R", "cut.off")))], digits = rounding)                                                                
            }
            spWeights <- as.data.frame(spWeights)
            spWeights
          }
          )

setMethod("rWeights",
          signature(occData = "matrix"),
          function(occData, Qmax = apply(occData, 2, max), Qmin = apply(occData, 2, min), wMethods = "W", rCutoff, normalised = T, assemblages, extended = F, rounding = 3)
          {
            if(any(is.na(Qmax)))
            {
              Qmax = apply(occData, 2, max, na.rm = T)
            }
            if(any(is.na(Qmin)))
            {
              Qmin = apply(occData, 2, min, na.rm = T)
            }
            
            if(is.numeric(rCutoff))
            {
              if(length(rCutoff) != ncol(occData))
              {
                stop("Specify a vector of cut-offs of same length than the number of occurrence datasets")
              }
              if(any(rCutoff > 1)|any(round(rCutoff * Qmax, 10) < Qmin))
              {
                stop("Cutoff values must be between (Qmin/Qmax) and 1")
              }
              if(any(round(rCutoff * Qmax, 10) == Qmin))
              {
                rCutoff[which(round(rCutoff * Qmax, 10) == Qmin)] <- rCutoff[which(round(rCutoff * Qmax, 10) == Qmin)]  + .000000001
                warning("The cutoff cannot be equal to Qmin for mathematical reasons. Setting a cutoff slightly above Qmin (Qmin + .000000001).")
              }
              cat("Rarity cut-off points:\n", paste(colnames(occData), rCutoff, "/", rCutoff * Qmax, collapse = "\n ", sep = " "), "\n")
            } else if (rCutoff == "Leroy")
            {
              if(!is.matrix(assemblages))
              {
                if(!is.data.frame(assemblages))
                {
                  stop("Please supply a correct assemblage matrix or data.frame for the Leroy rarity cut-off selection; rownames = species names, colnames = assemblage labels")
                }
                assemblages <- as.matrix(assemblages)
              }
              if(any(apply(assemblages, 2, sum) < 5))
              {
                if(any(apply(assemblages, 2, sum) < 1))
                {
                  stop("Some assemblages do not have any species, please remove them")
                }
                cat("Warning: some assemblages have very low species richnesses, which might bias cutoff calculations\n")
              }
              if(nrow(assemblages) != nrow(occData))
              {
                cat("Remark: Number of species different between assemblages and occData, calculating cutoff using species provided in assemblages\n")
                if(any(!(rownames(assemblages) %in% rownames(occData))))
                {
                  stop("Row names (sp. names) of assemblages do not correspond to row names (sp. names) of occData")
                }
                occData2 <- occData[match(rownames(assemblages), rownames(occData)), ]
              } else
              {
                occData2 <- occData
              }
              assemblages2 <- array(dim = c(ncol(occData),
                                            dim(assemblages)))
              
              for(x1 in 1:ncol(occData2))
              {
                assemblages2[x1, , ] <- occData2[, x1] * assemblages[, 1:ncol(assemblages)]
              }
              assemblages2[which(assemblages2 == 0)] <- NA
              
              rCutoff <- apply(assemblages2, 1, function(x) mean(apply(x, 2, quantile, .25, na.rm=T))) / Qmax
              if(any(round(rCutoff * Qmax, 10) == Qmin))
              {
                rCutoff[which(round(rCutoff * Qmax, 10) == Qmin)] <- rCutoff[which(round(rCutoff * Qmax, 10) == Qmin)]  + .000000001
                warning("The cutoff cannot be equal to Qmin for mathematical reasons. Setting a cutoff slightly above Qmin (Qmin + .000000001).")
              }
              cat("Rarity cut-off points:\n", paste(colnames(occData), rCutoff, "/", rCutoff * Qmax, collapse = "\n ", sep = " "), "\n")
              
            } else if (rCutoff == "Gaston")
            {
              rCutoff <- apply(occData, 2, quantile, .25, na.rm = T) / Qmax
              if(any(round(rCutoff * Qmax, 10) == Qmin))
              {
                rCutoff[which(round(rCutoff * Qmax, 10) == Qmin)] <- rCutoff[which(round(rCutoff * Qmax, 10) == Qmin)]  + .000000001
                warning("The cutoff cannot be equal to Qmin for mathematical reasons. Setting a cutoff slightly above Qmin (Qmin + .000000001).")
              }
              cat("Rarity cut-off points:\n", paste(colnames(occData), rCutoff, "/", rCutoff * Qmax, collapse = "\n ", sep = " "), "\n")
            } else stop("Please specify a correct cutoff (either a value or a correct method)")

            
            spWeights <- occData
            colnames(spWeights) <- paste("Q", 1:ncol(occData), sep = "")
            for(x1 in 1:ncol(occData))
            {
              spWeights <- cbind(spWeights,
                                 occData[, x1] <= rCutoff[x1] * Qmax[x1])
              colnames(spWeights)[ncol(spWeights)] <- paste("R", x1, sep = "")
            }
            if ("W" %in% wMethods)
            {
              for(x1 in 1:ncol(occData))
              {
                spWeights <- cbind(spWeights,
                                   exp( -(((occData[, x1] - Qmin[x1])/(rCutoff[x1] * Qmax[x1] - Qmin[x1])) * 0.97 + 1.05)^2))

                if(normalised)
                {
                  spWeights[, ncol(spWeights)] <- spWeights[, ncol(spWeights)] / exp( -(((Qmin[x1] - Qmin[x1])/(rCutoff[x1] * Qmax[x1] - Qmin[x1])) * 0.97 + 1.05)^2)
                }
                colnames(spWeights)[ncol(spWeights)] <- paste("W", x1, sep = "")
              }
              spWeights <- cbind(spWeights,
                                 W = rowSums(spWeights[, (ncol(spWeights) - ncol(occData) + 1):ncol(spWeights)]))
            }
            if ("invQ" %in% wMethods)
            {
              for(x1 in 1:ncol(occData))
              {
                spWeights <- cbind(spWeights,
                                   invQ = 1 / occData[, x1])
                if(normalised == T)
                {
                  spWeights[, ncol(spWeights)] <- spWeights[, ncol(spWeights)] / (1 / Qmin[x1])
                }
                colnames(spWeights)[ncol(spWeights)] <- paste("invQ", x1, sep = "")
              }
              spWeights <- cbind(spWeights,
                                 invQ = rowSums(spWeights[, (ncol(spWeights) - ncol(occData) + 1):ncol(spWeights)]))
            }
            if ("oldW" %in% wMethods)
            {
              for(x1 in 1:ncol(occData))
              {
                n <- .nFind(Qmin = Qmin[x1], Qmax = Qmax[x1], r = rCutoff[x1])
                spWeights <- cbind(spWeights,
                                   oldW = exp( -(((occData[, x1] / Qmax[x1]) * n + 1)^2) ))
                if(normalised == T)
                {
                  spWeights[, ncol(spWeights)] <- spWeights[, ncol(spWeights)] / exp( -(((Qmin[x1] / Qmax[x1])) * n + 1)^2)
                }
                colnames(spWeights)[ncol(spWeights)] <- paste("oldW", x1, sep = "")
              }
              spWeights <- cbind(spWeights,
                                 oldW = rowSums(spWeights[, (ncol(spWeights) - ncol(occData) + 1):ncol(spWeights)]))
            }

            if(!extended)
            {
              spWeights <- cbind(spWeights[, 1:(ncol(occData)*2)],
                                 spWeights[, wMethods])
              if(length(wMethods) == 1)
              {
                colnames(spWeights)[ncol(spWeights)] <- wMethods
              }
            }
            if(any(c("oldW", "W") %in% wMethods))
            {
              for(x1 in 1:ncol(occData))
              {
                spWeights <- cbind(spWeights,
                                   rCutoff[x1])
                colnames(spWeights)[ncol(spWeights)] <- paste("cut.off", x1, sep = "")
              }
            }
            if(rounding)
            {
              if(!is.numeric(rounding))
              {
                stop("Either provide an integer value (if you want rounding) or FALSE (if you do not want rounding) for argument 'rounding'")
              }
              rounding <- as.integer(rounding)
              spWeights[, which(!(colnames(spWeights) %in% c("Q", "R", "cut.off", 
                                                             paste("Q", 1:1000, sep = ""), 
                                                             paste("R", 1:1000, sep = ""),  
                                                             paste("cut.off", 1:1000, sep = ""))))] <- 
                round(spWeights[, which(!(colnames(spWeights) %in% c("Q", "R", "cut.off", 
                                                                     paste("Q", 1:1000, sep = ""), 
                                                                     paste("R", 1:1000, sep = ""),  
                                                                     paste("cut.off", 1:1000, sep = ""))))], digits = rounding)
            }
            spWeights <- as.data.frame(spWeights)
            spWeights
          }
          )

setMethod("rWeights",
          signature(occData = "data.frame"),
          function(occData, Qmax = apply(occData, 2, max), Qmin = apply(occData, 2, min), wMethods = "W", rCutoff, normalised = T, assemblages, extended = F)
          {
            occData <- as.matrix(occData)
            rWeights(occData = occData, Qmax = Qmax, Qmin = Qmin, wMethods = wMethods, rCutoff = rCutoff, normalised = normalised, assemblages = assemblages, extended = extended)
          }
          )


######################################
######    Old n calculation    #######
######################################
# nFind <- function (Qmin, Qmax, r)
# {
#   epsil = 0.001
#   n0 = 1 / r
#   Qt = Qmax * r
#   # First and second derivatives 
#   dW <- function(n)
#   {
#     dW=(exp(-((n*Qmin)/Qmax+1)^2)-exp(-(n+1)^2))/(Qmax-1)-(2*n*((n*Qt)/Qmax+1)*exp(-((n*Qt)/Qmax+1)^2))/Qmax
#     dW
#   }  
#   ddW<-function(n)
#   {
#     ddW=(4*n*Qt*((n*Qt)/Qmax+1)^2*exp(-((n*Qt)/Qmax+1)^2))/Qmax^2-(2*((n*Qt)/Qmax+1)*exp(-((n*Qt)/Qmax+1)^2))/Qmax-(2*n*Qt*exp(-((n*Qt)/Qmax+1)^2))/Qmax^2+(2*(n+1)*exp(-(n+1)^2)-(2*Qmin*((n*Qmin)/Qmax+1)*exp(-((n*Qmin)/Qmax+1)^2))/Qmax)/(Qmax-1)
#     ddW
#   }
#   # Newton's method for approximating n
#   n=n0	
#   delta=1
#   while (abs(delta)>epsil) 
#   {	
#     delta=dW(n)/ddW(n)
#     n= (n-delta)
#   }
#   n	
#   #n0 = round(n0, 2)
# }

.nFind <- function (Qmin, Qmax, r)
{
  epsil = 0.001
  n0 = 1 / r
  Qt = Qmax * r
  # First and second derivatives 
  dW <- function(n)
  {
    dW=(exp(-((n*Qmin)/Qmax+1)^2)-exp(-(n+1)^2))/(Qmax-Qmin)-(2*n*((n*Qt)/Qmax+1)*exp(-((n*Qt)/Qmax+1)^2))/Qmax
    dW
  }  
  ddW<-function(n)
  {
    ddW=(4*n*Qt*((n*Qt)/Qmax+1)^2*exp(-((n*Qt)/Qmax+1)^2))/Qmax^2 - 
      (2*((n*Qt)/Qmax+1)*exp(-((n*Qt)/Qmax+1)^2))/Qmax -
      (2*n*Qt*exp(-((n*Qt)/Qmax+1)^2))/Qmax^2 +
      (2*(n+1)*exp(-(n+1)^2)-(2*Qmin*((n*Qmin)/Qmax+1)*exp(-((n*Qmin)/Qmax+1)^2))/Qmax)/(Qmax-Qmin)
    ddW
  }
  # Newton's method for approximating n
  n=n0  
  delta=1
  while (abs(delta)>epsil) 
  {	
    delta=dW(n)/ddW(n)
    n= (n-delta)
  }
  n	
  #n0 = round(n0, 2)
}

######################################
######    Simple W function    #######
######################################
.W <- function(Q, Qmin = min(Q), Qmax = max(Q), r)
{
  exp(-(((Q-Qmin)/(r*Qmax - Qmin)) * 0.97 + 1.05)^2)
}


