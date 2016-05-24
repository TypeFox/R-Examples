##########################################################################################################
#
# pwrRasch: Statistical Power Simulation for Testing the Rasch Model
#
# Internal function: Simulation 
#
# Authors: Takuya Yanagida <takuya.yanagida@univie.ac.at>
#  		     Jan Steinfeld <jan.steinfeld@univie.ac.at>
#
##########################################################################################################

pwr.rasch.internal <- function(b = b, ipar = ipar, ppar = ppar, 
                               runs = runs, H0 = H0, sig.level = sig.level, 
                               method = method, output = output) {
  
  #-----------------------------------------------------------------------------------------------------#
  # Numer of observations b and number of items c
  
  # b (number of observation per group)
  if (length(eval(parse(text = ppar[[1]]))) != length(eval(parse(text = ppar[[2]])))) {
    
    stop("Different numbers of persons specified in both groups")
    
  }
  
  # c (number of items)
  c <- unique(unlist(lapply(ipar, length)))
  
  #-----------------------------------------------------------------------------------------------------#
  # Simulation
  
  method <- ifelse(all(c("loop", "vectorized") %in% method), "vectorized", method)
  
  if (!method %in% c("loop", "vectorized")) {
    
    stop(paste0("Option \"", method, "\" for argument method unknown, use \"loop\" or \"vectorized\""))
    
  } 
  
  gc()
  
  ### for-loop
  if (method == "loop") {
    
    x <- switch(as.character(nchar(format(runs, scientific = FALSE))), 
                "1" = 9, "2" = 7, "3" = 5, "4" = 3, "5" = 1, "6" = 1, "7" = 1, "8" = 1, "9" = 1, "10" = 1)
    
    group <- c(rep(0, times = b), rep(1, times = b))
    
    # H0 condition
    if (H0 == TRUE) {
      
      cat(paste0("  Conducting H0 simulation...    ", t1 <- Sys.time(), "\n"))
      
      H0.AC.p <- NULL 
      
      for (i in 1:runs) {
        
        cat("\r", paste0(" run ", formatC(i, digits = nchar(runs) - 1, format = "d", flag = 0)," of ",runs,"..."))
        
        # Simulate 0/1 data
        dat <- rbind(simul.rasch(eval(parse(text = ppar[[1]])), ipar[[1]]), 
                     simul.rasch(eval(parse(text = ppar[[2]])), ipar[[1]]))
        
        # Three-way analysis of variance with mixed classification
        H0.AC.p <- c(H0.AC.p, aov.rasch.sim(data = reshape.rasch(dat, group = group)))    
        
      }
      
      cat(paste0(" finished", paste(rep(" ", times = x), collapse = ""), Sys.time(),"\n"))

    }
    
    # H1 condition
    H1.AC.p <- NULL 
    
    cat(paste0("  Conducting H1 simulation...    ", if (H0 == FALSE) t1 <- Sys.time() else Sys.time(), "\n"))   
    
    for (i in 1:runs) {
      
      cat("\r", paste0(" run ", formatC(i, digits = nchar(runs) - 1, format = "d", flag = 0)," of ", runs,"..."))
      
      # Simulate 0/1 data
      dat <- rbind(simul.rasch(eval(parse(text = ppar[[1]])), ipar[[1]]),
                   simul.rasch(eval(parse(text = ppar[[2]])), ipar[[2]]))
      
      # Three-way analysis of variance with mixed classification
      H1.AC.p <- c(H1.AC.p, aov.rasch.sim(data = reshape.rasch(dat, group = group))) 
      
    }
    
    cat(paste0(" finished", paste(rep(" ", times = x), collapse = ""), t2 <- Sys.time(),"\n"))                            

  }
  
  ### Vectorized computing  
  if (method == "vectorized") {  
    
    group <- c(rep(0, times = b), rep(1, times = b))
    
    # H0 condition
    if (H0 == TRUE) {
      
      cat(paste0("  Conducting H0 simulation...    ", t1 <- Sys.time(), "\n"))          
      
      dat.mx <- rbind(simul.rasch(eval(parse(text = ppar[[1]])), ipar[[1]]), 
                      simul.rasch(eval(parse(text = ppar[[2]])), ipar[[1]]))
      dat.mx <- reshape.rasch(dat.mx, group = group)
      
      dat.mx1 <- simul.rasch(eval(parse(text = sub("b", "(runs-1)*b", ppar[[1]]))), ipar[[1]]) 
      dat.mx2 <- simul.rasch(eval(parse(text = sub("b", "(runs-1)*b", ppar[[2]]))), ipar[[1]])
      
      ###
      
      dat.mx1 <- array(dat.mx1, dim = c(b, runs - 1, ncol(dat.mx1)))
      dat.mx1 <- aperm(dat.mx1, c(1, 3, 2))
      dat.mx1 <- matrix(dat.mx1, nrow = b*ncol(dat.mx1))
      
      dat.mx2 <- array(dat.mx2, dim = c(b, runs - 1, ncol(dat.mx2)))
      dat.mx2 <- aperm(dat.mx2, c(1, 3, 2))
      dat.mx2 <- matrix(dat.mx2, nrow = b*ncol(dat.mx2))
      
      dat.mx <- cbind(dat.mx, rbind(dat.mx1, dat.mx2))
      
      ###
      
      H0.AC.p <- unname(aov.rasch.sim.vec(data = dat.mx))
      
    }
    
    # H1 condition    
    cat(paste0("  Conducting H1 simulation...    ", if (H0 == FALSE) t1 <- Sys.time() else Sys.time(), "\n"))   
    
    dat.mx <- rbind(simul.rasch(eval(parse(text = ppar[[1]])), ipar[[1]]),
                    simul.rasch(eval(parse(text = ppar[[2]])), ipar[[2]]))
    dat.mx <- reshape.rasch(dat.mx, group = group)
    
    dat.mx1 <- simul.rasch(eval(parse(text = sub("b", "(runs-1)*b", ppar[[1]]))), ipar[[1]]) 
    dat.mx2 <- simul.rasch(eval(parse(text = sub("b", "(runs-1)*b", ppar[[2]]))), ipar[[2]])
    
    ###
    
    dat.mx1 <- array(dat.mx1, dim = c(b, runs - 1, ncol(dat.mx1)))
    dat.mx1 <- aperm(dat.mx1, c(1, 3, 2))
    dat.mx1 <- matrix(dat.mx1, nrow = b*ncol(dat.mx1))
    
    dat.mx2 <- array(dat.mx2, dim = c(b, runs - 1, ncol(dat.mx2)))
    dat.mx2 <- aperm(dat.mx2, c(1, 3, 2))
    dat.mx2 <- matrix(dat.mx2, nrow = b*ncol(dat.mx2))
    
    dat.mx <- cbind(dat.mx, rbind(dat.mx1, dat.mx2))
    
    ###
    
    H1.AC.p <- unname(aov.rasch.sim.vec(data = dat.mx))
    
  }  
  
  if (method != "loop") cat(paste0("  Finished simulation...         ", t2 <- Sys.time(), "\n"))  
  
  cat("--------------------------------------------------------\n")
  
  print(round(t2 - t1, digits = 2))  
  
  if (H0 == TRUE) {
  
    return(list(b = b, ipar = ipar, c = c, ppar = ppar, runs = runs, sig.level = sig.level,
                H0.AC.p = H0.AC.p, H1.AC.p = H1.AC.p,
                power = sum(H1.AC.p < sig.level) / runs,
                type1 = sum(H0.AC.p < sig.level) / runs))
    
  } else {
    
    return(list(b = b, ipar = ipar, c = c, ppar = ppar, runs = runs, sig.level = sig.level,
                H1.AC.p = H1.AC.p,
                power = sum(H1.AC.p < sig.level) / runs))
    
  }

}