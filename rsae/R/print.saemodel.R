print.saemodel <-
function(x, ...){
   # distinguish synthetic (potentially contaminated) from ordinary data
   if (is.null(attr(x, "contam"))){
      cat("SAE MODEL TYPE: B (J.N.K. Rao's classification)\n")
      cat("--- \n")

      cat(paste("FIXED EFFECTS: ", attr(x, "yname"), " ~ ", paste(attr(x, "xnames"), collapse=" + "), sep=""), "\n")
      cat(paste("AREA-SPECIFIC RANDOM EFFECTS: ", attr(x, "areadef"), sep=""), "\n")

   }else{
      skeleton <- attr(x, "contam")$skeleton
      cat("SAE MODEL TYPE: B (J.N.K. Rao's classification)\n")
      cat("DATA: Synthetic, simulated data \n")
      cat("MODEL: \n")
      beta <- skeleton$beta
      intercept <- skeleton$intercept
      p <- length(beta)
      if (p > 1){
	 if (is.null(intercept)){
	    cat(paste("   y_ij = sum_k[ beta_k * x_kij] + v_i + e_ij", sep=""), "\n")
	    cat(" with \n")
	 }else{
	    cat(paste("   y_ij = intercept + sum_k[ beta_k * x_kij] + v_i + e_ij", sep=""), "\n")
	    cat(" with \n")
	 }
	 cat(paste("   each x_kij ~ N(0, 1), k=1,...,",p , sep=""), "\n")
      }else{
	 if (is.null(intercept)){
	    cat(paste("   y_ij = beta * x_ij + v_i + e_ij", sep=""), "\n")
	    cat(" with \n")
	 }else{
	    cat(paste("   y_ij = intercept + beta * x_ij + v_i + e_ij", sep=""), "\n")
	    cat(" with \n")
	 }
	 cat(paste("   x_ij ~ N(0, 1)", sep=""), "\n")

      }
      # random effects
       if(skeleton$vu.epsilon==0){
 	 cat("   v_i ~  N(0, ", skeleton$vu,  ")\n", sep="") 
       }
       else{
 	 cat("   v_i ~ (", 1-skeleton$vu.epsilon, ")*N(0, ", skeleton$vu, ") + ", skeleton$vu.epsilon, "*N(0, ", skeleton$vu.contam, ") \n", sep="") 
       }
       if(skeleton$ve.epsilon==0){
 	 cat("   e_ij ~ N(0, ", skeleton$ve,  ")\n", sep="")
       }else{
 	 cat("   e_ij ~ (", 1-skeleton$ve.epsilon, ")*N(0, ", skeleton$ve, ") + ", skeleton$ve.epsilon, "*N(0, ", skeleton$ve.contam, ") \n", sep="")
       }
   }   
 }

