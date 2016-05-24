############
#Generic function splits on sensitivity parameter type
############
calc.range = function(sensParam, grid.dim, zetaz.range, zetay.range, buffer, U.model, zero.loc, Xcoef.plot, Y, Z, X, Y.res, Z.res, v_Y, v_Z, theta, sgnTau0, control.fit, null.trt) {
	if(sensParam == "coef") 
		result = calc.range.coef(grid.dim, zetaz.range, zetay.range, buffer, U.model, zero.loc, Xcoef.plot, Y, Z, X, Y.res, Z.res, v_Y, v_Z, theta, sgnTau0, control.fit, null.trt)

	if(sensParam == "cor") 
		result = calc.range.cor(grid.dim, zetaz.range, zetay.range, U.model, zero.loc, Xcoef.plot, Y, Z, X, Y.res, Z.res, v_Y, v_Z, theta, sgnTau0, control.fit)

	return(result)
}


##############
#function to calculate vector of sensitivity paramters 
#using coefficients as SPs
##############
calc.range.coef = function(grid.dim, zetaz.range, zetay.range, buffer, U.model, zero.loc, Xcoef.plot, Y, Z, X, Y.res, Z.res, v_Y, v_Z, theta, sgnTau0, control.fit,null.trt) {
	extreme.coef = matrix(c(-sqrt((v_Y-buffer)/(1-buffer)), -sqrt(v_Z-buffer), sqrt((v_Y-buffer)/(1-buffer)), sqrt(v_Z-buffer)), nrow = 2) 
	  if(U.model == "binomial" & !is.binary(Z)){ 
	    extreme.coef = matrix(c(-sqrt(4*v_Y-buffer), -sqrt(v_Z/(theta*(1-theta))-buffer), sqrt(4*v_Y-buffer), sqrt(v_Z/(theta*(1-theta))-buffer)), nrow = 2) 
	  }
	  if(U.model == "binomial" & is.binary(Z)){ 
	    lp.quant = quantile(null.trt$linear.predictors, 0.25)
	    zetaz.min = max(-(2-lp.quant), -3) 
	    zetaz.max = min((2-lp.quant), 3)   
	    extreme.coef = matrix(c(-sqrt(4*v_Y-buffer), zetaz.min, sqrt(4*v_Y-buffer), zetaz.max), nrow = 2)
	  }


  if(!is.null(zetay.range) & !is.null(zetaz.range)){ #custom grid range.
    if(zetay.range[1] < extreme.coef[1,1] | zetay.range[2] > extreme.coef[1,2]){
      zetay.range[1] = max(zetay.range[1], extreme.coef[1,1])
      zetay.range[2] = min(zetay.range[2], extreme.coef[1,2])
      warning("Sensitivity parameter range for Y inconsistent with possible values given data.  Range restricted.")
    }
    if(zetaz.range[1] < extreme.coef[2,1] | zetaz.range[2] > extreme.coef[2,2]){
      zetaz.range[1] = max(zetaz.range[1], extreme.coef[2,1])
      zetaz.range[2] = min(zetaz.range[2], extreme.coef[2,2])
      warning("Sensitivity parameter range for Z inconsistent with possible values given data.  Range restricted.")
    }
    
    if(sign(zetaz.range[1])==sign(zetaz.range[2])|any(zetaz.range ==0)){ #one quadrant.
	      #define the vector of sens.parm for treatment
	      zetaZ <- seq(zetaz.range[1], zetaz.range[2], length.out = grid.dim[1]) 
	    }else{ #two quadrants.
	      
	        #change z-dimension to odd number
	        grid.dim[1]=ifelse(grid.dim[1]%%2==1,grid.dim[1],grid.dim[1]+1)
	        #create temporary seq to find border.
	        zetaZ <- seq(zetaz.range[1], zetaz.range[2], length.out = grid.dim[1])
	        #number of cells left and right of vertical axis.
	        dim.left = which.min(abs(zetaZ[which(zetaZ<0)]))
	        dim.right = grid.dim[1] - dim.left - 1
	        #define the vector of sens.parm for treatment
	        zetaZ <- c(seq(zetaz.range[1], zetaz.range[1]/(dim.left*3), length.out=dim.left), 0,
	                   seq(zetaz.range[2]/(dim.right*3), zetaz.range[2], length.out=dim.right))
	      
	      }
	    
    
	    if(sign(zetay.range[1])==sign(zetay.range[2])|any(zetay.range == 0)){ #one quadrant.
	      #define vector of sens.parm for treatment
	      zetaY <- seq(zetay.range[1], zetay.range[2], length.out = grid.dim[2])

	    }else{#two quadrants in vertical direction.
	      
	        #change y-dimension to odd number
	        grid.dim[2]=ifelse(grid.dim[2]%%2==1,grid.dim[2],grid.dim[2]+1)
	        #create temporary seq to find border.
	        zetaY <- seq(zetay.range[1], zetay.range[2], length.out = grid.dim[2])
	        #number of cells below and above horizontal axis.
	        dim.down = which.min(abs(zetaY[which(zetaY<0)]))
	        dim.up = grid.dim[2] - dim.down - 1
	        #define the vector of sens.parm for response
	        zetaY <- c(seq(zetay.range[1], zetay.range[1]/(dim.down*3), length.out=dim.down), 0,
	                   seq(zetay.range[2]/(dim.up*3), zetay.range[2], length.out=dim.up))        
	      
	    }
    
	  }else if(zero.loc == "full"){
	    #change z-dimension to odd number
	    grid.dim[1]=ifelse(grid.dim[1]%%2==1,grid.dim[1],grid.dim[1]+1)
	    #number of cells left and right of vertical axis.
	    dim.left = dim.right = (grid.dim[1]-1)/2
	    #define the vectors of sens.parms
	    zetaZ <- c(seq(extreme.coef[2,1]*.95, extreme.coef[2,1]*.95/(dim.left*3),	length.out=dim.left), 0,
	               seq(extreme.coef[2,2]*.95/(dim.right*3), extreme.coef[2,2]*.95, length.out=dim.right))
	    zetaY <- seq(0, extreme.coef[1,2]*.95, length.out = grid.dim[2])
	  }else{
	    #find ranges for final grid
	    cat("Finding grid range...\n")
	    grid.range = grid.search(extreme.coef, zero.loc, Xcoef.plot, Y, Z, X=X, 
	                             Y.res, Z.res, v_Y, v_Z, theta, sgnTau0, 
	                             control.fit = control.fit, sensParam = "coef")
    
	    zetaY <- seq(grid.range[1,1], grid.range[1,2], length.out = grid.dim[2])
	    zetaZ <- seq(grid.range[2,1], grid.range[2,2], length.out = grid.dim[1])
      
      zetaY[which.min(abs(zetaY))] = 0
      zetaZ[which.min(abs(zetaZ))] = 0
	    
	  }
	return(list(zetaZ = zetaZ, zetaY = zetaY))
}



##############
#function to calculate vector of sensitivity paramters 
#using correlations as SPs
##############
calc.range.cor = function(grid.dim, zetaz.range, zetay.range, U.model, zero.loc, Xcoef.plot, Y, Z, X, 
	                             Y.res, Z.res, v_Y, v_Z, theta, sgnTau0, 
	                             control.fit) {
  
	extreme.cors = maxCor(Y.res, Z.res)
	if(U.model == "binomial") {
 		extreme.cors = 2*dnorm(0)*extreme.cors
	}

	  if(!is.null(zetay.range) & !is.null(zetaz.range)){ #custom grid range.
	    if(zetay.range[1] < -extreme.cors[1,1] | zetay.range[2] > extreme.cors[1,2]){
	      zetay.range[1] = max(zetay.range[1], -extreme.cors[1,1])
	      zetay.range[2] = min(zetay.range[2], extreme.cors[1,2])
	      warning("Sensitivity parameter range for Y inconsistent with possible values given data.  Range restricted.")
	    }
	    if(zetaz.range[1] < extreme.cors[2,1] | zetaz.range[2] > extreme.cors[2,2]){
	      zetaz.range[1] = max(zetaz.range[1], extreme.cors[2,1])
	      zetaz.range[2] = min(zetaz.range[2], extreme.cors[2,2])
	      warning("Sensitivity parameter range for Z inconsistent with possible values given data.  Range restricted.")
	    }
	    
	    if(sign(zetaz.range[1])==sign(zetaz.range[2])|any(zetaz.range == 0)){ #one quadrant.
	      #define the vector of sens.parm for treatment
	      zetaZ <- seq(zetaz.range[1], zetaz.range[2], length.out = grid.dim[1])
 
	    }else{ #two quadrants.
	      
	        #change z-dimension to odd number
	        grid.dim[1]=ifelse(grid.dim[1]%%2==1,grid.dim[1],grid.dim[1]+1)
	        #create temporary seq to find border.
	        zetaZ <- seq(zetaz.range[1], zetaz.range[2], length.out = grid.dim[1])
	        #number of cells left and right of vertical axis.
	        dim.left = which.min(abs(zetaZ[which(zetaZ<0)]))
	        dim.right = grid.dim[1] - dim.left - 1
	        #define the vector of sens.parm for treatment
	        zetaZ <- c(seq(zetaz.range[1], zetaz.range[1]/(dim.left*3), length.out=dim.left), 0,
	                   seq(zetaz.range[2]/(dim.right*3), zetaz.range[2], length.out=dim.right))
	      
	      
	    }
    
	    if(sign(zetay.range[1])==sign(zetay.range[2])|any(zetay.range==0)){ #one quadrant.
	      #define vector of sens.parm for treatment
	      zetaY <- seq(zetay.range[1], zetay.range[2], length.out = grid.dim[2])

	    }else{#two quadrants in vertical direction.
	      
	        #change y-dimension to odd number
	        grid.dim[2]=ifelse(grid.dim[2]%%2==1,grid.dim[2],grid.dim[2]+1)
          #create temporary seq to find border.
	        zetaY <- seq(zetay.range[1], zetay.range[2], length.out = grid.dim[2])
	        #number of cells below and above horizontal axis.
	        dim.down = which.min(abs(zetaY[which(zetaY<0)]))
	        dim.up = grid.dim[2] - dim.down - 1
	        #define the vector of sens.parm for response
	        zetaY <- c(seq(zetay.range[1], zetay.range[1]/(dim.down*3), length.out=dim.down), 0,
	                   seq(zetay.range[2]/(dim.up*3), zetay.range[2], length.out=dim.up))        
	      
	    }
    
	  }else if(zero.loc == "full"){
	    #change z-dimension to odd number
	    grid.dim[1]=ifelse(grid.dim[1]%%2==1,grid.dim[1],grid.dim[1]+1)
	    #number of cells left and right of vertical axis.
	    dim.left = dim.right = (grid.dim[1]-1)/2
	    #define the vectors of sens.parms
	    zetaZ <- c(seq(extreme.cors[2,1]*.95, extreme.cors[2,1]*.95/(dim.left*3),	length.out=dim.left),0,
	               seq(extreme.cors[2,2]*.95/(dim.right*3), extreme.cors[2,2]*.95, length.out=dim.right))
	    zetaY <- seq(0, extreme.cors[1,2]*.95, length.out = grid.dim[2])
	  }else{
	    #find ranges for final grid
	    cat("Finding grid range...\n")
	    grid.range = grid.search(extreme.cors, zero.loc, Xcoef.plot, Y, Z, X, 
	                             Y.res, Z.res, v_Y, v_Z, theta, sgnTau0, 
	                             control.fit = control.fit, sensParam = "cor")
    
	    zetaY <- seq(grid.range[1,1], grid.range[1,2], length.out = grid.dim[2])
	    zetaZ <- seq(grid.range[2,1], grid.range[2,2], length.out = grid.dim[1])
	    
      zetaY[which.min(abs(zetaY))] = 0
	    zetaZ[which.min(abs(zetaZ))] = 0
	    
	  }
	return(list(zetaZ = zetaZ, zetaY = zetaY))
}
