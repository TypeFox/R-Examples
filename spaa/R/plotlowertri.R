plotlowertri <-
function(input, valuename = "r", pchlist = c(19, 17, 15, 1, 5, 2, 7),
         interval = 6, cex = 1, ncex = 1, int =1.2, add.number = TRUE, 
         size = FALSE, add.text = FALSE, show.legend = TRUE, digits = 2){
   if (!is.matrix(input)){ 
       input = as.matrix(input)
           if (any(!is.numeric(input))){
		       stop("Non numeric values found in the input matrix.\n")
			   }
      }
   if (!nrow(input)==ncol(input)){
       stop("The input data must be a square matrix.\n")
	   }
   if (!(interval+1)==length(pchlist)){
       stop("Pchlist must be equal to interval + 1.\n")
	   }
 
   if (nrow(input) < 10 ){
       cat(
       "Warning: Too few rows, please adjust the legend using \"cex\" and \"int\".\n", 
       "Some of the points may be missing due to the distribution of data.\n",  
       "You may also have to choose an appropriate number of intervals \n", 
       "using \"interval\" and \"pchlist\".\n")
      }
   
    if(nrow(input) > 40){
        cat("Warning: Too many rows, please adjust the legend using cex and int.\n")
    }
    if (is.null(pchlist)){
	    pchlist = c(19, 17, 15, 1, 5, 2, 7) 
	}
    else{ 
	    pchlist = pchlist 
	}
    interval = interval
    d <- as.dist(input)
    seqn  <- (seq(0, interval^2, by = interval))/(interval^2)
    limit0 <- round((min(d) + (max(d)-min(d))*seqn), digits)
    limit <- limit0[c(-1, -length(limit0))]
    limit <- sort(limit, decreasing = TRUE)

    matrix <- lower.tri(input)
    plot(x = 0:nrow(matrix),y = 0:nrow(matrix), type = "n", axes = FALSE, xlab="", ylab="") 
    xleft = 0
    ybottom = 0
    for (i in 2:nrow(matrix))
      for (j in 1:length(matrix[i,][matrix[i,]])){
           xleft <-  j - 1
           ybottom <- nrow(matrix) - i
           rect(xleft = xleft, ybottom = ybottom, xright = xleft + 1, ytop = ybottom + 1) 
           x <- xleft   + 0.5
           y <- ybottom + 0.5
           if (input[i,j] >= limit[1]){ 
			  pch = pchlist[1]
			  }
           if (input[i,j] < limit[length(limit)]){ 
			  pch = pchlist[(length(limit)+1)]
			  }
           else{
               for (n in 2:(length(limit))){
                    if ( input[i,j] >= limit[n]&input[i,j] < limit[n-1] ){ 
					   pch = pchlist[n] 
					   }
                }
            }
           if(size){ 
		   points(x,y, pch = pch, cex = cex*abs(input[i,j])) 
		   }
           else { points(x,y, pch = pch, cex = cex) }
           if(add.text) { 
             text(x,y+0.25, pch)
             text(x,y-0.25, round(input[i,j], digits))
           }
       }
      if (add.number){
            x1 = - 0.5
            for (textloc in 1:nrow(matrix)){
                x1 <- x1 + 1
                y1 <- nrow(matrix) - textloc + 0.5
                text(x1, y1, as.character(textloc), cex = ncex)
            }
       }
      if (show.legend){
    xlegend <- ncol(input)
        for (k in 0:(length(limit)))  { 
             ylegend <- nrow(input) - k * int
             points( ncol(input)*0.6, ylegend, pch = pchlist[k+1] )
             if (k == 0){
              text(xlegend*0.68, ylegend, paste(formatC(sprintf("%.2f",limit[k + 1]), width = 5)))
              text(xlegend*0.76, ylegend, expression(""<=""))
              text(xlegend*0.8,  ylegend, valuename)
             }
             if (! k == 0 & !k == length(limit)){ 
              text(xlegend*0.88, ylegend, paste("<", formatC(sprintf("%.2f", limit[k]), width = 5)))
              text(xlegend*0.8,  ylegend, valuename)
              text(xlegend*0.68, ylegend, paste(formatC(sprintf("%.2f", limit[k+1]), width = 5)))
              text(xlegend*0.76, ylegend, expression(""<=""))
              }
             if (k == length(limit)){
              text(xlegend*0.88, ylegend , paste("<", formatC(sprintf("%.2f",limit[k]), width = 5)))
              text(xlegend*0.8, ylegend, valuename)
              }
          }
       }

}

