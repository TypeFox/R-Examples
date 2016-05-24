# function to plot estimation results  2015-3-12 Thong Pham
plot.PAFit <-
function(x,data,true_f = NULL, plot = c("A","f","true_f"), plot_bin = TRUE,
         line = TRUE, confidence = TRUE, high_deg = NULL, shade = 0.5,
         ...) {
  if (plot_bin == TRUE) {
      x$k <- x$center_k
      x$A <- x$theta
      x$upper_A <- x$upper_bin
      x$lower_A <- x$lower_bin
  } 
  if ("A" == plot[1]) {
      non_zero <- which(x$A > 10^-20 & x$k >= x$deg_threshold)
      limit <- c(min(x$lower_A[non_zero]),max(x$upper_A[non_zero]))
      xlim  <- c(min(x$k[non_zero] + 1),max(x$k[non_zero] + 1))
      plot(x$k[non_zero][1] + 1,x$A[non_zero][1],xlab = "k + 1",
           ylab = "A_k",xlim = xlim, ylim = limit,log = "xy",...)
      if (TRUE == confidence) {
          arrows(x0 = x$k + 1, y0 = x$lower_A, x1 = x$k + 1, y1 = x$upper_A, code = 3,angle = 90, length = 0,col = rgb(0,0,0,shade))
      }
          points(x$k[non_zero] + 1,x$A[non_zero],pch = 1,...)
      if (TRUE == line) {
          alpha <- x$alpha
          beta <- x$linear_fit$coefficients[1]
          lines(x$k[non_zero] + 1,exp(beta)*(x$k[non_zero] + 1)^alpha,lwd= 2)
      }
  }
  else if ("f" == plot[1]) {
      if (FALSE == is.null(high_deg))
          non_zero <- x$lower_f > 10^-20 & data$increase >= high_deg
      else
          non_zero <- x$lower_f > 10^-20 & data$increase > 0
      if (length(non_zero) <= 0)
        stop("There is no data. Please decrease high_deg") 
      if (TRUE == confidence)
          lim_y = c(min(x$lower_f[non_zero]), 
                    max(x$upper_f[non_zero]))
      else lim_y = c(min(x$f[non_zero]),max(x$f[non_zero]))
      xlim <- c(min(data$increase[non_zero])+1,max(data$increase[non_zero]))
      plot(data$increase[non_zero][1],x$f[non_zero][1],log="xy",ylab = "Estimated fitness",xlim = xlim, 
           ylim = lim_y, xlab = "Number of edges acquired",...)
      if (TRUE == confidence)
          arrows(x0 = data$increase[non_zero], y0 = x$lower_f[non_zero], x1 = data$increase[non_zero], 
                 y1 = x$upper_f[non_zero], code = 3,angle = 90, length = 0,col = rgb(0,0,0,shade))
      points(data$increase[non_zero],x$f[non_zero],pch = 1,col = rgb(0,0,0,shade),...)
      abline(h = 1)
  }
  else if ("true_f" == plot[1]) {
          true_f1   <- length(true_f[data$node_id])*true_f[data$node_id]/sum(true_f[data$node_id])
          if (FALSE == is.null(high_deg)) {
              non_zero <- x$lower_f > 10^-20 & true_f1 > 10^-20 & data$increase > high_deg
          } else
              non_zero <- x$lower_f > 10^-20 & true_f1 > 10^-20 
          if (length(non_zero) <= 0)
             stop("There is no data. Please decrease high_deg")  
          b        <- lm(true_f1[non_zero] ~ 0 + x$f[non_zero], weights = data$increase[non_zero])$coefficients[1]
          upper_f <-  exp(log(b*x$f[non_zero]) + 2 * sqrt(x$var_f[non_zero] / x$f[non_zero] ^ 2))
          lower_f <-  exp(log(b*x$f[non_zero]) - 2 * sqrt(x$var_f[non_zero] / x$f[non_zero] ^ 2))
          xlim <- c(min(c(lower_f,true_f1[non_zero])), max(c(upper_f,true_f1)))
          ylim <- c(min(c(lower_f,true_f1[non_zero])), max(c(upper_f,true_f1)))
          plot(b*x$f[non_zero][1],true_f1[non_zero][1], xlim= xlim, ylim = ylim,
                xlab = "Estimated fitness", ylab= "True fitness",log = "xy", pch ="",...)
          if (TRUE == confidence)
              arrows(x0 = lower_f, y0 = true_f1[non_zero], x1 = upper_f, y1 = true_f1[non_zero], code = 3,angle = 90, length = 0,
                     col = rgb(0,0,0,shade))
          points(b*x$f[non_zero],true_f1[non_zero],pch ="",...)
          abline(a=0,b = 1)
          text(b*x$f[non_zero],true_f1[non_zero],data$increase[non_zero],col = rgb(0,0,0,0.5),...) 
      }
}
