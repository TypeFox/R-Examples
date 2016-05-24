
plot.gSlc <-
function(x, pages = 0,responseScale = FALSE, xlab = NULL, ylab = NULL, main = NULL, bty = NULL,...)
{
   gSlcObj <- x

#Get the information from gSlcObj
   grid.size <- 101 
   formula.infor <- gSlcObj$formulaInfor
   data <- gSlcObj$scaleData
   sigsqu <- gSlcObj$sigmaSquared
   nu <- gSlcObj$nu
   num.knots <- gSlcObj$nBasis
   random.factor <- gSlcObj$randomFactor
   family <- gSlcObj$family
   X.min <- gSlcObj$Xmin
   X.range <- gSlcObj$Xrange

   num.smooth.funcs <- length(formula.infor$svar)
   num.linear.funcs <- length(formula.infor$lvar)
   num.vars <- num.smooth.funcs + num.linear.funcs 

   linkfunction <- function(x,family) {
     if (family == "poisson") { return(exp(x))}
     if (family == "binomial") {return(exp(x)/(1+exp(x))) }
   }

   if (num.smooth.funcs < 1) {
   options(show.error.messages = FALSE)
   cat("\n \n There are no smooth functions in the model to plot.\n")
   cat(" The summary() function can be used to visualise\n")
   cat(" output for parametric components.\n\n\n")
   stop("")
   }
   nc.X <- num.vars + 1     
   nr.nu <- nrow(nu)
   nc.nu <- ncol(nu)
 
   upper.perc <- function(x,cred.level)
       return(quantile(x,(1+cred.level)/2))

   lower.perc <- function(x,cred.level)
       return(quantile(x,(1-cred.level)/2))
  
# Do slice plots of first function estimate and truth.
   
   X.grid.01 <- matrix(0,grid.size,nc.X)
   for (i in 1:nc.X) {
   X.grid.01[,i] <- seq(0,1,length=grid.size)
   }  

   X.grid.mean <- matrix(0,grid.size, nc.X)
   row.X.grid.mean <- c(1,apply(as.matrix(data[,formula.infor$varlist]),2,mean))
   for(j in 1:grid.size) {
   X.grid.mean[j,] <- row.X.grid.mean
   }
    
   # Find the matrix sqrt of Omega and knots
   # matrix.sqrt.Omega <- array(dim = c(num.knots,0))

   Zspline.grid.01 <- array(dim = c(grid.size,0))
   Zspline.grid.mean <- array(dim = c(grid.size,0))
   
   if (num.smooth.funcs > 0) {
     
     for (i in 1 : num.smooth.funcs) {
        
       knots.i <- seq(0,1,length=(num.knots[i]+2))[-c(1,(num.knots[i]+2))]
       knots.i <- quantile(unique(data[,formula.infor$svar[i]]),knots.i)

       svd.Omega.i <- svd(abs(outer(knots.i,knots.i,"-"))^3)
       matrix.sqrt.Omega.i <- t(svd.Omega.i$v%*%(t(svd.Omega.i$u)*sqrt(svd.Omega.i$d)))

       Zspline.i.grid.01 <- t(solve(matrix.sqrt.Omega.i,
                    t(abs(outer(X.grid.01[,i+1+num.linear.funcs],knots.i,"-")^3))))
       Zspline.i.grid.mean <- t(solve(matrix.sqrt.Omega.i,
                    t(abs(outer(X.grid.mean[,i+1+num.linear.funcs],knots.i,"-")^3))))
       Zspline.grid.01 <- cbind(Zspline.grid.01,Zspline.i.grid.01)

       Zspline.grid.mean <- cbind(Zspline.grid.mean,Zspline.i.grid.mean)
     }
    
     up.num.knots.i <- 0
     low.num.knots.i <- 0
     
       
     num.plots <- num.smooth.funcs
            

     if (pages < 0) pages <- 0
     if (pages > num.plots) pages <- num.plots
       
     if (pages !=0) {  
           ppp <- (num.plots - 1) %/% pages + 1 
           num.row <- round(sqrt(ppp))
           num.col <- (ppp-1) %/% num.row + 1  
           
           par(mfrow = c(num.row, num.col))
     }
       
     if (pages == 0) {
           if (prod(par("mfrow")) < 1) {ppp <- 1}
           else  {ppp <- prod(par("mfrow")) }
     }
       
     for (i in 1: num.smooth.funcs) {
        
         X.grid <- X.grid.mean
         X.grid[,1+ num.linear.funcs +i] <- X.grid.01[,1+ num.linear.funcs +i] 
   
         Z.grid <- Zspline.grid.mean

         up.num.knots.i <- up.num.knots.i + num.knots[i]
         low.num.knots.i <- up.num.knots.i - num.knots[i] + 1

         Z.grid[,low.num.knots.i:up.num.knots.i] <- Zspline.grid.01[,low.num.knots.i:up.num.knots.i]
   
         C.grid <- cbind(X.grid,Z.grid)
         etahat.grid <- t(C.grid%*%nu)
         if (responseScale) {
             etahat.grid <- linkfunction(etahat.grid,family)
         }
         mean.etahat.grid <- apply(etahat.grid,2,mean)

         lower.etahat.grid <- apply(etahat.grid,2,lower.perc,0.95)
         upper.etahat.grid <- apply(etahat.grid,2,upper.perc,0.95)
         ylim.val <- range(c(lower.etahat.grid,upper.etahat.grid)) 
         par(mai = c(1.02,0.82,0.42,0.42))

         X.grid.01[,1+num.linear.funcs+i] <- X.grid.01[,1+num.linear.funcs+i] * X.range[num.linear.funcs + i] + X.min[num.linear.funcs + i]
         
         if (is.null(xlab[i]) || is.na(xlab[i]) || (xlab[i] == "")) {xlabi <- formula.infor$svar[i]
            } else{xlabi <- xlab[i]}

         if (is.null(ylab[i]) || is.na(ylab[i])) {ylabi <- ""
            } else {ylabi <- ylab[i]}

         if (is.null(main[i]) || is.na(main[i])) {maini <- ""
            } else{maini <- main[i]}

         btylist <- c("o","l","7","c", "u", "]") 
         if (is.null(bty[i]) || is.na(bty[i]) || (!any(btylist == bty[i])) ) {btyi <- "o"
            } else {btyi <- bty[i]}

         plot(X.grid.01[,1+num.linear.funcs+i],mean.etahat.grid,type="n",ylim=ylim.val ,xlab = xlabi, ylab = ylabi, main = maini, bty = btyi,...)
         polygon(c(X.grid.01[,1+num.linear.funcs+i],rev(X.grid.01[,1+num.linear.funcs + i])),c(lower.etahat.grid,
              rev(upper.etahat.grid)),col="palegreen",border=F)
  
         lines(X.grid.01[,1+num.linear.funcs+i],mean.etahat.grid,col="forestgreen",lwd=3)
         
         if (i%%ppp ==0) {par(ask = TRUE)}
     }
    par(ask = FALSE)
  }
  
}
