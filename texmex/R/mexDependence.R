`mexDependence` <-
function (x, which, dqu, margins = "laplace", constrain=TRUE, v = 10, maxit=1000000, start=c(.01, .01), marTransform="mixture", nOptim = 1,
          PlotLikDo=FALSE, PlotLikRange=list(a=c(-1,1),b=c(-3,1)), PlotLikTitle=NULL)
{
   theCall <- match.call()
   if (class(x) != "migpd")
       stop("you need to use an object created by migpd")

   x$margins <-  casefold(margins)
   x <- mexTransform(x, margins = casefold(margins),method = marTransform)

   if (margins == "gumbel" & constrain){
     warning("With Gumbel margins, you can't constrain, setting constrain=FALSE")
     constrain <- FALSE
   }

   if (missing(which)) {
       warning("Missing 'which'. Conditioning on", dimnames(x$transformed)[[2]][1], "\n")
       which <- 1
   }
   else if (length(which) > 1)
       stop("which must be of length 1")
   else if (is.character(which))
       which <- match(which, dimnames(x$transformed)[[2]])

   if (missing(dqu)) {
       warning("Assuming same quantile for dependence thesholding as was used\n     to fit corresponding marginal model...\n")
       dqu <- x$mqu[which]
   }
   dth <- quantile(x$transformed[, which], dqu)

   dependent <- (1:(dim(x$data)[[2]]))[-which]
   if (length(dqu) < length(dependent))
       dqu <- rep(dqu, length = length(dependent))

   # Allowable range of 'a' depends on marginal distributions
   aLow <- ifelse(x$margins == "gumbel", 10^(-10),-1 + 10^(-10))

   if (missing(start)){
     start <- c(.01, .01)
   } else if(class(start) == "mex"){
     start <- start$dependence$coefficients[1:2,]
   }

   if( length(start) == 2 ){
     start <- matrix(rep(start,length(dependent)),nrow=2)
   }

   if( length(start) != 2*length(dependent)){
     stop("start should be of type 'mex' or be a vector of length 2, or be a matrix with 2 rows and ncol equal to the number of dependence models to be estimated")
   }

   if( ! missing(PlotLikRange) ){
     PlotLikDo <- TRUE
   }

   qfun <- function(X, yex, wh, aLow, margins, constrain, v, maxit, start){
     Qpos <- function(param, yex, ydep, constrain, v, aLow) {

  	   a <- param[1]
       b <- param[2]

       res <- PosGumb.Laplace.negProfileLogLik(yex, ydep, a, b, constrain, v, aLow) # defined in file mexDependenceLowLevelFunctions
       res$profLik
     } # Close Qpos <- function

     o <- try(optim(par=start, fn=Qpos,
              control=list(maxit=maxit),
              yex = yex[wh], ydep = X[wh], constrain=constrain, v=v, aLow=aLow), silent=TRUE)

     if (class(o) == "try-error"){
        warning("Error in optim call from mexDependence")
        o <- as.list(o)
        o$par <- rep(NA, 6)
        o$value <- NA
     } else if (o$convergence != 0) {
        warning("Non-convergence in mexDependence")
        o <- as.list(o)
        o$par <- rep(NA, 6)

     } else if(nOptim > 1) {

        for( i in 2:nOptim ){
           o <- try(optim(par=o$par, fn=Qpos,
                    control=list(maxit=maxit),
                    yex = yex[wh], ydep = X[wh], constrain=constrain, v=v, aLow=aLow), silent=TRUE)
           if (class(o) == "try-error"){
             warning("Error in optim call from mexDependence")
             o <- as.list(o)
             o$par <- rep(NA, 6)
             o$value <- NA
             break()
           } else if (o$convergence != 0) {
             warning("Non-convergence in mexDependence")
             o <- as.list(o)
             o$par <- rep(NA, 6)
             break()
           }
        }
     }

     if ( PlotLikDo ){# plot profile likelihood for (a,b)
       nGridPlotLik <- 50
       a.grid <- seq(PlotLikRange$a[1],PlotLikRange$a[2],length=nGridPlotLik)
       b.grid <- seq(PlotLikRange$b[1],PlotLikRange$b[2],length=nGridPlotLik)
       NegProfLik <- matrix(0,nrow=nGridPlotLik,ncol=nGridPlotLik)
       for(i in 1:nGridPlotLik){
         for(j in 1:nGridPlotLik){
           NegProfLik[i,j] <- PosGumb.Laplace.negProfileLogLik(yex=yex[wh], ydep=X[wh],
                                  a = a.grid[i],b=b.grid[j], constrain=constrain,v=v,aLow=aLow)$profLik
         }
       }
       NegProfLik[NegProfLik > 10^10] <- NA
       if(sum(!is.na(NegProfLik))){
          filled.contour(a.grid,b.grid,-NegProfLik,main=paste("Profile likelihood",PlotLikTitle),color.palette = terrain.colors,
                         xlab="a",ylab="b",plot.axes={ axis(1); axis(2); points(o$par[1],o$par[2]) })
       }
     }

     if (!is.na(o$par[1])) { # gumbel margins and negative dependence
        if (margins == "gumbel" & o$par[1] <= 10^(-5) & o$par[2] < 0) {
           lo <- c(10^(-10), -Inf, -Inf, 10^(-10), -Inf, 10^(-10))
           Qneg <- function(yex, ydep, param) {
               param <- param[-1]
               b <- param[1]
               cee <- param[2]
               d <- param[3]
               m <- param[4]
               s <- param[5]

               obj <- function(yex, ydep, b, cee, d, m, s) {
                      mu <- cee - d * log(yex) + m * yex^b
                      sig <- s * yex^b
                      log(sig) + 0.5 * ((ydep - mu)/sig)^2
                      }
               res <- sum(obj(yex, ydep, b, cee, d, m, s))
               res
          }
          o <- try(optim(c(0, 0, 0, 0, 0, 1), Qneg, method = "L-BFGS-B", lower=lo,
                   upper=c(1, 1-10^(-10), Inf, 1-10^(-10), Inf, Inf),
                   yex = yex[wh], ydep = X[wh]), silent=TRUE)

          if (class(o) == "try-error" || o$convergence != 0) {
             warning("Non-convergence in mexDependence")
             o <- as.list(o)
             o$par <- rep(NA, 6)
          }
        } else { # end if gumbel margins and neg dependence
          Z <- (X[wh] - yex[wh] * o$par[1]) / (yex[wh]^o$par[2])
          o$par <- c(o$par[1:2], 0, 0, mean(Z),sd(Z))
        }
    }
    c(o$par[1:6], o$value) # Parameters and negative loglik
   } # Close qfun <- function(

   yex <- c(x$transformed[, which])
   wh <- yex > unique(dth)

   res <- sapply(1:length(dependent),
                 function(X,dat,yex,wh,aLow,margins,constrain,v,maxit,start)qfun(dat[,X],yex,wh,aLow,margins,constrain,v,maxit,start[,X]),
                 dat=as.matrix(x$transformed[, dependent]), yex=yex, wh=wh, aLow=aLow, margins=margins, 
                 constrain=constrain, v=v, maxit=maxit, start=start)

   loglik <- -res[7,]
   res <- matrix(res[1:6,], nrow=6)

   dimnames(res)[[1]] <- c(letters[1:4],"m","s")
   dimnames(res)[[2]] <- dimnames(x$transformed)[[2]][dependent]
   gdata <- as.matrix(x$transformed[wh, -which])
   tfun <- function(i, data, yex, a, b, cee, d) {
       data <- data[, i]
       a <- a[i]
       b <- b[i]
       cee <- cee[i]
       d <- d[i]
       if (is.na(a))
           rep(NA, length(data))
       else {
           if (a < 10^(-5) & b < 0)
               a <- cee - d * log(yex)
           else a <- a * yex
           (data - a)/(yex^b)
       }
   }
   z <- try(sapply(1:(dim(gdata)[[2]]), tfun, data = gdata,
       yex = yex[wh], a = res[1, ], b = res[2, ], cee = res[3, ], d = res[4, ]))
   if (class(z) %in% c("Error", "try-error")) {
       z <- matrix(nrow = 0, ncol = dim(x$data)[[2]] - 1)
   }
   else if (is.R()) {
       if (!is.array(z)) {
           z <- matrix(nrow = 0, ncol = dim(x$data)[[2]] - 1)
       }
   }
   dimnames(z) <- list(NULL,dimnames(x$transformed)[[2]][dependent])
   res2 <- list(coefficients = res, Z = z, dth = unique(dth),
               dqu = unique(dqu), which = which, conditioningVariable= colnames(x$data)[which],
	             loglik=loglik, margins=margins, constrain=constrain, v=v)
   oldClass(res2) <- "mexDependence"

   output <- list(margins=x, dependence=res2, call=theCall)
   oldClass(output) <- "mex"
   output
}

