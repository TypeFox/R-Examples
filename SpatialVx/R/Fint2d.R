Fint2d <- function(X, Ws, s, method = c("round", "bilinear", "bicubic"), derivs=FALSE, ...) {
    ##
    ## Function to extract value of Forcast at locations and apply bi-linear interpolation
    ## where necessary.
    ##
    ## 'X' matrix of forecast values.
    ## 'Ws' 'nm X 2' matrix giving the warped coordinates for the entire domain.
    ## 's' 'nm X 2' matrix giving the original forecast coordinates.
    ## 'method' character giving the interpolation method to use.  Default is to take the nearest value (round).
    ##    Alternative is to use bi-linear interpolation.
    ## 'derivs' logical, should the gradient components be calculated also?
    ##
    ## Value: numeric vector of length mn giving he deformed forecast field.
    ##
 
    method <- tolower(method)
    method <- match.arg(method)
 
    dimout <- dim(X)
  
    if( method=="round") {
  
          minx <- miny <- 1
          maxx <- dimout[1]
          maxy <- dimout[2]
  
    } else if( method=="bilinear") {
  
          minx <- miny <- 1
          maxx <- dimout[1] - 1
          maxy <- dimout[2] - 1
  
    } else if( method=="bicubic") {
  
          minx <- miny <- 2
          maxx <- dimout[1] - 2
          maxy <- dimout[2] - 2
  
    } # end of if else 'method' stmts.
  
    u <- Ws[,1]
    v <- Ws[,2]

    u[ u > maxx ] <- maxx
    u[ u < minx ] <- minx
    v[ v > maxy ] <- maxy
    v[ v < miny ] <- miny
  
    n <- length( u)
  
    out <- matrix(NA, nrow=dimout[1], ncol=dimout[2])
  
    if( derivs) out.x <- out.y <- out

    if(method == "round") {

       x <- floor(u + 0.5)
       x[ which(x > maxx) ] <- maxx
       x[ which(x < minx) ] <- minx

       y <- floor(v + 0.5)
       y[ which(y > maxy) ] <- maxy
       y[ which(y < miny) ] <- miny

       # for( i in 1:n) out[s[i,1],s[i,2]] <- X[x[i],y[i]]
       out[ s ] <- X[ cbind(x,y) ]

       if( derivs) {

          out.x <- cbind(out[, 2:dimout[2]] - out[, 1:(dimout[2] - 1)], 0)
          out.y <- rbind(out[2:dimout[1], ] - out[1:(dimout[1]-1), ], 0)

       } # end of if 'derivs' stmt.

    } else if(is.element(method, c("bilinear","bicubic"))) {

	fu <- floor( u+0.5)
	fu[ which(fu>maxx) ] <- maxx
	fu[ which(fu<minx) ] <- minx

	fv <- floor(v + 0.5)
	fv[ which(fv > maxy) ] <- maxy
	fv[ which(fv < miny) ] <- miny

	ufrac <- u - fu
	vfrac <- v - fv

	fuv <- cbind(fu, fv)
	fuv1 <- fuv
	fuv1[,2] <- fuv1[,2]+1
	fu1v <- fuv
	fu1v[,1] <- fu1v[,1]+1
	fu1v1 <- cbind( fu1v[,1], fuv1[,2])
 

	if(method == "bilinear") {

	    out[ s ] <- (1-ufrac)*(1-vfrac)*X[ fuv] + (1-ufrac)*vfrac*X[ fuv1] + ufrac*(1-vfrac)*X[ fu1v] + ufrac*vfrac*X[ fu1v1]

            if(derivs) {

              out.x[ s] <- (1-vfrac)*(X[ fu1v] - X[ fuv])  + vfrac*(X[ fu1v1] - X[ fuv1])
              out.y[ s] <- (1-ufrac)*(X[ fuv1] - X[ fuv]) + ufrac*(X[ fu1v1] - X[ fu1v])

            } # end of if 'derivs' stmt.

        } else if( method == "bicubic") {

	    u.bneg1 <- (2*ufrac^2 - ufrac^3 - ufrac)/2
	    v.bneg1 <- (2*vfrac^2 - vfrac^3 - vfrac)/2
	    u.b0    <- (3*ufrac^3 - 5*ufrac^2 + 2)/2
	    v.b0    <- (3*vfrac^3 - 5*vfrac^2 + 2)/2
	    u.b1    <- (4*ufrac^2 - 3*ufrac^3 + ufrac)/2
	    v.b1    <- (4*vfrac^2 - 3*vfrac^3 + vfrac)/2
	    u.b2    <- ((ufrac - 1)*ufrac^2)/2
	    v.b2    <- ((vfrac - 1)*vfrac^2)/2
	    
	    fun1vn1 <- fuv - 1
	    fun1v <- cbind( fun1vn1[,1], fuv[,2])
	    fun1v1 <- cbind( fun1v[,1], fuv[,2]+1)
	    fun1v2 <- cbind( fun1v[,1], fuv[,2]+2)
	    fuvn1 <- cbind( fuv[,1], fun1vn1[,2])
	    fuv1 <- cbind( fuv[,1], fun1v1[,2])
	    fuv2 <- cbind( fuv[,1], fun1v2[,2])
	    fu1vn1 <- cbind( fuv[,1]+1, fun1vn1[,2])
	    fu1v <- cbind( fu1vn1[,1], fuv[,2])
	    fu1v1 <- cbind( fu1vn1[,1], fuv1[,2])
	    fu1v2 <- cbind( fu1vn1[,1], fun1v2[,2])
	    fu2vn1 <- cbind( fuv[,1]+2, fun1vn1[,2])
	    fu2v <- cbind( fu2vn1[,1], fun1v[,2])
	    fu2v1 <- cbind( fu2vn1[,1], fun1v1[,2])
	    fu2v2 <- cbind( fu2vn1[,1], fun1v2[,2])

	    if( derivs) {

		du.bneg1 <- (4*ufrac - 3*ufrac^2 - 1)/2
		dv.bneg1 <- (4*vfrac - 3*vfrac^2 - 1)/2
		du.b0 <- (9*ufrac^2 - 10*ufrac)/2
		dv.b0 <- (9*vfrac^2 - 10*vfrac)/2
		du.b1 <- (8*ufrac - 9*ufrac^2 + 1)/2
		dv.b1 <- (8*vfrac - 9*vfrac^2 + 1)/2
		du.b2 <- (3*ufrac^2 - 2*ufrac)/2
		dv.b2 <- (3*vfrac^2 - 2*vfrac)/2

            } # end of if 'derivs' stmt.

	    out[ s ] <- u.bneg1*(v.bneg1*X[ fun1vn1] + v.b0*X[ fun1v] + v.b1*X[ fun1v1] + v.b2*X[ fun1v2]) +
                      u.b0*( v.bneg1*X[ fuvn1] + v.b0*X[ fuv] + v.b1*X[ fuv1] + v.b2*X[ fuv2]) +
                      u.b1*( v.bneg1*X[ fu1vn1] + v.b0*X[ fu1v] + v.b1*X[ fu1v1] + v.b2*X[ fu1v2]) +
                      u.b2*( v.bneg1*X[ fu2vn1] + v.b0*X[ fu2v] + v.b1*X[ fu2v1] + v.b2*X[ fu2v2])

	    if( derivs) {

          	out.x[ s ] <- du.bneg1*( v.bneg1*X[ fun1vn1] + v.b0*X[ fun1v] + v.b1*X[ fun1v1] + v.b2*X[ fun1v2]) +
                           du.b0*( v.bneg1*X[ fuvn1] + v.b0*X[ fuv] + v.b1*X[ fuv1] + v.b2*X[ fuv2]) +
                           du.b1*( v.bneg1*X[ fu1vn1] + v.b0*X[ fu1v] + v.b1*X[ fu1v1] + v.b2*X[ fu1v2]) +
                           du.b2*( v.bneg1*X[ fu2vn1] + v.b0*X[ fu2v] + v.b1*X[ fu2v1] + v.b2*X[ fu2v2])

          	out.y[ s ] <- dv.bneg1*( u.bneg1*X[ fun1vn1] + u.b0*X[ fuvn1] + u.b1*X[ fu1vn1] + u.b2*X[ fu2vn1]) +
                           dv.b0*( u.bneg1*X[ fun1v] + u.b0*X[ fuv] + u.b1*X[ fu1v] + u.b2*X[ fu2v]) +
                           dv.b1*( u.bneg1*X[ fun1v1] + u.b0*X[ fuv1] + u.b1*X[ fu1v1] + u.b2*X[ fu2v1]) +
                           dv.b2*( u.bneg1*X[ fun1v2] + u.b0*X[ fuv2] + u.b1*X[ fu1v2] + u.b2*X[ fu2v2])

       	    } # end of if 'derivs' stmt.

	} # end of if method is bicubic stmts.

    } else stop("method must be one of round, bilinear or bicubic")

      # image( out, col=c("grey", tim.colors(256)))
      out <- zapsmall(out)
      if(derivs) return(list( xy=out, dx=out.x, dy=out.y))
      else return(out)

 } # end of 'Fint2d' function.
  
