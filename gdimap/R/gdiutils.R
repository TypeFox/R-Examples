##
## show gfas as 3D surface
##

gfasurf3d <-
function(im, zfactor=5, ...)
{
	z <- im*zfactor
	zlim <- range(z)
	zlen <- zlim[2] - zlim[1] + 1
	# colorlut <- colorspace::terrain.colors(zlen) 
	# colorlut <- colorspace::heat_hcl(zlen, h = c(0, -100), c = c(40, 80), l = c(75, 40), power = 1)
	colorlut <- colorspace::terrain_hcl(zlen, h = c(0, -100), c. = c(40, 80), l = c(75, 40), power = 1)
	col <- colorlut[ floor(z-zlim[1]+1) ] # assign colors to heights for each point
	d <- dim(im)
	surface3d(1:d[1], 1:d[2], z,  color=col, back="lines", add=TRUE, ...)
}

#----------------------------------------------
# norm01 <- function(x) { return ((x - min(x))/diff(range(x))) }
norm01 <- function(x) {
	if(diff(range(x)) != 0) { return ((x - min(x))/diff(range(x))) }
	else { return (x - min(x)) }
}

#----------------------------------------------
maxone <- function(x) { return(which(x == max(x))[1]) }

#----------------------------------------------
fnorm <- function(x) { sqrt(sum(x^2)) }

#----------------------------------------------
anisofn <- function(x, aniso=NULL) {
  if(!is.null(aniso)) {
    stopifnot(aniso >= 0 & aniso < 1)
    mx <- max(x)  
    mn <- min(x)  
    x <- pmax(x-mn-aniso*(mx-mn),0)
  }
  x <- norm01(x)
}
#----------------------------------------------
# GQI-standard
gqifn <-
function(odfvert, btable, lambda=NULL)
{
  if(is.null(lambda)) lambda <- 1.24
	lvalues <- sqrt(btable[,1]*0.01506)
	bvector <- btable[,2:4]*repmat(lvalues,1,3)
	q0 <- (odfvert %*% t(bvector) * lambda) / pi 
	invisible(gsl::sinc(q0))
}

# GQI2
gqifn2 <-
function(odfvert, btable, lambda=NULL)
{
	sqradial <- function(x, tol=0.01) {
    result <- (2 * x * cos(x) + (x * x - 2) * sin(x)) / (x ** 3)
    x.nearzero = (x < tol) & (x > -tol)
    ifelse(x.nearzero, 1./3, result)
	}
  if(is.null(lambda)) lambda <- 3
	lvalues <- sqrt(btable[,1]*0.01506)
	bvector <- btable[,2:4]*repmat(lvalues,1,3)
	q0 <- (odfvert %*% t(bvector) * lambda) / pi 
	invisible(sqradial(q0))
}

#----------------------------------------------
genfa <- 
function(odf)
{
	n <- length(odf)
	m1 <- sum(odf)
	m2 <- sum(odf^2)
  m1 <- m1^2;
  m1 <- m1 / n
  if (m2 == 0.0)
  	gfa <- 0
	else
	  gfa <- sqrt(n / (n-1)*(m2-m1)/m2)
}

#----------------------------------------------
repmat <-
function (a, m, n)
{
	nr <- nrow(a); nc <- ncol(a)
	b <- NULL; b1 <- NULL;
	for(i in seq(1,n)) { b1 <- cbind(b1,a) }
	for(i in seq(1,m)) { b <- rbind(b,b1) }
	dimnames(b) <- NULL
	invisible(b)
}

#----------------------------------------------
#
# Find ODF peaks 
#
findpeak <-
function(odf, odf.vertices, odf.faces)
{
  # cat("dim of odf.vertices:", dim(odf.vertices),"\n")
  # cat("dim of odf.faces:", dim(odf.faces),"\n")
  ispeak <- odf;
  odf.faces <- odf.faces - (odf.faces > length(odf)) * length(odf);
  ispeak[
    odf.faces[1,
       odf[odf.faces[2,]] >= odf[odf.faces[1,]] |
      odf[odf.faces[3,]] >= odf[odf.faces[1,]]
     ]
  ] <- 0
  ispeak[
      odf.faces[2,
        odf[odf.faces[1,]] >= odf[odf.faces[2,]] |
        odf[odf.faces[3,]] >= odf[odf.faces[2,]]
      ]
  ] <- 0
  ispeak[
      odf.faces[3,
        odf[odf.faces[2,]] >= odf[odf.faces[3,]] |
        odf[odf.faces[1,]] >= odf[odf.faces[3,]]
      ]
  ] <- 0;
  ##
  ox <- sort(-ispeak, index.return=TRUE );
  values <- ox$x
  ordering <- ox$ix
  p <- ordering[values < 0] # indices
  pcoords <- odf.vertices[,p] # coords of peak vertices
  pcoords <- matrix(pcoords, nrow=3, ncol=length(p))
  invisible(list(peaks=p, pcoords=pcoords, np=length(p)))
}

