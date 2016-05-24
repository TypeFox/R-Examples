
## Authors 
## Martin Schlather, schlather@math.uni-mannheim.de
##
##
## Copyright (C) 2015 Martin Schlather
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 3
## of the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  




RMcauchytbm <- function(alpha, beta, gamma, var, scale, Aniso, proj) {
  return(RMtbm(fulldim=gamma,
               RMgencauchy(alpha, beta, var, scale, Aniso, proj)))
}
  
RMcardinalsine <- function(var, scale, Aniso, proj) {
  return(RMwave(var, scale, Aniso, proj))
}
  
RMgneitingdiff <- function(nu, taper.scale, scale, var, Aniso, proj){
  return(RMmult(RMgengneiting(kappa=3, mu=1.5, scale=taper.scale) *
                RMwhittle(nu=nu, scale=scale),
                var=var, Aniso=Aniso, proj=proj))
}

RMparswmX <- function(nudiag, rho, var, scale, Aniso, proj) {
  return(RMschur(M=rho, RMparswm(nudiag, var, scale, Aniso, proj)))
}

RMpoweredexp <- function(alpha, var, scale, Aniso, proj) {
  return(RMstable(alpha, var, scale, Aniso, proj))
}

RMtent <- function(var, scale, Aniso, proj) {
  return(RMaskey(alpha=1.0, var, scale, Aniso, proj))
}

RMwendland <- function(kappa, mu, var, scale, Aniso, proj) {
  return(RMgengneiting(kappa, mu, var, scale, Aniso, proj))
}


RMcovariate <- function(c, x, y=NULL, z=NULL, T=NULL, grid,
                        var, scale, Aniso, proj, raw, norm, addNA) {
  Call <- iRMcovariate
  if (missing(x) && length(T)==0) {
    if (length(y)!=0 || length(T)!=0 || !missing(grid))
      stop("y, z, T, grid may only be given if 'x' is given")
    Call(norm=norm, c=c, scale=scale, Aniso=Aniso, proj=proj, var=var, raw=raw,
         addNA=addNA)
  } else {
    new <- C_CheckXT(x, y, z, T, grid, printlevel=0)
    Call(norm=norm, c=c, x=new, scale=scale, Aniso=Aniso, proj=proj, var=var,
         raw=raw, addNA=addNA)
  }
}
  
 
RMfixcov <- function(M, x, y=NULL, z=NULL, T=NULL, grid,
                     var, scale, Aniso, proj, raw, norm) {
  Call <- iRMfixcov
  if (missing(x) && length(T)==0) {
    if (length(y)!=0 || length(T)!=0 || !missing(grid))
      stop("y, z, T, grid may only be given if 'x' is given")
    Call(norm=norm, M=M, scale=scale, Aniso=Aniso, proj=proj, var=var, raw=raw)
  } else {
    new <- C_CheckXT(x, y, z, T, grid, printlevel=0)
    Call(norm=norm, M=M, x=new, scale=scale, Aniso=Aniso, proj=proj, var=var,
         raw=raw)
  }
}
 

R.lon <- function() R.p(1, "spherical system")
R.lat <- function() R.p(2, "spherical system")


RMchoquet <- function(b) stop("not implemented yet")

RMpolynome <- function(degree, dim, value=NA, varnames = c("x", "y", "z", "T"),
                       proj=1:4) {
  if (degree < 0 || degree > 5) stop("the degree is out of range")
  if (dim < 0  || dim > 4) stop("the dimension is out of range")
  x <- as.matrix(do.call("expand.grid", rep(list(0:degree), dim)))
  sums <- rowSums(x)
  y <- NULL
  for (i in 0:degree) {
    idx <- sums == i
    y <- rbind(y, x[idx, ])
  }
  n <- nrow(y)
#  y <- as.vector(y)
  z <- paste(rep(paste(" ", varnames[1:dim], sep=""), each=n),
             ifelse(y>1, "^", ""),
             ifelse(y>1, y, ""), sep="")
  z[ y == 0] <- ""
  dim(z) <- dim(y)
  z <- apply(z, 1, paste, collapse="", sep="")
  m <- length(z) - length(value)
  if (m > 0) value <- c(value, rep(NA, m)) else
  if (m < 0) value <- value[1:length(z)]
  cat( paste( value, z, collapse = " + ", sep=""), "\n" )
  
  z <- paste(rep(paste("R.p(", proj[1:dim], ")", sep=""), each=n),
             ifelse(y>1, "^", ""),
             ifelse(y>1, y, ""), sep="")
  z[ y == 0 ] <- ""
  if (length(z) > 100) stop("maximum is ", MAXSUB, "^2 terms")
  dim(z) <- dim(y)
  z <- apply(z, 1, function(x)  paste(x[x!=""] , collapse="*", sep=""))
  value <- paste(ZF_SYMBOLS_CONST, "(", value, ")", sep="")
  z <- paste(value , z, sep="*")
  z[1] <- value[1]
  idx <- as.integer(length(z) / MAXSUB) * MAXSUB  
  if (idx > 0) {
    zz <- z[ 1:idx ]
    dim(zz) <- c(MAXSUB, length(zz) / MAXSUB)
    zz <- apply(zz, 2, function(x) paste("RMplus(", paste(x, collapse=", "), ")" ))
  } else zz <- NULL
  if (idx < length(z)) {
#    Print(zz, idx, idx + 1 == length(z))
    zz[length(zz) + 1] <-  if (idx + 1 == length(z)) z[length(z)] else 
       paste("RMplus(", paste(z[(idx+1) : length(z)], collapse=", "), ")" )
      
#    Print("A", idx, zz, (idx + 1) : length(z))
  }

  if (length(zz) > 1)
    zz <- paste("RMplus(", paste(zz, collapse=", "), ")")

 #  Print( zz)
  ##invisible
  return(eval(parse(text = zz)))
}


