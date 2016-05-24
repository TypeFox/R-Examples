
##### Copyright 1986-2009 David W. Scott
#####
##### This program is free software; you can redistribute it and/or
##### modify it under the terms of the GNU General Public License as 
##### published by the Free Software Foundation; either version 2 of 
##### the License, or (at your option) any later version.
#####
##### This program is distributed in the hope that it will be useful,
##### but WITHOUT ANY WARRANTY; without even the implied warranty of
##### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##### See the GNU General Public License for more details.
#####
##### You should have received a copy of the GNU General Public
##### License along with this program; if not, write to the Free 
##### Software Foundation, Inc.,
##### 51 Franklin St, Fifth Floor, 
##### Boston, MA  02110-1301 USA
#####
##### On Debian GNU/Linux systems, the complete text of the GNU
##### General Public License can be found in
##### /usr/share/common-licenses/GPL-2.

bin1 <- function(x,ab=nicerange(x),nbin=50) {
  n <- length(x)
  if(ab[1]>=ab[2])
    stop("Interval vector has negative orientation")
  if(nbin<=0)
    stop("Number of bin intervals nonpositive")
  r<-.Fortran("bin1",
           as.double(x),
           as.integer(n),
           as.double(ab),
           as.integer(nbin),
           nc=integer(nbin),
           nskip=integer(1),
           PACKAGE="ash")
  list(nc=r$nc,ab=ab,nskip=r$nskip)
}

ash1 <- function(bins,m=5,kopt=c(2,2)){
  nc <- bins$nc
  ab <- bins$ab
  nbin <- length(nc)
  r <- .Fortran("ash1",
                as.integer(m),
                as.integer(nc),
                as.integer(nbin),
                as.double(ab),
                as.integer(kopt),
                t=double(nbin),
                f=double(nbin),
                double(m),
                ier=integer(1),
                PACKAGE="ash")
  if(r$ier==1) print("ash estimate nonzero outside interval ab")
  list(x=r$t,y=r$f,m=m,ab=ab,kopt=kopt,ier=r$ier)
}

nicerange <- function(x, beta = .1) {
 ab <- range(x)	# interval surrounding data
 del <- ((ab[2] - ab[1]) * beta)/2
 return(c(ab + c( - del, del)))
}

bin2 <- function(x,ab,nbin=c(20,20)){
  if(missing(ab)){
    ab <- t(array(c(nicerange(x[,1]),nicerange(x[,2])),c(2,2)))
  }
  n <- nrow(x)
  r <- .Fortran("bin2",
           as.double(x),
           as.integer(n),
           as.double(ab),
           as.integer(nbin[1]),
           as.integer(nbin[2]),
           nc=integer(nbin[1]*nbin[2]),
           nskip=integer(1),
           PACKAGE="ash")
  list(nc=matrix(r$nc,nbin[1],nbin[2]),ab=ab,nskip=r$nskip)
}

ash2 <- function(bins,m=c(5,5),kopt=c(2,2)){
  nc <- bins$nc;  if(!is.matrix(nc)) stop("bins does not contain bin count matrix")
  ab <- bins$ab;  if(!is.matrix(ab)) stop("ab not a matrix - should be 2 by 2")
  nbin <- dim(nc)
  r <- .Fortran("ash2",
                as.integer(m[1]),
                as.integer(m[2]),
                as.integer(nc),
                as.integer(nbin[1]),
                as.integer(nbin[2]),
                as.double(ab),
                as.integer(kopt),
                f=double(nbin[1]*nbin[2]),
                double(m[1]*m[2]),
		        ier=double(1),
                PACKAGE="ash")
  if(r$ier==1) print(" estimate nonzero outside ab rectangle")
  list(z=matrix(r$f,nbin[1],nbin[2]),
 x=center(ab[1,],nbin[1])[[1]],y=center(ab[2,],nbin[2])[[1]],
	ab=ab,m=m,kopt=kopt,ier=r$ier)
}

center <- function(ab,k) {
	h <- diff(ab)/k
 list( seq(ab[1]+h/2,by=h,length=k) ) }


