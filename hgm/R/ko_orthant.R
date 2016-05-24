# $OpenXM: OpenXM/src/R/r-packages/hgm/R/ko_orthant.R,v 1.6 2014/04/09 08:31:06 tkoyama Exp $
#dyn.load("hgm_ko_orthant.so")
if (!is.loaded("hgm")) library.dynam("hgm",package="hgm",lib.loc=NULL);

hgm.ko.ncorthant <- function(x,y,rk_step_size=1e-3){
      if (class(x) != "matrix") {
	 print("Error: x is not a square matrix.");
	 return(NULL);      
      }

      if (nrow(x) != ncol(x)) {
	 print("Error: x is not a square matrix.");
	 return(NULL);
      }

      dim <- nrow(x);
      if (!identical(x,t(x))) {
	 print("Error: x is not a symmetric matrix.");
	 return(NULL);
      }

      if (!all(eigen(x)$values > 0)){
	 print("Error: x is not positive definite.");
	 return(NULL);
      }

      if (dim != length(y)) {
         print("Error: The dimensions of x and y differ.");
	 return(NULL);
      }

      .C("hgm_ko_orthant",
	as.integer(dim),
	as.double(t(x)), 
	as.double(y),
	result=double(1))$result
}

HgmKoNcorthant <- hgm.ko.ncorthant
