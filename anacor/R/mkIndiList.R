mkIndiList<-function(data, type = rep("C",dim(data)[2]), knots, ord) 
{
# mkIndiList takes a data frame, a vector of types, a list of knot vectors, and a vector
# of orders. It returns a list of codings for the variable. If type is "C" it returns
# a crisp indicator, if type is "A" it returns a numerical version of the variable, and
# if type is "F" is returns the b-spline basis as a fuzzy indicator. In the last case,
# a knot sequence and an order must also be defined for the variable.

m <- dim(data)[2]; n<-dim(data)[1]; fz<-list()
for (j in 1:m) 
   {
   if (type[j]=="C") 
	 {
     fz<-c(fz,list(mkCrisp(data[[j]])))
     names(fz)[[j]] <- colnames(data)[j]
     if (!is.null(rownames(data))) rownames(fz[[j]]) <- rownames(data) 
     if (!is.null(levels(data[,j]))) colnames(fz[[j]]) <- levels(data[,j])
   }
   if (type[j]=="A") {
	  fz <- c(fz,list(as.double(data[[j]])))
    names(fz)[[j]] <- colnames(data)[j]
   }
   if (type[j]=="F") {
	  fz <- c(fz,list(bsplineS(data[[j]],knots[[j]],ord[j])))
	  names(fz)[[j]] <- colnames(data)[j]
   }
 }

return(fz)
}