"x2u"<-
function(x, labels = seq(along = x))
{
#  DATE WRITTEN:  21 January 1994       LAST REVISED:  21 January 1994
#  AUTHOR:  Scott Chasalow
#
#  DESCRIPTION:
#        Convert an x-encoded simplex-lattice point to a u-encoded
#        simplex-lattice point  (equivalently,  "untabulate" bin counts)
#
#  USAGE:
#        x2u(x)
#
#  ARGUMENTS:
#  x:    A numeric vector.  x[i] is interpreted as the count in bin i.
#  labels:  A vector.  Interpreted as the bin labels;  default value is
#        seq(along = x), which causes return of a u-encoded simplex-lattice 
#        point.  Other values of labels cause return of the result of 
#        subscripting labels with the u-encoded simplex-lattice point that 
#        would have been obtained if the default value of labels were used.
#
#        Arguments x and labels must be of equal length.
#
#  VALUE:
#        rep(labels, x), a vector of length sum(x).  If labels = seq(along = x)
#        (the default),  value is the u-encoded translation of the simplex 
#        lattice point, x.  Equivalently,  value gives the bin numbers, 
#        in lexicographic order,  for the objects represented by the counts in 
#        x.  For other values of argument "labels", value gives the bin labels 
#        for the objects represented by the counts in x (equivalent to 
#        labels[x2u(x)]).
#
#  SEE ALSO:
#        tabulate,  rep
#
	if(length(labels) != length(x)) stop(
			"Arguments x and labels not of equal length")
	rep(labels, x)
}

