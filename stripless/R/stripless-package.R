#'Trellis Displays Without Strips For Lattice Graphics
#'
#'Stripless provides a simple interface to make trellis plots without strip 
#'labels using lattice graphics. Instead of having strip labels identify each 
#'panel's settings of the conditioning variables, their position in a structured
#'row/column layout encodes this information, and a legend decodes them. This
#'avoids cluttering the plot area with redundant layers of strip labels when
#'there are many conditioning variables, as is often the case for factorial
#'designs, especially the (fractions of) two level factorials widely used in
#'industry.
#'
#'The package consists of a generic \code{strucplot} function and methods, an
#'associated \code{print} method for objects inheriting from class
#'\code{"structured"}, and some (currently 1) additional trellis
#'panel functions. This functionality is implemented as a wrapper
#'to lattice's \code{xyplot} function, so all other plot specifications
#'(panel functions, colors, titles, etc.) are given through the usual
#'\code{\link[lattice]{xyplot}} arguments.
#'
#'See \code{\link{strucplot}} for a full explanation of the structured plotting 
#'paradigm, how \code{strucplot} implements it, and examples illustrating how
#'it works.
#'
#'@note The author has made a considerable effort to provide clear, complete
#'documentation of the package's functionality. He would therefore appreciate
#'receiving any reports of errors, inconsistencies, or infelicities in the Help
#'docs, as well as suggestions for improvement in the docs or underlying 
#'package functionality.
#'  
#'@name stripless-package
#'@docType package
#'  
NULL