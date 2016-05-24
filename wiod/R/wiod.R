#' @name wiod
#' @docType package
#' @title World Input Output Database 1995-2011
#' @seealso http://wwww.wiod.org/ http://qua.st/wiod
#' @references {Timmer, Marcel P. (ed) (2012), "The World Input-Output Database (WIOD): Contents Sources and Methods", WIOD Working Paper Number 10, downloadable at http://www.wiod.org/publications/papers/wiod10.pdf }
#' @import decompr
#' @import gvc
#' @examples
#' data(wiod95)
#'
#' library(decompr)
#'
#' w95 <- decomp(inter95,
#'               final95,
#'               countries,
#'               industries,
#'               output95,
#'               method="leontief")
#'
#' library(gvc)
#'
#' i2e95 <- i2e(w95)
#'
NULL
#' @name countries
#' @docType data
#' @title WIOD countries
#' @description the names of the countries
NULL
#' @name industries
#' @docType data
#' @title WIOD industries
#' @description the names of the industries
NULL
#' @name final95
#' @docType data
#' @title WIOD 1995 final
#' @description WIOD 1995 final demand data
NULL
#' @name inter95
#' @docType data
#' @title WIOD 1995 inter
#' @description WIOD 1995 intermediate demand data
NULL
#' @name output95
#' @docType data
#' @title WIOD 1995 output
#' @description WIOD 1995 final output
NULL
#' @name final96
#' @docType data
#' @title WIOD 1996 final
#' @description WIOD 1996 final demand data
NULL
#' @name inter96
#' @docType data
#' @title WIOD 1996 inter
#' @description WIOD 1996 intermediate demand data
NULL
#' @name output96
#' @docType data
#' @title WIOD 1996 output
#' @description WIOD 1996 final output
NULL
#' @name final97
#' @docType data
#' @title WIOD 1997 final
#' @description WIOD 1997 final demand data
NULL
#' @name inter97
#' @docType data
#' @title WIOD 1997 inter
#' @description WIOD 1997 intermediate demand data
NULL
#' @name output97
#' @docType data
#' @title WIOD 1997 output
#' @description WIOD 1997 final output
NULL
#' @name final98
#' @docType data
#' @title WIOD 1998 final
#' @description WIOD 1998 final demand data
NULL
#' @name inter98
#' @docType data
#' @title WIOD 1998 inter
#' @description WIOD 1998 intermediate demand data
NULL
#' @name output98
#' @docType data
#' @title WIOD 1998 output
#' @description WIOD 1998 final output
NULL
#' @name final99
#' @docType data
#' @title WIOD 1999 final
#' @description WIOD 1999 final demand data
NULL
#' @name inter99
#' @docType data
#' @title WIOD 1999 inter
#' @description WIOD 1999 intermediate demand data
NULL
#' @name output99
#' @docType data
#' @title WIOD 1999 output
#' @description WIOD 1999 final output
NULL
#' @name final00
#' @docType data
#' @title WIOD 2000 final
#' @description WIOD 2000 final demand data
NULL
#' @name inter00
#' @docType data
#' @title WIOD 2000 inter
#' @description WIOD 2000 intermediate demand data
NULL
#' @name output00
#' @docType data
#' @title WIOD 2000 output
#' @description WIOD 2000 final output
NULL
#' @name final01
#' @docType data
#' @title WIOD 2001 final
#' @description WIOD 2001 final demand data
NULL
#' @name inter01
#' @docType data
#' @title WIOD 2001 inter
#' @description WIOD 2001 intermediate demand data
NULL
#' @name output01
#' @docType data
#' @title WIOD 2001 output
#' @description WIOD 2001 final output
NULL
#' @name final02
#' @docType data
#' @title WIOD 2002 final
#' @description WIOD 2002 final demand data
NULL
#' @name inter02
#' @docType data
#' @title WIOD 2002 inter
#' @description WIOD 2002 intermediate demand data
NULL
#' @name output02
#' @docType data
#' @title WIOD 2002 output
#' @description WIOD 2002 final output
NULL
#' @name final03
#' @docType data
#' @title WIOD 2003 final
#' @description WIOD 2003 final demand data
NULL
#' @name inter03
#' @docType data
#' @title WIOD 2003 inter
#' @description WIOD 2003 intermediate demand data
NULL
#' @name output03
#' @docType data
#' @title WIOD 2003 output
#' @description WIOD 2003 final output
NULL
#' @name final04
#' @docType data
#' @title WIOD 2004 final
#' @description WIOD 2004 final demand data
NULL
#' @name inter04
#' @docType data
#' @title WIOD 2004 inter
#' @description WIOD 2004 intermediate demand data
NULL
#' @name output04
#' @docType data
#' @title WIOD 2004 output
#' @description WIOD 2004 final output
NULL
#' @name final05
#' @docType data
#' @title WIOD 2005 final
#' @description WIOD 2005 final demand data
NULL
#' @name inter05
#' @docType data
#' @title WIOD 2005 inter
#' @description WIOD 2005 intermediate demand data
NULL
#' @name output05
#' @docType data
#' @title WIOD 2005 output
#' @description WIOD 2005 final output
NULL
#' @name final06
#' @docType data
#' @title WIOD 2006 final
#' @description WIOD 2006 final demand data
NULL
#' @name inter06
#' @docType data
#' @title WIOD 2006 inter
#' @description WIOD 2006 intermediate demand data
NULL
#' @name output06
#' @docType data
#' @title WIOD 2006 output
#' @description WIOD 2006 final output
NULL
#' @name final07
#' @docType data
#' @title WIOD 2007 final
#' @description WIOD 2007 final demand data
NULL
#' @name inter07
#' @docType data
#' @title WIOD 2007 inter
#' @description WIOD 2007 intermediate demand data
NULL
#' @name output07
#' @docType data
#' @title WIOD 2007 output
#' @description WIOD 2007 final output
NULL
#' @name final08
#' @docType data
#' @title WIOD 2008 final
#' @description WIOD 2008 final demand data
NULL
#' @name inter08
#' @docType data
#' @title WIOD 2008 inter
#' @description WIOD 2008 intermediate demand data
NULL
#' @name output08
#' @docType data
#' @title WIOD 2008 output
#' @description WIOD 2008 final output
NULL
#' @name final09
#' @docType data
#' @title WIOD 2009 final
#' @description WIOD 2009 final demand data
NULL
#' @name inter09
#' @docType data
#' @title WIOD 2009 inter
#' @description WIOD 2009 intermediate demand data
NULL
#' @name output09
#' @docType data
#' @title WIOD 2009 output
#' @description WIOD 2009 final output
NULL
#' @name final10
#' @docType data
#' @title WIOD 2010 final
#' @description WIOD 2010 final demand data
NULL
#' @name inter10
#' @docType data
#' @title WIOD 2010 inter
#' @description WIOD 2010 intermediate demand data
NULL
#' @name output10
#' @docType data
#' @title WIOD 2010 output
#' @description WIOD 2010 final output
NULL
#' @name final11
#' @docType data
#' @title WIOD 2011 final
#' @description WIOD 2011 final demand data
NULL
#' @name inter11
#' @docType data
#' @title WIOD 2011 inter
#' @description WIOD 2011 intermediate demand data
NULL
#' @name output11
#' @docType data
#' @title WIOD 2011 output
#' @description WIOD 2011 final output
NULL
