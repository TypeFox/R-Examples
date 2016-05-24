#' ViSiElse: A visual tool for behaviour analysis 
#' 
#' VisiElse is \strong{a graphical tool} designed to visualize and to \strong{give an overview of behavioural observations}
#' realized on individuals or groups. For example, ViSiElse allows \strong{visualization of raw data}
#' during experimental observations of the \strong{realization of a procedure} like a medical algorithm.
#' It graphically presents an overview of individuals and group \strong{actions usually acquired from timestamps 
#' during video recorded sessions}. Options of the package allow adding graphical information 
#' as statistical indicators (mean, standard deviation, quantiles or statistical tests) but also for
#' each action green or black zones providing visual information about the accuracy of the realized actions.
#'
#' @section Principal:
#' ViSiElse concerns \strong{behavioural processes} that, like a simulated healthcare procedure, 
#' can be \strong{decomposed in actions}. We define two different types of actions in ViSiElse: 
#' \strong{punctual and long}. The actions called punctual are brief actions defined as a \strong{time points}. 
#' They \strong{does not last enough to be measured on the chosen time scale}. The actions called long are the 
#' ones defined by \strong{duration}. They are characterized by both \strong{a beginning punctual action and an
#'  ending one}. In order to \strong{model a procedure}, there is a need to sort actions in the way they 
#'  are supposed to be realized as defined for example by medical consensus and/or algorithms for a 
#'  medical procedure. This \strong{structure is stored in a S4-class object ViSibook}, it can be modified, 
#'  printed or plotted.
#'   
#' @section S4-methods:
#' This Package contains two S4class of object: \code{\linkS4class{ViSibook}} and ViSigrid. 
#' Basically a \code{\linkS4class{ViSibook}} object will store information on a process 
#' like a health care procedure, then a \code{\linkS4class{ViSigrid}} object is built with 
#' observations of this procedure and the procedure \code{\linkS4class{ViSibook}}, and finally the graphic 
#' is optained by plotting this \code{\linkS4class{ViSigrid}} object.
#' 
#' 
#' @section ViSibook-class:
#' The \strong{minimum stucture} for a ViSibook must give for each action its name (without special characters but"_"), 
#' its label, and its type (punctual or long). For a long action, in addition it is necessary to provide the two specific
#'  actions that defined its starting and ending.
#'The \strong{order} by which actions are supposed to happen is also required. It is possible to not attribute a rank order at
#' some actions, but not all. 
#' This is useful in the case of a procedure focuses on a long action but 
#' does not on punctual actions defining its end or beginning thus those two punctual actions can
#'  have an order not set.
#' 
#' 
#' Beyond the minimum structure green zones or/and black zones can be defined to help to visualize if a behaviour 
#' is realized on time.
#' Black zones can be defined on both punctual and long action but green zones can only be defined
#'  on punctual actions. To get the details see
#' the definition of the class \code{\linkS4class{ViSibook}}.
#' 

#' @section ViSiGrid-class:
#' A ViSigrid object is built using the function \code{ \link{buildViSiGrid}} with at least : \itemize{
#' \item{The procedure structure}{: It is stored in a ViSibook object.}
#' \item{Observations}{ : Dataset which stores individuals realisation times for each punctal action.}
#' }
#' Options can be add like group, indicator,... see \code{ \link{buildViSiGrid}} and \code{\linkS4class{ViSigrid}} 
#' to get more details explanations.
#' 
#' 
#' At the end the plot of a \code{ViSigrid} object provides a graphics of the realisation of a process. 
#' @docType package
#' @name ViSiElse
NULL

#' @import methods
#' @import grid
#' @import chron
#' @import Matrix
#' @import colorspace
#' @import stringr
#' @importFrom stats median var quantile wilcox.test mood.test
NULL
########################################################################################################################################
#' Class \code{ViSibook} defines the structure of the process to be plotted.
#' @title Class \code{ViSiBook}
#' @name ViSibook-class
#' @rdname ViSibook-class
#' @slot vars a vector storing names of actions.
#' @slot label a vector storing brief description of actions.
#' @slot typeA a vector storing type of actions, "l" for long ( which have a stating time and an ending time ), "p" for punctual.
#' @slot showorder a vector storing order in which actions will be plotted. When an actions is not to be plot
#'  \code{showorder} should be \code{NA}.
#' @slot deb a vector. \itemize{ 
#' \item{Long actions}{ \code{deb} stores the punctual action names that corresponds to long actions beginning.}
#' \item{Punctual action}{ \code{NA} .} 
#' }
#' @slot fin a vector. \itemize{ 
#' \item{Long actions}{ \code{fin} stores the punctual action names that corresponds to long actions ending.}
#' \item{Punctual actions}{ \code{NA} .} 
#' }
#' @slot GZDeb a vector, optional, \code{GZdeb} stores punctual actions green zone starting time.
#' @slot GZFin a vector, optional, \code{GZFin} stores punctual actions green zone ending time.
#' @slot Repetition optional a vector, optional, When a green zone is defined, 
#'  \code{Repetition} stores the length of the time interval between green zones.
#' @slot BZBeforeDeb a vector, optional, \code{BZBeforeDeb} a vector storing punctual black zone 1 starting time.
#' @slot BZBeforeFin a vector, optional, \code{BZBeforeFin} storing punctual black zone 1 ending time.
#' @slot BZAfterDeb a vector, optional, \code{BZAfterDeb} stores punctual black zone 2 starting time .
#' @slot BZAfterFin a vector, optional, \code{BZAfterFin} stores punctual black zone 2 ending time. 
#' @slot BZLong a vector, optional, \code{BZLong} stores the long action black zome time. 
#' @slot BZLtype a vector, optional, \code{BZLtype} stores the type of the black zone,
#'  "time" if the action should be finish at a time, "span" if the action should be finish in a time.
#' @slot NAMES a vector storing names of slots that are to be defined.
#' @seealso \code{\link{plot-ViSibook-method}},
#' \code{\link{changeShoworder-ViSibook-method}}, \code{\link{ConvertFromViSibook-ViSibook-method}}, 
#' \code{\link{ConvertoViSibook}}, \code{\link{intubation}}.
#' and see \code{\link{plot-ViSigrid-method}} for examples.
#' @exportClass ViSibook 
ViSibook <- setClass( "ViSibook",
                   slots = c( vars = "vector" , 
                            label = "vector" , 
                            typeA = "vector" , 
                            showorder = "vector" ,
                            deb = "vector" , 
                            fin = "vector" ,
                            GZDeb = "vector" ,
                            GZFin = "vector" ,
                            Repetition = "vector" ,
                            BZBeforeDeb = "vector" ,
                            BZBeforeFin = "vector" ,
                            BZAfterDeb = "vector" ,
                            BZAfterFin = "vector" ,
                            BZLong = "vector" ,
                            BZLtype = "vector" ,
                            NAMES = "vector"
                   ))
########################################################################################################################################
#' Method get for ViSibook object.
#' @rdname get-ViSibook-methods
#' @aliases [,ViSibook,numeric,missing-method
#' @param x a ViSibook object.
#' @param i a numeric.
#' @param j a numeric.
#' @param drop = TRUE.
#' @return obj.
#' @docType methods
#' @seealso \code{\linkS4class{ViSibook}}.
#' @exportMethod [ 
#'  
setMethod("[", signature = c( "ViSibook","numeric","missing"), function(x , i ,  j ,drop=TRUE) { 
  obj <- c() 
  if (missing( i )) {i <- seq_along( methods::slot( x , "vars") )}
  for (y in methods::slot( x, "NAMES" )[ j ] ) {
     obj <- cbind( obj , methods::slot( x , y )[ i ] ) 
  }  
  colnames( obj ) <- methods::slot( x, "NAMES" )[ j ]
  return(obj  ) 
} )
########################
#' @rdname get-ViSibook-methods
#' @aliases [,ViSibook,missing,numeric-method
#' @docType methods 
#' @exportMethod [  
setMethod("[", c( "ViSibook","missing","numeric"), function(x , i ,  j ,drop=TRUE) { 
  obj <- c() 
  if ( missing(i)) { i <- seq_along( methods::slot( x , "vars") )}
  for (y in methods::slot( x, "NAMES" )[ j ] ) {
    obj <- cbind( obj , methods::slot( x , y )[ i ] ) 
  }  
  colnames( obj ) <- methods::slot( x, "NAMES" )[ j ]
  return(obj  ) 
} )
###################
#' @rdname get-ViSibook-methods
#' @aliases [,ViSibook,numeric,numeric-method
#' @docType methods 
#' @exportMethod [  
setMethod("[", c( "ViSibook","numeric","numeric"), function(x , i ,  j ,drop=TRUE) { 
  obj <- c() 
  if (missing(i)) { i <- seq_along( methods::slot( x , "vars") )}
  for (y in methods::slot( x, "NAMES" )[ j ] ) {
    obj <- cbind( obj , methods::slot( x , y )[ i ] ) 
  }  
  colnames( obj ) <- methods::slot( x, "NAMES" )[ j ]
  return(obj  ) 
} )
########################################################################################################################################
#' Method set for ViSibook object.
#' @name set-ViSibook-method
#' @rdname set-ViSibook-methods
#' @aliases [<-,ViSibook,numeric,numeric,ANY-method
#' @param x a ViSibook object.
#' @param i a numeric.
#' @param j a numeric.
#' @param value object to allocate.
#' @return a ViSibook object.
#' @docType methods
#' @seealso \code{\linkS4class{ViSibook}}
#' @exportMethod [<-  
setReplaceMethod( f = "[" , 
                  signature = c( x = "ViSibook", i = "numeric", j = "numeric" , value = "ANY") , 
                  definition = function(x , i , j , value ){ 
  if (length( j ) > 1 ) { stop( " length( j ) > 1 : \n For ViSibook objects allocation is not defined for several slots at the same time \n ")}
  methods::slot( x, methods::slot( x, "NAMES" )[ j ] )[ i ] <- value; return( x )})
################################
#' @rdname set-ViSibook-methods
#' @aliases [<-,ViSibook,missing,numeric,ANY-method
#' @exportMethod [<-
setReplaceMethod( f = "[" , 
                  signature = c( x = "ViSibook", i = "missing", j = "numeric" , value = "ANY") , 
                  definition = function(x , i , j , value ){ 
                    if (length( j ) > 1 ) { stop( " length( j ) > 1 : \n For ViSibook objects allocation is not defined for several slots at the same time \n ")}
                    methods::slot( x, methods::slot( x, "NAMES" )[ j ] )[ i ] <- value; return( x )})
#################################
#' @rdname set-ViSibook-methods
#' @aliases [<-,ViSibook,numeric,missing,ANY-method
#' @exportMethod [<-
setReplaceMethod( f = "[" , 
                  signature = c( x = "ViSibook", i = "numeric", j = "missing" , value = "ANY") , 
                  definition = function(x , i , j , value ){ 
                    if (length( j ) > 1 ) { stop( " length( j ) > 1 : \n For ViSibook objects allocation is not defined for several slots at the same time \n ")}
                    methods::slot( x, methods::slot( x, "NAMES" )[ j ] )[ i ] <- value; return( x )})
########################################################################################################################################
#' Method Dim for ViSibook object.
#' @name dim-ViSibook-method
#' @title Method \code{dim-ViSibook}
#' @rdname dim-ViSibook-methods
#' @docType methods
#' @aliases dim,ViSibook-method
#' @param x a ViSibook object.
#' @return Vector \itemize{
#' \item{[1]}{ The number of actions defined in x.}
#' \item{[2]}{ The number of characteristics defined in x, 
#' its minimum value is 6 and its maximum is 15.}
#' }
#' @seealso \code{\linkS4class{ViSibook}}
#' @exportMethod dim
setMethod( "dim" , "ViSibook" , function(x ){ return( c( length( methods::slot( x , "vars" ) ) , length( methods::slot( x , "NAMES" ) ) ) ) } )
########################################################################################################################################
#' Method show for ViSibook object.
#' @name show-ViSibook-method
#' @title Method \code{show-ViSibook}
#' @rdname show-ViSibook-methods
#' @aliases show,ViSibook-method
#' @exportMethod show
#' @docType methods
#' @seealso \code{\linkS4class{ViSibook}}.
#' @param object a ViSibook .
setMethod( "show" , "ViSibook" , function(object ){ 
  print( object[ ,  seq( 1 ,length( slot(object , "NAMES" ) ) , 1 ) ] ) }
)
########################################################################################################################################
#' The method \code{changeShoworder} allows the slot "showorder" of a ViSibook to be modified.
#' @name changeShoworder-ViSibook-method
#' @title Method \code{changeShoworder}
#' @rdname changeShoworder-ViSibook-methods
#' @seealso  \code{\linkS4class{ViSibook}} and see \code{\link{plot-ViSigrid-method}} for examples.
setGeneric( "changeShoworder" , function(book , v) standardGeneric( "changeShoworder" ) )
##############
#' @exportMethod changeShoworder
#' @rdname changeShoworder-ViSibook-methods
#' @aliases changeShoworder,ViSibook,numeric-method 
#' changeShoworder
#' @param book a ViSibook object.
#' @param v a vector containing the number from the slot "showorder" of book.
#' @return book a ViSibook object.
setMethod( "changeShoworder" , c("ViSibook","numeric") , function(book,v ) { 
  temp <- c()
  for (i in seq_along(book[ , 4])) { # i = 1
    temp[i] <- switch( as.character( book[ , 4][i] %in% v),
                       "FALSE" = 0,
                       "TRUE" = match( book[ , 4][i] ,v))
  }
 if ( sum(temp == 0) > 0) { temp[ which(temp == 0) ] <- rep( NA , sum(temp == 0)  ) }
  methods::slot(book, "showorder") <- temp
  return( book )
})
########################################################################################################
#' The method \code{ConvertFromViSibook} converts a ViSibook in a data.frame object.
#' @name ConvertFromViSibook-ViSibook-method
#' @title Method \code{ConvertFromViSibook-ViSibook}
#' @rdname ConvertFromViSibook-ViSibook-methods
#' @seealso \code{\linkS4class{ViSibook}} and see \code{\link{plot-ViSigrid-method}} for examples.
setGeneric( "ConvertFromViSibook" , function(x ) standardGeneric( "ConvertFromViSibook" ) )
#' @rdname ConvertFromViSibook-ViSibook-methods
#' @aliases ConvertFromViSibook,ViSibook-method
#' ConvertFromViSibook
#' @exportMethod ConvertFromViSibook
#' @docType methods
#' @param x a ViSibook object.
#' @return a data.frame.
setMethod( "ConvertFromViSibook" , "ViSibook" , function(x ){ 
  ret <- as.data.frame( x[ seq( 1 , length( methods::slot( x , "vars" ) ), 1 ) , 
                           seq( 1 , length( methods::slot( x , "NAMES" ) ) , 1 ) ], colnames = methods::slot( x , "NAMES" ) , stringsAsFactors = FALSE 
  )
  colnames( ret ) <- methods::slot( x , "NAMES" )
  return(ret)
})
