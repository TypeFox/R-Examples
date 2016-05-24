#' The function \code{buildViSiGrid} build an object class
#' \code{ViSigrid} with at least an object class \code{ViSibook} and data of individual realisation 
#' times for each punctual action defined in the \code{ViSibook}.
#' @name buildViSiGrid
#' @title Function \code{buildViSiGrid}
#' @rdname buildViSiGrid
#' @aliases buildViSiGrid
#' @export buildViSiGrid
#' @param X  A \code{data.frame} or \code{matrix}. \code{X} stores punctual action realization times. 
#' The actions are defined in \code{book}, and X columns names should correspond to the slot "vars" of \code{book}.
#' X must also have a column to identify individuals.
#' @param book  A \code{ViSiBook} object. 
#' It stores the process structure. 
#' @param group A \code{factor} with two \code{levels}. 
#' \code{group} indicates the group attributed to the individuals, 
#' it has same the length as the number of rows of \code{X}.
#' @param Xsup  A \code{data.frame} or \code{matrix} storing supplementary time data,  \strong{all indivuals in}  \code{Xsup}  \strong{must be in } \code{X}.
#' @param method   In \{ \code{"global"} ,  \code{"cut"} ,  \code{"join"} , \code{"within"} \}.
#' \code{method} specifies the plotting method, see \code{details}. If \code{group} is \code{NULL},
#'  \code{method} is set to \code{"global"}.
#' @param grwithin  A level of \code{group}. 
#' If \code{method} is set to \code{within}, \code{grwithin} specifies the group to consider. 
#' @param quantity  
#' In \{ "\code{N}" , "\code{dens}" \}. \code{quantity} allows choosing the quantity represented for punctual action
#' When \code{quantity} is set to "N" the number of individuals is considered. Otherwise when it 
#' is set to "dens" proportion of individuals is considered instead. If \code{group} is defined and 
#' \code{method} set to "cut" or "within", this proportion is calculated regarding each represented group.
#' @param informer  In \{ "\code{NULL}" , "\code{median}" , "\code{mean}" \}. 
#' If \code{informer} is set 
#' to "median" the median and quartiles are computed, 
#' if it is set to "mean" the  mean and standard deviation are. If \code{informer} is \code{NULL} no indicators are computed.
#' @param tests  A boolean. 
#' When \code{informer} is not NULL and \code{group} is defined, if \code{tests} is \code{TRUE}, tests are computed to compare groups. 
#' If the parameter \code{informer} is set to "\code{mean}", 
#' the function \code{wilcox.test()} is used, if \code{informer} is set to "\code{median}" the function \code{mood.test()} is used.
#' @param threshold.test  A numeric between 0 and  1. 
#' \code{threshold.test} is the value of the p-value under which the H0 hypothesis
#'  of the test is rejected when \code{tests} is \code{TRUE}.
#' @param max_tps  A numeric, \code{>0}. \code{max_tps} is the maximum time used to build the grid in the plot.
#' \code{max_tps} is useful when \code{Xsup} is given. If \code{max_tps} is \code{NULL} it is automatically computed.
#' @param pixel An integer. 
#' It is the number of unit of time under which individuals are aggregated in the plot.
#' @param t_0  either 0, either a value of the slot "\code{vars}" in \code{book}, 
#' \code{t_0} indicates the starting time to plot. 
#' @param decrgr2  A boolean. When sorted.line is TRUE and decrgr2 is TRUE, long actions of the second group are plotted in decreasing order by starting times.
#' @param sorted.line  A boolean. 
#' When \code{sorted.line} is \code{TRUE}, it allows long actions to be sorted by starting time.
#' @param times  A boolean. If \code{times} is \code{TRUE}, it incidicates that \code{X} contains data in a time format.
#' @param timeformat  time format.  If \code{times} is \code{TRUE}.
#' @param idsubject  An integer betweem 1 and \code{dim(X)[2]}.  \code{idsubject} indicates the 
#' number of the column of X that contains individuals id numbers.
#' @param colvect   A \code{matrix} containing colors. 
#'  Colors are automatically computed if \code{colvect} is \code{NULL}. 
#'  If \code{group} is not \code{NULL} colvect should have two rows otherwise one. 
#' @param ncolvect  A \code{numeric}. 
#'  \code{ncolvect} indicates the number of columns of \code{colvect}. Its default setting is \code{dim(X)[1]}. 
#'  \code{ncolvect} is considered only if \code{colvect} is \code{NULL}. 
#' @return a ViSigrid object.
#' @details 
#' \itemize{
#' \item{ \code{method} }{   \itemize{
#'  \item{ \code{global} }{:  The plot of the ViSigrid object returned will not consider the parameter \code{group} and plot indistinctly all individuals. }
#'  \item{ \code{cut} }{: In the plot of the ViSigrid object returned each group will be plotted apart within each action line. }
#'  \item{  \code{join} }{: In the plot of the ViSigrid object returned groups will be plotted gathered within each action line. }
#'  \item{\code{within}}{ : In the plot of the ViSigrid object returned , within each action line, there will be two lines, 
#' as for the method \code{cut}, the difference is that the first line will plot all individuals and 
#' the second one individuals belonging to the group specified in \code{grwithin}. }
#'  } } 
#' \item{ \code{informer}}{ 
#' 
#' The parameter \code{informer} allows choosing an indicator. \code{informer} can take three values:
#' \itemize{
#' \item{ \code{median}:}{ Median and quartiles are calculated for each action, using the function 
#' quantile from the package stats. This is the default value.}
#' \item{\code{mean}:}{ Mean and standard deviation are calculated for each action, using the 
#' functions mean and var from the package stats.}
#' \item{\code{NULL}:}{ no indicators are computed.}
#' }
#' When a group is declared indicators are calculated by group if the method cut or 
#' within is chosen. 
#' 
#' When plotting the \code{\linkS4class{ViSigrid}} object, indicators for a punctual action are represented by 
#' white circles linked by a line. For long action, only a black line is plotted from 
#' the median (or mean) of the punctal action staring it. The line length represents 
#' the median (or mean) of the long action duration.
#' Informers are computed directly on the given matrix for punctual action. 
#' And for a long action it is calculated on the difference between the beginning punctual action and the ending one. 
#' }
#' \item{ \code{tests} and \code{threshold.test}}{ 
#' 
#' As for the parameter informer, tests are computed on the given 
#' matrix or data.frame X for a punctual action. And for a long action it is calculated on its difference between its beginning and ending punctual actions.
#' In  \code{\link{plot-ViSigrid-method}}, results of the tests are represented by a star only when the
#'  resulted p-value is bellow or equal to the parameter threshold.test.
#' }
#' \item{ \code{pixel}}{
#' 
#' The parameter pixel represents the number of unit of time under which individuals 
#' are aggregated for punctual action in the plot. When the parameter pixel is too 
#' small the information represented will be too much aggregated to allow interpretation.
#' 
#' For punctual actions data are aggregated in a matrix \eqn{M} . The number of row of \eqn{M} is the 
#' number of action and its number of columns is  \eqn{[ ( max(X)-t_{0} )/pixel]}.
#' 
#' \eqn{M_{i,j}} contains the number of observations of the \eqn{i}-th punctual action (by the order 
#' of the ViSibook object) between \eqn{t_0 + (j-1)pixel} included and 
#' \eqn{t_0 + j*pixel} excluded.
#' }
#' \item{ \code{t_0}}{
#' 
#' The origin of the graphic can be set using the parameter t_0. There is two ways to define it:
#' \itemize{
#' \item{A number:}{  set to 0__. It can be change at convenience, but for long actions black zones will not
#' be drawn, and for punctual actions black and green zones will not be translated.}
#' \item{The name of a punctual action:}{ To set the origin of the graphic to the moment
#'  when the action was done for each individual. Black and green zones will not be translated as well.}
#' }
#' }
#' 
#' }
#' 
#' @details x can also has the colunms : GZDebn,  GZFin, Repetition, BZBeforeDeb, BZBeforeFin, BZAfterDeb, BZAfterFin, BZLong , BZLtype 
#' @seealso Dataset \code{\link{intubation}}. Classes \code{\linkS4class{ViSigrid}} and \code{\linkS4class{ViSibook}}.
#'  The method plot for ViSigrid object \code{\link{plot-ViSigrid-method}} for examples.
#'
buildViSiGrid <- function(X , book , Xsup = NULL, 
                       method = "global" , group = NULL  , grwithin = NULL , 
                       informer = "median" , tests = TRUE , threshold.test = 0.01 ,
                       quantity = "N" , pixel = 20 , t_0 = 0, sorted.line = TRUE ,  
                       decrgr2 = FALSE, max_tps = NULL, colvect = NULL , ncolvect = dim(X)[1] ,
                       times = FALSE , timeformat = c('hh:mm:ss') ,	idsubject = 1	) {
  # A data.frame containing the times performances for punctuals actions, colnames( X ) must corresponds to names given in the book in the slot vars
  ### Verification group if method !="global"
  if (method != "global" ) {
    if (is.null( group ) ) {
      warning( " Group is NULL, method set to \"global\" \n " )
      method <- "global"
    }else{
      group <- factor(group)
      if (nlevels( factor( group ) ) > 2) {
        warning( " Incorrect number of groups (2 max), method set to \"global\" \n " )
        method <- "global"
      }
    }
  } 
  #...................................................................................................................................................................................................................
  ### contruction of the grid matrix
  if (method == "global" ) {
    temp <- MATgrid( X , book , pixel = pixel , times = times , timeformat = timeformat , idsubject = idsubject , retX = TRUE , t_0 = t_0 , max_tps = max_tps )
    MATp <- temp$MATGrid
    X <- temp$X
    vect_tps <- temp$vect_tps 
    t_0 <- temp$t_0
  }else{
    group <- as.factor(group)
    if (method == "join" || method == "cut" || method == "within" ) {
      ####For within method coying of the data 
      if (method == "within" ) {
        # if the level to considered is not given the function arbitraly select the first level
        grwithin <- switch( as.character( is.null( grwithin ) ) , "TRUE" = levels( group )[ 1 ], "FALSE" = grwithin )
        temp  <-  dim( X )[ 1 ]
        # Coping time data 
        X <- Matrix::rBind( X , X[ which( group == grwithin ), ] )
        # Adapting the group vector to the new time matrix 
        group <- gl( 2 , temp , labels = c( "All" , as.character( grwithin ) ) ) [ seq( 1 , temp + sum( group == grwithin ) , 1 ) ]
      }
      vect_tps <- MATgrid( X , book , pixel = pixel , times = times , timeformat = timeformat , idsubject = idsubject , onlyvect_tps = TRUE , t_0 = t_0 , max_tps = max_tps)
      MATp <- Matrix::Matrix( rep( 0 , sum( methods::slot( book , "typeA") == "p" ) * 2 *  length( vect_tps ) ) , nrow =  sum( methods::slot( book , "typeA") == "p" ) * 2 , length( vect_tps ) , sparse = TRUE )
      temp <- MATgrid( X[ which( group == levels( group )[ 1 ] ) , ] , book , pixel = pixel , times = times , timeformat = timeformat , idsubject = idsubject ,  vect_tps = vect_tps , retX = TRUE)
      MATp[ seq( 1 , sum( methods::slot( book , "typeA") == "p" ) , 1 ), ] <- temp$MATGrid
      tempX <- temp$X
      temp <- MATgrid( X[ which( group == levels( group )[ 2 ] ) , ] , book , pixel = pixel , times = times , timeformat = timeformat ,  vect_tps = vect_tps , idsubject = idsubject , retX = TRUE)
      MATp[ seq( sum( methods::slot( book , "typeA") == "p" ) + 1 , sum( methods::slot( book , "typeA") == "p" ) *2 , 1 ), ] <-  temp$MATGrid 
      t_0 <- temp$t_0
      X <- Matrix::rBind( tempX , temp$X )
      rm( tempX )
    }
  }
  #...................................................................................................................................................................................................................
  #### Density or count  & colors grid matrice
  if (quantity == "dens" ) {
    if (method != "global" & method != "join" ) {
      MATp <- Matrix::rBind( round( MATp[ seq( 1 , dim( MATp )[ 1 ] / 2 , 1 ), ] / sum( group == levels( group )[ 1 ] ) * 10 ^ nchar( sum( group == levels( group )[ 1 ] ) ) , 0 ) , round( MATp[ seq( dim( MATp )[ 1 ] / 2 + 1 , dim( MATp )[ 1 ] , 1 ) , ] / sum( group == levels( group )[ 2 ] ) * 10 ^ nchar( sum( group == levels( group )[ 2 ] ) ) , 0 ) )
    }else{
      MATp <- round( MATp / (dim( X )[ 1 ] ) * 10 ^ (nchar( dim( X )[ 1 ] ) ) , 0 )
    }
    ncolvect <- max(MATp) + 1
  }
  if (method == "global" ) {
    if (is.null( colvect ) ) {
      colvect <- matrix( colorspace::sequential_hcl( n = ncolvect , 
                                         h = 264 , 
                                         c. = c( 80, 85 ), 
                                         l = c( 30, 95 ), 
                                         power = 0.7 
                                        )[ seq( ncolvect , 1 , -1 ) ] , nrow = 1 )
    }else{
      if (length( colvect ) < max( MATp ) ) {
        warning( paste( "  length( colvect ) ", max( MATp ) , "\n" ) )
        MATp <- round( MATp / (dim( X ) [ 1 ] ) * (length( colvect ) ) , 0 )
      }
      colvect <- matrix( colvect , nrow = 1 )
    }
  }else{
    if (is.null(colvect) ) {
      colvect <- as.matrix( Matrix::rBind( colorspace::sequential_hcl( n = ncolvect , 
                                                   h = 264, 
                                                   c. = c( 80 , 85 ) , 
                                                   l = c( 30 , 95 ) , 
                                                   power = 0.7 
                                                   )[ seq( ncolvect , 1 , -1 ) ] , 
                                   colorspace::sequential_hcl( n = ncolvect , 
                                                   h = -32 , 
                                                   c. = c( 80 , 85 ) , 
                                                   l = c( 30 , 95 ) , 
                                                   power = 0.7 
                                                   )[ seq( ncolvect , 1 , -1 ) ] ) )
    }else{
      colvect <- as.matrix( colvect )
      if ((dim( colvect )[ 1 ] ) < max( MATp ) ) {
        warning( paste( " length(colvect)< ", max( MATp ) , "\n" ) )
        MATp <- round( MATp / (dim( X )[ 1 ] ) * (dim( colvect )[ 1 ] ) , 0 )
      }
    }
  }
  #...................................................................................................................................................................................................................
  ### informers Punctuals action
  if (is.null( informer ) == FALSE ) {
    if (method == "join" || method == "global" ) {
      if (informer == "mean") {
        informers <- apply( X[ , -idsubject ] , MARGIN = 2 , FUN = function(x )(c( mean( x , na.rm = TRUE ) - sqrt( stats::var( x , na.rm = TRUE ) ) , mean( x , na.rm = TRUE ) , mean( x , na.rm = TRUE ) + sqrt( var( x , na.rm = TRUE ) ) ) ) )
      }else{
        if (informer == "median") {
          informers <- apply( X[ , -idsubject ] , MARGIN = 2 , FUN = function(x )(c( stats::quantile( x , na.rm = TRUE )[ 2 ]  , stats::median( x , na.rm = TRUE )  , stats::quantile( x , na.rm = TRUE )[ 4 ] ) ) )
        }else{
          warning( " informer value not recognized, it must be \"mean\", \" median\" or NULL - informer set to NULL \n " )
          informer <- NULL
        }	
      }
    }else{
      if ( method == "cut" || method == "within" ) {
        if (informer == "mean") {
          informers <- Matrix::rBind(
            apply( X[ which( group == levels( group ) [ 1 ] ) , -idsubject ] , MARGIN = 2 , FUN = function(x )(c( mean( x , na.rm = TRUE ) - sqrt( stats::var( x , na.rm = TRUE ) ) , mean( x , na.rm = TRUE ) , mean( x , na.rm = TRUE ) + sqrt( stats::var( x , na.rm = TRUE ) ) ) ) ),
            apply( X[ which( group == levels( group ) [ 2 ] ) , -idsubject ] , MARGIN = 2 , FUN = function(x )(c( mean( x , na.rm = TRUE ) - sqrt( stats::var( x , na.rm = TRUE ) ) , mean( x , na.rm = TRUE ) , mean( x , na.rm = TRUE ) + sqrt( stats::var( x , na.rm = TRUE ) ) ) ) )
          )
        }else{
          if (informer == "median") {
            informers <- Matrix::rBind( apply( X[ which( group == levels( group )[ 1 ] ) , -idsubject ] , MARGIN = 2 , FUN = function(x )(c( stats::quantile( x , na.rm = TRUE )[ 2 ]  , stats::median( x , na.rm = TRUE )  , stats::quantile( x , na.rm = TRUE )[ 4 ] ) ) ) ,
                                apply( X[ which( group == levels( group )[ 2 ] ) , -idsubject ] , MARGIN = 2 , FUN = function(x )(c( stats::quantile( x , na.rm = TRUE )[ 2 ]  , stats::median( x , na.rm = TRUE )  , stats::quantile( x , na.rm = TRUE )[ 4 ] ) ) ) )
          }else{
            warning( " informer value not recognized, it must be \"mean\", \" median\" or NULL - informer set to NULL \n " )
            informer <- NULL
          }	
        }
      }
    }
  }else(
    informers <- new( "matrix" )
  )
  #...................................................................................................................................................................................................................
  ################# Tests Punctuals Actions
  if (tests == TRUE & method != "global" & is.null( informer ) == FALSE  ) {
    if (informer == "mean") {
      testsP <- apply( X[ , -idsubject ] , MARGIN = 2 , FUN = function(x ) { stats::wilcox.test( x[ which( group == levels( group )[ 1 ] ) ] , x[ which( group == levels( group )[ 2 ] ) ] )$p.value } ) < threshold.test
    }else{
      if (informer == "median") {
        testsP <- apply( X[ , -idsubject ] , MARGIN = 2 , function(x ) { stats::mood.test( x[ which( group == levels( group )[ 1 ] ) ] , x[ which( group == levels( group )[ 2 ] ) ] )$p.value } ) < threshold.test
      }	
    }
  }else{
    if (tests == TRUE & method == "global" ) {
      tests <- FALSE
      testsP <- NULL
    }
    if (tests == TRUE & is.null( informer ) ) {
      warning("if test=TRUE, informer should not be NULL, no tests are computed")
    }
    tests <- FALSE
    testsP <- NULL
  }
  #...................................................................................................................................................................................................................
  ################### Long Actions 
  if (any(methods::slot( book , "typeA") == "l" &  is.na( book[ , 4] ) == FALSE ) ) {
  sia <- which( methods::slot( book , "typeA") == "l" &  is.na( book[ , 4] ) == FALSE )
  ia <- sia[ which.min(book[ , 4 ][ sia ]) ]
  temp <- buildL( X , book , ia , group , decrgr2 , sorted.line , method , vect_tps)
  idsort <- temp$idsort 
  L <- temp$L 
  BZL <- new( "dgCMatrix" )
  BZL <- Matrix::Matrix( rep( 0 , dim( X )[ 1 ] * sum( methods::slot( book , "typeA") == "l" &  is.na( book[ , 4] ) == FALSE ) ) , 
                 nrow =  dim( X )[ 1 ] , 
                 ncol = sum( methods::slot( book , "typeA") == "l" &  is.na( book[ , 4] ) == FALSE ) , sparse = TRUE , doDiag = FALSE )
  BZL <- temp$BZL 
  for (ia in sia[ order(book[ , 4][ sia ]) ][ -1 ] ) {    
    temp <- buildL( X , book , ia , group , decrgr2 , sorted.line , method , vect_tps)
    idsort <- cbind( idsort , temp$idsort )
    L <- cbind( L , temp$L )
    BZL <- cbind( BZL , temp$BZL )
  }
  #...................................................................................................................................................................................................................
  ################### Long Actions Informers
  if (is.null( informer ) == FALSE ) {
    if (method == "join" || method == "global" ) {
      if (informer == "mean") {
        for (ia in seq( 1 , dim( L )[ 2 ] , 2 ) ) {  # ia = 1
          informers <- cbind( informers , 
                              c( mean( L[ , ia + 1 ] - L[ , ia ] , na.rm = TRUE ) - sqrt( stats::var( L[ , ia + 1 ] - L[ , ia ] , na.rm = TRUE ) ),
                                 mean( L[ , ia + 1 ] - L[ , ia ] , na.rm = TRUE ) ,
                                 mean( L[ , ia + 1 ] - L[ , ia ] , na.rm = TRUE ) + sqrt( stats::var( L[ , ia + 1 ] - L[ , ia ] , na.rm = TRUE ) ) ) )
        }
      }else{
        for (ia in seq( 1 , dim( L )[ 2 ] , 2 ) ) {
          if (informer == "median" ) {
            informers <- cbind( informers , 
                                c( stats::quantile( L[ , ia + 1 ] - L[ , ia ] , na.rm = TRUE )[ 2 ] ,
                                   stats::quantile( L[ , ia + 1 ] - L[ , ia ] , na.rm = TRUE )[ 3 ],
                                   stats::quantile( L[ , ia + 1 ] - L[ , ia ] , na.rm = TRUE )[ 4 ] ) )	
          }
        }
      }
    }else{
      if (informer == "mean" ) {
        for (ia in seq( 1 , dim( L )[ 2 ] , 2 ) ) {
          informers <- cbind( informers , 
                              c( mean( L[ which( group == levels( group ) [ 1 ] ), ia + 1 ] - L[ which( group == levels( group ) [ 1 ] ) , ia ] , na.rm = TRUE ) - sqrt( stats::var( L[ which( group == levels( group ) [ 1 ] ) , ia + 1 ] - L[ which( group == levels( group ) [ 1 ] ) , ia ] , na.rm = TRUE ) ),
                                 mean( L[ which( group == levels( group ) [ 1 ] ) , ia + 1 ] - L[ which( group == levels( group ) [ 1 ] ) , ia ] , na.rm = TRUE ),
                                 mean( L[ which( group == levels( group ) [ 1 ] ), ia + 1 ] - L[ which( group == levels( group ) [ 1 ] ) , ia ] , na.rm = TRUE ) + sqrt( stats::var( L[ which( group == levels( group ) [ 1 ] ) , ia + 1 ] - L[ which( group == levels( group ) [ 1 ] ) , ia ] , na.rm = TRUE ) )  ,
                                 mean( L[ which( group == levels( group ) [ 2 ] ), ia + 1 ] - L[ which( group == levels( group ) [ 2 ] ) , ia ] , na.rm = TRUE ) - sqrt( stats::var( L[ which( group == levels( group ) [ 2 ] ) , ia + 1 ] - L[ which( group == levels( group ) [ 2 ] ) , ia ] , na.rm = TRUE ) ),
                                 mean( L[ which( group == levels( group ) [ 2 ] ) , ia + 1 ] - L[ which( group == levels( group ) [ 2 ] ) , ia ] , na.rm = TRUE ),
                                 mean( L[ which( group == levels( group ) [ 2 ] ), ia + 1 ] - L[ which( group == levels( group ) [ 2 ] ) , ia ] , na.rm = TRUE ) + sqrt( stats::var( L[ which( group == levels( group ) [ 2 ] ) , ia + 1 ] - L[ which( group == levels( group ) [ 2 ] ) , ia ] , na.rm = TRUE ) ) ) )
        }
      }else{
        if (informer == "median" ) {
          for (ia in seq( 1 , dim( L )[ 2 ] , 2 ) ) {
            informers <- cbind( informers , 
                                c( stats::quantile( L[ which( group == levels( group ) [ 1 ] ) , ia + 1 ] - L[ which( group == levels( group ) [ 1 ] ) , ia ] , na.rm = TRUE )[ 2 ] ,
                                   stats::quantile( L[ which( group == levels( group ) [ 1 ] ), ia + 1 ] - L[ which( group == levels( group ) [ 1 ] ) , ia ] , na.rm = TRUE )[ 3 ],
                                   stats::quantile( L[ which( group == levels( group ) [ 1 ] ), ia + 1 ] - L[ which( group == levels( group ) [ 1 ] ) , ia ] , na.rm = TRUE )[ 4 ],
                                   stats::quantile( L[ which( group == levels( group ) [ 2 ] ) , ia + 1 ] - L[ which( group == levels( group ) [ 2 ] ) , ia ] , na.rm = TRUE )[ 2 ] ,
                                   stats::quantile( L[ which( group == levels( group ) [ 2 ] ), ia + 1 ] - L[ which( group == levels( group ) [ 2 ] ) , ia ] , na.rm = TRUE )[ 3 ],
                                   stats::quantile( L[ which( group == levels( group ) [ 2 ] ), ia + 1 ] - L[ which( group == levels( group ) [ 2 ] ) , ia ] , na.rm = TRUE )[ 4 ] 
                                ))	
          } 
        }
      }
    }
  index <- unlist( lapply(methods::slot(book ,"deb")[ sia[ sort(book[ , 4][ sia ], index.return = TRUE)$ix ]],
                          function(x )(which( colnames(informers)[ seq(1 , dim(informers)[ 2 ] - dim(L)[ 2 ] / 2 , 1)] == x  )))
  )
  informers <- cbind( informers , informers[ ,seq(dim(informers)[ 2 ] - dim(L)[ 2 ] / 2 + 1 , dim(informers)[ 2 ] ,1) ]  + informers[ , index ] )
    colnames(informers)[seq( sum( methods::slot( book , "typeA") == "p") + 1 ,sum( methods::slot( book , "typeA") == "p" ) + dim( L )[ 2 ]  , 1 )] = c(
      unlist( lapply( methods::slot( book , "vars")[ which( methods::slot( book , "typeA") == "l" &  is.na( book[ , 4] ) == FALSE ) ], function(x)(paste0("span_" , x)))) ,
      unlist( lapply( methods::slot( book , "vars")[ which( methods::slot( book , "typeA") == "l" &  is.na( book[ , 4] ) == FALSE ) ], function(x)(paste0("plot_" , x)))) )
  }
  if (tests == TRUE & method != "global" & is.null( informer ) == FALSE  ) {
    temp <- matrix( unlist( lapply( seq( 1 , dim( L )[ 2 ] , 2 ) ,  function(x ) (L[ , x + 1 ] - L[ , x ] ) )  ) , nrow = dim( L )[ 1 ] , byrow = FALSE )
    if (informer == "mean") {
      testsP <- c( testsP , apply( temp , MARGIN = 2 , FUN = function(x ) { stats::wilcox.test( x[ which( group == levels( group )[ 1 ] ) ] , x[ which( group == levels( group )[ 2 ] ) ] )$p.value } ) < threshold.test )
    }else{
      if (informer == "median") {
        testsP <- c( testsP , apply( temp , MARGIN = 2 , function(x ) { stats::mood.test( x[ which( group == levels( group )[ 1 ] ) ] , x[ which( group == levels( group )[ 2 ] ) ] )$p.value } ) < threshold.test )
      }	
    }
  }
  #...................................................................................................................................................................................................................
  }else{
  L <- NULL
  BZL <- NULL
  idsort <- NULL
  }
  if ( is.matrix(idsort) == FALSE & is.null(idsort) == FALSE ) {
    idsort <- matrix(idsort)
  }
   ################### Supplementary individuals 
  if (is.null( Xsup ) == FALSE ) {
    if ( method == "global" || is.null( group ) || nlevels( factor( group ) ) > 2 ) {
      temp 	<- MATgrid( Xsup , book , pixel = pixel , times = times , timeformat = timeformat , idsubject = idsubject ,  vect_tps = vect_tps , retX = TRUE )
      MATpsup <- temp$MATGrid
      Xsup 	<- temp$X
    }else{
      if (method == "join" || method == "cut" || method == "within" ) {
        idsup 		<- Xsup[ , idsubject ]
        groupsup	<- group[ idsup ]				
        if (method == "within" ) { 					
          temp 	<- dim( Xsup )[ 1 ]
          temp1 	<- which( groupsup == rep( grwithin , length( groupsup ) ) )
          if (any( groupsup == rep( grwithin , length( groupsup ) ) ) ) {
            Xsup 		<- Matrix::rBind( Xsup , Xsup[ temp1 , ] )						
            groupsup	<- gl( n = 2 , k = temp , labels = c( "All" , grwithin ) )[ seq( 1 , temp + length( temp1 ) , 1) ]
            idsup 		<- c( idsup , idsup[ temp1 ] )
          }
        }
        groupsup 	<- factor(groupsup)
        MATpsup 	<- Matrix::Matrix( rep( 0 , sum( methods::slot( book , "typeA") == "p" ) * 2 *  length( vect_tps ) ) , nrow =  sum( methods::slot( book , "typeA") == "p" ) * 2 , length( vect_tps ) , sparse = TRUE )
        temp <- MATgrid( Xsup[ which( groupsup == levels( group )[ 1 ] ) , ] , book , pixel = pixel , times = times , timeformat = timeformat , idsubject = idsubject ,  vect_tps = vect_tps , retX = TRUE )
        MATpsup[ seq( 1 , sum( methods::slot( book , "typeA") == "p" ) , 1), ] <- temp$MATGrid
        tempX <- temp$X
        temp <- MATgrid( Xsup[ which( groupsup == levels( group )[ 2 ] ) , ] , book , pixel = pixel , times = times , timeformat = timeformat , idsubject = idsubject , vect_tps = vect_tps , retX = TRUE ) 
        MATpsup[ seq( 1 + sum( methods::slot( book , "typeA") == "p" ), 2 * sum( methods::slot( book , "typeA") == "p" ) , 1 ) , ] <- temp$MATGrid
        Xsup <- Matrix::rBind( tempX , temp$X )
      }
    }
    if (quantity == "dens" ) {
      if (method != "global" & method != "join" ) {
        MATpsup <- Matrix::rBind( round( MATpsup[ seq( 1 , dim( MATpsup )[ 1 ] / 2 , 1 ), ] / sum( group == (levels( group )[ 1 ] ) ) * 10 ^ (nchar( sum( group == (levels( group )[ 1 ] ) ) ) ) , 0 ) ,
                          round( MATpsup[ seq( dim( MATpsup )[ 1 ] / 2 + 1 , dim( MATpsup )[ 1 ] , 1 ),] / sum( group == (levels( group )[ 2 ] ) ) * 10 ^ (nchar( sum( group == (levels( group )[ 2 ] ) ) ) ) , 0 ) )
      }else{
        MATpsup <- round( MATpsup / (dim( X )[ 1 ] ) * 10 ^ (nchar( dim( X ) [ 1 ] ) ) , 0 )
      }
    }
    if (any(methods::slot( book , "typeA") == "l") ) {
      ia <- which( methods::slot( book , "typeA") == "l" &  is.na( book[ , 4] ) == FALSE )[ 1 ]
      temp <- buildL( Xsup , book , ia , group = groupsup , decrgr2 , sorted.line = FALSE , method , BZL = FALSE , vect_tps)
      idsortsup <- temp$idsort 
      Lsup <- temp$L 
      for (ia in which( methods::slot( book , "typeA") == "l" &  is.na( book[ , 4] ) == FALSE )[ -1 ] ) {    # ia =2
        temp <- buildL( Xsup , book , ia , group , decrgr2 , sorted.line = FALSE , method , BZL = FALSE , vect_tps)
        idsortsup <- cbind( idsortsup , temp$idsort )
        Lsup <- cbind( Lsup , temp$L )
      }
    }else{
      Lsup <- NULL
    }
  }
  ret <- ViSigrid( 	MATp = MATp ,
                    MATpsup = switch( as.character( is.null( Xsup ) ) , "TRUE" = new( "dgCMatrix" ) , "FALSE" = MATpsup ) ,
                    idsup = switch( as.character( is.null( Xsup ) ) , "TRUE" = vector() , "FALSE" = Xsup[ , idsubject] ) ,
                    colvect = colvect ,
                    L = switch( as.character( is.null( L ) ) , "TRUE" = new( "data.frame" ) , "FALSE" = L )  ,
                    idsort = switch( as.character( is.null(idsort ) ) , "TRUE" = new( "matrix" ) , "FALSE" = idsort )   ,
                    BZL = switch( as.character( is.null( BZL ) ) , "TRUE" = new( "dgCMatrix" ) , "FALSE" = BZL ) ,
                    Lsup = switch( as.character( is.null( Xsup ) ) , "TRUE" = new( "data.frame" ) , "FALSE" = Lsup ) ,
                    book = book ,
                    group = switch( as.character( is.null( group ) ) , "TRUE" = factor() , "FALSE" = group ) ,
                    vect_tps = vect_tps ,
                    informers = informers ,
                    testsP = switch( as.character( is.null( testsP ) ) , "TRUE" = vector(), "FALSE" = testsP ),
                    parameters = list( 
                      method = method , 
                      grwithin = switch( as.character( is.null( grwithin ) ) ,"TRUE" = NULL, "FALSE" = grwithin ) ,
                      quantity = quantity ,
                      informer = informer ,
                      tests = tests ,
                      threshold.test = threshold.test ,
                      pixel = pixel ,
                      t_0 = t_0
                    )
  )
  return(ret)
}
