#' Method plot for ViSigrid object. This method provides a graphic
#' of raw data during experimental observations of the realization of a procedure
#' like a medical algorithm. It graphically presents an overview of individuals
#' and group actions usually acquired from timestamps during video recorded sessions. 
#' @title Method \code{plot-ViSigrid}
#' @name plot-ViSigrid-method
#' @rdname plot-ViSigrid-methods
#' @aliases plot,ViSigrid-method
#' @exportMethod plot
#' @docType methods
#' @param x A \code{ViSigrid} object built using the \code{\link{buildViSiGrid}} function.
#' @param scal.unit.tps Unity of time for the grey grid legend.
#' @param unit.tps Unit of time (s,min,..).
#' @param Fontsize.label.Action Fontsize of labels of plotted actions.
#' @param Fontsize.label.Time Fontsize of the time axe.
#' @param Fontsize.label.color  Fontsize of legends.
#' @param Fontsize.title Fontsize of the title.
#' @param colgreenzone Color of the green zones.
#' @param colblackzone Color of black zones.
#' @param alphaZones Alpha of green and black zones.
#' @param linA Height of the plotting area in each actions lines, < 1.
#' @param alphainf Alpha of informers circles.
#' @param alphasup Alpha of supplementary times.
#' @param rcircle circle radius of informers circles.
#' @param lwdline line width of lines linking the 3 informers circles.
#' @param main Title.
#' @param size.main Title size.
#' @param col.main Title color.
#' @param col.grid Color of the legend box.
#' @param lwd.grid Lines width of the legend box.
#' @param lty.grid Lines type of the legend box.
#' @param ncharlabel Maximum number of plotted characters for labels of actions.
#' @param vp0h Height of the main plot window, <1.
#' @param vp0w Width  of the main plot window, <1.
#' @seealso \code{\linkS4class{ViSigrid}}, \code{\linkS4class{ViSibook}},
#'  \code{\link{buildViSiGrid}}, \code{\link{changeShoworder-ViSibook-method}}.
#' @examples 
#' intubation
#' ## Construction of the ViSiBook for the data intubation
#' vars <- c("time_in_intub","time_insert_probe","time_out_intub","delay_intub_prob") 
#' label <- c( "the blade is in the mouth", "Insertion of the tube into the mouth",
#'            "The blade is out of the mouth","Blade in the month - Tube is inserted")
#' typeA <- c( "p","p","p","l")
#' showorder <- c( 1,2,4,3)
#' deb   <- c( NA, NA,NA,"time_in_intub")
#' fin <- c( NA,NA,NA,"time_insert_probe")
#' bookdataframe <- data.frame(vars,label,typeA,showorder,deb,fin)
#' bookintubation <- ConvertoViSibook(bookdataframe)
#' plot(bookintubation)
#' ##
#' x <- buildViSiGrid( X = intubation, book = bookintubation,pixel = 5 )
#' plot( x )
setMethod( f = "plot",
           signature = "ViSigrid", 
           definition = function(x , scal.unit.tps = 10 , unit.tps = "s" , main = " " , 
                                 ncharlabel=30 , size.main = 12 ,  Fontsize.title = 11 ,
                                 Fontsize.label.Action = 11, Fontsize.label.Time = 11 , Fontsize.label.color = 9,
                                 col.main = "black" , col.grid = "grey" , colgreenzone = "green" , colblackzone = "black" ,
                                 alphainf = 0.8 , alphasup = 0.6, alphaZones = 0.2  ,
                                 vp0h	= 0.6, vp0w = 0.6, linA = 0.7 , rcircle = 15 , lwdline = 2 , lwd.grid = 1 , lty.grid = 1 
                               ) {	
             book	<- methods::slot( x , "book" )	
             sortindex <- sort( book[ ,4]  , index.return = TRUE)$ix
             if ( any( is.na( book[ , 4] )  ) ) {
               for (i in seq( 1 , sum( is.na( book[ , 4]  ) ), 1) ) {  # i =2
                 sortindex <- mapply( FUN = function(x , y )(return( if (y >= x ) { return( y + 1) }else{return( y ) } ) ) , 
                                      x = which( is.na( book[ , 4] )  )[ i ] , 
                                      y = sortindex )  
               }
             }
             ################################################
             # Definition of constants 
             inftps 	<- max( methods::slot( x , "vect_tps" ) ) 
             lgv	 	<- length( sortindex )  # Number of actions
             lgH	 	<- length( methods::slot( x , "vect_tps" ) ) - 1# Number of times
             newx	<- sapply( 	c( -0.06 , 0 , 0.06 ) , function(x ) { x * cos( seq( -pi , pi , 2 * pi / 8 ) ) } 
             ) - sapply(	c( 0 , 0.3 , 0 ) , function(x ) { -x * sin( seq( -pi , pi , 2 * pi / 8 ) ) } ) + 0.5
             newy    <- sapply(	c( -0.06 , 0 , 0.06 ) , function(x ) { x * sin( seq( -pi , pi , 2 * pi / 8 ) ) } 
             ) - sapply( c( 0 , 0.3 , 0 ) , function(x ) { x * cos( seq( -pi , pi , 2 * pi / 8 ) ) } ) + 0.5
             # Defintion of viewport and layout	
             vp0  <- grid::viewport( 	x = grid::unit( (1 - vp0w ) * 2/3 , "npc" ) , 
                                   y = grid::unit( (1 - vp0h ) * 2/3 , "npc" ) , 	
                                   width = vp0w , 
                                   height = vp0h , 
                                   just = c( 0 , 0 ) ,
                                   name = "vp0" ) # cadre
             layoutAction 	<- grid::viewport( layout = grid.layout( lgv , 1 , widths = 1 , heights = 1 ) , name = "layoutAction" )
             vplayoutA <- lapply(  sortindex , function(x )( 
               grid::viewport( layout.pos.row = which( sortindex == x ) , layout.pos.col = 1 , name = paste0( "vp" , methods::slot( book , "vars")[ x ] ) ) ))
             names(vplayoutA) = paste0( rep("vp" , lgv ) , methods::slot( book , "vars")[ seq( 1 , lgv , 1 ) ] )
             # PLoting ................................................................................................................................................
             grid.newpage()
              grid::pushViewport( vp0 )
              grid::pushViewport( layoutAction ) 
             for (ia in sortindex ) { #		ia=2
               ### Enter into the action lines
                grid::pushViewport( vplayoutA[[ which( sortindex == ia ) ]] )	
               # Punctuals treatment __________________________________________________________________________________________________________________________________________________________________	
               if (methods::slot( book , "typeA")[  ia  ] == "p" ) {
                 #.............................................................................................................................................................
                 # Green Zone 
                 if (	length( methods::slot( book , "GZDeb") ) > 0 ) {   #### If some green zone are defined in book
                   if (	is.na( methods::slot( book , "GZDeb")[ ia ] ) == FALSE ) {   ###If the action had a grren zone define
                     temp <- c( methods::slot( book , "GZDeb")[ ia ] , methods::slot( book , "GZFin")[ ia ])
                     if (temp[ 2 ] == Inf ) { 
                       temp[ 2 ] <- inftps
                     }
                     grid::pushViewport( grid::viewport( x = grid::unit( temp[ 1 ] / inftps , "npc" ) , 
                                             y = grid::unit( 0 , "npc" ) , 
                                             width = grid::unit( (temp[ 2 ] - temp[ 1 ] ) / inftps , "npc" ) , 
                                             height = grid::unit( 1 , "npc" ) , just = c( 0 , 0 ) ) )
                     grid::grid.rect( gp = gpar( col = FALSE , fill = colgreenzone , alpha = alphaZones ) )
                      grid::upViewport()
                     if (	length( methods::slot( book , "Repetition") ) > 0 ) {
                       if (	is.na( methods::slot( book , "Repetition")[ia] ) == FALSE ) { 
                       while (temp[2] + methods::slot( book , "Repetition")[ ia ] < inftps ) {
                           temp <- temp + rep( methods::slot( book , "Repetition")[ ia ] ,2 )
                            grid::pushViewport( grid::viewport( x = grid::unit( temp[ 1 ] / inftps , "npc" ) ,
                                                   y = grid::unit( 0 , "npc" ) ,
                                                   width = grid::unit( (temp[ 2 ] - temp[ 1 ] ) / inftps , "npc" ) , 
                                                   height = grid::unit( 1 , "npc" ) , just = c( 0 , 0 ) ) )
                           grid::grid.rect( gp = gpar( col = FALSE , fill = colgreenzone , alpha = alphaZones ) )
                            grid::upViewport()
                       }
                     }
                      }
                   }
                 }
                 #.............................................................................................................................................................
                 # Fin Green Zone 
                 #.............................................................................................................................................................
                 # Black Zone Zone
                 if (	length( methods::slot( book , "BZBeforeDeb" ) ) > 0 ) {   #### If some Black 1 zone are defined in book
                   if (is.na( methods::slot( book , "BZBeforeDeb" )[ ia ] ) == FALSE  ) {
                     temp <- c( methods::slot( book , "BZBeforeDeb" )[ ia ] , methods::slot( book , "BZBeforeFin" )[ ia ] )  
                     if (temp[ 2 ] == Inf ) {
                       temp[ 2 ] <- inftps
                     }
                      grid::pushViewport(grid::viewport(	x = grid::unit( temp[ 1 ] / max( inftps ) , "npc" ) , 
                                            y = grid::unit( 0, "npc" ) ,
                                            width = grid::unit( (temp[ 2 ] - temp[ 1 ] ) / inftps , "npc" ) ,
                                            height = grid::unit( 1 , "npc" ) , 
                                            just = c( 0 , 0 )))
                     grid::grid.rect( gp = gpar( col = FALSE , fill = colblackzone , alpha = alphaZones ) )
                      grid::upViewport()											
                   }	
                 }
                 if (	length( methods::slot( book , "BZAfterDeb" ) ) > 0 ) {   #### If some Black 1 zone are defined in book
                   if (is.na( methods::slot( book , "BZAfterDeb" )[ ia ] ) == FALSE  ) {
                     temp <- c( methods::slot( book , "BZAfterDeb" )[ ia ] , methods::slot( book , "BZAfterFin" )[ ia ] )  
                     if (temp[ 2 ] == Inf ) {
                       temp[ 2 ] <- inftps
                     }
                      grid::pushViewport(grid::viewport(	x = grid::unit( temp[ 1 ] / max( inftps ) , "npc" ) , 
                                            y = grid::unit( 0, "npc" ) ,
                                            width = grid::unit( (temp[ 2 ] - temp[ 1 ] ) / inftps , "npc" ) ,
                                            height = grid::unit( 1 , "npc" ) , 
                                            just = c( 0 , 0 )
                     ))
                     grid::grid.rect( gp = gpar( col = FALSE , fill = colblackzone , alpha = alphaZones ) )
                      grid::upViewport()											
                   }	
                 }
                 #  Fin Black Zone Zone 
                 #.............................................................................................................................................................
                 # Plot the punction action
                 plotpunctual( 	mat = methods::slot( x , "MATp" ) , 
                                iip = which( sortindex[ which( methods::slot( book , "typeA" )[sortindex] == "p" )] == ia ),
                                book = book, colvect = methods::slot( x , "colvect" ) , 
                                lgH = lgH , 
                                method = methods::slot( x , "parameters")$method , 
                                linA = linA 
                                )
                 # END Plot the punction action
                 #.............................................................................................................................................................
                 # Informers + Tests
                 if (is.null( methods::slot( x , "parameters" )$informer ) == FALSE) {
                   if (any( c( "global" , "join") == methods::slot( x , "parameters")$method ) ) {	
                     plotinformers( informers = methods::slot( x , "informers" ) ,  
                                     inftps = inftps , 
                                     iip = which(methods::slot( book , "vars" )[ which( methods::slot( book , "typeA" ) == "p" ) ] == methods::slot( book , "vars" )[ ia ]), 
                                     t_0 = methods::slot( x , "parameters")$t_0,
                                     alphainf = alphainf, 
                                     lwdline = lwdline , 
                                     rcircle = rcircle  )
                   }else{
                     ### Plot group 1 
                      grid::pushViewport( grid::viewport( x = grid::unit( 0 , "npc" ) , 
                                             y = grid::unit( 1/2 , "npc" ) ,
                                             width = grid::unit( 1 , "npc" ) , 
                                             height = grid::unit( linA/2 , "npc" ) ,
                                             just = c( 0 , 0 )))
                      plotinformers( informers = methods::slot( x , "informers" )[ c( 1 , 2 , 3 ) , ] ,  
                                     inftps = inftps , 
                                     iip = which(methods::slot( book , "vars" )[ which( methods::slot( book , "typeA" ) == "p" ) ] == methods::slot( book , "vars" )[ ia ]) , 
                                     t_0 = methods::slot( x , "parameters")$t_0,
                                     alphainf = alphainf, 
                                     lwdline = lwdline  , 
                                     rcircle = rcircle  )
                     grid::upViewport()
                     ### Plot group 2 
                     grid::pushViewport( grid::viewport( x = grid::unit( 0 , "npc" ) , 
                                             y = grid::unit( (1 - linA ) / 2 , "npc" ) ,
                                             width = grid::unit( 1 , "npc" ) , 
                                             height = grid::unit( linA/2 , "npc" ) ,
                                             just = c( 0 , 0 )))
                     plotinformers( informers = methods::slot( x , "informers" )[ c( 4 , 5 ,6 ) , ] ,  
                                     inftps = inftps , 
                                     iip = which(methods::slot( book , "vars" )[ which( methods::slot( book , "typeA" ) == "p" ) ] == methods::slot( book , "vars" )[ ia ]) , 
                                     t_0 = methods::slot( x , "parameters")$t_0,
                                     alphainf = alphainf, 
                                     lwdline = lwdline  , 
                                     rcircle = rcircle  )
                     grid::upViewport()
                   }
                   # Test stars
                   if ( any( c( "cut", "within" , "join") == methods::slot( x , "parameters")$method ) & methods::slot( x , "parameters" )$test == TRUE ) {	
                     if (methods::slot( x , "testsP" )[ which(methods::slot( book , "vars" )[ which( methods::slot( book , "typeA" ) == "p" ) ] == methods::slot( book , "vars" )[ ia ]) ] == TRUE ) {
                        grid::pushViewport( grid::viewport( x = grid::unit( 1 , "npc" ) , 
                                               y = grid::unit( 0.3 , "npc" ) ,
                                               width = grid::unit( 0.035 , "npc" ) , 
                                               height = grid::unit(  0.4 , "npc" ) ,
                                               just = c( 0 , 0 ) , 
                                               clip = TRUE))
                       grid.polygon(	x = t( newx ) , 
                                     y = t( newy ) ,
                                     id = NULL , 
                                     id.lengths = rep( 3 , dim( newx )[ 1 ] ) ,
                                     default.units = "npc" , 
                                     gp = gpar( col = FALSE , fill = "black" ), 
                                     draw = TRUE, 
                                     vp = NULL)
                       popViewport()
                     }
                   }
                 }
                 # END Informers + Tests
                 #.............................................................................................................................................................
                 #.............................................................................................................................................................
                 # Supplementary times points 
                 if (length( methods::slot( x , "MATpsup" ) ) > 0 ) {     
                   plotpunctualsup(X = methods::slot( x , "MATpsup" ) , 
                                   idsup = methods::slot( x , "idsup" ) , 
                                   iip = which( sortindex[ which( methods::slot( book , "typeA" )[sortindex] == "p" )] == ia ) ,
                                   book = book,
                                   method = methods::slot( x , "parameters")$method ,
                                   linA = linA , 
                                   lgH = lgH ,
                                   colvect = methods::slot( x , "colvect" ) ,
                                   alphasup = alphasup)
                 }	 
               }else{
                iipp = which(methods::slot( book , "vars")[order(book[  , 4 ] )][which(methods::slot( book , "typeA")[order(book[  , 4 ] )] == "l")] == methods::slot( book , "vars" )[ia] ) 
                 plotL(	L = methods::slot( x , "L" )[ , c( 2 * sum( methods::slot( book , "typeA" )[sortindex][ seq( 1 ,  which( sortindex == ia ) , 1) ] == "l" ) - 1 ,
                                                  2 * sum( methods::slot( book , "typeA" )[sortindex][ seq( 1 ,  which( sortindex == ia ) , 1) ] == "l" ) ) ] ,
                        idsort = methods::slot( x , "idsort" )[ ,iipp ] ,
                        inftps = inftps ,
                        group =  methods::slot( x , "group" ) ,
                        BZL = methods::slot( x , "BZL" ) ,
                        Lsup = methods::slot( x , "Lsup" ) ,
                        idsup = methods::slot( x , "idsup" ) ,
                        iip = iipp,
                        t_0 = methods::slot( x , "parameters")$t_0,
                        cols = methods::slot( x , "colvect" ) ,
                        linA = linA ,
                        alphaZones = alphaZones ,
                        alphasup = alphasup ,
                        colblackzone = colblackzone)
                 if (is.null( methods::slot( x , "parameters" )$informer ) == FALSE ) {
                   iip  <- 	sum(methods::slot( book , "typeA") == "p") + sum(methods::slot( book , "typeA")[ sortindex ] == "l" ) + iipp
                   iip2 <-  which( colnames(methods::slot( x , "informers")) == methods::slot( book , "deb")[ ia ] )
                   if (any(is.na(c(methods::slot( x , "informers" )[ , iip ] , methods::slot( x , "informers" )[ , iip2 ]))) == FALSE ) {   
                     if ( is.character( methods::slot( x , "parameters")$t_0)==FALSE & methods::slot( x , "parameters")$t_0 != 0 ){
                       methods::slot( x , "informers" )[ , iip ] <- methods::slot( x , "informers" )[ , iip ] - rep(methods::slot( x , "parameters")$t_0, length(methods::slot( x , "informers" )[ , iip ]))
                       methods::slot( x , "informers" )[ , iip2 ] <- methods::slot( x , "informers" )[ , iip2 ] - rep(methods::slot( x , "parameters")$t_0, length(methods::slot( x , "informers" )[ , iip2 ]))
                      }
                   if (any( c( "global",  "join") == methods::slot( x , "parameters")$method ) ) {
                      grid::pushViewport( viewport( x = grid::unit(  methods::slot( x , "informers" )[ 2 , iip2 ] / inftps , "npc" ) , 
                                             y = grid::unit( 0 , "npc" ) ,
                                             width = grid::unit(  (methods::slot( x , "informers" )[ 2 , iip ] -  methods::slot( x , "informers" )[ 2 , iip2 ] ) / inftps , "npc" ) , 
                                             height = grid::unit( 1 , "npc" ) , 
                                             just = c( 0 , 0)))
                      grid::grid.lines( x = grid::unit( c( 0 , 1 ) , "npc" ) , y = grid::unit( c( 1/2 , 1/2 ) , "npc" ) , default.units = "npc" , arrow = NULL , gp = gpar( col = "black" , lwd = lwdline))
                      grid::upViewport()	
                   }else{
                     ### Plot group 1 
                      grid::pushViewport(grid::viewport( x = grid::unit( 0 , "npc" ) , 
                                             y = grid::unit( 1/2 , "npc" ) ,
                                             width = grid::unit( 1 , "npc" ) , 
                                             height = grid::unit( linA/2 , "npc" ) ,
                                             just = c( 0 , 0 )))
                        grid::pushViewport(grid::viewport( x = grid::unit(  methods::slot( x , "informers" )[ 2 , iip2 ] / inftps , "npc" ) , 
                                             y = grid::unit( 0 , "npc" ) ,
                                             width = grid::unit( (methods::slot( x , "informers" )[ 2 , iip ] -  methods::slot( x , "informers" )[ 2 , iip2 ] ) / inftps , "npc" ) , 
                                             height = grid::unit( 1 , "npc" ) , 
                                             just = c( 0 , 0) ))
                          grid::grid.lines( x = grid::unit( c( 0 , 1 ) , "npc" ) , y = grid::unit( c( 1/2 , 1/2 ) , "npc" ) , default.units = "npc" , arrow = NULL , gp = gpar( col = "black" , lwd = lwdline))
                        grid::upViewport()	
                      grid::upViewport()
                     ### Plot group 2 
                      grid::pushViewport(grid::viewport( x = grid::unit( 0 , "npc" ) , 
                                             y = grid::unit( (1 - linA ) / 2 , "npc" ) ,
                                             width = grid::unit( 1 , "npc" ) , 
                                             height = grid::unit( linA/2 , "npc" ) ,
                                             just = c( 0 , 0 )))
                        grid::pushViewport(grid::viewport( x = grid::unit(  methods::slot( x , "informers" )[ 5 , iip2 ] / inftps , "npc" ) , 
                                             y = grid::unit( 0 , "npc" ) ,
                                             width = grid::unit(  (methods::slot( x , "informers" )[ 5 , iip ] -  methods::slot( x , "informers" )[ 5 , iip2 ] ) / inftps , "npc" ) , 
                                             height = grid::unit( 1 , "npc" ) , 
                                             just = c( 0 , 0)))
                          grid::grid.lines( x = grid::unit( c( 0 , 1 ) , "npc" ) , y = grid::unit( c( 1/2 , 1/2 ) , "npc" ) , default.units = "npc" , arrow = NULL , gp = gpar( col = "black" , lwd = lwdline))
                        grid::upViewport()
                      grid::upViewport()
                   }
                   if ( any( c( "cut", "within" , "join") == methods::slot( x , "parameters")$method ) & methods::slot( x , "parameters" )$test == TRUE ) {	
                     if (methods::slot( x , "testsP" )[ sum( methods::slot( book , "typeA" ) == "p" ) + sum( methods::slot( book , "typeA" )[sortindex][ seq( 1 ,  which( sortindex == ia ) , 1) ] == "l" ) ] == TRUE ) {
                        grid::pushViewport(grid::viewport( x = grid::unit( 1 , "npc" ) , 
                                               y = grid::unit( 0.3 , "npc" ) ,
                                               width = grid::unit( 0.035 , "npc" ) , 
                                               height = grid::unit(  0.4 , "npc" ) ,
                                               just = c( 0 , 0 ) , 
                                               clip = TRUE))
                           grid::grid.polygon( x = t( newx ) , 
                                         y = t( newy ) ,
                                         id = NULL , 
                                         id.lengths = rep( 3 , dim( newx )[ 1 ] ) ,
                                         default.units = "npc" , 
                                         gp = gpar( col = FALSE , fill = "black" ), 
                                         draw = TRUE, 
                                         vp = NULL)
                        popViewport()
                     }
                   }
                   }
                 }
               }
                grid::upViewport() # Out of the line of action
             }
              grid::upViewport() # Out of the line of action
              grid::upViewport() # Out of the cadre
             #################################################### Legends ###################################################
             #Legend Axe Y
             #___________________________________________________________________________________________________________________
             # Legend Title graph
              grid::pushViewport(grid::viewport( x = grid::unit( (1 - vp0w ) * 2/3 , "npc" ) ,
                                     y = grid::unit( (1 - vp0h ) * 2/3 + vp0h , "npc") ,
                                     width = grid::unit( vp0w , "npc" ) ,
                                     height = grid::unit( (1 - vp0h ) * 0.5 , "npc" ) , 
                                     just = c( 0 , 0 ) ))
                grid::grid.text( main , gp = gpar( fontsize = size.main , col = col.main) )
              grid::upViewport()
             #___________________________________________________________________________________________________________________
             # Legend names action 
              grid::pushViewport(grid::viewport( 	x = grid::unit( 0 , "npc" ) , 
                                      y = grid::unit( (1 - vp0h ) * 2/3 , "npc" ) , 	
                                      width = (1 - vp0w ) * 2/3  , 
                                      height = vp0h , 
                                      just = c( 0 , 0 ) ,
                                      name = "vp0" ))
                grid::pushViewport( layoutAction )
                 for (ii in sortindex ) { 
                    grid::pushViewport( vplayoutA[[ which( sortindex == ii ) ]] )	
                      grid::grid.text( 	substr( as.character(methods::slot( book , "label" )[ ii ] ), 1 , ncharlabel ) ,
                               rot = 0.5 , 
                               gp = gpar( col = "black" , fontsize = Fontsize.label.Action , fontface = "plain" ))
                    grid::upViewport()
               }
                grid::upViewport()
              grid::upViewport()
             #___________________________________________________________________________________________________________________
             #___________________________________________________________________________________________________________________
             # Grid times
              grid::pushViewport(vp0)
                grid::grid.grill( h = grid::unit(	seq( 1 / (lgv ) , 1 - 1 / lgv , 1 / lgv ) ,"npc" ) ,
                         v = grid::unit( seq( 1 ,floor( inftps / scal.unit.tps ) , 1 ) /inftps* scal.unit.tps  , "npc" ),
                         gp = gpar( col = col.grid , lwd = lwd.grid , lty = lty.grid )
                )
              grid::grid.lines( x = grid::unit(c( 0 , 1 ), "npc"),
                         y = grid::unit( c( 0 ) , "npc" ) ,
                         default.units = "npc" ,
                         arrow = arrow( angle = 20 , length = grid::unit( 0.10 , "inches" ) ) ,
                         gp = gpar( col = "black" , lwd = 2 ) )
              grid::grid.lines( x = grid::unit( c( 0 , 0 ) , "npc" ) ,
                         y = grid::unit( c( 0 , 1 ) , "npc" ) ,
                         default.units = "npc" ,
                         gp = gpar( col = "black" , lwd = 2 ))
             grid::upViewport() # Out of the cadre
             if (vp0w / floor( inftps / scal.unit.tps ) < 0.03 ) {
                grid::pushViewport(grid::viewport( grid::unit( x = (1 - vp0w ) * 2 / 3 , "npc" ) , 
                                       y = grid::unit( (1 - vp0h ) * 2 / 3 , "npc" ) , 
                                       width = grid::unit( 0.03 , "npc" ) ,
                                       height = grid::unit( (1 - vp0h ) * 1 / 6 , "npc" ) , just = c( 1 , 1 )))
               grid::grid.text( 	x = grid::unit( 0.1 , "npc" ) ,
                           y = grid::unit( 0.5 , "npc" ) , paste( methods::slot( x , "vect_tps" )[ 1 ] , "s" ) ,
                           gp = gpar( col = "black" , fontsize = Fontsize.label.Time ) , just = c( 0 , 0 ))
               grid::grid.lines( x = grid::unit( c( 1 , 1 ) , "npc" ) , 
                           y = grid::unit( c( 0.5 , 1 ) , "npc" ) , 
                           default.units = "npc" ,
                           gp = gpar( col = "black" ) )
                grid::upViewport()   
             }else{
                grid::pushViewport(grid::viewport( x = grid::unit( (1 - vp0w ) * 2 / 3 , "npc" ) ,
                                       y = grid::unit( (1 - vp0h ) * 2 / 3 , "npc" ) ,
                                       width = grid::unit( 0.03	, "npc" ) ,
                                       height = grid::unit( (1 - vp0h ) * 1 / 6 , "npc" ) , 
                                       just = c( 0 , 1)))
               grid::grid.text( 	x = grid::unit( 0.1 , "npc" ) , 
                           y = grid::unit( 0.5 , "npc" ) , 
                           paste( methods::slot( x , "vect_tps" )[ 1 ] , unit.tps ) ,
                           gp = gpar( col = "black" , fontsize = Fontsize.label.Time ) ,
                           just = c( 0 , 0 ))
               grid::grid.lines( x = grid::unit( c( 0 , 0 ) , "npc" ) ,
                           y = grid::unit( c( 0.5 , 1 ) , "npc" ) ,
                           default.units = "npc" ,
                           gp = gpar( col = "black" ))
                grid::upViewport()
             }
             #___________________________________________________________________________________________________________________
              grid::pushViewport(grid::viewport( x = grid::unit( (1 - vp0w ) * 2 / 3 + vp0w /floor( inftps / scal.unit.tps ) , "npc" ) ,
                                     y = grid::unit( (1 - vp0h ) * 2 / 3 , "npc" ) ,
                                     width = grid::unit( 0.03 + 0.01 * nchar( scal.unit.tps ) , "npc" ) ,
                                     height = grid::unit( (1 - vp0h ) / 6 , "npc" ) ,
                                     just = c( 0 , 1 ) 
             )
             )
                 grid::grid.text( 	x = grid::unit( 0.1 , "npc" ) ,
                             y = grid::unit( 0.5 , "npc" ) ,
                             paste( methods::slot( x , "vect_tps" )[ 1 ] + scal.unit.tps , unit.tps ) ,
                             gp = gpar( fontsize = Fontsize.label.Time ) ,
                             just = c( 0 , 0)
                 )
                 grid::grid.lines( x = grid::unit( c( 0 , 0 ) , "npc" ) ,
                             y = grid::unit( c( 0.5 , 1 ) , "npc" ) ,
                             default.units = "npc" ,
                             gp = gpar( col = "black" ) )
              grid::upViewport()
             #___________________________________________________________________________________________________________________
              grid::pushViewport(grid::viewport( x = grid::unit( vp0w * (floor( inftps / scal.unit.tps ) - 1 ) / floor( inftps / scal.unit.tps ) + (1 - vp0w ) * 2 / 3 , "npc" ) ,
                                     y = grid::unit( (1 - vp0h ) * 2 / 3 , "npc" ) ,
                                     width = grid::unit( vp0w / floor( inftps / scal.unit.tps ) , "npc" ) ,
                                     height = grid::unit( (1 - vp0h ) / 6 , "npc" ) , 
                                     just = c( 0 , 1 )))
               grid::grid.text( paste0( "Times (" , unit.tps , ")" ) ,
                          gp = gpar( fontsize = Fontsize.label.Time ) , just = c( 1 , 0 ))
              grid::upViewport()
             #________________________________________________________________________________________________________________
             #___________________________________________________________________________________________________________________
             ####Legend bottom
             ###  &  any( methods::slot( book , "typeA" )[ sortindex ] == "l")
             lcol <- switch( as.character( length( methods::slot( x , "colvect" )[ 1 , ] ) == 1), "TRUE" = 1 ,"FALSE" = length( methods::slot( x , "colvect" )[ 1 , ] ) %/% 2)
             if (any( methods::slot( book , "typeA" )[ sortindex ] == "p" )) {
               part <- switch( 	as.character( any( c( length( methods::slot( book , "BZBeforeDeb" ) ) , length( methods::slot( book , "BZAfterDeb" ) ) ) > 0 ) ) , "TRUE" = 1 , "FALSE" = 0 ) + switch( as.character( length( methods::slot( book , "GZDeb" ) ) > 0 ) , "TRUE" = 1 , "FALSE" = 0 ) + switch( as.character( is.null(methods::slot( x , "parameters" )$informer )) ,"FALSE" = 1 , "TRUE" = 0 )
               if (part == 0 ) {
                  grid::pushViewport(grid::viewport( y = grid::unit( 0.005 , "npc"),
                                         x = grid::unit( 0.04 , "npc" ) , 
                                         height = grid::unit( (1 - vp0h ) * 1 / 3 + 0.02, "npc" ) , 
                                         width = grid::unit( 0.9 , "npc" ) ,
                                         just = c( 0 , 0 )))
                  grid::grid.rect( gp = gpar( col = "black" , fill = "grey",alpha = 0.3 ))	
                  grid::upViewport()		
                  grid::pushViewport(grid::viewport( y = grid::unit( 0.02 , "npc") , 
                                         x = grid::unit( 0.05 , "npc" ) ,
                                         height = grid::unit( (1 - vp0h ) * 1 / 3 , "npc" ) ,
                                         width = grid::unit( 0.4 , "npc" ),
                                         just = c( 0 , 0 ) ) )
                   grid::pushViewport(grid::viewport( x = grid::unit( 0 , "npc" ) , 
                                          y = grid::unit( 5 / 6 , "npc" ) , 
                                          width = grid::unit( 1 , "npc" ) , 
                                          height = grid::unit( 1/6 , "npc" ) , 
                                          just = c( 0 , 0 ) ) )
                    grid::grid.text( x = grid::unit( 0 , "npc" ) , y = grid::unit( 1, "npc" ), "Punctual Action", gp = gpar( col = "black" , fontsize = Fontsize.label.color ) , just = c( 0 , 1 ) )
                   grid::upViewport()
                   grid::pushViewport(grid::viewport( x = grid::unit( 0 , "npc" ) , 
                                          y = grid::unit( 0 , "npc" ) , 
                                          width = grid::unit( 1 , "npc" ) , 
                                          height = grid::unit( 5 / 6 , "npc" ) , 
                                          just = c( 0 , 0 ) ) )
                  if ( methods::slot( x , "parameters" )$method == "global" ) {
                    grid::pushViewport(grid::viewport(	x = grid::unit( 0 , "npc" ) ,  
                                           y = grid::unit( 8/12 , "npc" ) , 
                                           width = grid::unit( 1/3 , "npc" ) ,
                                           height = grid::unit( 1/3 , "npc" ) , 
                                           just = c( 0 , 0 ) ) )
                    grid::grid.text(	x = grid::unit( 0 , "npc" ) , y = grid::unit( 1/2  , "npc" ) , paste( switch( methods::slot( x , "parameters" )$quantity , "ind" = "N : 1 ->" , "dens" = " % : 1 -> ") , length(   methods::slot( x , "colvect")[ 1 , ] ) ) , gp = gpar( col = "black" , fontsize = Fontsize.label.color ) , just = c( 0 , 1/2 ) )
                    grid::upViewport()
                    grid::pushViewport(grid::viewport( x = grid::unit( 7/12 , "npc" ) ,  
                                           y = grid::unit( 9/12 , "npc" ) , 
                                           width = grid::unit( 1/3 , "npc" ) , 
                                           height = grid::unit( 2/12 , "npc" ) ,
                                           just = c( 0 , 0 )))
                     grid::pushViewport(grid::viewport( layout = grid.layout( 1 , length(   methods::slot( x , "colvect")[ 1 , ] ) , widths = 1 , heights = 1 ) ) )
                     for (ii in seq( 1 , length( methods::slot( x , "colvect")[ 1 , ] ) , 1 ) ) {
                        grid::pushViewport(grid::viewport( layout.pos.col = ii , layout.pos.row = 1 ) )
                        grid::grid.rect( gp = gpar( col = FALSE , fill = methods::slot( x , "colvect")[ 1 , ii ] ) )
                        grid::upViewport()
                     }
                     grid::upViewport()					
                    grid::upViewport()
                  }else{                           	
                    grid::pushViewport(grid::viewport(	x = grid::unit( 0 , "npc" ) ,  
                                           y = grid::unit(8/12 , "npc" ) , 
                                           width = grid::unit( 1/3 , "npc" ) , 
                                           height = grid::unit( 1/3 , "npc" ) ,
                                           just = c( 0 , 0 )  )  )
                    grid::grid.text( x = grid::unit( 0 , "npc" ) , y = grid::unit( 1 / 2 , "npc" ) , paste( switch( methods::slot( x , "parameters" )$quantity , "ind" = "N : 1 ->" , "dens" = " % : 1 -> ") , length(   methods::slot( x , "colvect")[ 1 , ] ) ) , gp = gpar( col = "black" , fontsize = Fontsize.label.color ) ,  just = c( 0 , 1/2 ) )
                    grid::upViewport()
                    grid::pushViewport(grid::viewport( x = grid::unit( 0 , "npc" ) , 
                                           y = grid::unit(4/12 , "npc" ) , 
                                           width = grid::unit( 1/3 , "npc" ) , 
                                           height = grid::unit( 1/3 , "npc" ), 
                                           just = c( 0 , 0 )))
                    grid::grid.text( 	x = grid::unit( 0 , "npc" ) , y = grid::unit( 1 / 2 , "npc" ) ,  paste( "  Group" , levels( methods::slot( x , "group" ) )[ 1 ] ) , gp = gpar( col = "black" , fontsize = Fontsize.label.color ) , just = c( 0 , 1/2 ) )
                    grid::upViewport()
                    grid::pushViewport(grid::viewport( x = grid::unit( 0 , "npc" ) ,  
                                           y = grid::unit(0/12 , "npc" ) , 
                                           width = grid::unit( 1/3 , "npc" ) , 
                                           height = grid::unit( 1/3 , "npc" ) ,
                                           just = c( 0 , 0 )))
                    grid::grid.text( x = grid::unit( 0 , "npc" ) , y = grid::unit( 1 / 2 , "npc" ) ,  paste( "  Group" , levels( methods::slot( x , "group" ) )[ 2 ] ) , gp = gpar( col = "black" , fontsize = Fontsize.label.color ) , just = c( 0 , 1/2 )  )
                    grid::upViewport()
                    grid::pushViewport(grid::viewport( x = grid::unit( 7/12 , "npc" ) , 
                                           y = grid::unit( 5/12 , "npc" ) ,  
                                           width = grid::unit( 1/3 , "npc" ) , 
                                           height = grid::unit( 2/12 , "npc" ) , 
                                           just = c( 0 , 0 )))
                     grid::pushViewport(grid::viewport( layout = grid.layout( 1 , length(   methods::slot( x , "colvect")[ 1 , ] ) , widths = 1 , heights = 1 ) ) )
                      for (ii in seq( 1 , length(   methods::slot( x , "colvect")[ 1 , ] ) , 1 ) ) {
                         grid::pushViewport(grid::viewport( layout.pos.col = ii , layout.pos.row = 1 ) )
                          grid::grid.rect( gp = gpar( col = FALSE , fill = methods::slot( x , "colvect")[ 1 , ii ] ) )
                         grid::upViewport()
                      }
                     grid::upViewport()					
                    grid::upViewport()
                    grid::pushViewport(grid::viewport( x = grid::unit( 7/12 , "npc" ) ,
                                           y = grid::unit( 1/12 , "npc" ) , 
                                           width = grid::unit( 1/3 , "npc" ) ,
                                           height = grid::unit( 2/12 , "npc" ) ,
                                           just = c( 0 , 0 )))
                     grid::pushViewport(grid::viewport( layout = grid.layout( 1 , length(   methods::slot( x , "colvect")[ 1 , ] ) , widths = 1 , heights = 1 ) ) )
                      for (ii in seq( 1 , length(   methods::slot( x , "colvect")[ 1 , ] ) , 1 ) ) {
                         grid::pushViewport(grid::viewport( layout.pos.col = ii , layout.pos.row = 1 ) )
                          grid::grid.rect( gp = gpar( col = FALSE , fill = methods::slot( x , "colvect")[ 2 , ii ] ) )
                         grid::upViewport()
                      }
                     grid::upViewport()					
                    grid::upViewport()
                 }
                  grid::upViewport()
                  grid::upViewport()
               }else{
                 part = 3
                  grid::pushViewport(grid::viewport( y = grid::unit( 0.005 , "npc"),
                                         x = grid::unit( 0.04 , "npc" ) , 
                                         height = grid::unit( (1 - vp0h ) * 1 / 3 + 0.02, "npc" ) , 
                                         width = grid::unit( 0.9 , "npc" ) ,
                                         just = c( 0 , 0 )))
                 grid::grid.rect( gp = gpar( col = "black" , fill = "grey",alpha = 0.3 ))	
                  grid::upViewport()		
                  grid::pushViewport(grid::viewport( y = grid::unit( 0.02 , "npc") , x = grid::unit( 0.05 , "npc" ) ,  height = grid::unit( (1 - vp0h ) * 1 / 3 , "npc" ) ,  width = grid::unit( 0.4 , "npc" ) , just = c( 0 , 0 ) ) )
                  grid::pushViewport(grid::viewport( x = grid::unit( 0 , "npc" ) , y = grid::unit( 5 / 6 , "npc" ) , width = grid::unit( 1 , "npc" ) , height = grid::unit( 1/6 , "npc" ) , just = c( 0 , 0 ) ) )
                 grid::grid.text( 	x = grid::unit( 0 , "npc" ) , y = grid::unit( 1, "npc" ), "Punctual Action", gp = gpar( col = "black" , fontsize = Fontsize.label.color ) , just = c( 0 , 1 ) )
                  grid::upViewport()
                  grid::pushViewport(grid::viewport( x = grid::unit( 0 , "npc" ) , y = grid::unit( 0 , "npc" ) , width = grid::unit( 1 , "npc" ) , height = grid::unit( 5 / 6 , "npc" ) , just = c( 0 , 0 ) ) )
                  grid::pushViewport(grid::viewport( layout = grid.layout( 1 , 2 ) , width = 1 , height = 1 ) )
                  grid::pushViewport(grid::viewport( layout.pos.col = 2 , layout.pos.row = 1 ) )
                 if (length( methods::slot( book , "GZDeb" ) ) > 0 ) {
                    grid::pushViewport(grid::viewport( x = grid::unit( 0 , "npc" ) , y = grid::unit( (part - 1 ) / part , "npc" ) , width = grid::unit( 1 , "npc" ) , height = grid::unit( 1 / part , "npc" ) , just = c( 0 , 0 ) ) )
                   grid::grid.text(  x = grid::unit( 0 , "npc" ) , y = grid::unit( 0.5 , "npc" ) , "Green Zone" , gp = gpar( col = "black" , fontsize = Fontsize.label.color ) , just = c( 0 , 1 / 2 ))
                   grid::grid.circle( 	x = grid::unit( 5 / 6 , "npc" ) ,   y = grid::unit( 0.4 , "npc" ) , r = grid::unit( 0.4 , "npc" ) , gp = gpar( col = FALSE , fill = colgreenzone , alpha = alphaZones ) )
                    grid::upViewport() 
                 }
                 if (any( c( length( methods::slot( book , "BZBeforeDeb" ) ) , length( methods::slot( book , "BZAfterDeb" ) ) ) > 0 ) ) {
                    grid::pushViewport(grid::viewport( x = grid::unit( 0 , "npc" ) , y = grid::unit( (part - 1 - switch( as.character(length( methods::slot( book , "GZDeb" ) ) > 0 ) , "TRUE" = 1 , "FALSE" = 0 ) ) / part , "npc" ) , width = grid::unit( 1 , "npc" ) , height = grid::unit( 1 / part , "npc" ) , just = c( 0 , 0 )  ) )
                   grid::grid.text( 	x = grid::unit( 0 , "npc" ) , y = grid::unit( 0.5 , "npc" ) , "Black Zone" , gp = gpar( col = "black" , fontsize = Fontsize.label.color ) , just = c( 0 , 1 / 2 ) )
                   grid::grid.circle( 	x = grid::unit( 5 / 6 , "npc" ) , y = grid::unit( 0.4 , "npc" ) , r = grid::unit( 0.4 , "npc" ) , gp = gpar( col = FALSE , fill = colblackzone , alpha = alphaZones ) )
                    grid::upViewport()
                 }
                 if (is.null( methods::slot( x , "parameters" )$informer ) == FALSE ) {
                    grid::pushViewport(grid::viewport( x = grid::unit( 0 , "npc" ) , y = grid::unit( (part - 1 - 
                                                                                 switch( 	as.character( any( c( length( methods::slot( book , "BZBeforeDeb" ) ) ,
                                                                                                                length( methods::slot( book , "BZAfterDeb" ) ) ) > 0 ) ) , 
                                                                                          "TRUE" = 1 , "FALSE" = 0 ) - switch( as.character(length( methods::slot( book , "GZDeb" ) ) > 0 ) , 
                                                                                                                               "TRUE" = 1 ,"FALSE" = 0 ) ) / part , "npc" ) ,
                                           width = grid::unit( 1 , "npc" ) , height = grid::unit( 1 / part , "npc" ) , just = c( 0 , 0 )  ) )
                   grid::grid.text(	x = grid::unit( 0 , "npc" ) , y = grid::unit( 0.5 , "npc" ) , paste( switch( methods::slot( x , "parameters" )$informer , "mean" = "Mean +/- sd" , "median" = "Median Q1-Q3" ) ) ,  gp = gpar( col = "Black" , fontsize = Fontsize.label.color ) ,  just = c( 0 , 1/2 ) )
                   grid::grid.lines( x = grid::unit( c( 4/6 , 6/6 ) , "npc" ) , y = grid::unit( c( 0.4 , 0.4 ) , "npc" ), default.units = "npc" ,  arrow = NULL ,   gp = gpar( col = "black" , lwd = 1 ) )
                   grid::grid.circle( x = grid::unit( 4 / 6 , "npc" ) , y = grid::unit( 0.4 , "npc" ) , r = grid::unit( 0.2 , "npc" ) , gp = gpar( col = "black" , fill = "white" ) )
                   grid::grid.circle( x = grid::unit( 5 / 6 , "npc" ) , y = grid::unit( 0.4 , "npc" ) , r = grid::unit( 0.2 , "npc" ) , gp = gpar( col = "black" , fill = "white" ) )
                   grid::grid.circle( x = grid::unit( 6 / 6 , "npc" ) , y = grid::unit( 0.4 , "npc" ) , r = grid::unit( 0.2 , "npc" ) , gp = gpar( col = "black" , fill = "white" ) )
                    grid::upViewport()
                 }
                  grid::upViewport()
                  grid::pushViewport(grid::viewport( layout.pos.col = 1 , layout.pos.row = 1 ))
                 if ( methods::slot( x , "parameters" )$method == "global" ) {
                    grid::pushViewport(grid::viewport(	x = grid::unit( 0 , "npc" ) ,  y = grid::unit( 8/12 , "npc" ) , width = grid::unit( 1/3 , "npc" ) ,height = grid::unit( 1/3 , "npc" ) , just = c( 0 , 0 ) ) )
                   grid::grid.text( 	x = grid::unit( 0 , "npc" ) , y = grid::unit( 1/2  , "npc" ) , paste( switch( methods::slot( x , "parameters" )$quantity , "ind" = "N : 1 ->" , "dens" = " % : 1 -> ") , length(   methods::slot( x , "colvect")[ 1 , ] ) ) , gp = gpar( col = "black" , fontsize = Fontsize.label.color ) , just = c( 0 , 1/2 ) )
                    grid::upViewport() 
                    grid::pushViewport(grid::viewport( x = grid::unit( 7/12 , "npc" ) ,  y = grid::unit( 9/12 , "npc" ) , width = grid::unit( 1/3 , "npc" ) , height = grid::unit( 2/12 , "npc" ) ,  just = c( 0 , 0 )  )  )
                    grid::pushViewport(grid::viewport( layout = grid.layout( 1 , length(   methods::slot( x , "colvect")[ 1 , ] ) , widths = 1 , heights = 1 ) ) )
                   for (ii in seq( 1 , length(   methods::slot( x , "colvect")[ 1 , ] ) , 1 )) {
                      grid::pushViewport(grid::viewport( layout.pos.col = ii , layout.pos.row = 1 ) )
                     grid::grid.rect( gp = gpar( col = FALSE , fill = methods::slot( x , "colvect")[ 1 , ii ] ) )
                      grid::upViewport()
                   }
                    grid::upViewport()					
                    grid::upViewport()
                }else{  
                    grid::pushViewport(grid::viewport(	x = grid::unit( 0 , "npc" ) ,  y = grid::unit(8/12 , "npc" ) , width = grid::unit( 1/3 , "npc" ) , height = grid::unit( 1/3 , "npc" ) ,just = c( 0 , 0 )  )  )
                   grid::grid.text( 	x = grid::unit( 0 , "npc" ) , y = grid::unit( 1 / 2 , "npc" ) , paste( switch( methods::slot( x , "parameters" )$quantity , "ind" = "N : 1 ->" , "dens" = " % : 1 -> ") , length(   methods::slot( x , "colvect")[ 1 , ] ) ) , gp = gpar( col = "black" , fontsize = Fontsize.label.color ) ,  just = c( 0 , 1/2 ) )
                    grid::upViewport()
                    grid::pushViewport(grid::viewport(	x = grid::unit( 0 , "npc" ) , y = grid::unit(4/12 , "npc" ) , width = grid::unit( 1/3 , "npc" ) , height = grid::unit( 1/3 , "npc" ) , just = c( 0 , 0 ) ) )
                   grid::grid.text( 	x = grid::unit( 0 , "npc" ) , y = grid::unit( 1 / 2 , "npc" ) ,  paste( "  Group" , levels( methods::slot( x , "group" ) )[ 1 ] ) , gp = gpar( col = "black" , fontsize = Fontsize.label.color ) , just = c( 0 , 1/2 ) )
                    grid::upViewport()
                    grid::pushViewport(grid::viewport(	x = grid::unit( 0 , "npc" ) ,  y = grid::unit(0/12 , "npc" ) , width = grid::unit( 1/3 , "npc" ) , height = grid::unit( 1/3 , "npc" ) ,  just = c( 0 , 0 ) )  )
                   grid::grid.text( 	x = grid::unit( 0 , "npc" ) , y = grid::unit( 1 / 2 , "npc" ) ,  paste( "  Group" , levels( methods::slot( x , "group" ) )[ 2 ] ) , gp = gpar( col = "black" , fontsize = Fontsize.label.color ) , just = c( 0 , 1/2 )  )
                    grid::upViewport()
                    grid::pushViewport(grid::viewport( x = grid::unit( 7/12 , "npc" ) , y = grid::unit( 5/12 , "npc" ) ,  width = grid::unit( 1/3 , "npc" ) , height = grid::unit( 2/12 , "npc" ) , just = c( 0 , 0 )  ) )
                    grid::pushViewport(grid::viewport( layout = grid.layout( 1 , length(   methods::slot( x , "colvect")[ 1 , ] ) , widths = 1 , heights = 1 ) ) )
                   for (ii in seq( 1 , length(   methods::slot( x , "colvect")[ 1 , ] ) , 1 ) ) {
                      grid::pushViewport(grid::viewport( layout.pos.col = ii , layout.pos.row = 1 ) )
                     grid::grid.rect( gp = gpar( col = FALSE , fill = methods::slot( x , "colvect")[ 1 , ii ] ))
                      grid::upViewport()
                   }
                    grid::upViewport()					
                    grid::upViewport()
                    grid::pushViewport( grid::viewport( x = grid::unit( 7/12 , "npc" ) ,  y = grid::unit( 1/12 , "npc" ) , width = grid::unit( 1/3 , "npc" ) ,  height = grid::unit( 2/12 , "npc" ) , just = c( 0 , 0 ) ) )
                    grid::pushViewport( grid::viewport( layout = grid.layout( 1 , length(   methods::slot( x , "colvect")[ 1 , ] ) , widths = 1 , heights = 1 ) ) )
                   for (ii in seq( 1 , length(   methods::slot( x , "colvect")[ 1 , ] ) , 1 ) ) {
                      grid::pushViewport( grid::viewport( layout.pos.col = ii , layout.pos.row = 1 ) )
                     grid::grid.rect( gp = gpar( col = FALSE , fill = methods::slot( x , "colvect")[ 2 , ii ] ) )
                      grid::upViewport()
                   }
                    grid::upViewport()					
                    grid::upViewport()
                 }
                  grid::upViewport()
                  grid::upViewport()
                  grid::upViewport()
                 grid::upViewport()
               }
             ###### Long action..................................................................................................................................
              grid::pushViewport( grid::viewport( y = grid::unit( 0.02 , "npc") ,   x = grid::unit( 0.5 , "npc" ) ,  height = grid::unit( (1 - vp0h ) * 1/3 , "npc" ) , width = grid::unit( 0.4 , "npc" ) ,  just = c( 0 , 0 ) ) )
              ##########
               grid::pushViewport( grid::viewport( x = grid::unit( 0 , "npc" ) ,  y = grid::unit( 5/6 , "npc" ) ,  width = grid::unit( 1 , "npc" ) , height = grid::unit( 1/6 , "npc" ) ,  just = c( 0 , 0 )  ))
               grid::grid.text(  x = grid::unit( 0 , "npc" ) ,  y = grid::unit( 1 , "npc" ) ,  "Long Action" ,  gp = gpar( col = "black" , fontsize = Fontsize.label.color ) , just = c( 0 , 1 )  )
                grid::upViewport()
                grid::pushViewport( grid::viewport( x = grid::unit( 0 , "npc" ) ,  y = grid::unit( 0 , "npc" ) ,  width = grid::unit( 1 , "npc" ) ,  height = grid::unit( 5/6 , "npc" ) , just = c( 0 , 0 )  )  )
               part = switch( methods::slot( x , "parameters" )$method , "global" = 1 , 2 ) + switch( as.character( is.null( methods::slot( x , "parameters" )$informer ) == FALSE ) , "TRUE" = 1 , "FALSE" = 0 )
               if (is.null( methods::slot( x , "parameters" )$informer ) == FALSE ) {
                  grid::pushViewport( grid::viewport( x = grid::unit( 0 , "npc" ) ,   y = grid::unit( 0 , "npc" ) , width = grid::unit( 1 , "npc" ) ,  height = grid::unit( 1/part , "npc" ) , just = c( 0 , 0 ) )   )
                 grid::grid.text( x = grid::unit( 0 , "npc" ) , paste( switch( 	methods::slot( x , "parameters" )$informer , "mean" = "Mean" ,  "median" = "Median" ), "of the span" ) , gp = gpar( col = "black" , fontsize = Fontsize.label.color ) ,  just = c( 0 , 1/2 )  )
                 grid::grid.lines( x = grid::unit( c( 0.5 , 0.9 ) , "npc" ) , y = grid::unit( c( 0.4 , 0.4 ) , "npc" ) ,  default.units = "npc" ,  arrow = NULL,   gp = gpar( col = "black" , lwd = 2 )  )
                  grid::upViewport()
               }
                grid::pushViewport( grid::viewport( x = grid::unit( 0 , "npc" ) , y = grid::unit( switch( as.character( is.null( methods::slot( x , "parameters" )$informer ) == FALSE ) , "TRUE" = 1 , "FALSE" = 0 ) / part , "npc" ), width = grid::unit( 1 , "npc" ) ,   height = grid::unit( 1/part , "npc" ) , just = c( 0 , 0 )  ) )
               if (methods::slot( x , "parameters"  )$method != "global"  ) {		
                  grid::pushViewport( grid::viewport( layout = grid.layout( 1 , 3 ) , width = 1 , height = 1 ) )
                 if ( length( methods::slot( book , "BZLong" ) )  > 0 ) {
                    grid::pushViewport( grid::viewport( layout.pos.col = 3 , layout.pos.row = 1 ) )
                   grid::grid.text( x = grid::unit( 0 , "npc" ) , "Not in time" , gp = gpar( col = "black" , fontsize = Fontsize.label.color ) , just = c( 0 , 1/2 ) )
                   grid::grid.circle( 	x = grid::unit( 0.8 , "npc" ) ,   r = grid::unit( 0.4 , "npc" ) ,  gp = gpar( 	col = FALSE , fill = methods::slot( x , "colvect" )[ 2 , lcol ]    ) )
                   grid::grid.circle( 	x = grid::unit( 0.8 , "npc" ) ,  r = grid::unit( 0.4 , "npc" ) ,  gp = gpar( 	col = FALSE , fill = colblackzone , alpha = alphaZones ) )
                    grid::upViewport()
                 }
                  grid::pushViewport( grid::viewport( layout.pos.col = 2 , layout.pos.row = 1 ) )
                 grid::grid.text( 	x = grid::unit( 0 , "npc" ) ,  "Done" , gp = gpar( col = "black" , fontsize = Fontsize.label.color ) , just = c( 0 , 1/2 ) )
                 grid::grid.circle( 	x = grid::unit( 0.6 , "npc" ) , r = grid::unit( 0.4 , "npc" ) , gp = gpar( col = FALSE , fill = methods::slot( x , "colvect" )[ 2 , lcol] ) )
                  grid::upViewport()
                  grid::pushViewport( grid::viewport( layout.pos.col = 1 , layout.pos.row = 1 ) )
                 grid::grid.text(	x = grid::unit( 0 , "npc" ) ,  paste( "Group" , levels( methods::slot( x , "group") )[ 2 ] ) ,  gp = gpar( col = "black" , fontsize = Fontsize.label.color ) , just = c( 0 , 1/2 ) )
                  grid::upViewport()
                  grid::upViewport()
                  grid::upViewport()
                 grid::pushViewport( grid::viewport( x = grid::unit( 0 , "npc" ) , 
                                         y = grid::unit( (switch(  as.character( is.null( methods::slot( x , "parameters" )$informer ) == FALSE )  , "TRUE" = 1 , "FALSE" = 0 ) + 1 ) / part , "npc" ),  
                                         width = grid::unit( 1 , "npc" ) , 
                                         height = grid::unit( 1/part , "npc" ) , just = c( 0 , 0 )  ) )
                  grid::pushViewport( grid::viewport( layout = grid.layout( 1 , 3 ) , width = 1 , height = 1 ) )
                 if ( length( methods::slot( book , "BZLong" ) )  > 0 ) {
                    grid::pushViewport( grid::viewport( layout.pos.col = 3 , layout.pos.row = 1 ) )
                   grid::grid.text( x = grid::unit( 0 , "npc" ) , "Not in time" , gp = gpar( col = "black" , fontsize = Fontsize.label.color ) , just = c( 0 , 1/2 ) )
                   grid::grid.circle( 	x = grid::unit( 0.8 , "npc" ) ,   r = grid::unit( 0.4 , "npc" ) , gp = gpar( 	col = FALSE ,  fill = methods::slot( x , "colvect" )[ 1 , lcol ]   ) )
                   grid::grid.circle( 	x = grid::unit( 0.8 , "npc" ) ,   r = grid::unit( 0.4 , "npc" ) , gp = gpar( 	col = FALSE , fill = colblackzone , alpha = alphaZones ) )
                    grid::upViewport()
                 }
                  grid::pushViewport( grid::viewport( layout.pos.col = 2 , layout.pos.row = 1 ) )
                 grid::grid.text( 	x = grid::unit( 0 , "npc" ) , "Done" , gp = gpar( col = "black" , fontsize = Fontsize.label.color ) ,   just = c( 0 , 1/2 ))
                 grid::grid.circle( 	x = grid::unit( 0.6 , "npc" ) , r = grid::unit( 0.4 , "npc" ) , gp = gpar( col = FALSE , fill = methods::slot( x , "colvect" )[ 1 ,  lcol ] ) )
                  grid::upViewport()
                  grid::pushViewport( grid::viewport( layout.pos.col = 1 , layout.pos.row = 1 ) )
                 grid::grid.text(	x = grid::unit( 0 , "npc" ) ,  paste( "Group" , levels( methods::slot( x , "group") )[ 1 ] ) ,   gp = gpar( col = "black" , fontsize = Fontsize.label.color ) , just = c( 0 , 1/2 )  )
                  grid::upViewport()
                  grid::upViewport()
                  grid::upViewport()
               }else{
                  grid::pushViewport( grid::viewport( layout = grid.layout( 1 , 2 ) , width = 1 , height = 1 ) )
                 if ( length( methods::slot( book , "BZLong" ) )  > 0 ) {
                    grid::pushViewport( grid::viewport( layout.pos.col = 2, layout.pos.row = 1 ) )
                   grid::grid.text( x = grid::unit( 0 , "npc" ) , "Not in time" , gp = gpar( col = "black" , fontsize = Fontsize.label.color ) , just = c( 0 , 1/2 ) )
                   grid::grid.circle( 	x = grid::unit( 0.8 , "npc" ) ,   r = grid::unit( 0.4 , "npc" ) , gp = gpar( 	col = FALSE , fill = methods::slot( x , "colvect" )[ 1 , lcol ]   ) )
                   grid::grid.circle( 	x = grid::unit( 0.8 , "npc" ) ,  r = grid::unit( 0.4 , "npc" ) ,  gp = gpar( 	col = FALSE , fill = colblackzone , alpha = alphaZones )  )
                    grid::upViewport()
                 }
                  grid::pushViewport( grid::viewport( layout.pos.col = 1 , layout.pos.row = 1 ) )
                 grid::grid.text( 	x = grid::unit( 0 , "npc" ) , "Done" ,  gp = gpar( col = "black" , fontsize = Fontsize.label.color ) , just = c( 0 , 1/2 )  )
                 grid::grid.circle( 	x = grid::unit( 0.6 , "npc" ) , r = grid::unit( 0.4 , "npc" ) ,   gp = gpar( col = FALSE , fill = methods::slot( x , "colvect" )[ 1 , lcol ] ) )
                  grid::upViewport()
                  grid::upViewport()
                  grid::upViewport()
                }
                grid::upViewport()
                grid::upViewport()			
             }
           if (any( methods::slot( book , "typeA" )[ sortindex ] == "p" ) == FALSE &  any( methods::slot( book , "typeA" )[ sortindex ] == "l")) {
                grid::pushViewport( grid::viewport( y = grid::unit( 0.01 , "npc"),
                                       x = grid::unit( (1 - vp0w ) * 2 / 3  , "npc" ) , 
                                       height = grid::unit( (1 - vp0h ) * 1 / 3 + 0.02, "npc" ) , 
                                       width = grid::unit( vp0w , "npc" ) ,
                                       just = c( 0 , 0 )))
               grid::grid.rect( gp = gpar( col = "black" , fill = "grey",alpha = 0.3 ))	
               ###### Long action..................................................................................................................................
               ###### Long action..................................................................................................................................
                grid::pushViewport( grid::viewport( x = grid::unit( 0 , "npc" ) ,  y = grid::unit( 5/6 , "npc" ) ,  width = grid::unit( 1 , "npc" ) , height = grid::unit( 1/6 , "npc" ) ,  just = c( 0 , 0 )  ))
               grid::grid.text(  x = grid::unit( 0 , "npc" ) ,  y = grid::unit( 1 , "npc" ) ,  paste("Long Action  -  N = ",  dim(methods::slot(x,"L"))[1] ),  gp = gpar( col = "black" , fontsize = Fontsize.label.color ) , just = c( 0 , 1 )  )
                grid::upViewport()
                grid::pushViewport( grid::viewport( x = grid::unit( 0 , "npc" ) ,  y = grid::unit( 0 , "npc" ) ,  width = grid::unit( 1 , "npc" ) ,  height = grid::unit( 5/6 , "npc" ) , just = c( 0 , 0 )  )  )
               part = switch( methods::slot( x , "parameters" )$method , "global" = 1 , 2 ) + switch( as.character( is.null( methods::slot( x , "parameters" )$informer ) == FALSE ) , "TRUE" = 1 , "FALSE" = 0 )
               if (is.null( methods::slot( x , "parameters" )$informer ) == FALSE ) {
                  grid::pushViewport( grid::viewport( x = grid::unit( 0 , "npc" ) ,   y = grid::unit( 0 , "npc" ) , width = grid::unit( 1 , "npc" ) ,  height = grid::unit( 1/part , "npc" ) , just = c( 0 , 0 ) )   )
                 grid::grid.text( x = grid::unit( 0 , "npc" ) , paste( switch( 	methods::slot( x , "parameters" )$informer , "mean" = "Mean" ,  "median" = "Median" ), "of the span" ) , gp = gpar( col = "black" , fontsize = Fontsize.label.color ) ,  just = c( 0 , 1/2 )  )
                 grid::grid.lines( x = grid::unit( c( 0.6 , 0.9 ) , "npc" ) , y = grid::unit( c( 0.4 , 0.4 ) , "npc" ) ,  default.units = "npc" ,  arrow = NULL,   gp = gpar( col = "black" , lwd = 2 )  )
                  grid::upViewport()
               }
                grid::pushViewport( grid::viewport( x = grid::unit( 0 , "npc" ) , y = grid::unit( switch( as.character( is.null( methods::slot( x , "parameters" )$informer ) == FALSE ) , "TRUE" = 1 , "FALSE" = 0 ) / part , "npc" ), 
                                       width = grid::unit( 1 , "npc" ) ,
                                       height = grid::unit( 1/part , "npc" ) , 
                                       just = c( 0 , 0 )  ) )
               if (methods::slot( x , "parameters"  )$method != "global"  ) {		
                  grid::pushViewport( grid::viewport( layout = grid.layout( 1 , 3 ) , width = 1 , height = 1 ) )
                 if ( length( methods::slot( book , "BZLong" ) )  > 0 ) {
                    grid::pushViewport( grid::viewport( layout.pos.col = 3 , layout.pos.row = 1 ) )
                   grid::grid.text( x = grid::unit( 0 , "npc" ) , "Not in time" , gp = gpar( col = "black" , fontsize = Fontsize.label.color ) , just = c( 0 , 1/2 ) )
                   grid::grid.circle( 	x = grid::unit( 0.8 , "npc" ) ,   r = grid::unit( 0.4 , "npc" ) ,  gp = gpar( 	col = FALSE , fill = methods::slot( x , "colvect" )[ 2 , lcol ]    ) )
                   grid::grid.circle( 	x = grid::unit( 0.8 , "npc" ) ,  r = grid::unit( 0.4 , "npc" ) ,  gp = gpar( 	col = FALSE , fill = colblackzone , alpha = alphaZones ) )
                    grid::upViewport()
                 }   #    grid::grid.rect( gp = gpar( col = "black" , fill= FALSE ))
                  grid::pushViewport( grid::viewport( layout.pos.col = 2 , layout.pos.row = 1 ) )
                 grid::grid.text( 	x = grid::unit( 0 , "npc" ) ,  "Done" , gp = gpar( col = "black" , fontsize = Fontsize.label.color ) , just = c( 0 , 1/2 ) )
                 grid::grid.circle( 	x = grid::unit( 0.6 , "npc" ) , r = grid::unit( 0.4 , "npc" ) , gp = gpar( col = FALSE , fill = methods::slot( x , "colvect" )[ 2 , lcol ] ) )
                  grid::upViewport()
                  grid::pushViewport( grid::viewport( layout.pos.col = 1 , layout.pos.row = 1 ) )
                  grid::grid.text(	x = grid::unit( 0 , "npc" ) ,  paste( "Group" , levels( methods::slot( x , "group") )[ 2 ] ) ,  gp = gpar( col = "black" , fontsize = Fontsize.label.color ) , just = c( 0 , 1/2 ) )
                  grid::upViewport()
                  grid::upViewport()
                  grid::upViewport()
                 grid::pushViewport( grid::viewport( x = grid::unit( 0 , "npc" ) , 
                                         y = grid::unit( (switch(  as.character( is.null( methods::slot( x , "parameters" )$informer ) == FALSE )  , "TRUE" = 1 , "FALSE" = 0 ) + 1 ) / part , "npc" ),  
                                         width = grid::unit( 1 , "npc" ), 
                                         height = grid::unit( 1/part , "npc" ) , just = c( 0 , 0 )  ) )
                 grid::pushViewport( grid::viewport( layout = grid.layout( 1 , 3 ) , width = 1 , height = 1 ) )
                 if ( length( methods::slot( book , "BZLong" ) )  > 0 ) {
                   grid::pushViewport( grid::viewport( layout.pos.col = 3 , layout.pos.row = 1 ) )
                   grid::grid.text( x = grid::unit( 0 , "npc" ) , "Not in time" , gp = gpar( col = "black" , fontsize = Fontsize.label.color ) , just = c( 0 , 1/2 ) )
                   grid::grid.circle( 	x = grid::unit( 0.8 , "npc" ) ,   r = grid::unit( 0.4 , "npc" ) , gp = gpar( 	col = FALSE ,  fill = methods::slot( x , "colvect" )[ 1 , lcol ]   ) )
                   grid::grid.circle( 	x = grid::unit( 0.8 , "npc" ) ,   r = grid::unit( 0.4 , "npc" ) , gp = gpar( 	col = FALSE , fill = colblackzone , alpha = alphaZones ) )
                   grid::upViewport()
                 }
                  grid::pushViewport( grid::viewport( layout.pos.col = 2 , layout.pos.row = 1 ) )
                    grid::grid.text( 	x = grid::unit( 0 , "npc" ) , "Done" , gp = gpar( col = "black" , fontsize = Fontsize.label.color ) ,   just = c( 0 , 1/2 ))
                    grid::grid.circle( 	x = grid::unit( 0.6 , "npc" ) , r = grid::unit( 0.4 , "npc" ) , gp = gpar( col = FALSE , fill = methods::slot( x , "colvect" )[ 1 , lcol ] ) )
                  grid::upViewport()
                  grid::pushViewport( grid::viewport( layout.pos.col = 1 , layout.pos.row = 1 ) )
                    grid::grid.text(	x = grid::unit( 0 , "npc" ) ,  paste( "Group" , levels( methods::slot( x , "group") )[ 1 ] ) ,   gp = gpar( col = "black" , fontsize = Fontsize.label.color ) , just = c( 0 , 1/2 )  )
                  grid::upViewport()
                  grid::upViewport()
                  grid::upViewport()
               }else{
                  grid::pushViewport( grid::viewport( layout = grid.layout( 1 , 2 ) , width = 1 , height = 1 ) )
                 if ( length( methods::slot( book , "BZLong" ) )  > 0 ) {
                    grid::pushViewport( grid::viewport( layout.pos.col = 2, layout.pos.row = 1 ) )
                       grid::grid.text( x = grid::unit( 0 , "npc" ) , "Not in time" , gp = gpar( col = "black" , fontsize = Fontsize.label.color ) , just = c( 0 , 1/2 ) )
                       grid::grid.circle( 	x = grid::unit( 0.8 , "npc" ) ,   r = grid::unit( 0.4 , "npc" ) , gp = gpar( 	col = FALSE , fill = methods::slot( x , "colvect" )[ 1 , lcol ]   ) )
                       grid::grid.circle( 	x = grid::unit( 0.8 , "npc" ) ,  r = grid::unit( 0.4 , "npc" ) ,  gp = gpar( 	col = FALSE , fill = colblackzone , alpha = alphaZones )  )
                    grid::upViewport()
                 }
                  grid::pushViewport( grid::viewport( layout.pos.col = 1 , layout.pos.row = 1 ) )
                    grid::grid.text( 	x = grid::unit( 0 , "npc" ) , "Done" ,  gp = gpar( col = "black" , fontsize = Fontsize.label.color ) , just = c( 0 , 1/2 )  )
                    grid::grid.circle( 	x = grid::unit( 0.6 , "npc" ) , r = grid::unit( 0.4 , "npc" ) ,   gp = gpar( col = FALSE , fill = methods::slot( x , "colvect" )[ 1 , lcol ] ) )
                  grid::upViewport()
                  grid::upViewport()
                  grid::upViewport()
               }
                grid::upViewport()
                grid::upViewport()			
            }
}
)
