#' Method plot for ViSibook object.
#' @name plot-ViSibook-method
#' @title Method \code{plot-ViSibook}
#' @rdname plot-ViSibook-methods
#' @aliases plot,ViSibook-method
#' @exportMethod plot
#' @docType methods
#' @param x a ViSibook object.
#' @param ncharmax is the number maximum of plotted character for the labels of punctual actions, set to 10.
#' @param ncharmaxdelay number maximum of  plotted character for the labels of long actions, set to 50.
#' @seealso \code{\linkS4class{ViSibook}}, \code{\link{buildViSiGrid}}, 
#' \code{\link{changeShoworder-ViSibook-method}} and see \code{\link{plot-ViSigrid-method}} for examples.
setMethod( f = "plot",
           signature = "ViSibook",
           definition = function(x, ncharmax = 10, ncharmaxdelay  =50 ) {
           grid::grid.newpage()
           grid::pushViewport(grid::viewport( 	x = grid::unit(1/2 , "npc" ) , 
                                   y = grid::unit( 1/2 , "npc" ) , 	
                                   width = grid::unit( 0.95, "npc" ) , 
                                   height = grid::unit( 0.95, "npc" ) , 
                                   just = c( 1/2 , 1/2 ) 
           ))
          
             index_p <- which(  x[ , 3 ] == "p" )
             index_l <- which( is.na(x[ , 4 ] ) == FALSE & x[ , 3 ] == "l" )
             index_l <- index_l[sort( x[ , 4 ][index_l] , index.return = TRUE)$ix]
             grid::pushViewport(grid::viewport( layout = grid::grid.layout( length( index_l ) + 3, length( index_p ) + 1 , widths = 1 , heights = 1 ) ))
               for (i in seq_along(index_p)) { # i =1
                 grid::pushViewport( grid::viewport( layout.pos.row = 3 , layout.pos.col = i  ) ) 
                  grid::grid.lines( x = grid::unit(c(0,1),"npc"),y = grid::unit(0.5,"npc") ,gp = grid::gpar(col = "black") )      
                  grid::grid.lines( x = grid::unit(1,"npc"),y = grid::unit(c(0.5,0.7),"npc") ,gp = grid::gpar(col = "black") )
                 grid::upViewport()
                 grid::pushViewport( grid::viewport( layout.pos.row = 1 , layout.pos.col = i  ) ) 
                  temp <- .strInsert_cha( x[ , 2 ][index_p[i]], strInsert = " \n " , strWhere = " " ,lchar = ncharmax  )
                  grid::grid.text( label = temp, x = grid::unit(1,"npc"), y = grid::unit(1/2,"npc") , rot = 0,gp = grid::gpar( col = "black",fontsize = 10))
                 grid::upViewport()      
               }
               grid::pushViewport( grid::viewport( layout.pos.row = 3 , layout.pos.col = length( index_p ) +1 ) ) 
               grid::grid.lines( x = grid::unit(c(0,1),"npc"),y = grid::unit(0.5,"npc") ,gp = grid::gpar(col = "black") )      
               grid::upViewport()
               index_ps <- index_p[ which( is.na(x[ , 4 ][index_p] ) == FALSE  )]
    
               for (i in seq_along(index_ps)) { # i =1
                 grid::pushViewport( grid::viewport( layout.pos.row = 3 , layout.pos.col = which(index_p == index_ps[i] ) ) ) 
                  grid::grid.circle( x = grid::unit(1,"npc"),y = grid::unit(1,"npc"), r = grid::unit(0.4,"npc"), gp = grid::gpar(col = FALSE,fill = "blue",alpha = 0.6))
                  grid::grid.text( label = paste( x[ , 4 ][index_ps[i]] ), x = grid::unit(1,"npc"), y = grid::unit(1,"npc"), gp = grid::gpar( col = "white", fontsize = 11))
                 grid::upViewport()
               }
               for (j in seq_along( index_l )) {# j = 1
                 ipdeb <- which(x[ , 1 ][which(x[ , 3 ] == "p"  )] == x[ , 5 ][index_l[j]] )
                 
                 ipfin <- which(x[ , 1 ][which(x[ , 3 ] == "p"  )] == x[ , 6 ][index_l[j]] )
                 temp <- .strInsert_cha( x[ , 2 ][index_l[j]], strInsert = " \n ", strWhere = " ", lchar = ncharmaxdelay)
                 
                 for (ii in seq( ipdeb + 1, ipfin)) { #ii=2
                   grid::pushViewport( grid::viewport( layout.pos.row = j + 3 , layout.pos.col = ii  ) ) 
                    grid::grid.rect( x = grid::unit(0.5, "npc"), y = grid::unit(0.5, "npc"), gp = grid::gpar(col = FALSE, fill = "navy", alpha = 0.1))
                    grid::grid.lines(x = grid::unit(c(0,1),"npc"),y = grid::unit(1,"npc") , gp = grid::gpar(col = "black"))
                   grid::upViewport() 
                 }
                 grid::pushViewport( grid::viewport( layout.pos.row = j + 3 , layout.pos.col = ipdeb  ) ) 
                   grid::grid.circle( x = grid::unit(0.5,"npc"), y = grid::unit(1,"npc"), r = grid::unit(0.4,"npc"), gp = grid::gpar(col = FALSE, fill = "blue", alpha = 0.6))
                   grid::grid.text( label = paste( x[ , 4 ][index_l[j]] ), x = grid::unit(0.5,"npc"), y = grid::unit(1,"npc"), gp = grid::gpar( col = "white",fontsize = 11))
                   grid::grid.circle( x = grid::unit(1,"npc"), y = grid::unit(1,"npc"), r = grid::unit(0.1,"npc"), gp = grid::gpar(col = "black", fill = "black"))
                   grid::grid.text( label = temp, x = grid::unit(1.5,"npc"), y = grid::unit(1/2,"npc"), rot = 0, gp = grid::gpar(col = "black", fontsize = 10))
                 grid::upViewport() 
                 grid::pushViewport( grid::viewport( layout.pos.row = j + 3 , layout.pos.col = ipfin  ) ) 
                  grid::grid.circle( x = grid::unit(1,"npc"), y = grid::unit(1,"npc"), r = grid::unit(0.1,"npc"), gp = grid::gpar( col = "black", fill = "black"))
                 grid::upViewport()
                  } 
             grid::upViewport()
           grid::upViewport()
           })
##
.strInsert_cha <- function(str, strInsert = " \n ", strWhere = " ", lchar = 20 ){
  # str : string to change
  temp <-  as.vector( unlist( gregexpr(strWhere, str) ))
  if ( any( temp > lchar ) ) {
    temp2 <- list()
    while (any( temp > lchar ) ) {
      temp2[[ length( temp2 ) + 1 ]] <- stringr::str_sub( str, 
                                                start = 1, 
                                                end =  temp[ which.min( abs( temp - lchar ) ) ] - 1 )
      str <- stringr::str_sub(str, 
                     start =  temp[ which.min( abs( temp - lchar ) ) ]  + 1, 
                     end =  nchar(str) 
      )
      temp <-  as.vector( unlist( gregexpr(strWhere, str) ) )
    }
    temp2[[ length( temp2 ) + 1 ]] <- str
    str <- temp2[[ 1 ]]
    for (i in seq( 2 , length( temp2) , 1 ) ) {
      str <- stringr::str_c( str, temp2[[ i ]] , sep = strInsert )
    }
  }
  return( str )
}
