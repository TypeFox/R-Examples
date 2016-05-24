plotpunctual <- function(mat , iip , book , colvect , lgH , method , linA , alpha = 1) {    
  if (method == "global" ) {
    grid::pushViewport( grid::viewport( x = unit(  1/2 , "npc" ) ,
                            y = unit( 1/2 , "npc" ) ,
                            width = unit( 1 , "npc" ) ,
                            height = unit( linA , "npc" ) ,
                            just = c(  1/2 , 1/2 ) ))
    grid::pushViewport( grid::viewport( layout = grid.layout( 1 , lgH , widths = 1 , heights = 1  , just = c( 0 , 0 )) ) )
    temp 	<- which( mat[ iip, ] != 0 )
    if (length( temp ) > 0 ) { 
      for (jj in temp ) { 
        grid::pushViewport( grid::viewport( layout.pos.row = 1, layout.pos.col = jj , just = c( 0 , 0 ) ) )
        grid::grid.rect( gp = gpar( 	col    = colvect[ 1 , mat[ iip , jj ] ] ,
                               fill   = colvect[ 1 , mat[ iip , jj ] ] ))
        grid::upViewport()
      }
    }
    grid::upViewport()
    grid::upViewport()
  }
  if (method == "join" ) {   # Division per group only at last two person of the qsame group are on the same pixel of time
    grid::pushViewport( grid::viewport( x = unit( 0 , "npc" ) ,
                            y = unit( 1/2 , "npc" ) ,
                            width = unit( 1 , "npc" ) , 
                            height = unit( linA , "npc" ) ,
                            just = c( 0 , 1/2 ) ) )
    grid::pushViewport( grid::viewport( layout = grid.layout( 1 , lgH , widths = 1 , heights = 1 ) ) )
    temp 	<- which( mat[ iip, ] != 0 & mat[ iip + (dim( mat )[ 1 ] / 2 ) , ] != 0)    # Retrieving group 1 gridMatrice
    if (length( temp ) > 0 ) {
      for (jj in temp ) { 
        grid::pushViewport( grid::viewport( layout.pos.row = 1 , layout.pos.col = jj ) )
        grid::pushViewport( grid::viewport( x = unit( 0 , "npc" ) ,
                                y = unit( 1/2 , "npc" ) ,
                                width = unit( 1 , "npc" ) , 
                                height = unit( 0.8 , "npc" ) ,
                                just = c( 0 , 1/2 ) ) ) 
        grid::pushViewport( grid::viewport( layout = grid.layout( 2 , 1 , widths = 1 , heights = 1 ) ) )
        grid::pushViewport( grid::viewport( layout.pos.row = 1 , layout.pos.col = 1 ) )
        grid::grid.rect( gp = gpar( 	col = colvect[ 1 , ][ mat[ iip , jj ] ] , 
                               fill = colvect[ 1 , ][ mat[ iip , jj ] ] ,
                               alpha = alpha ) )
        grid::upViewport()
        grid::pushViewport( grid::viewport( layout.pos.row = 2 , layout.pos.col = 1 ) )
        grid::grid.rect( 	gp = gpar( 	col = colvect[ 2 , ][ mat[ iip + (dim( mat )[ 1 ] / 2 ) , jj ] ] , 
                                fill = colvect[ 2 , ][ mat[ iip + (dim( mat )[ 1 ] / 2 ) , jj ] ] ,
                                alpha = alpha ) )
        grid::upViewport()
        grid::upViewport()
        grid::upViewport()
        grid::upViewport()
      }
    }
    temp  <- which( mat[ iip, ] != 0 & mat[ iip + (dim( mat )[ 1 ] / 2 ) , ] == 0 )
    if (length( temp ) > 0 ) {
      for (jj in temp ) {
        grid::pushViewport( grid::viewport( layout.pos.row = 1 , layout.pos.col = jj ) )
        grid::grid.rect( gp = gpar( 	col = colvect[ 1 , mat[ iip , jj ] ] , 
                               fill = colvect[ 1 , mat[ iip , jj ] ] , 
                              alpha = alpha))
        grid::upViewport()
      }
    }
    temp  <- which( mat[ iip + (dim( mat )[ 1 ] / 2 ) , ] != 0 & mat[ iip , ] == 0 )
    if (length( temp ) > 0 ) {
      for (jj in temp ) {
        grid::pushViewport( grid::viewport( layout.pos.row = 1 , layout.pos.col = jj ) ) 
        grid::grid.rect( gp = gpar( 	col = colvect[ 2 , mat[ iip + (dim( mat )[ 1 ] / 2 ) , jj ] ] , 
                               fill = colvect[ 2 , mat[ iip + (dim( mat ) [ 1 ] / 2 ) ,jj ] ] , 
                   alpha = alpha ))
        grid::upViewport()
      }
    }
    grid::upViewport()
    grid::upViewport()
  }
  if (method == "within" | method == "cut" ) {
    grid::pushViewport( grid::viewport( x = unit( 0 , "npc" ) ,
                            y = unit( 1/2 , "npc" ) , 
                            width = unit( 1 , "npc" ) ,
                            height = unit( linA , "npc" ) , 
                            just = c( 0 , 1/2 ) ) )
    grid::pushViewport( grid::viewport( layout = grid.layout( 1 , lgH ) ) )   ### Grid time
    temp <- which( mat[ iip , ] != 0 )
    if (length( temp ) > 0 ) {
      for (jj in temp ) {
        grid::pushViewport( grid::viewport( layout.pos.row = 1 , layout.pos.col = jj ) )  						# push into the grid case
        grid::pushViewport( grid::viewport( layout = grid.layout( 2 , 1 ) ) )   						# division of the grid in two 
        grid::pushViewport( grid::viewport( layout.pos.row = 1 , layout.pos.col = 1 ) )  				# go in the top part
        grid::grid.rect( gp = gpar( 	col = colvect[1 , mat[ iip , jj ] ] , 
                               fill = colvect[ 1 , mat[ iip , jj ] ] ,
                               alpha = alpha))
        grid::upViewport()
        grid::upViewport()
        grid::upViewport()
      }
    }
    temp <- which( mat[ iip + (dim( mat )[ 1 ] / 2 ) , ] != 0 )
    if (length( temp ) > 0 ) {
      for (jj in temp ) {
        grid::pushViewport( grid::viewport( layout.pos.row = 1 , layout.pos.col = jj ) )						# push into the grid case
        grid::pushViewport( grid::viewport( layout = grid.layout( 2 , 1 , widths = 1 , heights = 1 ) ) )			# division of the grid in two 
        grid::pushViewport( grid::viewport( layout.pos.row = 2 , layout.pos.col = 1 ) )					# # go in the bellow part
        grid::grid.rect( gp = gpar( 	col = colvect[ 2 , mat[ iip + (dim( mat )[ 1 ] / 2 ) , jj ] ] ,
                               fill = colvect[ 2 , mat[ iip + (dim( mat )[ 1 ] / 2 ) , jj ] ] ,
                               alpha = alpha))
        grid::upViewport()
        grid::upViewport()
        grid::upViewport()
      }
    }
    grid::upViewport()
    grid::upViewport()
  }
}
