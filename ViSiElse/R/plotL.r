plotL <- function(L , idsort , inftps , group , BZL , Lsup , idsup , iip ,t_0, cols , linA , alphaZones , alphasup , colblackzone){
  ngl <- dim( L )[ 1 ]
  if( is.character( t_0)==FALSE & t_0 != 0 ){
    L<- L- matrix( rep(t_0, dim(L)[1]*dim(L)[2] ),nrow=dim(L)[1])
  }
  # Retrieval of non -1 line
  r = which( L[idsort , 1] < 0 | L[idsort , 2] <0 ) 
  # Retrieval of non NA line
  r = unique( c(r, which( is.na(L[idsort , 1]) | is.na(L[idsort , 2]) ) 	) )
  if (length(r) > 0) {
    idsort <- idsort[-r]
  }
  idsort <- as.vector(idsort)
 if ( length(idsort) > 0 ) {
  if ( dim( L )[ 1 ] > 0 ) {  																	# If there is some individuals having done the action
    grid::pushViewport( grid::viewport( x = grid::unit( 0 , "npc" ) ,
                            y = grid::unit( 1/2 , "npc" ) ,
                            width = grid::unit( 1 , "npc" ) , 
                            height = grid::unit( linA , "npc" ) ,
                            just = c( 0 , 1/2 )))
    grid::pushViewport( grid::viewport( layout = grid.layout( ngl , 1 ) ) )					# layout individuals
   if ( length(idsort) > 1 ) {
     for (ind in idsort ) {
       grid::pushViewport( grid::viewport( layout.pos.row = which( idsort == ind ) , 				# goto layout corresponding to ViSigrid object sorting
                              layout.pos.col = 1 ))
       grid::pushViewport( grid::viewport( x = grid::unit(  L[ ind , 1 ]  / inftps , "npc" ) ,
                              y = grid::unit( 0 , "npc" ) , 
                              width = grid::unit( (L[ ind , 2 ] -  L[ ind , 1 ] ) / inftps , "npc" ) ,
                              height = grid::unit( 1 , "npc" ) , 
                              just = c( 0 , 0 ) ))
      if (length( group ) > 0 ) {
        grid::grid.rect( 	gp = gpar(	col = FALSE , 
                               fill = cols[ which( levels( group ) == group[ ind ] ) , dim( cols )[ 2 ] %/% 2 ] ,
                               alpha = 1))
      }else{
        grid::grid.rect( 	gp = gpar(	col = FALSE , 
                               fill = cols[ , dim( cols )[ 2 ] %/% 2 ] ,
                               alpha = 1))
      }
       grid::upViewport()
       grid::upViewport()
    }				
   }else{
     ind = idsort
     grid::pushViewport( grid::viewport( x = grid::unit(  L[ ind , 1 ]  / inftps , "npc" ) ,
                             y = grid::unit( 0 , "npc" ) , 
                             width = grid::unit( (L[ ind , 2 ] -  L[ ind , 1 ] ) / inftps , "npc" ) ,
                             height = grid::unit( 1 , "npc" ) , 
                             just = c( 0 , 0 )))
     if (length( group ) > 0 ) { 
       grid::grid.rect( 	gp = gpar(	col = cols[ which( levels( group ) == group[ ind ] ) , dim( cols )[ 2 ]] , 
                              fill = cols[ which( levels( group ) == group[ ind ] ) , dim( cols )[ 2 ]  ] ,
                              alpha = 1 ) )
     }else{
       grid::grid.rect( 	gp = gpar(	col = cols[ , dim( cols )[ 2 ]  ] , 
                              fill = cols[ , dim( cols )[ 2 ]  ] ,
                              alpha = 1))
     }
     grid::upViewport()
   } 
    if ( length( BZL )  > 0 & ( is.character( t_0) | t_0 == 0)) {															# If a black zone was defined in the book
     if (length(r) > 0) {
        BZL[ r , iip ] <- rep(0,  length(r) )
     }
      if ( any( BZL[ , iip ]  > 0 ) ) {													# If a black zone was defined in the book for this action 
          for (ind in  which( BZL[ , iip ] > 0)  ) {	 						# for the individuals being in the black zone
         grid::pushViewport( grid::viewport( layout.pos.row = which( idsort == ind ), layout.pos.col = 1 )  )		# goto layout corresponding to ViSigrid object sorting
         grid::pushViewport( grid::viewport( x = grid::unit(  BZL[ ind , iip ]  / inftps , "npc" ) ,
                                  y = grid::unit( 0 , "npc" ) , 
                                  width = grid::unit( (L[ ind , 2 ] -  BZL[ ind , iip ] ) / inftps , "npc" ) ,
                                  height = grid::unit( 1 , "npc" ) , 
                                  just = c( 0 , 0 )))
         grid::grid.rect( gp = gpar(	col = FALSE , 
                                fill = colblackzone ,
                                alpha = alphaZones))
         grid::upViewport()
         grid::upViewport()
        }
      }
    }
    if ( dim( Lsup )[ 1 ]  > 0 ) {															              # If Lsup was add
      if (any(  Lsup[ , iip*2 - 1 ] > 0 | Lsup[ , iip * 2 ]  > 0 ) ) {						# If Lsup was addd for this action 
        for (ind in unique( idsup[ which( Lsup[ , iip*2 - 1 ] > 0 | Lsup[ , iip * 2 ]  > 0 ) ] ) ) {# for the individuals has sup value for the action		
          grid::pushViewport( grid::viewport( layout.pos.row = which( idsort == ind ), layout.pos.col = 1)) # goto layout corresponding to ViSigrid object sorting	
          for (i in which( idsup == ind ) ) {										# for all an individual repetitions					
            grid::pushViewport( grid::viewport( x = grid::unit(  Lsup[  i , iip * 2 - 1 ]  / inftps , "npc" ) ,
                                    y = grid::unit( 0 , "npc" ) , 
                                    width = grid::unit( (Lsup[ i , iip * 2 ] -  Lsup[ i , iip * 2 - 1 ] ) / inftps , "npc" ) ,
                                    height = grid::unit( 1 , "npc" ) , 
                                    just = c( 0 , 0 ) ) )
            if (length( group ) > 0 ) {
              grid::grid.rect( gp = gpar(	col = cols[ which( levels( group ) == group[ ind ] )  , dim( cols )[ 2 ] %/% 2 ] , 
                                    fill = cols[ which( levels( group ) == group[ ind ] )  , dim( cols )[ 2 ] %/% 2 ] ,
                                    alpha = alphasup))
            }else{
              grid::grid.rect( gp = gpar(	col = cols[ , dim( cols )[ 2 ] %/% 2 ] , 
                                    fill = cols[ , dim( cols )[ 2 ] %/% 2 ] ,
                                    alpha = alphasup))
            }
            grid::upViewport()
          }
          grid::upViewport()
        }
      }
    }
    grid::upViewport()			
    grid::upViewport()
  }
}
}
