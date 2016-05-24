plotpunctualsup <- function(X , idsup , iip , book , method , linA , lgH , colvect , alphasup ) {
  if ( any( c( "join" , "cut" , "within" ) == method ) ) {
    grid::pushViewport( grid::viewport( x = unit( 0 , "npc" ) ,											# top 
                            y = unit( linA + (1 - linA ) / 2 , "npc" ) ,
                            width = unit( 1 , "npc" ) ,
                            height = unit( (1 - linA ) / 2 , "npc" ) ,
                            just = c( 0 , 0 ) ) ) 
    plotpunctual(  mat = X[ seq( 1 , dim( X )[ 1 ] / 2 , 1 ) , ] , 
                   iip = iip , 
                   book = book , 
                   colvect = matrix(  colvect[ 1 , ] , 
                   nrow = 1 , 
                   ncol = length(  colvect[ 1 , ] ) ) , 
                   lgH = lgH , 
                   method = "global" , 
                   linA = 0.1 ,
                   alpha = alphasup)
    grid::upViewport()
    grid::pushViewport( grid::viewport( x = unit( 0 , "npc" ) ,										
                            y = unit( 0 , "npc" ) ,
                            width = unit( 1 , "npc" ) ,
                            height = unit( (1 - linA ) / 2 , "npc" ) ,
                            just = c( 0 , 0 ) ) )
    plotpunctual( 	mat = X[ seq( dim( X )[ 1 ] / 2 + 1 , dim( X )[ 1 ] , 1 ) , ] , 
                   iip = iip , 
                   book = book, matrix(  colvect[ 2 , ] , nrow = 1 , ncol = length(  colvect[ 2 , ] ) )  , 
                   lgH = lgH , 
                   method = "global" , 
                   linA = 0.1 ,
                   alpha = alphasup )
    grid::upViewport()
  }else{
    grid::pushViewport( grid::viewport( x = unit( 0 , "npc" ) ,											# top 
                            y = unit( linA + (1 - linA ) / 2 , "npc" ) ,
                            width = unit( 1 , "npc" ) ,
                            height = unit( (1 - linA ) / 2 , "npc" ) ,
                            just = c( 0 , 0 ) ))
    plotpunctual( 	mat = X , 
                   iip = iip , 
                   book = book, colvect =  colvect , 
                   lgH = lgH , 
                   method = "global" , 
                   linA = 0.1 ,
                   alpha = alphasup)
    grid::upViewport()
  }
}
