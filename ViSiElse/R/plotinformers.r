plotinformers <- function(informers ,  inftps , iip ,t_0, alphainf , lwdline , rcircle = 1 ){
  # informer matrix or data.frame with 3 row, 
  # inftps : max (vect temps)
  # iip : column indextation in informers thats will be used 
  # alphainf : alpha for the color of circles 
  # hline : size of the line set to 1
  # unit=0,just1=c(0,1)
  ## Plot line 
  ###  
  if (any(is.na(informers[ ,iip ])) == FALSE ) {
    if( is.character( t_0)==FALSE & t_0 != 0 ){
      informers<- informers- matrix( rep(t_0, dim(informers)[1]*dim(informers)[2] ),nrow=dim(informers)[1])
    }
    grid::pushViewport( grid::viewport( x = grid::unit(  informers[ 1 , iip ] / inftps , "npc" ) , 
                          y = grid::unit( 0 , "npc" ) ,
                          width = grid::unit(  (informers[ 3 , iip ] -  informers[ 1 , iip ] ) / inftps , "npc" ) , 
                          height = grid::unit( 1 , "npc" ) , 
                          just = c( 0 , 0)))
      grid::grid.lines( x = grid::unit( c( 0 , 1 ) , "npc" ) , y = grid::unit( c( 1/2 , 1/2 ) , "npc" ) , default.units = "npc" , arrow = NULL , gp = gpar( col = "black" , lwd = lwdline))
    grid::upViewport()
  #### Plot circle min
    grid::pushViewport( grid::viewport( x = grid::unit( informers[ 1 , iip ] / inftps , "npc" ) , 
                          y = grid::unit( 0.2 , "npc" ) ,
                          width = grid::unit( rcircle , "points" ) , 
                          height = grid::unit( 0.6 , "npc" ) ,
                          just = c( 1/2 , 0 )))
    grid::grid.circle( gp = gpar( col = "black" , fill = "white" , alpha = alphainf ) )
    grid::upViewport()
  ##### Plot circle middle 
    grid::pushViewport( grid::viewport( x = grid::unit( informers[ 2 , iip ] / inftps , "npc" ) , 
                          y = grid::unit( 0.2 , "npc" ) ,
                          width = grid::unit( rcircle , "points" ) , 
                          height = grid::unit( 0.6 , "npc" ) ,
                          just = c( 1/2 , 0 ) ))
    grid::grid.circle( gp = gpar( col = "black" , fill = "white" , alpha = alphainf ) )
  grid::upViewport()
  ##### Plot max circle
  grid::pushViewport(viewport( 	x = grid::unit( informers[ 3 , iip ] / inftps , "npc" ) ,
                          y = grid::unit( 0.2 , "npc" ) ,
                          width = grid::unit( rcircle , "points" ) , 
                          height = grid::unit( 0.6 , "npc" ) , just = c( 1/2 , 0 ) ))
    grid::grid.circle( gp = gpar( col = "black" , fill = "white" , alpha = alphainf ) )
  grid::upViewport()
  } 
}	
