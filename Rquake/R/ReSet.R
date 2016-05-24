ReSet <- function(nh, g)
  {
####   RSEIS::swig button - reset the station/channel choices for grep in RSEIS::Mine.seis
  
    dev.set( dev.next() )
    s1 = SELstaDB(nh$IDB, sel=1, newdev=TRUE, STAY=FALSE)

    g$RESET = list(gs=g$IDB$sta[s1$sta]  , gc=g$IDB$comp[s1$comp] )

    
    dev.set( g$MAINdev)
    g$zloc = list(x=NULL, y=NULL)
    g$action="donothing"
    invisible(list(global.vars=g))
  }
