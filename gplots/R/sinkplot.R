# $Id: sinkplot.R 1673 2013-06-27 20:20:46Z warnes $
sinkplotEnv <-new.env()
  
sinkplot <- function(operation=c("start","plot","cancel"),...)
  {
    operation <- match.arg(operation)

    if( operation=="start" )
      {
        if (exists(".sinkplot.conn", envir=sinkplotEnv) &&
            get(".sinkplot.conn", envir=sinkplotEnv) )
          stop("sinkplot already in force")


        .sinkplot.conn <- textConnection(".sinkplot.data", "w", local=FALSE)
        assign(x=".sinkplot.conn", value=.sinkplot.conn, envir=sinkplotEnv)

        on.exit(sink())
        sink(.sinkplot.conn)
        on.exit()
      }
    else
      {
        if (!exists(".sinkplot.conn", envir=sinkplotEnv) ||
            !get(".sinkplot.conn", envir=sinkplotEnv) )
          stop("No sinkplot currently in force")

        sink()

        data <- get(".sinkplot.data", envir=sinkplotEnv)

        if( operation=="plot" )
            textplot( paste( data, collapse="\n"), ... )

        close(get(".sinkplot.conn", envir=sinkplotEnv))

        if(exists(".sinkplot.data", envir=sinkplotEnv, inherits=FALSE))
          rm(list=".sinkplot.data", pos=sinkplotEnv, inherits=FALSE)

        if(exists(".sinkplot.conn", envir=sinkplotEnv, inherits=FALSE))
          rm(list=".sinkplot.conn", pos=sinkplotEnv, inherits=FALSE)

        invisible(data)
      }
  }
