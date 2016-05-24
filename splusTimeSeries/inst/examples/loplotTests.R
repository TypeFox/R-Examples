library("splusTimeSeries", lib.loc="~/rPackages/tmpLib")

{
  ## source compare.ps.output() that uses ImageMagick to compare files:
    source(file.path(getenv("STEST"), "loop", "testfuns", "comparePS.q"))
    ps.options(reset=TRUE)
    TRUE
}
{
  if (existsFunction("use.legacy.graphics")) {
      while(length(dev.list()))
           dev.off() # ensure no devices are open before call to 
	             # use.legacy.graphics()
      use.legacy.graphics(TRUE)
  }
  ### SET NEWGEN to TRUE in the environment if you want to generate 
  ###  files for a new platform
  if (getenv("NEWGEN") != "") file.generate.new <- TRUE else 
	  file.generate.new <- FALSE
  plat <- platform()
  if (plat == "WIN386" || plat == "WINX64") {
	compare.loc <- paste(getenv("STEST"), "\\loop\\timeseries\\",
		plat, "_", version$major, 
#		".", version$minor, 
		"\\", sep="")
  } else 
	compare.loc <- paste(getenv("STEST"), "/loop/timeseries/", 
		plat, "_", version$major, 
#		".", version$minor, 
		"/", sep="")

  # get a temporary file name and make a device
  tmpname <- tempfile()
#  compare.loc <- paste( "/homes/jhodgdon/work/slocal/timeser/stest/",
#	plat, ".", version$major, "_", version$minor, "/", sep="")
  goodfile <- paste( compare.loc, "lowplot.ps", sep="")

  if( file.generate.new )
    postscript(file=goodfile,
		setfont=ps.setfont.std, bullet=ps.bullet.std)
  else
    postscript(file=tmpname,
		setfont=ps.setfont.std, bullet=ps.bullet.std)

  if( !file.generate.new )
    cat( "Comparing file for \"", plat, "\" platform\n", sep="" )
  else
    cat( "Generating new file for \"", plat, "\" platform\n", sep="" )

  nchar(plat) > 0
}

{
  # generate a plot from low-level routines

  library("splusTimeSeries")
  user.x <- 1:10
  user.y <- 1:10
  y.tck <- 1:10
  scale.x <- list( c( 0, 11 ), c( 0, 1 ))
  scale.y <- NULL

  x <- c( 0, 1 )
  y <- c( 0, 11 )

  plot( x, y, col = 1, axes = FALSE, type = "n" )

  out1 <- axis.render( 1, ticks = list( list( at = user.x, lwd = 4, col=3 ), 
			              list( at = user.x[-10] + .5, 
						lwd = 2, col=3 )),
		      labels = list( at = user.x, text = user.x, col=3 ),
		      breaks = list( list( at = 2.5, style = "g", col=3 ),
			list( at = 4.5, style = "b", col = 2 ),
			list( at = 7.5, style = "d", width=.02, col=3 )),
		      scale = scale.x)

  out2 <- axis.render( 2, ticks = list( at = user.y, col=4 ),
		      labels = list( at = user.y, text = user.y,
			 srt = 90, col=4 ),
		      scale = scale.y )
  out3 <- axis.render( 3, ticks = list(at = user.x, lwd = 4, col=5 ),
		      breaks =list( list( at = 2.5, style = "g", col=5 ),
			list( at = 4.5, style = "b", col = 5 ),
			list( at = 7.5, style = "d", width=.02, col=5 )),
		      scale = scale.x )
  out4 <- axis.render( 4, ticks = list( at = user.y, col=6), 
	scale = scale.y )

  out5 <- grid.render( list( x = user.x, y = user.y, lty = 2, col=7 ), 
	x.scale = scale.x )

  hloc.stuff <- data.frame( high = user.y + .5,
	   low = user.y - .5,
	   open = c( 1, 1.9, 3.5, 4.2, 4.8, 6.4, 7.0, 8, 8.5, 10.4 ),
	   close = c( 1, 2.1, 2.7, 4, 5, 6, 7, 7.6, 9, 10 ) )

  hloc.render( user.x[1:5], hloc.stuff[1:5,,drop=FALSE], 
	       x.scale = scale.x, y.scale = scale.y, 
	       style = "l", col = 8 )

  hloc.render( user.x[6:10], hloc.stuff[6:10,,drop=FALSE], 
	       x.scale = scale.x, y.scale = scale.y, 
	       style = "c", col = 9 )

  lines.render( x=user.x, hloc.stuff[,3:4,drop=FALSE], 
	x.scale = scale.x, y.scale = scale.y,
	col = 10 )

  stackbar.render( x = user.x, y = hloc.stuff[,c("low", "high"),drop=FALSE] / 2, 
	x.scale = scale.x, y.scale = scale.y,
	line.col = 11, fill.col=c(12,13), cum=TRUE, 
	x.width=(user.x[2]-user.x[1]), density = 20, angle=c(45,-45))

  dev.off()
  TRUE
}

{
  # check returned output of plotting functions
  ( all( out1$ticks[[1]]$at == user.x ) &&
    all( out1$ticks[[2]]$at == user.x[-10] + .5 ) &&
    all( out1$labels$at == user.x ) &&
    all( out1$labels$text == user.x ) &&
    ( length( out1$breaks ) == 3 ))
}

{
  # check output file
  out6 <- NULL
  tmp2 <- NULL

  if( file.generate.new )
    TRUE
  else
    compare.ps.output( goodfile, tmpname )
  
## don't use diff to compare files
## tests started failing 3/29/07 due to fix for bug 44998
#  {
#    tmp2 <- paste( tmpname, ".2", sep = "" )
#    shell( paste( "grep -v \"%\"", tmpname, ">", tmp2 ))
#    out6 <- shell( paste( "grep -v \"%\"", goodfile, "| diff - ", tmp2 ))
#    unlink( tmp2 )
#    if( length( out6 ) == 0 )
#    {
#      unlink( tmpname)
#      TRUE
#    } else
#    {
#      cat( "****Files do not match -- compare", goodfile, "with", tmpname, "\n" )
#      FALSE
#    }
#  }
 
}	

{
  # cleanup
  rm( tmpname, tmp2, goodfile, user.x, user.y, y.tck, scale.x, scale.y, x, y, 
      out1, out2, out3, out4, out5, out6, hloc.stuff, file.generate.new, plat, compare.loc )
  rm( FILE.GENERATE.NEW, TOLERANCE, VERBOSE, COMPARE.LOC,
      makeName, openFile, graphTest, compare.ps.output )

  TRUE
}
