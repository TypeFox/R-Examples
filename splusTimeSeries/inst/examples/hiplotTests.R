library("splusTimeSeries")

{
  ### SET NEWGEN to TRUE in the environment if you want to generate 
  ###  files for a new platform
  if (getenv("NEWGEN") != "") file.generate.new <- TRUE
  else file.generate.new <- FALSE
  plat <- platform()
  if (plat == "WIN386" || plat == "WINX64") {
    compare.loc <- paste(getenv("STEST"), "\\loop\\timeseries\\",
                         plat, "_", version$major, 
                         ##		".", version$minor, 
                         "\\", sep="")
  } else 
	compare.loc <- paste(getenv("STEST"), "/loop/timeseries/", 
		plat, "_", version$major, 
                             ##		".", version$minor, 
                             "/", sep="")
  ## setup
  if( !file.generate.new )
    cat( "Comparing files for \"", plat, "\" platform\n", sep="" )
  else
    cat( "Generating new files for \"", plat, "\" platform\n", sep="" )

  ##  compare.loc <- paste( "/homes/jhodgdon/work/slocal/timeser/stest/",
  ##	plat, ".", version$major, "_", version$minor, "/", sep="")

  ## source compare.ps.output() that uses ImageMagick to compare files:
  source(file.path(getenv("STEST"), "loop", "testfuns", "comparePS.q"))
  ps.options(reset=TRUE)

  nchar(plat) > 0
}

{
  # generate hloc plot from Dow Jones data
  data("djia", package="splusTimeSeries")
  
  dow <- djia[positions(djia)>=timeDate("09/01/87") &
     positions(djia)<=timeDate("11/01/87"), 1:4]

  tmpname <- tempfile()
  goodfile <- paste( compare.loc, "hipl1.ps", sep="" )

  if( file.generate.new )
    trellis.device( "postscript", file=goodfile,
		setfont=ps.setfont.std, bullet=ps.bullet.std)
  else
    trellis.device( "postscript", file=tmpname,
		setfont=ps.setfont.std, bullet=ps.bullet.std)
  plot(dow, plot.type="hloc")
  dev.off()

  # check postscript file output
  if( file.generate.new )
    TRUE
  else
    compare.ps.output( goodfile, tmpname )
}

{

  # generate hloc plot from tbond data
  library("splusTimeSeries")
  data("tbond", package="splusTimeSeries")
  tb <- tbond[positions(tbond)>=timeDate("02/01/94") & 
     positions(tbond)<=timeDate("02/13/94"),]

  tmpname <- tempfile()
  goodfile <- paste( compare.loc, "hipl2.ps", sep="" )
  if( file.generate.new )
    trellis.device( "postscript", file=goodfile,
		setfont=ps.setfont.std, bullet=ps.bullet.std)
  else
    trellis.device( "postscript", file=tmpname,
		setfont=ps.setfont.std, bullet=ps.bullet.std)
  plot(tb, plot.type="hloc")
  dev.off()

  # check postscript file output
  if( file.generate.new )
    TRUE
  else
    compare.ps.output( goodfile, tmpname )
}


{

  # generate line plot from tbill auction data
  library("splusTimeSeries")
  data("tbauc.3m", package="splusTimeSeries")
  tb3m <- tbauc.3m[positions(tbauc.3m)>=timeDate("01/01/96") &
     positions(tbauc.3m)<=timeDate("07/01/97"),]
  data("tbauc.6m", package="splusTimeSeries")
  tb6m <- tbauc.6m[positions(tbauc.6m)>=timeDate("01/01/96") & 
    positions(tbauc.6m)<=timeDate("07/01/97"),]
  data("tbauc.1y", package="splusTimeSeries")
  tb1y <- tbauc.1y[positions(tbauc.1y)>=timeDate("01/01/96") & 
     positions(tbauc.1y)<=timeDate("07/01/97"),]

  tmpname <- tempfile()
  goodfile <- paste( compare.loc, "hipl3.ps", sep="" )

  if( file.generate.new )
    trellis.device( "postscript", file=goodfile,
		setfont=ps.setfont.std, bullet=ps.bullet.std)
  else
    trellis.device( "postscript", file=tmpname,
		setfont=ps.setfont.std, bullet=ps.bullet.std)

  plot(tb3m, tb6m, tb1y)
  dev.off()

  # check postscript file output

  if( file.generate.new )
    TRUE
  else
    compare.ps.output( goodfile, tmpname )
}

{
  # generate point plot from tbill auction data

  goodfile <- paste( compare.loc, "hipl4.ps", sep="" )
  tmpname <- tempfile()
  if( file.generate.new )
    trellis.device( "postscript", file=goodfile,
		setfont=ps.setfont.std, bullet=ps.bullet.std)
  else
    trellis.device( "postscript", file=tmpname,
		setfont=ps.setfont.std, bullet=ps.bullet.std)

  plot(tb3m, tb6m, tb1y, merge.args=list(pos="union"),
       plot.args=list(type="p"))
  dev.off()

  # check postscript file output
  if( file.generate.new )
    TRUE
  else
    compare.ps.output( goodfile, tmpname )
}

{
  # generate line plot from network packet data
  library("splusTimeSeries")
  data("net.packet", package="splusTimeSeries")
  np <- net.packet[net.packet[,1]=="TCP",2]
  np <- np[1:1000,]

  tmpname <- tempfile()
  goodfile <- paste( compare.loc, "hipl5.ps", sep="" )
  if( file.generate.new )
    trellis.device( "postscript", file=goodfile,
		setfont=ps.setfont.std, bullet=ps.bullet.std)
  else
    trellis.device( "postscript", file=tmpname,
		setfont=ps.setfont.std, bullet=ps.bullet.std)

  plot(np)
  dev.off()

  # check postscript file output

  if( file.generate.new )
    TRUE
  else
    compare.ps.output( goodfile, tmpname )
}

{
  # generate moving avg plot
  tb.hloc <- aggregateSeries( tbond, pos=timeSeq( timeDate("1/7/1994" ),
					 timeDate( "2/4/1995" ), by = "days" ),
			   colnames = c( "high", "low", "open", "close" ),
			   FUN = hloc, together=TRUE )
  tb.avg <- aggregate( tb.hloc[,"close"], moving=20, FUN=mean, adj=1 )

  tmpname <- tempfile()
  goodfile <- paste( compare.loc, "hipl6.ps", sep="" )

  if( file.generate.new )
    trellis.device( "postscript", file=goodfile,
		setfont=ps.setfont.std, bullet=ps.bullet.std)
  else
    trellis.device( "postscript", file=tmpname,
		setfont=ps.setfont.std, bullet=ps.bullet.std)
  plout <- plot( tb.hloc, plot.type="hloc", main="T-Bonds")
  lines.render( positions(tb.avg), seriesData(tb.avg), 
	      x.scale = plout$scale, col=3 )
  dev.off()

  # check postscript file output

  if( file.generate.new )
    TRUE
  else
    compare.ps.output( goodfile, tmpname )
}

{
  # generate signal plot from say.wavelet
  library("splusTimeSeries")
  data("say.wavelet", package="splusTimeSeries")

  tmpname <- tempfile()
  goodfile <- paste( compare.loc, "hiplsig1.ps", sep="" )

  if( file.generate.new )
    trellis.device( "postscript", file=goodfile,
		setfont=ps.setfont.std, bullet=ps.bullet.std)
  else
    trellis.device( "postscript", file=tmpname,
		setfont=ps.setfont.std, bullet=ps.bullet.std)
  plot( say.wavelet )
  dev.off()

  # check postscript file output

  if( file.generate.new )
    TRUE
  else
    compare.ps.output( goodfile, tmpname )
}

{
  # generate log plot from say.wavelet
  tmpname <- tempfile()
  goodfile <- paste( compare.loc, "hiplsig2.ps", sep="" )

  if( file.generate.new )
    trellis.device( "postscript", file=goodfile,
		setfont=ps.setfont.std, bullet=ps.bullet.std)
  else
    trellis.device( "postscript", file=tmpname,
		setfont=ps.setfont.std, bullet=ps.bullet.std)
  plot( say.wavelet[-1,] + 300, log.axes="xy" )
  dev.off()

  # check postscript file output

  if( file.generate.new )
    TRUE
  else
    compare.ps.output( goodfile, tmpname )
}

{
  # generate db plot from say.wavelet
  tmpname <- tempfile()
  goodfile <- paste( compare.loc, "hiplsig3.ps", sep="" )

  if( file.generate.new )
    trellis.device( "postscript", file=goodfile,
		setfont=ps.setfont.std, bullet=ps.bullet.std)
  else
    trellis.device( "postscript", file=tmpname,
		setfont=ps.setfont.std, bullet=ps.bullet.std)
  plot( say.wavelet + 300, dB=TRUE )
  dev.off()

  # check postscript file output

  if( file.generate.new )
    TRUE
  else
    compare.ps.output( goodfile, tmpname )
}

{
  # cleanup
  rm( tmpname, goodfile, plat, file.generate.new, 
      compare.ps.output, compare.loc )
  rm( FILE.GENERATE.NEW, TOLERANCE, VERBOSE, COMPARE.LOC,
      makeName, openFile, graphTest )
  rm( dow, tb, tb3m, tb6m, tb1y, np, tb.hloc, tb.avg, plout )
  TRUE
}




