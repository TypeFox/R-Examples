BA.plot <-
function( y1, y2, meth.names = NULL,
                     wh.comp = 1:2,
                     pl.type = "BA",
                    dif.type = "const",
                     sd.type = "const",
                       model = if( inherits(y1,"Meth") & has.repl(y1) ) "exch"
                               else NULL,
                        eqax = FALSE,
                       axlim = if( is.data.frame(y1) ) range(y1$y) else range(c(y1,y2)),
                      diflim = NULL,
                        grid = TRUE,
                      N.grid = 10,
                    col.grid = grey(0.9),
                      points = TRUE,
                  col.points = "black",
                  cex.points = 1,
                  pch.points = 16,
                         lwd = c(3,1,1),
                   col.lines = "blue",
                   repl.conn = FALSE,
                    col.conn = "gray",
                    lwd.conn = 1,
                        xlab = NULL,
                        ylab = NULL,
                         eqn = FALSE,
                     col.eqn = col.lines,
                    font.eqn = 2,
                      digits = 2,
                   Transform = if( mult ) "log" else NULL,
                        mult = FALSE,
                       alpha = NULL,
                         ... )
{
# Allow sloppy definition of arguments
 pl.type <- ifelse( tolower( substr( pl.type,1,1) ) == "c", "conv" , "BA"  )
dif.type <- ifelse( tolower( substr(dif.type,1,1) ) == "c", "const", "lin" )
 sd.type <- ifelse( tolower( substr( sd.type,1,1) ) == "c", "const", "lin" )
 if( !is.null(model) )
   model <- ifelse( tolower( substr(   model,1,1) ) == "l", "linked", "exch" )

if( is.vector( y1 ) )
  # If we have a vector, check if it has the same length as the second argument
  if( length(y1)!=length(y2) )
    stop( "Arguments y1 and y2 must have same length, but",
          "length(y1)=",  length(y1),
        ", length(y2)=",  length(y2) )
  else
  {
  # And if they are of same length, make a Meth object out of it using
  # supplied names if givem
  tmp <- data.frame(y1,y2)
  if( is.character(meth.names) ) names(tmp) <- meth.names
  y1 <- Meth( tmp, y=1:2, print=FALSE )
  }

if( is.data.frame( y1 ) )
  {
  # If the dataframe is not a Meth object make it, Meth will take
  # care of the possible errors if the right columns are not there.
  if( !inherits(y1,"Meth") )
    y1 <- Meth( y1, print=FALSE )
  # Select the two methods to compare and subset the Meth object to
  # the two methods that we plot, and make sure wh.comp hold their names
  if( is.numeric(wh.comp) ) wh.comp <- levels(y1$meth)[wh.comp]
  obj <- Meth( y1[y1$meth %in% wh.comp,], print=FALSE )
  }

else stop("Wrong data structrue for y1 supplied: str(y1):", str(y1) )

# So we turn this into a MethComp object
if( is.null(model) )
M.obj <- DA.reg( obj, DA.slope = dif.type=="lin",
                      Transform = Transform )
else
M.obj <- BA.est( obj, linked = model=="linked",
                      Transform = Transform )

# Then we compute various parameters for the plotting

# axlim:
if( is.null(axlim) ) axlim <- range( obj$y )

# diflim:
if( is.null(diflim) ) diflim <- axlim - mean(axlim)

# And then we can use the default machinery to plot this:
plot.MethComp( M.obj,
             wh.comp = wh.comp,
             pl.type = pl.type,
            dif.type = dif.type,
             sd.type = sd.type,
               axlim = axlim,
              diflim = diflim,
              points = points,
                grid = grid,
              N.grid = N.grid,
            col.grid = col.grid,
                 lwd = lwd,
           col.lines = col.lines,
          col.points = col.points,
          pch.points = pch.points,
           repl.conn = repl.conn,
            col.conn = col.conn,
            lwd.conn = lwd.conn,
                 eqn = eqn,
             col.eqn = col.eqn,
            font.eqn = font.eqn,
              digits = digits,
               alpha = alpha,
                mult = mult,
                 ... )

attr( M.obj, "pl.type" ) <- c( pl.type =  pl.type,
                              dif.type = dif.type,
                               sd.type =  sd.type,
                                 model = model )
invisible( M.obj )
}

######################################################################
# Utility to connect points to means.

connect2mean <-
function( obj, wh.comp,
               pl.type = "conv",
              col.conn = "gray",
              lwd.conn = 1,
                   ... )
{
wob <- to.wide( obj )
              # The points
wob <- cbind( wob[,wh.comp[2]],
              wob[,wh.comp[1]],
              # - and the item averages
          ave(wob[,wh.comp[2]],wob[,"item"]),
          ave(wob[,wh.comp[1]],wob[,"item"]) )
# Convert to D-A coordinates if required
if( pl.type == "BA" )
  {
  ADmat <- rbind(c(0.5,-1),
                 c(0.5, 1))
  wob <- wob %*% rbind( cbind( ADmat, ADmat*0 ),
                        cbind( ADmat*0, ADmat ) )
  }
segments( wob[,1], wob[,2], wob[,3], wob[,4],
          col = col.conn, lwd=lwd.conn, ... )
}
