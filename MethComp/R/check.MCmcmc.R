# The new method functions (pirs is already a method
# NOTE: "trace" is a debugging function from the base package
#       which will be masked by this one
trace <- function (obj, ...) UseMethod("trace")
post  <- function (obj, ...) UseMethod("post")

trace.MCmcmc <-
function( obj, what="sd",
       scales = c("same","free"),
       layout = "col",
       aspect = "fill",
          ... )
{
require( coda )
# Expand scales so a single value is sufficient
scales <- c(scales,scales)
# Decide which function to invoke
res <-
if( tolower(what) %in% c("sd","var") )
  trace.sd( obj,
         scales = scales,
         layout = layout,
         aspect = aspect,
            ... )
else
if( tolower(what) %in% c("beta","slope") )
  trace.mean( obj,
           scales = scales,
           layout = layout,
           aspect = aspect,
         par.type = "beta",
              ... )
else
if( tolower(what) %in% c("alpha","int") )
  trace.mean( obj,
           scales = scales,
           layout = layout,
           aspect = aspect,
         par.type = "alpha",
              ... )
# print( res )
return( res )
}

find.vars <-
function( obj,
       layout = "col" )
{
# Find where all the variance estimates are
wh <- grep("sigma",Nam <- varnames(obj))
# Get the names of these
nam <- varnames(obj)[wh]
# Build up the subset indices by method
m.nam <- attributes(obj)$methods
Nm <- length( m.nam )
sb <- numeric(0)
for( i in m.nam ) sb <- c( sb, wh[grep(i,nam)] )
# Make the layout of the traceplot
if( is.character(layout) )
  {
  if( layout=="col" )
    {
    layout <- c( Nm, length(sb)/Nm )
    sb <- as.vector( matrix( sb, nrow=Nm, byrow=TRUE ) )
    }
  else
    {
    layout <- c( length(sb)/Nm, Nm )
    sb <- as.vector( t( matrix( sb, nrow=Nm, byrow=TRUE ) ) )
    }
  }
return( list( sb=sb, layout=layout ) )
}

find.mean <-
function( obj,
       layout = "col",
     par.type = "beta" )
{
# Build up the subset indices by method in the right order
m.nam <- attributes(obj)$methods
Nm <- length( m.nam )
sb <- numeric(0)
for( ir in 1:Nm ) for( ic in 1:Nm ) if( ir != ic )
{
sb <- c(sb, grep( paste(par.type,"\\[",m.nam[ir],".",m.nam[ic],sep=""),
                  varnames(obj) ) )
}
# Create a layout of panels
if( is.character(layout) )
  {
  if( layout=="col" ) layout <- c(Nm-1,Nm)
  else
    {
    layout <- c(Nm,Nm-1)
    sb <- as.vector( t( matrix( sb, nrow=Nm-1, byrow=FALSE ) ) )
    }
  }
return( list( sb=sb, layout=layout ) )
}

trace.sd <-
function( obj,
       scales = c("free","free"),
       layout = "col",
       aspect = "fill",
          ... )
{
fv <- find.vars( obj, layout )
lattice::xyplot( subset( obj, fv$sb ),
                   scales = list(x=list(relation=scales[1]),
                                 y=list(relation=scales[2])),
                   layout = fv$layout,
                 as.table = TRUE,
                   aspect = aspect,
             par.settings = list(strip.background=list(col=gray(0.95))),
           par.strip.text = list(font=2),
                      ... )
}

trace.mean <-
function( obj,
       scales = c("free","free"),
       layout = "col",
       aspect = "fill",
     par.type = "beta",
          ... )
{
fm <- find.mean( obj, layout=layout, par.type=par.type )
lattice::xyplot( subset.MCmcmc( obj, fm$sb ),
                   scales = list(x=list(relation=scales[1]),
                                 y=list(relation=scales[2])),
                   layout = fm$layout,
                   aspect = aspect,
                 as.table = TRUE,
             par.settings = list(strip.background=list(col=gray(0.95))),
           par.strip.text = list(font=2),
                      ...)
}

post <- function( obj, ... ) UseMethod("post")

post.MCmcmc <-
function( obj, what="sd",
        check = TRUE,
       scales = "same",
       layout = "row",
          lwd = 2,
          col,
  plot.points = FALSE,
       aspect = "fill",
          ... )
{
require( coda )
scales <- c(scales,scales)
# Decide on coloring
col <- if( check ) rainbow(length(obj)) # number of chains
       else "black"
# Decide which function to invoke
if( tolower(what) %in% c("sd","var","vc","sigma") )
res <- post.sd( obj,
              check = check,
             scales = scales,
             layout = layout,
                lwd = lwd,
                col = col,
        plot.points = plot.points,
             aspect = aspect,
                ... )

sel <- ( tolower(what) %in% c("alpha","int") ) +
       ( tolower(what) %in% c("beta","slope") )*2
if( sel > 0 )
res <- post.mean( obj,
                check = check,
               scales = scales,
               layout = layout,
                  lwd = lwd,
                  col = col,
          plot.points = plot.points,
               aspect = aspect,
             par.type = c("alpha","beta")[sel],
                  ... )
return( res )
}

post.sd <-
function( obj,
        check = TRUE,
       scales = c("free","free"),
       layout = "row",
          lwd = 2,
          col,
  plot.points = FALSE,
       aspect = "fill",
          ... )
{
fv <- find.vars( obj )
obj <- subset.MCmcmc( obj, fv$sb )
if( !check ) obj <- coda::as.mcmc(as.matrix(obj))

lattice::densityplot( obj,
          layout = fv$layout,
             lwd = lwd,
             col = col,
     plot.points = plot.points,
          aspect = aspect,
  default.scales = list(x=list(relation=scales[1]),
                        y=list(relation=scales[2])),
    par.settings = list(strip.background=list(col=gray(0.95))),
  par.strip.text = list(font=2),
             ... )
}

post.mean <-
function( obj,
        check = TRUE,
       scales = c("free","free"),
       layout = "row",
          lwd = 2,
          col,
  plot.points = FALSE,
       aspect = "fill",
     par.type = "beta",
          ... )
{
fm <- find.mean( obj, layout=layout, par.type=par.type )
obj <- subset.MCmcmc( obj, fm$sb )
# If we believe in convergence
if( !check ) obj <- coda::as.mcmc(as.matrix(obj))

lattice::densityplot( obj,
          layout = fm$layout,
             lwd = lwd,
             col = col,
     plot.points = plot.points,
          aspect = aspect,
  default.scales = list(x=list(relation=scales[1]),
                        y=list(relation=scales[2])),
    par.settings = list(strip.background=list(col=gray(0.95))),
  par.strip.text = list(font=2),
             ... )
}

pairs.MCmcmc <-
function( x, what = "sd",
           subset = NULL,
              col = NULL,
              pch = 16,
              cex = 0.2,
           scales = "free",
              ... )
{
# Select colunms from posterior based on what=
sbset <- NULL
if( any( what %in% c("sd","sigma") ) ) sbset <- c(sbset,"mi","ir","res")
if( any( what %in% c("all")        ) ) sbset <- c(sbset,"tot")
if( any( what %in% c("alpha")      ) ) sbset <- c(sbset,"alpha")
if( any( what %in% c("beta")       ) ) sbset <- c(sbset,"beta")
if( any( what %in% c("mn","mean")  ) ) sbset <- c(sbset,"alpha","beta")
# Select columns from posterior based on subset=
if( is.character(subset) )
  {
  sbset <- NULL
  for( i in 1:length(subset) )
     sbset <- c( sbset, grep( subset[i], colnames(x[[1]]) ) )
  }
if( is.numeric(subset) ) sbset <- subset
sobj <- subset.MCmcmc( x, subset=sbset )
sobj <- as.matrix( sobj )
sobj <- sobj[,order(colnames(sobj))]
if( !is.null(col) )
  {
  col <- if( length(col)==nrow(sobj) ) col
         else
         if( length(col)==coda::nchain(x) ) rep(col,each=coda::niter(x))
         else
         rep( col[1], nrow(sobj) )
  }
else col <- rep(rainbow(coda::nchain(x)),coda::niter(x))
if( toupper(scales)=="SAME" )
  {
  rg <- range( sobj )
  sobj <- rbind( sobj, rg[1], rg[2] )
  col <- c(col,rep("transparent",2))
  }

pairs( sobj, gap=0, pch=pch, cex=cex, col=col, ... )
invisible( varnames( sobj ) )
}
