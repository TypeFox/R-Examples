# First the utility functions


######################################################################
lint <-
function( ci, tt, u )
{
# Makes a linear interpolation, but does not crash if all ci values are
# identical, but requires that both ci and tt are non-decreasing.
# ci plays the role of cumulative intensity, tt of time
if( any( diff(ci)<0 ) | any( diff(tt)<0 ) ) stop("Non-icreasing arguments")
c.u <- min( c( ci[ci>u], max(ci) ) )
c.l <- max( c( ci[ci<u], min(ci) ) )
t.u <- min( c( tt[ci>u], max(tt) ) )
t.l <- max( c( tt[ci<u], min(tt) ) )
# c.u==c.l if u is outside the range of ci
ifelse( c.u==c.l, t.l, t.l + (u-c.l)/(c.u-c.l)*(t.u-t.l) )
}


######################################################################
sim1 <-
function( rt, time.pts )
{
# Simulates a single transition time and state based on the data frame
# rt with columns lex.id and timescales. It is assumed that the coumns
# in in rt are the id, followed by the set of estimated transition
# rates to the different states reachable from the current one.
ci <- apply( rbind(0,rt[,-1,drop=FALSE]), 2, cumsum )[1:nrow(rt),,drop=FALSE]
tt <- uu <- -log( runif(ncol(ci)) )
for( i in 1:ncol(ci) ) tt[i] <- lint( ci[,i], time.pts, uu[i] )
# Note this resulting data frame has 1 row and is NOT a Lexis object
data.frame( lex.id  = rt[1,1],
            lex.dur = min(tt,na.rm=TRUE),
            lex.Xst = factor( if( min(tt) < max(time.pts) )
                                  colnames(ci)[tt==min(tt)]
                              else NA ) )
}


######################################################################
simX <-
function( nd, Tr, time.pts, tS )
{
# Simulation is done from the data frame nd, in chunks of starting
# state, lex.Cst. This is necessary because different states have
# different (sets of) exit rates. Therefore, this function simulates
# for a set of persons from the same starting state.
np <- length( time.pts )
nr <- nrow( nd )
if( nr==0 ) return( NULL )

# The 'as.character' below is necessary because indexing by a factor
# by default is by the number of the level, and we are not indexing by
# this, but by components of Tr which just happens to have names that
# are a subset of the levels of lex.Cst.
cst <- as.character( unique(nd$lex.Cst) )
if( length(cst)>1 ) stop( "More than one lex.Cst present:\n", cst, "\n" )

# Expand each person by the time points
prfrm <- nd[rep(1:nr,each=np),]
prfrm[,tS] <- prfrm[,tS] + rep(time.pts,nr)
prfrm$lex.dur <- il <- min( diff(time.pts) )
# Poisson-models should use the estimated rate at the midpoint of the
# intervals:
prfrp <- prfrm
prfrp[,tS] <- prfrp[,tS]+il/2

# Make a data frame with predicted rates for each of the transitions
# out of this state for these times
rt <- data.frame( lex.id = prfrm$lex.id )
for( i in 1:length(Tr[[cst]]) )
   {
   if( inherits( Tr[[cst]][[i]], "glm" ) )
   rt <- cbind( rt, predict( Tr[[cst]][[i]],
                             type="response",
                             newdata=prfrp ) )
   else
   if( inherits( Tr[[cst]][[i]], "coxph" ) )
   rt <- cbind( rt, predict( Tr[[cst]][[i]],
                             type="expected",
                             newdata=prfrm ) )
   else
   if( is.function( Tr[[cst]][[i]] ) )
   rt <- cbind( rt, Tr[[cst]][[i]](prfrm) )
   else
   stop( "Invalid object supplied as transition, elements of the list must be either:\n",
         "- a glm(poisson) object fitted to a Lexis object\n",
         "- a coxph object fitted to a Lexis object\n",
         "- a function that takes a Lexis object as argument and returns\n",
         "  average rates for each record in the same units as lex.dur.")
   }
names( rt )[-1] <- names( Tr[[cst]] )

# Then find the transition time and exit state for each person:
xx <- match( c("lex.dur","lex.Xst"), names(nd) )
if( any(!is.na(xx)) ) nd <- nd[,-xx[!is.na(xx)]]
merge( nd,
       do.call( rbind,
                lapply( split( rt,
                               rt$lex.id ),
                        sim1,
                        time.pts ) ),
       by="lex.id" )
}


######################################################################
get.next <-
function( sf, tr.st, tS, tF )
{
# Produces an initial Lexis object for the next simulation for those
# who have ended up in a transient state.
# Note that this exploits the existence of the "time.since" attribute
# for Lexis objects and assumes that a character vector naming the
# transient states is supplied as argument.
if( nrow(sf)==0 ) return( sf )
nxt <- sf[sf$lex.Xst %in% tr.st,]
if( nrow(nxt) == 0 ) return( nxt )
nxt[,tS] <- nxt[,tS] + nxt$lex.dur
wh <- tF
for( i in 1:length(wh) )
   if( wh[i] != "" ) nxt[nxt$lex.Xst==wh[i],tS[i]] <- 0
nxt$lex.Cst <- nxt$lex.Xst
return( nxt )
}


######################################################################
chop.lex <-
function( obj, tS, cens )
{
# A function that chops off all follow-up beyond cens since entry for
# each individual
# Entry times on 1st timescale
zz <- entry( obj, 1, by.id=TRUE )
# Merge with the revised exit times on this timescale
ww <- merge( obj, data.frame( lex.id = as.numeric(names(zz)),
                                cens = zz+cens ) )
# Only retain records with an entry time prior to the revised exit time
ww <- ww[ww[,tS[1]] < ww$cens,]
# Revise the duration according the the revised exit time
x.dur <- pmin( ww$lex.dur, ww[,"cens"]-ww[,tS[1]] )
# Change lex.Xst to lex.Cst for those with shortened follow-up
ww$lex.Xst[x.dur<ww$lex.dur] <- ww$lex.Cst[x.dur<ww$lex.dur]
# Insert the updated follow-yp time
ww$lex.dur <- pmin( x.dur, ww$lex.dur )
ww
}


######################################################################
simLexis <-
function( Tr, # List of lists of transition objects
        init, # Lexis object of persons to simulate.
           N = 1, # No. persons simulated per line in init
      lex.id,
     t.range = 20, # Range for rate computation in the simulation
       n.int = 101, # length of time intervals
    time.pts = seq(0,t.range,length.out=n.int)
         )
{
# Expand the input data frame using N and put in lex.id
if( time.pts[1] !=0 )
    stop( "First time point must be 0, time.pts[1:3]= ",
          time.pts[1:3] )

# Expand init
if( !missing(N) )
  {
  if( length(N) == 1 )
       init <- init[rep(1:nrow(init),each=N),]
  else init <- init[rep(1:nrow(init),     N),]
  }
# and update lex.id if necessary
if( !missing(lex.id) )
  {
  if( length(lex.id)==nrow(init) )
       init$lex.id <- lex.id
  else init$lex.id <- 1:nrow(init)
  }
else   init$lex.id <- 1:nrow(init)

# Check/fix attributes
if( is.null( tS <- attr(init,"time.scales") ) )
  stop( "No time.scales attribute for init" )
if( is.null( tF <- attr(init,"time.since") ) )
  {
  attr(init,"time.since") <- tF <- rep( "", tS )
  warning( "'time.since' attribute set to blanks" )
  }

# Convenience constants
np <- length( time.pts )
tr.st <- names( Tr )

# Set up a NULL object to hold the follow-up records
sf <- NULL

# Take as initiators only those who start in a transient state
nxt <- init[init$lex.Cst %in% tr.st,]

# If some are not in a transient state then say so
if( nrow(nxt) < nrow(init) )
  {
  tt <- table(init$lex.Cst)
  tt <- tt[tt>0]
  nt <- length(tt)
  warning("\nSome initiators start in a absorbing state\n",
          "Initiator states represented are: ",
          paste( rbind( names(tt), rep(":",nt),
                        paste(tt), rep(" ",nt) ), collapse="" ), "\n",
          "Transient states are: ", paste( names( Tr ), coll=" " ) )
  if( nrow(nxt)==0 ) stop( "\nNo initiators in transient states!" )
  }

# Then we update those who are in a transient states and keep on doing
# that till all are in absorbing states or censored
while( nrow(nxt) > 0 )
{
nx <- do.call( rbind.data.frame,
               lapply( split( nxt,
                              nxt$lex.Cst ),
                       simX,
                       Tr, time.pts, tS ) )
sf <- rbind.data.frame( sf, nx )
nxt <- get.next( nx, tr.st, tS, tF )
}

# Doctor lex.Xst levels, fix values for the censored
sf$lex.Xst <- factor( sf$lex.Xst, levels=levels(sf$lex.Cst) )
sf$lex.Xst[is.na(sf$lex.Xst)] <- sf$lex.Cst[is.na(sf$lex.Xst)]

# Nicely order the output by persons, then times and states
nord <- match( c( "lex.id", tS,
                  "lex.dur",
                  "lex.Cst",
                  "lex.Xst" ), names(sf) )
noth <- setdiff( 1:ncol(sf), nord )
sf <- sf[order(sf$lex.id,sf[,tS[1]]),c(nord,noth)]
rownames(sf) <- NULL
# Finally, supply attributes - note we do not supply the "breaks"
# attribute as this is irrelevant for simulated objects
attr( sf, "time.scales" ) <- tS
attr( sf, "time.since"  ) <- tF
chop.lex( sf, tS, max(time.pts) )
}


######################################################################
nState <-
function ( obj,
            at,
          from,
     time.scale = 1 )
{
# Counts the number of persons in each state of the Lexis object 'obj'
# at the times 'at' from the time 'from' in the time scale
# 'time.scale'

# Determine timescales and absorbing and transient states
tS <- check.time.scale(obj,time.scale)
TT <- tmat(obj)
absorb <- rownames(TT)[apply(!is.na(TT),1,sum)==0]
transient <- setdiff( rownames(TT), absorb )

# Expand each record length(at) times
tab.frm <- obj[rep(1:nrow(obj),each=length(at)),
               c(tS,"lex.dur","lex.Cst","lex.Xst")]

# Stick in the corresponding times on the chosen time scale
tab.frm$when <- rep( at, nrow(obj) ) + from

# For transient states keep records that includes these points in time
tab.tr <- tab.frm[tab.frm[,tS]                 <= tab.frm$when &
                  tab.frm[,tS]+tab.frm$lex.dur >  tab.frm$when,]
tab.tr$State <- tab.tr$lex.Cst

# For absorbing states keep records where follow-up ended before
tab.ab <- tab.frm[tab.frm[,tS]+tab.frm$lex.dur <= tab.frm$when &
                  tab.frm$lex.Xst %in% absorb,]
tab.ab$State <- tab.ab$lex.Xst

# Make a table using the combination of those in transient and
# absorbing states.
with( rbind( tab.ab, tab.tr ), table( when, State ) )
}


######################################################################
pState <-
function( nSt, perm=1:ncol(nSt) )
{
# Compute cumulative proportions of persons across states in order
# designate by 'perm'
tt <- t( apply( nSt[,perm], 1, cumsum ) )
tt <- sweep( tt, 1, tt[,ncol(tt)], "/" )
class( tt ) <- c("pState","matrix")
tt
}


######################################################################
plot.pState <-
function( x,
        col = rainbow(ncol(x)),
     border = "transparent",
       xlab = "Time",
       ylim = 0:1,
       ylab = "Probability", ... )
{
# Function to plot cumulative probabilities along the time scale.
matplot( as.numeric(rownames(x)), x, type="n",
         ylim=ylim, yaxs="i", xaxs="i",
         xlab=xlab, ylab=ylab, ... )
lines.pState( x,
        col = col,
     border = border, ... )
}


######################################################################
lines.pState <-
function( x,
        col = rainbow(ncol(x)),
     border = "transparent", ... )
{
# Function to plot cumulative probabilities along the time scale.

# Fixing the colors:
nc <- ncol(x)
col    <- rep( col   , nc )[1:nc]
border <- rep( border, nc )[1:nc]

# Just for coding convenience when plotting polygons
pSt <- cbind( 0, x )
for( i in 2:ncol(pSt) )
   {
   polygon( c(    as.numeric(rownames(pSt)) ,
              rev(as.numeric(rownames(pSt))) ),
            c(    pSt[,i  ],
              rev(pSt[,i-1]) ),
            col=col[i-1], border=border[i-1], ... )
   }
}
