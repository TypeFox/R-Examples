do.rose.diag <-
function( Sector, freqs, shrink=1/2, col='salmon', prop=1,
rotation=NULL)
{
### works with output from kf.sort
### constructs rose.diagram for frequencies per sector
###
# warning("This function only runs after kf.sort has been run")
n.sectors<- max( Sector)
angle1<- 360/n.sectors
circ.freqs<- table( Sector, freqs)[,2]
angles<- angle1*(0:n.sectors)[-(n.sectors+1)]
vals.angles<- circular(rad( rotation+rep( angles, circ.freqs)), 
units='radians',
                       template='none')
### radians transformed to circular object
rose.diag(vals.angles, bins=n.sectors,shrink=shrink, col=col, prop=prop)
invisible( circ.freqs)
}
