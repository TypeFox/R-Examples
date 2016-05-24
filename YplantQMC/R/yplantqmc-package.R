

#'Two example hemiphotos in YplantQMC format
#'
#'A canopy with a small gap, and a canopy with a large gap. Both are objects of
#'class \code{yphemi}, which is normally constructed with the function
#'\code{\link{setHemi}}.
#'
#'
#'@name examplehemi
#'@aliases smallgap largegap examplehemi
#'@docType data
#'@source Thanks to Daniel Falster and Bob Pearcy.
#'@keywords datasets
#'@examples
#'
#'
#'\dontrun{
#'# Small gap
#'plot(smallgap)
#'
#'# Large gap, with solar path in June in the south of france.
#'southfrance <- setLocation(44)
#'juneday <- setMet(southfrance, month=6, day=21)
#'plot(largegap, juneday)
#'}
#'
NULL





#'Three example 3D plants
#'
#'Three virtual plant, in \code{YplantQMC} format. Objects of this kind would
#'normally be produced with \code{\link{constructplant}}.
#'
#'
#'@name plantexamples
#'@aliases plantexamples toona sugarmaple pilularis
#'@docType data
#'@source Plants provided by Jeff Kelly and Kerrie Sendall.
#'@keywords datasets
#'@examples
#'
#'
#'\dontrun{
#'# A sugar maple (Acer saccharum)
#'plot(sugarmaple)
#'
#'# A blackbutt (Eucalyptus pilularis)
#'plot(pilularis)
#'
#'# Australian red cedar (Toona australis)
#'plot(toona)
#'}
#'
#'
NULL





#'A turtle sky with 58 points
#'
#'These are the angles used in \code{\link{STARbar}} when \code{integration =
#'"Turtlesky"}.
#'
#'
#'@name turtle
#'@docType data
#'@format A data frame with 59 observations on the following 2 variables.
#'\describe{ \item{altitude}{a numeric vector} 
#'\item{azimuth}{a numeric vector} }
#'@keywords datasets
NULL





#'A turtle sky with 244 points
#'
#'Not currently used.
#'
#'
#'@name turtle244
#'@docType data
#'@format A data frame with 244 observations on the following 2 variables.
#'\describe{ 
#'\item{altitude}{a numeric vector} 
#'\item{azimuth}{a numeric vector} }
#'@keywords datasets
NULL





#'A turtle sky with 482 points
#'
#'A set of points that are uniformly distributed across the hemisphere, so that
#'each point represents (approx.) the same solid angle. Used in the diffuse
#'radiation calculations in \code{\link{YplantDay}}.
#'
#'
#'@name turtle482
#'@docType data
#'@format A data frame with 482 observations on the following 2 variables.
#'\describe{ 
#'\item{altitude}{a numeric vector} 
#'\item{azimuth}{a numeric vector} }
#'@keywords datasets
NULL





#'Altitude and azimuth angles
#'
#'The 160 angles used in the original Yplant (Pearcy & Yang 1996), in radians.
#'Used in \code{\link{STARbar}} when \code{integration = "Yplant"}.
#'
#'
#'@name yplantaltaz
#'@docType data
#'@format A data frame with 160 observations on the following 2 variables.
#'\describe{ 
#'\item{azimuth}{Azimuth angle (radians)}
#'\item{altitude}{Altitude (radians)} }
#'@keywords datasets
NULL



