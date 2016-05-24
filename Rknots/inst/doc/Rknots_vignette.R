## ----eval=TRUE-----------------------------------------------------------
library(Rknots)
str( head(Rolfsen.table, 5) )
head( names(Rolfsen.table) )

## ----eval=TRUE-----------------------------------------------------------
str( Rolfsen.table['3.1'] )

## ----eval=TRUE-----------------------------------------------------------
str( head(link.table, 5) )
head( names(link.table) )

## ----eval=TRUE-----------------------------------------------------------
set.seed(123)

knot <- makeExampleKnot(k = TRUE) #a knot
link <- makeExampleKnot(k = FALSE) #a link

## ----eval=TRUE-----------------------------------------------------------
knot.cls <- new('Knot', points3D = knot)
link.cls <- new('Knot', points3D = link$points3D, ends = link$ends)

( knot.cls <- newKnot(points3D = knot) )
( link.cls <- newKnot(points3D = link$points3D, ends = link$ends) )

## ----eval=TRUE-----------------------------------------------------------
head( getCoordinates(knot.cls), 5 )
getEnds(knot.cls)

head( link.cls['points3D'], 5)
link.cls['ends']

## ----eval=TRUE-----------------------------------------------------------
knot.bu <- knot.cls #copy
new.coordinates <- matrix( runif(60), ncol = 3 ) #generate random coordinates
setCoordinates(knot.cls) <- new.coordinates #replace
knot.cls
knot.cls <- knot.bu #restore

link.bu <- link.cls #copy
setEnds(link.cls) <- c(10, 50, 90) #replace separators
getEnds(link.cls)
link.cls <- link.bu #restore

## ----eval=TRUE-----------------------------------------------------------
knot.AB <- AlexanderBriggs(points3D = knot, ends = c())
str(knot.AB)
knot.msr <- msr(points3D = knot, ends = c())
str(knot.msr)

link.AB <- AlexanderBriggs(points3D = link$points3D, ends = link$ends)
str(link.AB)
link.msr <- msr(points3D = link$points3D, ends = link$ends)
str(link.msr)

## ----eval=TRUE-----------------------------------------------------------
knot.cls.AB <- reduceStructure(knot.cls, algorithm = 'AB' )
knot.cls.MSR <- reduceStructure(knot.cls, algorithm = 'MSR' )
link.cls.AB <- reduceStructure(link.cls, algorithm = 'AB' )
link.cls.MSR <- reduceStructure(link.cls, algorithm = 'MSR' )
link.cls.AB

## ----eval=TRUE-----------------------------------------------------------
par( mfrow=c(1,2) )

#plotDiagram(knot.AB$points3D, knot.AB$ends, lwd = 2, main = 'Alexander-Briggs')
#plotDiagram(knot.msr$points3D, knot.msr$ends, lwd = 2, main = 'MSR')
plotDiagram(link.AB$points3D, link.AB$ends, lwd = 2, main = 'Alexander-Briggs')
plotDiagram(link.msr$points3D, link.msr$ends, lwd = 2, main = 'MSR')

## ----eval=TRUE-----------------------------------------------------------
plot(link.cls.AB, lend = 2, lwd = 3, main = 'using par()')

## ----eval=TRUE-----------------------------------------------------------
( delta.k <- computeInvariant( knot.cls.AB, invariant = 'Alexander') )
jones.k <- computeInvariant( knot.cls.AB, invariant = 'Jones', skein.sign = -1)
homfly.k <- computeInvariant( knot.cls.AB, invariant = 'HOMFLY', skein.sign = -1)

( delta.l <- computeInvariant( link.cls.AB, invariant = 'Alexander') )
jones.l <- computeInvariant( link.cls.AB, invariant = 'Jones', skein.sign = -1)
homfly.l <- computeInvariant( link.cls.AB, invariant = 'HOMFLY', skein.sign = -1)

## ----eval=TRUE-----------------------------------------------------------
converted <-HOMFLY2Jones( homfly.k )
identical( converted, jones.k)

## ----eval=TRUE-----------------------------------------------------------
( computeInvariant( link.cls.AB, invariant = 'LK') )

## ----eval=FALSE----------------------------------------------------------
#  text <- names(Rolfsen.table)[1 : 6]
#  
#  par(mfrow = c(3,2))
#  for(i in 1 : 6) {
#    k <- Rolfsen.table[[i]]
#    k <- newKnot(k)
#    plot(k, lwd = 2, main = text[i],
#         sub = paste( computeInvariant(k, 'Alexander'),
#                      computeInvariant(k, 'Jones'),
#                      computeInvariant(k, 'HOMFLY'), sep = '\n'))
#  }

## ----eval=TRUE-----------------------------------------------------------
protein <- loadProtein(system.file("extdata/2K0A.pdb", package="Rknots"))
protein<- loadProtein('2K0A') #from the PDB
str(protein)

## ----eval=FALSE----------------------------------------------------------
#  ramp <- colorRamp(c('blue', 'white', 'red'))
#  pal <- rgb( ramp(seq(0, 1, length = 109)), max = 255)
#  
#  plotKnot3D(protein$A, colors = list( pal ),
#  	lwd = 8, radius = .4, showNC = TRUE, text = TRUE)

## ----eval=TRUE-----------------------------------------------------------
protein <- newKnot(protein$A)
protein

## ----eval=TRUE-----------------------------------------------------------
protein.cp <- closeAndProject( protein, w = 2 )

## ----eval=FALSE----------------------------------------------------------
#  par(mfrow = c(1,2))
#  plot(protein, main = 'original', lwd = 2)
#  plot(protein.cp, main = 'closed & projected', lwd = 2)

## ----eval=TRUE-----------------------------------------------------------
( delta.p <- computeInvariant( protein.cp, invariant = 'Alexander') )
( jones.p <- computeInvariant( protein.cp, invariant = 'Jones', skein.sign = -1) )
( homfly.p <- computeInvariant( protein.cp, invariant = 'HOMFLY', skein.sign = -1) )

## ----eval=TRUE-----------------------------------------------------------
getKnotType(polynomial = delta.p, invariant = 'Alexander')
getKnotType(polynomial = homfly.p, invariant = 'HOMFLY')
getKnotType(polynomial = homfly.p, invariant = 'HOMFLY', full.output = TRUE)

## ----eval=TRUE-----------------------------------------------------------
trefoil <- Rolfsen.table[[1]]
trefoil <- newKnot(trefoil)
( homfly.tr <- computeInvariant(trefoil, 'HOMFLY') )
( homfly.tl <- HOMFLY2mirror(homfly.tr) )
identical( homfly.p, homfly.tl )

## ----eval=TRUE-----------------------------------------------------------
processChain <- function(protein, i) {
	chain <- newKnot(protein[[i]])
	chain <- closeAndProject( chain )
	return( computeInvariant(chain, 'HOMFLY') )
}        

lengthChain <- function(protein, i) return( nrow(protein[[i]]))

protein <- loadProtein('1AJC', join.gaps = FALSE, cutoff = 7)
str(protein)
chains <- names(protein)        
polynomials <- lapply( 1: length(chains) , 
		function(i) {
			ifelse(lengthChain(protein, i) > 6, processChain(protein, i), 1) } )
cbind(chains, polynomials)     

## ----eval=TRUE-----------------------------------------------------------
protein <- loadProtein('1AJC', join.gaps = TRUE)
str(protein)
chains <- names(protein)        
polynomials <- lapply( 1: length(chains) , 
		function(i) {
			ifelse(lengthChain(protein, i) > 6, processChain(protein, i), 1) } )
cbind(chains, polynomials)        

## ----eval=TRUE-----------------------------------------------------------
sessionInfo()

