### R code from vignette source 'IBDLabels_tutorial.Rnw'

###################################################
### code chunk number 1: install (eval = FALSE)
###################################################
## install.packages( "IBDLabels", repos = "http://R-Forge.R-project.org")


###################################################
### code chunk number 2: library
###################################################
library("IBDLabels")


###################################################
### code chunk number 3: tutorial (eval = FALSE)
###################################################
## vignette("IBDLabels_tutorial") 


###################################################
### code chunk number 4: list
###################################################
## list everything in the package
ls( "package:IBDLabels", all=TRUE)
## list all functions and their arguments
lsf.str("package:IBDLabels")


###################################################
### code chunk number 5: vec
###################################################
## Vectors for all labels, note the invalid IBD state produced from 
## label 2 has NA.
allVec( ngam = 4 ) 
       
## Convert vector to label, the vectors are always renumbered with fgl2vec.
vec2label( c(1,1,1,3))
vec2label( c(1,1,1,2))

## Convert label to vector
label2vec( 1, ngam = 4 ) ## Coverts to valid vector
label2vec( 2, ngam = 4 ) ## Converts to invalid vector

## renumbering
fgl2vec( c( 1,1,1,3) )
fgl2vec( c(5,1,5,6) ) 

## maximum label for given number of gametes. 
maxlabel( ngam = 4 ) 
maxlabel( ngam = 5 )
maxlabel( ngam = 6 )


###################################################
### code chunk number 6: lex
###################################################

## Vector of all lexicograghic states with labels ( names of elements )
## labels are printed on the top line
allLex( ngam = 4 )

## Convert lex to label
lex2label( lex = 5 , ngam = 4 )

## Convert label to lex
label2lex( label = 3, ngam = 4 )

## Invalid labels do not convert to lex
label2lex( label = 2, ngam = 4)



###################################################
### code chunk number 7: jaq
###################################################

## list all jacquard states
allJaq()

## Convert label to jacquard
label2jaq( 3, phased = TRUE )
label2jaq( 3, phased = FALSE )

## Convert jacquard to label
jaq2label( 9, phased = TRUE )
jaq2label( 9, phased = FALSE )


###################################################
### code chunk number 8: IBDLabels_tutorial.Rnw:147-151
###################################################

## Lookup table of all states. Always puts NA for invalid.
allStates( ngam = 4 )



