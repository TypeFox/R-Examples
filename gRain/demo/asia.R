##
## gRain demo
## The Asia (chest clinic) example
##

library(gRain)

#' Build Bayesian network from conditional probability tables
#' ----------------------------------------------------------

#' Define conditional probability tables (CPTs)
#'
yn <- c("yes","no")
a    <- cptable(~asia, values=c(1,99), levels=yn)
t.a  <- cptable(~tub | asia, values=c(5,95,1,99), levels=yn)
s    <- cptable(~smoke, values=c(5,5), levels=yn)
l.s  <- cptable(~lung | smoke, values=c(1,9,1,99), levels=yn)
b.s  <- cptable(~bronc | smoke, values=c(6,4,3,7), levels=yn)
e.lt <- cptable(~either | lung:tub,values=c(1,0,1,0,1,0,0,1),
                levels=yn)
x.e  <- cptable(~xray | either, values=c(98,2,5,95), levels=yn)
d.be <- cptable(~dysp | bronc:either, values=c(9,1,7,3,8,2,1,9),
                levels=yn)

#' Collect CPTs, do some internal checks
cpt.list <- compileCPT(list(a, t.a, s, l.s, b.s, e.lt, x.e, d.be))

#' cpt.list is a list of multidimensional arrays
cpt.list
cpt.list$tub
cpt.list$either

#' Build network
bn <- grain( cpt.list )
bn

if (require(Rgraphviz))
    plot(bn) # plot dag

#' Compile network (moralize, triangulate, establish clique potential
#' representation)
bn <- compile(bn)
bn

#' Notice for a compiled network we can get these plots
if (require(Rgraphviz)){
    par(mfrow=c(1,2))
    plot(bn) # triangulated moral graph
    plot(bn, type="dag")
}

#' Find clique marginal representation.
bn <- propagate(bn); bn

#' Set evidence and query network
#' ------------------------------

qnodes <- c("lung","bronc","tub")

querygrain(bn, nodes=qnodes)
querygrain(bn, nodes=qnodes, type="joint")

#' Evidence (hard evidence) can be given in different forms:
#'

x1 <- setEvidence( bn, nodes=c("asia","dysp"), states=c("yes","yes"))
x1
getEvidence(x1); pEvidence(x1)
querygrain(x1, nodes=qnodes)

x2 <- setEvidence( bn, nslist=list(asia="yes", dysp="yes"))
x2
getEvidence(x2); pEvidence(x2)
querygrain(x2, nodes=qnodes)

#' Soft evidence (likelihood evidence) can be given as
#'

x3 <- setEvidence( bn, nslist=list(asia=c(.7, .1), dysp="yes"))
x3
getEvidence(x3); pEvidence(x3)
querygrain(x3, nodes=qnodes)

#' A small shortcut: A typical use of gRain involves first setting
#' evidence and then querying nodes. This can be achieved in one step
#'

querygrain(bn, nodes=qnodes, nslist=list(asia=c(.7, .1), dysp="yes"))

#' Incremental evidence
#'

x4  <- setEvidence2( bn, evidence=list(asia=c(1,0)))
x4
x5  <- setEvidence2( x4, evidence=list(dysp="yes"))
x5

#' No effect below; asia is already set
x6 <- setEvidence2(x5, evidence=list(asia="yes"))
getEvidence(x3); pEvidence(x3)

#' retractEvidence
#'

x7 <- retractEvidence(x6, "asia")
getEvidence(x7)


#' absorbEvidence
#' --------------
#'

#' When evidence is entered, the clique potentials are changed,
#' i.e. we obtain a new network It is sometimes useful to make this
#' new network act as if no evidence has been entered.

x66 <- absorbEvidence( x6 )

getEvidence( x6 )
getEvidence( x66 )

querygrain( x6,  nodes=qnodes )
querygrain( x66, nodes=qnodes )

querygrain( x6,  nodes=c("asia","dysp") )
querygrain( x66, nodes=c("asia","dysp") )



#' Multivariate evidence
#' ---------------------
#' This is a recent addition and the implementation is somewhat experimental

t.db <- parray(c("dysp","bronc"), list(yn,yn), values=c(.1,.2,.9,.8))
ev   <- list(asia=c(1,0), dysp="yes", t.db)
ev

#' NOTICE: setMEvidence is used here:
x8 <- setMEvidence(bn, evidence=ev, details=1)
getEvidence(x1); pEvidence(x1)

#' remove 3rd evidence element (NOTICE: retractMEvidence is used here):
#' NOTICE: We can not refer to a node here, must refer to the elements in
#' the evidence list.
x9 <- retractMEvidence(x8, 3)
x9
getEvidence(x9); pEvidence(x9)

x10 <- setMEvidence(x9, evidence=list(t.db))
x10
getEvidence(x10)

#' Evidence can just be added as pleased
x11 <- setMEvidence(x10, evidence=ev[1])
getEvidence(x11); pEvidence(x11)

#' Hence we can add evidence that leads to a zero probability
x12 <- setMEvidence(x11, evidence=list(asia=c(0,1)), details=1)
getEvidence(x12); pEvidence(x12)

#' Some details
#' ------------
#'

#' This may change, but at least this is the current implementation
#' (March 2015).
#'
#' Internally clique potentials are storede in the slot temppot and
#' the clique marginals in equipot

#' Before compilation neither of these exist

bn <- grain( cpt.list )
bn
getgrain(bn, "temppot")
getgrain(bn, "equipot")

#' After compilation, the former has the clique potentials and the latter has NAs
bn <- compile( bn )
getgrain(bn, "temppot")

#' These potentials are not probabilities!
lapply( getgrain(bn, "temppot"), sum )

#' After propagation, equipot has clique marginals
bn <- propagate( bn )
getgrain(bn, "temppot")
lapply( getgrain(bn, "temppot"), sum )
lapply( getgrain(bn, "equipot"), sum )

