aoc <-
function(veg,o.rgr,o.sgr){
#
# in this version of 10.10.12 o.rlab is calculated internally
  o.rlab<- rownames(veg)
  o.slab<- colnames(veg)
  o.rlab<- o.rlab[order(o.rgr)]
  o.slab<- o.slab[order(o.sgr)]
  o.rlab<- order(order(o.rlab))
  o.slab<- order(order(o.slab))

  nrel <- nrow(veg)
  nspec <- ncol(veg)
  nele <- nrel*nspec                                 # no. of elements in the entire table
  xall <- as.numeric(as.matrix(t(veg),nrow=nele))    # entire matrix in one dimension
  rowgroups <- rep(o.sgr,nrel)
  col <- c(rep(o.rgr,nspec)) 
  col <- matrix(col,nrow=nrel,ncol=nspec)
  col <- t(col)
  colgroups <- matrix(col,nrow=nele)
  F <- tapply(sign(xall),list(rowgroups, colgroups),sum)      # contingency table, counts
  Z <- tapply(sign(xall+1),list(rowgroups, colgroups),sum)    # contingency table, all elements
  A <- F/Z*sum(F)/sum(F/Z)                                    # adjusting to equal block size
  nrel<- max(as.integer(o.rgr))  
  nspec<- max(as.integer(o.sgr))
  rownames(A) <- seq(1,nspec,1)
  colnames(A) <- seq(1,nrel,1)
  aoc <- cca(A,diss=FALSE)                                    # CA, rows=spec.gr, col=rel.gr
# v are releve groups
# u are species groups 
# reorder table. This should yield the new releve and species order.
# But it does not!

  gr <- as.integer(o.rgr[o.rlab])
  grr<- gr
  k.r<- length(levels(o.rgr))       # no of releve groups
  k.s<- length(levels(o.sgr))       # no of species groups

  for(i in 1:k.r) gr[gr == i] <- aoc$CA$v[i,1]

  grr<-grr[order(gr)]
  gs <- as.integer(o.sgr[o.slab])
  gss<- gs

  for(i in 1:k.s) gs[gs == i] <- aoc$CA$u[i,1]

  gss<-gss[order(gs)]


# Next line would be for printint a table:
# out1 <- vegemite(veg,sp.ind=o.slab[order(gs)], site.ind=o.rlab[order(gr)], scale="Br", zero="-")
  m<- min(k.r,k.s)
  mscc<- aoc$CA$tot.chi/(aoc$CA$grand.total*(m-1))
  rorder<- o.rlab[order(gr)]               # final order of releves
  sorder<- o.slab[order(gs)]               # final order of species
# mean square contingency coefficient MSCC
  eigval<- aoc$CA$eig
  rank<- aoc$CA$rank
  gtot<- sum(F)
  cancorr<- eigval^0.5
  chisquare<- eigval*gtot
  mscc<- sum(chisquare)/(gtot*rank)
# list output
  outaoc<- list(rank=aoc$CA$rank,rgrscores=aoc$CA$v,sgrscores=aoc$CA$u,eigval=aoc$CA$eig,chisquare=chisquare,cancorr=cancorr,grand.total=gtot,MSCC=mscc,new.relorder=rorder,new.sporder=sorder,new.relgr=grr,new.spgr=gss,cont.table=A)
  }
