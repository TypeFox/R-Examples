recolor <- function(id.class, id.target) {
  ##
  ## This function colors id.target in accordance with the most likely candidate
  ## in id.class. It returns a list as id.trcl (which is a factor version of
  ## id.class) and id.prcl (which is a factored version of the colored id.target)
  ## written Ranjan Maitra, Ames, IA 50011-1210, 2015/10/17
  ##
  fac.cl <- as.factor(id.class)
  fac.tg <- as.factor(id.target)
  fac.new <- fac.tg
  trlevels <- levels(fac.cl)
  prlevels <- levels(fac.tg)
  x.tab <- table(fac.cl, fac.tg)
  for (i in 1:length(trlevels)) {
    newcl <- which.max(as.vector(x.tab[i,]))
    fac.new[fac.tg == prlevels[newcl]] <- trlevels[i]
    fac.new[fac.tg == prlevels[i]] <- trlevels[newcl]
    fac.tg <- fac.new
    prlevels <- levels(fac.tg)
    x.tab <- table(fac.cl, fac.tg)
  }
  return(list(id.trcl = fac.cl, id.prcl = fac.tg))
}

