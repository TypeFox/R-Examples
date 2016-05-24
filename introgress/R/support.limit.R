support.limit <-
function(h, geno, locustype, r, s, alleles, max.obj){
  -max.obj + like.h(geno, locustype, r, s, alleles, h) + 2
}

