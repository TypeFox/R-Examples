s.wrapper <-
function(a, b, max.obj, geno, locustype, r, s, alleles){
  f.a <- support.limit(h=a, geno, locustype, r, s, alleles, max.obj)
  f.b <- support.limit(h=b, geno, locustype, r, s, alleles, max.obj)
  f.a * f.b <= 0
}

