is.unknown <-
function(x,dag)
{ # internally used by brute.search();
  (dag$names[x]=='unknown' || dag$cov.types[x]==2);
}

