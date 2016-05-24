#
# map.functions.R
#
# Copyright (c) 2010-2012 GBIC: Danny Arends, Konrad Zych and Ritsert C. Jansen
# last modified Apr, 2013
# first written Aug, 2012
# Contains: avg_map_distance, map_distances
#

avg_map_distances <- function(m){
  if(any(class(m)=="cross")) m <- pull.map(m)
  avg_lengths <- unlist(lapply(m,function(x){(max(x) - min(x)) / length(x)}))
  names(avg_lengths) <- names(nmar(m))
  avg_lengths
}

map_distances <- function(m){
  if(any(class(m)=="cross")) m <- pull.map(m)
  chr_lengths <- unlist(lapply(m,function(x){(max(x) - min(x))}))
  names(chr_lengths) <- names(nmar(m))
  chr_lengths
}
