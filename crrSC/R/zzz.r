#.First.lib <- function(lib, pkg)
#{
#    library.dynam("crrSC", pkg, lib)
#}

# use .onLoad instead of .First.lib 
.onLoad <- function(lib, pkg) {
  library.dynam("crrSC", pkg, lib)
}
#end of .onLoad

