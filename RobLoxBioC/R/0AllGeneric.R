############# preparations ################
if(!isGeneric("robloxbioc")){
    setGeneric("robloxbioc", 
        function(x, ...) standardGeneric("robloxbioc"))
}

if(!isGeneric("KolmogorovMinDist")){
    setGeneric("KolmogorovMinDist", 
        function(x, D, ...) standardGeneric("KolmogorovMinDist"))
}
