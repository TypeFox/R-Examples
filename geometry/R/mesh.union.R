"mesh.union" <-
function (p, regionA, regionB, ...) 
matmin(regionA(p, ...), regionB(p, ...))
