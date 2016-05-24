"mesh.intersect" <-
function (p, regionA, regionB, ...) 
matmax(regionA(p, ...), regionB(p, ...))
