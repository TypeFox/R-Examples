#' @importFrom rgl rgl.triangles
#' @importFrom geometry surf.tri
#' @importFrom geometry delaunayn
plot3dtriangles <- function(m,...){
	tm <- try(t(surf.tri(m, delaunayn(m, options="Pp"))))
	if(!inherits(tm, "try-error"))
		rgl.triangles(m[tm,1], m[tm,2], m[tm,3],...)	
	else
		warning("Problem in surf.tri, tree skipped")
}
