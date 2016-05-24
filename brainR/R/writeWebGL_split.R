#' Write WebGL with split triangles
#'
#' Adapted \link[rgl]{writeWebGL} function that splits the triangles into 
#' 65535 vertices
#' 
#' @param ids - rgl ids (see \link[rgl]{rgl.ids})
#' @param writeIt - (logical) write the file out
#' @param verb - verbose output
#' @param ... - further arguments passed to \link[rgl]{writeWebGL}
#' @export
#' @examples
#' \dontrun{
#' #Brain Template from Copyright (C) 1993-2009 Louis Collins, 
#' #McConnell Brain Imaging Centre, 
#' #Montreal Neurological Institute, McGill University
#' #6th generation non-linear symmetric brain
#' template <- readNIfTI(system.file("MNI152_T1_2mm_brain.nii.gz", package="brainR")
#' , reorient=FALSE) 
#' dtemp <- dim(template)
#' ### 4500 - value that empirically value that presented a brain with gyri
#' ### lower values result in a smoother surface
#' brain <- contour3d(template, x=1:dtemp[1], y=1:dtemp[2], 
#' z=1:dtemp[3], level = 4500, alpha = 0.1, draw = FALSE)
#' drawScene.rgl(brain)
#' ### this would be the ``activation'' or surface you want to render - 
#' # hyper-intense white matter
#' contour3d(template, level = c(8200, 8250), 
#' alpha = c(0.5, 0.8), add = TRUE, color=c("yellow", "red"))
#' ### add text
#' text3d(x=dtemp[1]/2, y=dtemp[2]/2, z = dtemp[3]*0.98, text="Top")
#' text3d(x=dtemp[1]*0.98, y=dtemp[2]/2, z = dtemp[3]/2, text="Right")
#' fname <- "knitted_webGL.html"
#' writeWebGL_split(dir=getwd(), filename =fname, 
#' template = system.file("my_template.html", package="brainR"), width=500, 
#' writeIt=TRUE)
#' browseURL(fname)
#' }
#' @return if writeIt is TRUE, then returns the value from \link[rgl]{writeWebGL}.
#' Otherwise, returns the split triangles from the rgl objects


writeWebGL_split <- function(ids=rgl.ids()$id, writeIt= TRUE, verb=FALSE, ...){

	if (verb) print("Splitting Triangles")
	split_triangles <- function(ids = ids, maxsize=65535) {
	
		if (maxsize %% 3 != 0)
			stop("maxsize must be a multiple of 3")
	
		save <- par3d(skipRedraw=TRUE)
		on.exit(par3d(save))
	
		allids <- rgl.ids()
		ids <- with(allids, id[ id %in% ids & type == "triangles" ])
		for (id in ids) {
			count <- rgl.attrib.count(id, "vertices")
			if (count <= maxsize) next
			verts <- rgl.attrib(id, "vertices")
			norms <- rgl.attrib(id, "normals")
			cols <- rgl.attrib(id, "colors")
		
			rgl.pop(id=id)
			while (nrow(verts) > 0) {
				n <- min(nrow(verts), maxsize)
				triangles3d(verts[1:n,], normals=norms[1:n,], color=rgb(cols[1:n,1], cols[1:n,2], cols[1:n,3]), alpha=cols[1:n,4])
				verts <- verts[-(1:n),,drop=FALSE]
				norms <- norms[-(1:n),]
				cols <- cols[-(1:n),]
			}
		}
	}
	split_triangles(ids)
	if (writeIt) rgl::writeWebGL(...)
}

