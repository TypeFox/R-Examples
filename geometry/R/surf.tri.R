##' Find surface triangles from tetrahedra mesh
##' 
##' Find surface triangles from tetrahedron mesh typically obtained with
##' \code{\link{delaunayn}}.
##' 
##' \code{surf.tri} and \code{\link{convhulln}} serve a similar purpose in 3D,
##' but \code{surf.tri} also works for non-convex meshes obtained e.g. with
##' \code{\link{distmeshnd}}.  It also does not produce currently unavoidable
##' diagnostic output on the console as \code{convhulln} does at the Rterm
##' console--i.e., \code{surf.tri} is silent.
##' 
##' @param p An \code{n}-by-\code{3} matrix. The rows of \code{p} represent
##' \code{n} points in \code{dim}-dimensional space.
##' @param t Matrix with 4 columns, interpreted as output of
##' \code{\link{delaunayn}}.
##' @return An \code{m}-by-\code{3} index matrix of which each row defines a
##' triangle. The indices refer to the rows in \code{p}.
##' @note \code{surf.tri} was based on matlab code for mesh of Per-Olof Persson
##' (\url{http://persson.berkeley.edu/distmesh/}).
##' @author Raoul Grasman
##' @seealso \code{\link[tripack]{tri.mesh}}, \code{\link{convhulln}},
##' \code{\link{surf.tri}}, \code{\link{distmesh2d}}
##' @keywords math optimize dplot
##' @examples
##' 
##' \dontrun{
##' # more extensive example of surf.tri
##' 
##' # url's of publically available data:
##' data1.url = "http://neuroimage.usc.edu/USCPhantom/mesh_data.bin"
##' data2.url = "http://neuroimage.usc.edu/USCPhantom/CT_PCS_trans.bin"
##' 
##' meshdata = R.matlab::readMat(url(data1.url))
##' elec = R.matlab::readMat(url(data2.url))$eeg.ct2pcs/1000
##' brain = meshdata$mesh.brain[,c(1,3,2)]
##' scalp = meshdata$mesh.scalp[,c(1,3,2)]
##' skull = meshdata$mesh.skull[,c(1,3,2)]
##' tbr = t(surf.tri(brain, delaunayn(brain)))
##' tsk = t(surf.tri(skull, delaunayn(skull)))
##' tsc = t(surf.tri(scalp, delaunayn(scalp)))
##' rgl::rgl.triangles(brain[tbr,1], brain[tbr,2], brain[tbr,3],col="gray")
##' rgl::rgl.triangles(skull[tsk,1], skull[tsk,2], skull[tsk,3],col="white", alpha=0.3)
##' rgl::rgl.triangles(scalp[tsc,1], scalp[tsc,2], scalp[tsc,3],col="#a53900", alpha=0.6)
##' rgl::rgl.viewpoint(-40,30,.4,zoom=.03)
##' lx = c(-.025,.025); ly = -c(.02,.02);
##' rgl::rgl.spheres(elec[,1],elec[,3],elec[,2],radius=.0025,col='gray')
##' rgl::rgl.spheres( lx, ly,.11,radius=.015,col="white")
##' rgl::rgl.spheres( lx, ly,.116,radius=.015*.7,col="brown")
##' rgl::rgl.spheres( lx, ly,.124,radius=.015*.25,col="black")
##' }
##'
##' @export
"surf.tri" <-
function(p,t){
    # original by Per-Olof Persson (c) 2005 for MATLAB
    # ported to R and modified for efficiency by Raoul Grasman (c) 2005

    # construct all faces
    faces = rbind(t[,-4], t[,-3], t[,-2], t[,-1]);
    node4 = rbind(t[, 4], t[, 3], t[, 2], t[, 1]);

#    #original translated from MATLAB:
#    # select the faces that occur only once --> these are the surface boundary faces
#    faces = t(apply(faces,1,sort));                                         # sort each row
#    foo   = apply(faces,1,function(x) do.call("paste",as.list(x,sep=" "))); # makes a string from each row
#    vec   = table(foo);                                                     # tabulates the number of occurences of each string
#    ix    = sapply(names(vec[vec==1]),function(b) which(b==foo))            # obtain indices of faces with single occurence
#    tri   = faces[ix,];
#    node4 = node4[ix];


    # we wish to achieve
    #   > faces = t(apply(faces,1,sort));
    # but this is much too slow, we therefore use max.col and the fact
    # that there are only 3 columns in faces
    i.max = 3*(1:nrow(faces)-1) + max.col(faces)
    i.min = 3*(1:nrow(faces)-1) + max.col(-faces)
    faces = t(faces)
    faces = cbind(faces[i.min], faces[-c(i.max,i.min)], faces[i.max])
    ix = order(faces[,1], faces[,2], faces[,3])

    # Next, we wish to detect duplicated rows in faces, that is,
    #   > qx = duplicated(faces[ix,],MARGIN=1)              # logical indicating duplicates
    # but this is also much to slow, we therefore use the fact that
    # faces[ix,] has the duplicate rows ordered beneath each other
    # and the fact that each row occurs exactly once or twice
    fo = apply(faces[ix,],2,diff)
    dup = (abs(fo) %*% rep(1,3)) == 0        # a row of only zeros indicates duplicate
    dup = c(FALSE,dup)                       # first is never a duplicate
    qx = diff(dup)==0                        # only zero if two consecutive elems are not duplicates
    qx = c(qx, !dup[length(dup)])            # last row is either non-duplicate or should not be selected
    tri = faces[ix[qx],]                     # ix[qx] are indices of singly occuring faces
    node4 = node4[ix[qx]]

    # compute face orientations
    v1 = p[tri[,2],] - p[tri[,1],]; # edge vectors
    v2 = p[tri[,3],] - p[tri[,1],];
    v3 = p[node4,]   - p[tri[,1],];
    ix = which( apply(extprod3d(v1,v2) * v3, 1, sum) > 0 )
    tri[ix,c(2,3)] = tri[ix,c(3,2)]
    rownames(tri) = NULL
    tri
}
