
tol3d <- function(X,Y,TOL,RHO,levels=c(0.05),raise=0.1,col="white",lwd=2){
    gridx <- sort(rep(X,length(X)))
    gridy <- rep(Y,length(Y))
    levelseq <- 1:length(levels)
    
    contour.lines <- contourLines(X,Y,TOL,levels=levels)

    for(i in 1:length(contour.lines)){
        temp.tolcoords <- tol3d.component(data.frame(cbind(contour.lines[i][[1]]$x,contour.lines[i][[1]]$y)),gridx,gridy,RHO,raise)
        levelnum <- which(levels==contour.lines[i][[1]]$level)
        lines3d(contour.lines[i][[1]]$x,contour.lines[i][[1]]$y,temp.tolcoords,lwd=lwd,col=col)
    }
}
