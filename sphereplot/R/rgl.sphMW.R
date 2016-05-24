rgl.sphMW=function (radius = 1, col='purple', type='s', MWcenrad=0.02, addMWplane=TRUE) 
{
    if(addMWplane){rgl.sphcirc(CrossEq = -76.75, PeakDec = 62.6,col=col,radius=radius)}
    MWcenloc = cbind(266.42,-29,radius)
    plot3d(sph2car(MWcenloc,deg=TRUE),type=type,col=col,radius=MWcenrad,add=TRUE)
}

