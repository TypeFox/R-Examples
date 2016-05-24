`BOXarrows3D` <-
function(x1,y1,z1,x2,y2,z2, aglyph=NULL, Rview=ROTX(0), col=grey(.5),  border='black',
                      len = .7, basethick=.05, headlen=.3,headlip=.02)
  {
    if(missing(len)) { len = .7}
    if(missing(basethick)) {basethick=.05 }
    if(missing(headlip)) { headlip=.02}
    
    if(missing(headlen)) { headlen=.3 }

    if(missing(col)) {col=grey(.5) }
    if(missing(border)) { border='black'}
    
    if(missing(Rview)) { Rview=ROTX(0)}
    if(missing(aglyph)) {aglyph=NULL }
    
    up = par("usr")
    ui = par("pin")
    
    alpha = 65

    ratx = (up[2]-up[1])/ui[1]
    raty=  (up[4]-up[3])/ui[2]
    ##   if(length(siz)<length(x)) { siz=rep(siz,length=length(x)) }
     if(length(col)<length(x1)) { col=rep(col,length=length(x1)) }
    
    ##     thick=.4; headlength=4; headthick=1
    
    ##  x1=-12;y1=-10;x2=12; y2=3
    if(is.null(aglyph))
      aglyph = Z3Darrow(len = len , basethick =basethick , headlen =headlen , headlip=headlip )
    myzee = matrix(c(0,0,1, 1), nrow=1, ncol=4) 
    for(i in 1:length(x1))
      {
        T = TRANmat(x1[i],y1[i],z1[i] )
        thelen = sqrt((y2[i]-y1[i])^2+ (x2[i]-x1[i])^2+(z2[i]-z1[i] )^2)
        ##  T[4,4 ] = thelen
        SC = diag(1,nrow = 4)
        SC[3,3 ] = thelen
        gamma = 180*atan2(y2[i]-y1[i], x2[i]-x1[i] )/pi
        beta  = 180*atan2( sqrt( (x2[i]-x1[i])^2+(y2[i]-y1[i])^2),   z2[i] -  z1[i] )/pi
        gamma =  RPMG::fmod(gamma, 360)
        beta = RPMG::fmod(beta, 360)
        R3 = ROTZ(gamma)
        R2 = ROTY(beta)
        R1 = ROTZ(alpha)
        M =      SC %*%  R1  %*% R2  %*%  R3  %*% T  %*% Rview
        M2 =       R1  %*% R2  %*%  R3 %*% Rview
### print(paste("_________________________________________", i, j))
        pglyph3D(aglyph$aglyph, anorms=aglyph$anorm  , M=M, M2=M2, zee=myzee , col=col[i] )
### text(L$x1[i],L$y1[i], labels=i)
###  title(main=j)
      }
  }

