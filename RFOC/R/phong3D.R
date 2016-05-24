`phong3D` <-
function(aglyph,  M=diag(1, nrow=4), M2=diag(1, nrow=4), Light=c(45,45)   , anorms=list() , zee=c(0,0,1), col="white", border="black")
  {

    if(missing(M)) { M = diag(1, nrow=4)  }
    if(missing(M2)) { M2 =  M }
    
    if(missing(col)) { col=rep("white", length(aglyph)) }
    if(missing(border)) { border="black" }
    if(missing(zee)) { zee=c(0,0,1) }
    if(missing(Light)) { Light=c(45, 45) }

    ###########   light = direction of light vector at infinity, light1=azimuth, light2=nadir(zaxis)
    Light = Light*180/pi
    lightvec = -1*c(  cos(Light[1])*sin(Light[2]) ,   sin(Light[1])*sin(Light[2])  ,    cos(Light[2])          )

    
    if(missing(anorms))
      {
        anorms = list()
        for(i in 1:length(aglyph))
          {
            XX = RSEIS::xprod(aglyph[[i]][3,]-aglyph[[i]][2,], aglyph[[i]][2,]-aglyph[[i]][1,])
           ## print(paste(sep=' ', c(i, XX) ))
            anorms[[i]] = XX
          }
      }

    if(length(col)==1) { col=rep(col, length(aglyph))  } 

    bglyph = list()
    bvec =   vector()
    for( i in 1:length(aglyph))
      {
        Xt = cbind(aglyph[[i]], 1)  %*% M
####  get cross product of second vector versus first
#########   vectors should be oriented so this points out from object
###   XX = RSEIS::xprod( Xt[2,]-Xt[1,], Xt[3,]-Xt[2,])
        XX = c(anorms[[i]], 0)  %*% M2
        
        ##  print(XX[3] )
        ## if(XX[3]>0) polygon(Xt[,1], Xt[,2], col="white", border="black")
######   get depth of faces
        zd = mean(Xt[,3])
        bvec[i] = zd
        bglyph[[i]] = list(x=Xt[,1], y= Xt[,2], z= Xt[,3]    , xp=XX, zd=zd  )
      }
 ####   print('###################phonger##################')
    for( i in order(bvec))
      {
        zdot = (bglyph[[i]]$xp[1]*zee[1]+bglyph[[i]]$xp[2]*zee[2]+bglyph[[i]]$xp[3]*zee[3]   )
        
   ####   print(c( i, bvec[i], bglyph[[i]]$xp, zdot))
        if(zdot>=0)
          {
            phongdot = abs((bglyph[[i]]$xp[1]*lightvec[1]+bglyph[[i]]$xp[2]*lightvec[2]+bglyph[[i]]$xp[3]*lightvec[3]   ))
            
            
            mycol = col2rgb(col[i])
            myhsv = rgb2hsv(mycol)

            ## print(paste(sep=' ', i,phongdot))
            
            vcol = phongdot*myhsv[3,1]
            if(vcol > 0.95) vcol=0.95
            if(vcol < 0.1) vcol=0.1
            
            thecol = hsv(h=myhsv[1,1], s=myhsv[2,1] , v =vcol )
            ######## hsv(h = 1, s = 1, v = 1, gamma = 1, alpha)
            polygon(bglyph[[i]]$x,bglyph[[i]]$y, col=thecol, border=border, xpd=TRUE)

          }
        ###  polygon(bglyph[[i]]$x,bglyph[[i]]$y, col=NA, border=grey(0.8))
      }
  }

