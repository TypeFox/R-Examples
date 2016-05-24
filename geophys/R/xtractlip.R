
xtractlip<-function(AF)
  {

    ## require(cluster)



    jlipse<-function (x, y, r1, r2, phi, by = 5)
      {
        if (missing(by)) {
          by = 5
        }
      
        theta = seq(0, 360, by = by) * pi/180
        cosp = cos(phi * pi/180)
        sinp = sin(phi * pi/180)
        r = matrix(c(cosp, -sinp, sinp, cosp), ncol = 2)
        m = matrix(rep(0, 2 * length(theta)), ncol = 2)
        m[, 1] = r1 * cos(theta)
        m[, 2] = r2 * sin(theta)
        nm = m %*% r

        return(list(x=x + nm[, 1], y=y + nm[, 2]))
       
      }



    delphi = 6*pi/180
    
    angs = seq(from=-pi, to=pi, length=360/3)
    
    tang = atan2(AF$y-AF$my , AF$x-AF$mx)
    dis = sqrt(  (AF$x-AF$mx)^2 +  (AF$y-AF$my)^2 )

    wmins = vector()
    gx = vector()
    gy = vector()

    d = mean(dis)
    
######plot(AF$x-AF$mx, AF$y-AF$my, asp=1, pch=".", cex=2, col=grey(.8) )
######points(0,0)
    

    
    for(i in 1:(length(angs))  )
      {
        ######  use overlap
        alpha1 = angs[i]
        ## alpha2 = angs[i+2]
         alpha2 = angs[i]+delphi


        
      ######  segments(0, 0, d*cos(alpha1), d*sin(alpha1), col=i)
      ######  segments(0, 0, d*cos(alpha2), d*sin(alpha2), col=i)
        
      
        w = which( tang>alpha1 & tang<alpha2)

     ######   points(AF$x[w]-AF$mx, AF$y[w]-AF$my, col='green', pch=".", cex=2)
     

        
        kw = which.min(dis[w])

        
        wmins = c(wmins,w[kw] )
    
        ###### points(AF$x[w[kw]]-AF$mx, AF$y[w[kw]]-AF$my, col='red', pch=3, cex=.8)
     
        
      }


    wpicks = unique(wmins)
    
    gx = AF$x[wpicks]-AF$mx
    gy = AF$y[wpicks]-AF$my
    

    if(FALSE)
      {
        par(mfrow=c(1,1))

    plot(gx, gy, type='n', asp=1)
    
     points(gx, gy, pch=3, col='purple')
    points(0,0, pch=4, col='black')


  }
    ##text(gx, gy,labels=1:length(gx), pos=3)
    
dix = sqrt( (gx)^2 +  (gy)^2)

  w1 = which.min(dix)
  w2 = which.max(dix)
   ndix = length(dix)
      minax = dix[w1]
    maxax =  dix[w2]

    u1 = c(gx[w2], gy[w2])
    tol=0
    savew = vector()
    k = 0
    V1 = c(0, gx[w2], 0, gy[w2] )

    for(j in 1:ndix)
      {
        if(j==w2) next
        V2 = c(0, gx[j], 0, gy[j]  )
        PP =   perpproj( V1, V2  )
        
        kdis = sqrt( (V2[2] -PP$P1[1])^2 + (V2[4]- PP$P1[2])^2 )
        
        if(abs(kdis)< minax+tol)
          {
            k = k+1
            savew[k] = j
          }


      }

    fx = gx[savew]
    fy = gy[savew]
    
   ###  points(fx , fy, pch=2, col='brown')
  
 ######   points(c(gx[w2],gx[w1] ) , c(gy[w2], gy[w1]), col='blue', pch=6)
  ndix = length(dix)
    Sdix = order(dix)

    msmal = median(dix[Sdix[1:5]])
    mlarge = median(dix[Sdix[(ndix-5):ndix ]])
       rang = atan2(  (gy[w2]) ,     (gx[w2]))  

    JL = jlipse(AF$mx, AF$my, mlarge, msmal,  rang*180/pi, by = 5 )

  
if(FALSE)
  {
  ndix = length(dix)
    Sdix = order(dix)

    msmal = median(dix[Sdix[1:5]])
    mlarge = median(dix[Sdix[(ndix-5):ndix ]])
    
    
    minax = dix[w1]
    maxax =  dix[w2]

CIRC1 =  GEOmap::darc(rad =minax , ang1 = 0, ang2 = 360, x1 = 0, y1 = 0, n = 1)
lines(CIRC1)
    
    rang = atan2(  (gy[w2]) ,     (gx[w2]))  


        par(mfrow=c(1,1))

    plot(gx, gy, type='n', asp=1)
    
     points(gx, gy, pch=3, col='purple')
    points(0,0, pch=4, col='black')


    v = c(gx[w2], gy[w2])
 points(gx[w2], gy[w2])
        
        for(i in 1:length(gx))
          {
            points(gx[i], gy[i])
            VJ = vecproj(v, c(gx[i], gy[i]))
            
          }


        

JL = jlipse(AF$mx, AF$my, mlarge, msmal,  rang*180/pi, by = 5 )

     EX = seq(from=min(gx) ,to=max(gx), length=100  )
     WHY = seq(from=min(gy) ,to=max(gy), length=100  )
          mm = meshgrid(EX, WHY)
      desh(mm, add=TRUE, PTS=FALSE, colmesh=grey(.8) )
   
  }
######  JL = jlipse(0, 0, maxax, minax,  rang*180/pi, by = 5 )

######  lines(JL)


    exy <- cluster::ellipsoidhull(unname(cbind(fx+AF$mx, fy+AF$my)))
    ## lines(predict(exy))

    invisible(list(lip=list(x=fx+AF$mx, y=fy+AF$my), Elips=JL,    hull=exy))

  }
