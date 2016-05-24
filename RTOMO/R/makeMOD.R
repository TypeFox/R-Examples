`makeMOD` <-
function(xo, yo, ztop, x, y, z, r, v, bg )
  {
    #######   make a (synthetic three dimensional model
    #######  for tomographic inversion

    #######  model is desribed by balls of specfic radius

##########  xo and yo are the centers of the blocks

 ####### 
    
    MOD = list()

    dx = xo[2]-xo[1]
    dy = yo[2]-yo[1]
    
    M = meshgrid(xo, yo)

    N = length(ztop)

    
    
    for(j in 1:(N))
      {
        ZOD  = matrix(rep(bg[j], length(xo)*length(yo))  , ncol=length(xo), nrow=length(yo), byrow=TRUE)
        
        
        if(j<N)
          {
            zi = mean(c(ztop[j],  ztop[j+1] ) )
          }
        else
          {
            zi = ztop[j] + (ztop[j]-ztop[j-1])/2

          }

        if(length(x)>0)
          {
            
            for(i in 1:length(x))
              {

                dis = sqrt( (M$x-x[i])^2+ (M$y-y[i])^2+(zi-z[i])^2 )
                ZOD[dis<r[i]] = v[i]
              }

            image(xo, yo, t(ZOD) )


          }
        ##  locator()


        
        MOD[[j]]  = t(ZOD)
      }


    

    invisible( list(x=xo,y=yo, D=ztop, MOD=MOD) )
  }

