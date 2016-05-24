rot2Zplane<-function (vec, p) 
      {
        #####  given a vector (normal)  and a point
        ##### return the matrix that
        #####  rotates a 3D body to the
        ####   the x-y plan (z=0)
        
        v = vec/sqrt(sum(vec^2))
        r1 = GEOmap::roty4(v)
        r2 = GEOmap::rotx4(v)
        ##  r3 = rotdelta4(alpha)
        t1 =  GEOmap::trans4(c(-p[1], -p[2], -p[3]))
        M <- t1 %*% r2 %*% r1

        return(M)
        

      }
