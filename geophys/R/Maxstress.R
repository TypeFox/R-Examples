Maxstress<-function(NN, Stensor)
  {
    ###   NN is a normal vector to some plane
    
###   Stensor is an arbitrary stress tensor, must be symmetric

        if(is.vector(Stensor))
          {
            if(length(Stensor)==3)
              {
                Stensor = diag(Stensor)  
              }
            if(length(Stensor)==6)
              {
                nstress = Stensor[1:6]
                MAT = matrix(ncol=3, nrow=3)
                
                MAT[ lower.tri(MAT, diag=TRUE)] = nstress
                MAT[ upper.tri(MAT) ] = MAT[lower.tri(MAT)]
                
                Stensor = MAT
                
                
              }
            
            
          }

        ###  Normalize the normal vector,  just in case
        NN = NN/sqrt(sum(NN*NN))

 ES = eigen(Stensor)
###  E = diag(ES$values)

###   sigNORMmax = NN %*% E %*% NN
 
    sigNORMmax = NN[1]^2*ES$values[1] + NN[2]^2 * ES$values[2]  +NN[3]^2 *  ES$values[3]
    
    tauSHEARmax  = NN[1]^2*NN[2]^2*(ES$values[1]-ES$values[2])^2 +
      NN[2]^2 *NN[3]^2 * (ES$values[2]-ES$values[3])^2  +
        NN[3]^2 *NN[1]^2 *  (ES$values[1] -ES$values[3])^2

    tauSHEARmax = sqrt(tauSHEARmax)
    

    return(list(NN=NN, sigNORMmax=sigNORMmax, tauSHEARmax=tauSHEARmax)  )




  }
