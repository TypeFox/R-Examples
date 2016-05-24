flipGSVD <-
function(vs, d1=c(50, 50) , d2=c(48, 50)  )
{
##########  flip output of gsvd to match matlab!
 ##### d1 = dim(A)
 ##### d2 = dim(B)
  

####  diag(vs$C)

h1 = diag(vs$C[1:d1[2], 1:d1[2] ] )
  
oh1 = order(h1)
i1 = h1[oh1]

  temC = diag(i1)
  newC = temC
  if(d1[1]>d1[2])
    {
      
      newC = rbind( temC, vs$C[ (d1[2]+d1[1]-d1[2]):d1[1] ,  ] )
        
      }
  
##########
  newU =vs$U

  J1 = 1:d1[2]
  
  newU = vs$U
  newU[  , J1  ] = vs$U[ , J1[oh1] ]

  newX = vs$X[  , J1[oh1]] 
  
###  HH =  newU %*% newC %*% t(newX)

n = d2[2]
m = d2[1]

kel = (d2[2]-d2[1])

if(n<=m)
  {
    h2 = diag(vs$S[1:n, 1:n ] )
    
    oh2 = order(h2, decreasing = TRUE )
    i2 = h2[oh2]
    newS =diag(i2)

    if(n<m)
      {

        
        newS = rbind( newS, vs$S[ (d2[2]+d2[1]-d2[2]):d2[1] ,  ]   )

      }
    newV =vs$V

    J2 = 1:d1[2]
    
    newV[  , J2  ] = vs$V[ , J2[oh2] ]

  }
else
  {

    n2 = n-m+1
    zcols = n2:n
    
    h2 = diag(vs$S[ 1:m , zcols ] )
    
    oh2 = order(h2, decreasing = TRUE )
    i2 = h2[oh2]
    newS =diag(i2)

     newS = cbind(newS , vs$S[ 1:m , 1:(n2-1)  ])

    newV =vs$V[ ,oh2  ]

    
  }
  
  return(list(U=newU,V=newV, X=newX, C=newC,S=newS) )


}
