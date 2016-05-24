getEulers<-function(R)
  {
###  psi about X-axis
###  theta about Y axis
###  phi about Z-axis
    if(R[3,1]!=0)
      {
        theta1 = (-asin(R[3,1]))
        theta2 = pi - theta1
        psi1 = atan2(R[3,2]/cos(theta1) , R[3,3]/cos(theta1))
        psi2 = atan2(R[3,2]/cos(theta2) , R[3,3]/cos(theta2))
        phi1 = atan2(R[2,1]/cos(theta1) , R[1,1]/cos(theta1))
        phi2 = atan2(R[2,1]/cos(theta2) , R[1,1]/cos(theta2))
        SOL1 = c(phi1, theta1,  psi1  )
        SOL2 = c(phi2, theta2,  psi2  )


      }
    else
      {
        phi = 0
        delta = atan2(R[1,2], R[1,3] )
        if(R[3,1] == (-1))
          {
            theta = pi/2
            psi = phi + delta
          }
        else
          {
            theta = -pi/2
            psi = -phi + delta
          }
        
        SOL1 = c(phi, theta,  psi  )
        SOL2 = c(phi, theta,  psi )
        
      }
    return(SOL1) 
  }
