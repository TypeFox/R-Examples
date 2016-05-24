quadraticextrapolation <-
function(A1,A2,A3,lambda){
     #the fitted value of quadratic extrapolation at lambda=-1#
     reg1<-numeric()
   
     #the fitted value of quadratic extrapolation on diff at lambda=-1#
     reg2<-numeric() 
    
     #the fitted value of quadratic extrapolation on scale at lambda=-1#
     scalereg<-numeric()

     D=ncol(A1)

     #regression on estimates,tau#
     for(i in 1:D)
     {                      
             lambda2=lambda^2                                
             e1=coef(lm(A1[,i]~lambda + lambda2))
             a1 = e1[1] - e1[2] + e1[3]
             reg1 = c(reg1,a1)

             e2 = coef(lm(A2[,i]~lambda + lambda2))
             a2 = e2[1] - e2[2] + e2[3]
             reg2 = c(reg2, a2)                                     
     }

    
     #regression on scale#
     e3=coef(lm(A3[,1]~lambda+lambda2))
     a3= e3[1]-e3[2]+e3[3]
     scalereg= c(scalereg,a3)

     return(list("reg1"=reg1,"reg2"=reg2,"scalereg"=scalereg))
            
}
