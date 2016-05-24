linearextrapolation <-
function(A1,A2,A3,lambda){
  
     #the fitted value of linear extrapolation at lambda=-1#
     reg1<-numeric()
   
     #the fitted value of linear extrapolation on diff at lambda=-1#
     reg2<-numeric() 
    
     #the fitted value of linear extrapolation on scale at lambda=-1#
     scalereg<-numeric()
    
     D=ncol(A1)                         
                            
     #regression on estimates,tau#
     for(i in 1:D)
     {    
            e1=coef(lm(A1[,i]~lambda))
            a1= e1[1] - e1[2]
            reg1 = c(reg1,a1)
 
            e2 = coef(lm(A2[,i]~lambda))        
            a2 = e2[1] - e2[2]        
            reg2 = c(reg2, a2)
        
     }

     #regression on scale#
     e3 = coef(lm(A3[,1]~lambda))
     a3 =  e3[1] - e3[2]
     scalereg= c(scalereg,a3)

     return(list("reg1"=reg1,"reg2"=reg2,"scalereg"=scalereg))
}
