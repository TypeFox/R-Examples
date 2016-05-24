copgHsCont <- function(p1, p2, teta, teta.st, VC, Cont = FALSE){
     

der2c.derrho.derrho <- der2c.derp1.derp1 <- der2c.derp2.derp2 <- der2c.derp1.derp2 <- der2c.derp1.derrho <- der2c.derp2.derrho <- 1     
     
########################################################################################   
# Transformations on teta parameter
########################################################################################   

   
if(VC$BivD %in% c("N","AMH","FGM") ) {

derteta.derteta.st <- 1/cosh(teta.st)^2
der2teta.derteta.stteta.st <- -(2 * (sinh(teta.st) * cosh(teta.st))/(cosh(teta.st)^2)^2)
       
}  

if(VC$BivD %in% c("F") ) {

derteta.derteta.st <- 1
der2teta.derteta.stteta.st <- 0
       
} 
   

if(VC$BivD %in% c("C0", "C180","J0", "J180","G0", "G180") ) derteta.derteta.st <- der2teta.derteta.stteta.st <-  exp(teta.st) 
if(VC$BivD %in% c("C90","C270","J90","J270","G90","G270") ) derteta.derteta.st <- der2teta.derteta.stteta.st <- -exp(teta.st)  

  

########################################################################################
########################################################################################


########################################################################################
# Rotations
########################################################################################

if(VC$BivD %in% c("C90","J90","G90") ) {
p1 <- 1 - p1 
teta <- -teta
}  

if(VC$BivD %in% c("C180","J180","G180") ) {
p1 <- 1 - p1
p2 <- 1 - p2
}  

if(VC$BivD %in% c("C270","J270","G270") ) {
p2 <- 1 - p2 
teta <- -teta 
}   
   
########################################################################################   
########################################################################################
     
     
     
     
     
     
if(VC$BivD == "N"){


der2h.derp2teta <- -((1 + teta * ((qnorm(p1) - teta * qnorm(p2)) * (qnorm(p2) - 
    teta * (qnorm(p1) - teta * qnorm(p2))/(1 - teta^2)) + teta)/(1 - 
    teta^2)) * dnorm((qnorm(p1) - teta * qnorm(p2))/sqrt(1 - 
    teta^2)) * sqrt(2 * pi)/(exp(-(qnorm(p2)^2/2)) * sqrt(1 - 
    teta^2)))
    
    
der2h.derp2p2 <-  -(teta * dnorm((qnorm(p1) - teta * qnorm(p2))/sqrt(1 - 
    teta^2)) * (qnorm(p2) + teta * (qnorm(p1) - teta * qnorm(p2))/(1 - 
    teta^2))/(dnorm(qnorm(p2))^2 * sqrt(1 - teta^2)))

der2h.derteta.teta.st <- ((qnorm(p1) - teta * qnorm(p2)) * (qnorm(p2) - teta * 
    (qnorm(p1) - teta * qnorm(p2))/(1 - teta^2)) * (teta * (qnorm(p1) - 
    teta * qnorm(p2))/sqrt(1 - teta^2) - qnorm(p2) * sqrt(1 - 
    teta^2)) + (qnorm(p1) + teta * (2 * (teta * qnorm(p1)) - 
    3 * qnorm(p2)))/sqrt(1 - teta^2)) * dnorm((qnorm(p1) - teta * 
    qnorm(p2))/sqrt(1 - teta^2))/(1 - teta^2)^2
  
  
der2h.derp1p2 <-   -((exp(-(teta^2*(qnorm(p2)^2+qnorm(p1)^2)-(2*teta*qnorm(p2)*qnorm(p1)))/
                     (1 - teta^2)/2))/sqrt(1 - teta^2))*(teta*(sqrt(2*pi)/exp(-qnorm(p2)^2/2))/
                    (1 - teta^2))*(teta*qnorm(p2)-qnorm(p1)) 

der2h.derp1teta <- exp(-(teta * (teta * (qnorm(p1)^2 + qnorm(p2)^2) - 
    2 * (qnorm(p1) * qnorm(p2)))/(2 * (1 - teta^2)))) * (qnorm(p1) * 
    qnorm(p2) + teta * (1 - (qnorm(p1)^2 + qnorm(p2)^2 + teta * 
    (teta * (qnorm(p1)^2 + qnorm(p2)^2) - 2 * (qnorm(p1) * qnorm(p2)))/(1 - 
    teta^2))))/((1 - teta^2) * sqrt(1 - teta^2))
                                                            

 
der2h.derp1p1 <-   -((exp(-(teta^2*(qnorm(p1)^2+qnorm(p2)^2)-(2*teta*qnorm(p1)*qnorm(p2)))/
	             (1 - teta^2)/2))/sqrt(1 - teta^2))*(teta*(sqrt(2*pi)/exp(-qnorm(p1)^2/2))/
	            (1 - teta^2))*(teta*qnorm(p1)-qnorm(p2))    




if(Cont == TRUE){

c <- 1/sqrt(1 - teta^2)*exp(  - (teta^2*( qnorm(p1)^2 +  qnorm(p2)^2 ) - 2*teta*qnorm(p1)*qnorm(p2) ) / (2*(1 - teta^2)) ) 

			t6 = qnorm(p1) 
			t7 = qnorm(p2)
			t1 = t6*t7
			t2 = teta*teta
			t3 = 1-t2
			t4 = 4*t3*t3
			t5 = 1/t4
			t12 = t6*t6
			t13 = t7*t7
			t14 = 2*teta*t6*t7-t12-t13
			t21 = t14*t5
			t26 = 1/t3/2
			t29 = exp(t12/2+t13/2+t14*t26)
			t31 = sqrt(t3)
			t32 = 1/t31
			t38 = 2*t1*t26+4*t21*teta
			t39 = t38*t38
			t44 = 1/t31/t3
			t48 = t3*t3
der2c.derrho.derrho <- (16.0*t1*t5*teta+16.0*t14/t4/t3*t2+4.0*t21)*t29*t32+t39*t29*t32+2.0*t38*t29*t44*teta+3.0*t29/t31/t48*t2+t29*t44

			diffc = der2h.derp1p1 
			t1=qnorm(p1)
			t2=qnorm(p2)
			t3=dnorm(t1)
			# //t4=dnorm(t2,0.0,1.0,0);
			t5=1.0-teta*teta
			t6=teta*t1-t2
			t7=teta+t6*t1
			t8=t7/t3/t3
der2c.derp1.derp1 <- -teta/t5 * (diffc*t6/t3 + c*t8) 
			
			diffc = der2h.derp1p2
			t1=qnorm(p2)
			t2=qnorm(p1)
			t3=dnorm(t1)
			t5=1.0-teta*teta
			t6=teta*t1-t2
			t7=teta+t6*t1
			t8=t7/t3/t3
der2c.derp2.derp2 <- -teta/t5 * (diffc*t6/t3 + c*t8) 

			diffc = der2h.derp1p2

			t1=qnorm(p1)
			t2=qnorm(p2)
			t3=dnorm(t1)
			t4=dnorm(t2)
			t5=1.0-teta^2;
			t7=teta*t1-t2;
der2c.derp1.derp2 <- -teta/t3/t5 * (diffc*t7 - c/t4);


			diffc <- der2h.derp1teta # diffPDF(&p1,&p2,&k,param,copula,&diffc);

			t1=qnorm(p1)
			t2=qnorm(p2)
			t3=1/dnorm(t1)
			t4=teta*(teta*t1-t2);
			t5=1.0-teta^2
			t6=-t4/t5;
			t7=-2.0*teta*t1+t2+t2*teta*teta;
			t8=(-1+teta^2)^2
			t9=t7/t8;
der2c.derp1.derrho <- diffc*t6*t3 + c*t3*t9;

			diffc <- der2h.derp1teta 

			t1=qnorm(p2)
			t2=qnorm(p1)
			t3=1/dnorm(t1)
			t4=teta*(teta*t1-t2);
			t5=1.0-teta^2
			t6=-t4/t5;
			t7=-2.0*teta*t1+t2+t2*teta*teta;
			t8=(-1+teta^2)^2
			t9=t7/t8;
der2c.derp2.derrho <- diffc*t6*t3 + c*t3*t9;

}


}













if(VC$BivD == "FGM"){

der2h.derp2teta <-  -(2 * (p1 * (1 - p1)))
   
der2h.derp2p2 <- 0

der2h.derteta.teta.st <- 0

der2h.derp1p2 <-  -(2 * (teta * (1 - 2 * p1)))

der2h.derp1teta <-  (1 - 2 * p1) * (1 - 2 * p2)

der2h.derp1p1 <-  -(2 * (teta * (1 - 2 * p2)))





if(Cont == TRUE){


der2c.derrho.derrho <- 0


der2c.derp1.derp1 <- 0


der2c.derp2.derp2 <- 0




der2c.derp1.derp2 <- 4 * teta




der2c.derp1.derrho <- -(2 * (1 - 2 * p2))



der2c.derp2.derrho <- -(2 * (1 - 2 * p1))


}





}

















if(VC$BivD == "AMH"){

der2h.derp2teta <-  -(p1 * (1 - p1) * (2 + teta * (1 - p1) * (2 * (1 - p2 * (2 + 
    2 * (teta * (1 - p1) * (1 - p2)/(1 - teta * (1 - p1) * (1 - 
        p2))))) + 2 * (1 - p2) - 2 * (p2 * (1 + teta * (1 - p1) * 
    (1 - p2)/(1 - teta * (1 - p1) * (1 - p2)))))/(1 - teta * 
    (1 - p1) * (1 - p2)))/(1 - teta * (1 - p1) * (1 - p2))^2)
   
der2h.derp2p2 <- p1 * teta^2 * (1 - p1)^2 * (2 * (1 - p2 * teta * (1 - p1)/(1 - 
    teta * (1 - p1) * (1 - p2))) + 2 * (2 - 2 * (p2 * teta * 
    (1 - p1)/(1 - teta * (1 - p1) * (1 - p2)))))/(1 - teta * 
    (1 - p1) * (1 - p2))^3


der2h.derteta.teta.st <- -(p1 * (1 - p1)^2 * (1 - p2) * (2 * (p2 * (1 + teta * (1 - p1) * 
    (1 - p2)/(1 - teta * (1 - p1) * (1 - p2)))) - 2 * (1 - p2 * 
    (2 + 2 * (teta * (1 - p1) * (1 - p2)/(1 - teta * (1 - p1) * 
        (1 - p2))))))/(1 - teta * (1 - p1) * (1 - p2))^3)


der2h.derp1p2 <- -(teta * ((1 - 2 * (teta * (p1 * (1 - p2 * (2 + 2 * (teta * (1 - 
    p1) * (1 - p2)/(1 - teta * (1 - p1) * (1 - p2))))) + p2 * 
    (1 - p1))/(1 - teta * (1 - p1) * (1 - p2)))) * (1 - p1) + 
    1 - p1 * (3 + teta * (1 - p1) * (2 * (1 - p2) - 2 * (p2 * 
    (1 + teta * (1 - p1) * (1 - p2)/(1 - teta * (1 - p1) * (1 - 
        p2)))))/(1 - teta * (1 - p1) * (1 - p2))))/(1 - teta * 
    (1 - p1) * (1 - p2))^2)
 

der2h.derp1teta <-  -((((1 - p2) * (2 * (teta * (p1 * (1 - p2 * (2 + 2 * (teta * 
    (1 - p1) * (1 - p2)/(1 - teta * (1 - p1) * (1 - p2))))) + 
    p2 * (1 - p1))/(1 - teta * (1 - p1) * (1 - p2))) - 1) + p2 * 
    (1 - 2 * (p1 * teta * (1 - p2) * (1 + teta * (1 - p1) * (1 - 
        p2)/(1 - teta * (1 - p1) * (1 - p2)))/(1 - teta * (1 - 
        p1) * (1 - p2))))) * (1 - p1) + p1 * (1 - p2 * (2 + 2 * 
    (teta * (1 - p1) * (1 - p2)/(1 - teta * (1 - p1) * (1 - p2))))))/(1 - 
    teta * (1 - p1) * (1 - p2))^2)


der2h.derp1p1 <- -(teta * ((1 - 2 * (teta * (p1 * (1 - p2 * (2 + 2 * (teta * (1 - 
    p1) * (1 - p2)/(1 - teta * (1 - p1) * (1 - p2))))) + p2 * 
    (1 - p1))/(1 - teta * (1 - p1) * (1 - p2)))) * (1 - p2) + 
    1 + p2 * (teta * (1 - p2) * (2 * (p1 * (1 + teta * (1 - p1) * 
    (1 - p2)/(1 - teta * (1 - p1) * (1 - p2)))) - 2 * (1 - p1))/(1 - 
    teta * (1 - p1) * (1 - p2)) - 3))/(1 - teta * (1 - p1) * 
    (1 - p2))^2) 





if(Cont == TRUE){


der2c.derrho.derrho <- -((((1 - p2) * (teta * (p1 * (1 - p2 * (2 + 2 * (teta * (1 - 
    p1) * (1 - p2)/(1 - teta * (1 - p1) * (1 - p2))))) + p2 * 
    (1 - p1))/(1 - teta * (1 - p1) * (1 - p2)) - 2 * (1 - teta * 
    (p1 * (1 - p2 * (2 + 2 * (teta * (1 - p1) * (1 - p2)/(1 - 
        teta * (1 - p1) * (1 - p2))))) + p2 * (1 - p1))/(1 - 
    teta * (1 - p1) * (1 - p2)))) + 3 * (teta * (1 - p2) * (p1 * 
    (1 - p2 * (2 + 2 * (teta * (1 - p1) * (1 - p2)/(1 - teta * 
        (1 - p1) * (1 - p2))))) + p2 * (1 - p1))/(1 - teta * 
    (1 - p1) * (1 - p2))) + 4 * p2 * (1 - 2 * (p1 * teta * (1 - 
    p2) * (1 + teta * (1 - p1) * (1 - p2)/(1 - teta * (1 - p1) * 
    (1 - p2)))/(1 - teta * (1 - p1) * (1 - p2))))) * (1 - p1) + 
    p1 * (1 - p2 * (2 + 2 * (teta * (1 - p1) * (1 - p2)/(1 - 
        teta * (1 - p1) * (1 - p2))))) + p1 * (3 - p2 * (4 * 
    (1 + teta * (1 - p1) * (1 - p2)/(1 - teta * (1 - p1) * (1 - 
        p2))) + 6 + teta * (1 - p1) * (1 - p2) * (2 * (2 + 2 * 
    (teta * (1 - p1) * (1 - p2)/(1 - teta * (1 - p1) * (1 - p2)))) + 
    6)/(1 - teta * (1 - p1) * (1 - p2))))) * (1 - p1) * (1 - 
    p2)/(1 - teta * (1 - p1) * (1 - p2))^3)



der2c.derp1.derp1 <- -(teta^2 * ((1 - p2) * (4 * (teta * (p1 * (1 - p2 * (2 + 2 * 
    (teta * (1 - p1) * (1 - p2)/(1 - teta * (1 - p1) * (1 - p2))))) + 
    p2 * (1 - p1))/(1 - teta * (1 - p1) * (1 - p2))) - 2 * (1 - 
    teta * (p1 * (1 - p2 * (2 + 2 * (teta * (1 - p1) * (1 - p2)/(1 - 
        teta * (1 - p1) * (1 - p2))))) + p2 * (1 - p1))/(1 - 
        teta * (1 - p1) * (1 - p2)))) + p2 * (12 + 4 * (1 + teta * 
    (1 - p1) * (1 - p2)/(1 - teta * (1 - p1) * (1 - p2))) - teta * 
    (1 - p2) * (p1 * (2 * (2 + 2 * (teta * (1 - p1) * (1 - p2)/(1 - 
    teta * (1 - p1) * (1 - p2)))) + 8 * (1 + teta * (1 - p1) * 
    (1 - p2)/(1 - teta * (1 - p1) * (1 - p2)))) - 8 * (1 - p1))/(1 - 
    teta * (1 - p1) * (1 - p2))) - 4) * (1 - p2)/(1 - teta * 
    (1 - p1) * (1 - p2))^3)



der2c.derp2.derp2 <- -(teta^2 * ((1 - p1) * (4 * (teta * (p1 * (1 - p2 * (2 + 2 * 
    (teta * (1 - p1) * (1 - p2)/(1 - teta * (1 - p1) * (1 - p2))))) + 
    p2 * (1 - p1))/(1 - teta * (1 - p1) * (1 - p2))) - 2 * (1 - 
    teta * (p1 * (1 - p2 * (2 + 2 * (teta * (1 - p1) * (1 - p2)/(1 - 
        teta * (1 - p1) * (1 - p2))))) + p2 * (1 - p1))/(1 - 
        teta * (1 - p1) * (1 - p2)))) + p1 * (12 + 4 * (1 + teta * 
    (1 - p1) * (1 - p2)/(1 - teta * (1 - p1) * (1 - p2))) + teta * 
    (1 - p1) * (8 * (1 - p2) - p2 * (2 * (2 + 2 * (teta * (1 - 
    p1) * (1 - p2)/(1 - teta * (1 - p1) * (1 - p2)))) + 8 * (1 + 
    teta * (1 - p1) * (1 - p2)/(1 - teta * (1 - p1) * (1 - p2)))))/(1 - 
    teta * (1 - p1) * (1 - p2))) - 4) * (1 - p1)/(1 - teta * 
    (1 - p1) * (1 - p2))^3)





der2c.derp1.derp2 <- -(teta * (teta * (((1 - p2) * (3 * (teta * (p1 * (1 - p2 * (2 + 
    2 * (teta * (1 - p1) * (1 - p2)/(1 - teta * (1 - p1) * (1 - 
        p2))))) + p2 * (1 - p1))/(1 - teta * (1 - p1) * (1 - 
    p2))) - (2 + 2 * (1 - teta * (p1 * (1 - p2 * (2 + 2 * (teta * 
    (1 - p1) * (1 - p2)/(1 - teta * (1 - p1) * (1 - p2))))) + 
    p2 * (1 - p1))/(1 - teta * (1 - p1) * (1 - p2))))) + p2 * 
    (2 * (1 + teta * (1 - p1) * (1 - p2)/(1 - teta * (1 - p1) * 
        (1 - p2))) + 8 - 2 * (teta * (1 - p2) * (2 * (p1 * (1 + 
        teta * (1 - p1) * (1 - p2)/(1 - teta * (1 - p1) * (1 - 
            p2)))) - 2 * (1 - p1))/(1 - teta * (1 - p1) * (1 - 
        p2)))) - 2) * (1 - p1) + p1 * (2 + 2 * ((1 - p2) * (1 + 
    teta * (1 - p1) * (1 - p2)/(1 - teta * (1 - p1) * (1 - p2)))) - 
    p2 * (2 * (1 + teta * (1 - p1) * (1 - p2) * (2 * (teta * 
        (1 - p1) * (1 - p2)/(1 - teta * (1 - p1) * (1 - p2))) + 
        3)/(1 - teta * (1 - p1) * (1 - p2))) + 4 + 4 * (teta * 
        (1 - p1) * (1 - p2)/(1 - teta * (1 - p1) * (1 - p2))))) - 
    (1 - p2) * (2 - (2 * (p1 * (3 + teta * (1 - p1) * (2 * (1 - 
        p2) - 2 * (p2 * (1 + teta * (1 - p1) * (1 - p2)/(1 - 
        teta * (1 - p1) * (1 - p2)))))/(1 - teta * (1 - p1) * 
        (1 - p2)))) + teta * (1 - p1) * (p1 * (1 - p2 * (2 + 
        2 * (teta * (1 - p1) * (1 - p2)/(1 - teta * (1 - p1) * 
            (1 - p2))))) + p2 * (1 - p1))/(1 - teta * (1 - p1) * 
        (1 - p2)))))/(1 - teta * (1 - p1) * (1 - p2)) - 4)/(1 - 
    teta * (1 - p1) * (1 - p2))^2)





der2c.derp1.derrho <- -(((1 - p2) * (1 + teta * (((1 - p2) * (2 * (1 - teta * (p1 * 
    (1 - p2 * (2 + 2 * (teta * (1 - p1) * (1 - p2)/(1 - teta * 
        (1 - p1) * (1 - p2))))) + p2 * (1 - p1))/(1 - teta * 
    (1 - p1) * (1 - p2))) - 4 * (teta * (p1 * (1 - p2 * (2 + 
    2 * (teta * (1 - p1) * (1 - p2)/(1 - teta * (1 - p1) * (1 - 
        p2))))) + p2 * (1 - p1))/(1 - teta * (1 - p1) * (1 - 
    p2)))) + 2 + p2 * (teta * (1 - p2) * (8 * (p1 * (1 + teta * 
    (1 - p1) * (1 - p2)/(1 - teta * (1 - p1) * (1 - p2)))) - 
    4 * (1 - p1))/(1 - teta * (1 - p1) * (1 - p2)) - 10)) * (1 - 
    p1) - 4 * (p1 * (1 - p2 * (2 + 2 * (teta * (1 - p1) * (1 - 
    p2)/(1 - teta * (1 - p1) * (1 - p2)))))))/(1 - teta * (1 - 
    p1) * (1 - p2))) + 1 + p2 * (teta * (1 - p2) * (p1 * (2 * 
    (1 + teta * (1 - p1) * (1 - p2) * (2 * (teta * (1 - p1) * 
        (1 - p2)/(1 - teta * (1 - p1) * (1 - p2))) + 3)/(1 - 
        teta * (1 - p1) * (1 - p2))) + 2 * (1 + teta * (1 - p1) * 
    (1 - p2)/(1 - teta * (1 - p1) * (1 - p2)))) - (1 - p1) * 
    (2 + 2 * (1 + teta * (1 - p1) * (1 - p2)/(1 - teta * (1 - 
        p1) * (1 - p2)))))/(1 - teta * (1 - p1) * (1 - p2)) - 
    3))/(1 - teta * (1 - p1) * (1 - p2))^2)



der2c.derp2.derrho <- -(((1 - p1) * (1 + teta * (((1 - p2) * (2 * (1 - teta * (p1 * 
    (1 - p2 * (2 + 2 * (teta * (1 - p1) * (1 - p2)/(1 - teta * 
        (1 - p1) * (1 - p2))))) + p2 * (1 - p1))/(1 - teta * 
    (1 - p1) * (1 - p2))) - 3 * (teta * (p1 * (1 - p2 * (2 + 
    2 * (teta * (1 - p1) * (1 - p2)/(1 - teta * (1 - p1) * (1 - 
        p2))))) + p2 * (1 - p1))/(1 - teta * (1 - p1) * (1 - 
    p2)))) - 2 * (p2 * (2 - 2 * (p1 * teta * (1 - p2) * (1 + 
    teta * (1 - p1) * (1 - p2)/(1 - teta * (1 - p1) * (1 - p2)))/(1 - 
    teta * (1 - p1) * (1 - p2)))))) * (1 - p1) + (1 - p2) * (2 - 
    (2 * (p1 * (3 + teta * (1 - p1) * (2 * (1 - p2) - 2 * (p2 * 
        (1 + teta * (1 - p1) * (1 - p2)/(1 - teta * (1 - p1) * 
            (1 - p2)))))/(1 - teta * (1 - p1) * (1 - p2)))) + 
        teta * (1 - p1) * (p1 * (1 - p2 * (2 + 2 * (teta * (1 - 
            p1) * (1 - p2)/(1 - teta * (1 - p1) * (1 - p2))))) + 
            p2 * (1 - p1))/(1 - teta * (1 - p1) * (1 - p2)))) - 
    4 * (p1 * (1 - p2 * (2 + 2 * (teta * (1 - p1) * (1 - p2)/(1 - 
        teta * (1 - p1) * (1 - p2)))))))/(1 - teta * (1 - p1) * 
    (1 - p2))) + 1 - p1 * (3 + teta * ((1 - p2) * (2 + 2 * (1 + 
    teta * (1 - p1) * (1 - p2)/(1 - teta * (1 - p1) * (1 - p2)))) - 
    p2 * (2 * (1 + teta * (1 - p1) * (1 - p2) * (2 * (teta * 
        (1 - p1) * (1 - p2)/(1 - teta * (1 - p1) * (1 - p2))) + 
        3)/(1 - teta * (1 - p1) * (1 - p2))) + 2 * (1 + teta * 
        (1 - p1) * (1 - p2)/(1 - teta * (1 - p1) * (1 - p2))))) * 
    (1 - p1)/(1 - teta * (1 - p1) * (1 - p2))))/(1 - teta * (1 - 
    p1) * (1 - p2))^2)


}





}






        


if(VC$BivD %in% c("C0","C90","C180","C270")){

der2h.derp2teta <-  -(p1^teta * p2^(-2 + teta) * ((p1^teta - 1) * (teta * 
    (1 + teta) * log(p2) * (p1^teta * (1 + teta) + p2^teta * 
    teta * (p1^teta - 1)) - ((1 + teta) * log(1/p1^teta + 1/p2^teta - 
    1) + teta^2) * (p1^teta * (p2^teta - 1) - p2^teta)) + teta * 
    (1 + teta) * log(p1) * (p1^teta * teta + p2^teta * (1 + teta) * 
    (p1^teta - 1)))/(teta^2 * (1/p1^teta + 1/p2^teta - 1)^(1/teta) * 
    (p1^teta * (p2^teta - 1) - p2^teta)^3))
   
   
t1 = -teta-1
t2 = p2^t1
t3 = t1^2
t5 = p2^2
t6 = 1/t5
t7 = p1^-teta
t8 = p2^-teta
t9 = t7+t8-1
t11 = -1-1/teta
t12 = t9^t11
t13 = t6*t12
t16 = t2*t1*t13
t18 = 1/t9
t23 = t2*t12
t24 = t11*t11
t26 = t8*t8
t27 = teta^2
t29 = t9*t9
t32 = t26*t27*t6/t29
t34 = t23*t11
t36 = t6*t18

    der2h.derp2p2 <- t2*t3*t13-t16-2*t16*t11*t8*teta*t18+t23*t24*t32+t34*t8*t27*t36+t34*t8*teta*t36-t34*t32


t2 = p2^(-1.0*teta-1.0);
t3 = log(p2);
t4 = t3*t3;
t6 = p1^(-1.0*teta);
t7 = p2^(-1.0*teta);
t8 = t6+t7-1.0;
t10 = -1.0-1/teta;
t11 = t8^(1.0*t10);
t14 = teta^2;
t15 = 1/t14;
t16 = log(t8);
t18 = log(p1);
t21 = -t6*t18-t7*t3;
t23 = 1/t8;
t25 = t15*t16+t10*t21*t23;
t29 = t2*t11;
t32 = t25*t25;
t12 = t18*t18;
t9 = t21*t21;
t5 = t8*t8;

der2h.derteta.teta.st <- t2*t4*t11-2.0*t2*t3*t11*t25+t29*t32+t29*(-2.0/t14/teta*t16+2.0*t15*t21*t23+t10*(t6*t12+t7*t4)*t23-t10*t9/t5);




t1 = 1.0+teta;
t3 = (p2*p1)^(-1.0*t1);
t4 = t1*t3;
t5 = 1/p2;
t7 = p2^(-1.0*teta);
t8 = p1^(-1.0*teta);
t9 = t7+t8-1.0;
t11 = -2.0-1/teta;
t12 = t9^(1.0*t11);

der2h.derp1p2 <-    -t4*t1*t5*t12-t4*t12*t11*t7*teta*t5/t9;


t1 = p1*p2;
t2 = -teta-1.0;
t3 = t1^(1.0*t2);
t4 = p1^(-1.0*teta);
t5 = p2^(-1.0*teta);
t6 = t4+t5-1.0;
t7 = -2.0-1/teta;
t8 = t6^(1.0*t7);
t9 = -t2*t3;
t10 = log(t1);
t11 = teta^2;
t12 = log(t6);
t13 = log(p1);
t14 = log(p2);

der2h.derp1teta <-  t3*t8-t9*t10*t8+t9*t8*(1/t11*t12+t7*(-t4*t13-t5*t14)/t6);






t1 = 1.0+teta;
t3 = (p2*p1)^(-1.0*t1);
t4 = t1*t3;
t5 = 1/p1;
t7 = p1^(-1.0*teta);
t8 = p2^(-1.0*teta);
t9 = t7+t8-1.0;
t11 = -2.0-1/teta;
t12 = t9^(1.0*t11);

der2h.derp1p1 <-    -t4*t1*t5*t12-t4*t12*t11*t7*teta*t5/t9;





if(Cont == TRUE){



			t1 = p1*p2
			t2 = -teta-1
			t3 = t1^t2
			t4 = log(t1)
			t6 = p1^-teta
			t7 = p2^-teta
			t8 = t6+t7-1
			t10 = -2-1/teta
			t11 = t8^t10
			t15 = teta^2
			t16 = 1/t15
			t17 = log(t8)
			t19 = log(p1)
			t21 = log(p2)
			t24 = -t6*t19-t7*t21
			t26 = 1/t8
			t27 = t16*t17+t10*t24*t26
			t30 = -t2*t3
			t32 = t4*t4
			t14 = t27*t27
			t13 = t19*t19
			t12 = t21*t21
			t9 = t24*t24
			t5 = t8*t8
der2c.derrho.derrho <- -2*t3*t4*t11+2*t3*t11*t27+t30*t32*t11-2*t30*t4*t11*t27+t30*t11*t14+t30*t11*(-2/t15/teta*t17+2*t16*t24*t26+t10*(t6*t13+t7*t12)*t26-t10*t9/t5)


			t1 = 1+teta
			t3 = (p1*p2)^-t1
			t4 = t1*t3
			t5 = t1*t1
			t6 = p1^2
			t7 = 1/t6
			t9 = p1^-teta
			t10 = p2^-teta
			t11 = t9+t10-1
			t13 = -2-1/teta
			t14 = t11^t13
			t17 = -t1*t7
			t21 = t14*t13
			t22 = t9*teta
			t23 = 1/t11
			t28 = t13*t13
			t31 = t9*t9
			t32 = teta^2
			t34 = t11*t11
			t37 = t31*t32*t7/t34
			t39 = t4*t21
			t41 = t7*t23
der2c.derp1.derp1 <- t4*t5*t7*t14-t4*t17*t14-2*t4*t17*t21*t22*t23+t4*t14*t28*t37+t39*t9*t32*t41+t39*t22*t41-t39*t37








			t1 = 1+teta
			t3 = (p1*p2)^-t1
			t4 = t1*t3
			t5 = t1*t1
			t6 = p2^2
			t7 = 1/t6
			t9 = p2^-teta
			t10 = p1^-teta
			t11 = t9+t10-1
			t13 = -2-1/teta
			t14 = t11^t13
			t17 = -t1*t7
			t21 = t14*t13
			t22 = t9*teta
			t23 = 1/t11
			t28 = t13*t13
			t31 = t9*t9
			t32 = teta^2
			t34 = t11*t11
			t37 = t31*t32*t7/t34
			t39 = t4*t21
			t41 = t7*t23
der2c.derp2.derp2 <- t4*t5*t7*t14-t4*t17*t14-2*t4*t17*t21*t22*t23+t4*t14*t28*t37+t39*t9*t32*t41+t39*t22*t41-t39*t37





			t1 = 1+teta
			t3 = (p1*p2)^-t1
			t4 = t1*t3
			t5 = t1*t1
			t7 = 1/p2
			t8 = 1/p1
			t10 = p1^-teta
			t11 = p2^-teta
			t12 = t10+t11-1
			t14 = -2-1/teta
			t15 = t12^t14
			t23 = 1/t12
			t35 = t14*t14
			t39 = teta^2
			t41 = t12*t12
			t42 = 1/t41
der2c.derp1.derp2 <- t4*t5*t7*t8*t15+t4*t1*t8*t15*t14*t11*teta*t7*t23+t4*t1*t7*t15*t14*t10*teta*t8*t23+t4*t15*t35*t11*t39*t7*t42*t10*t8-t4*t15*t14*t10*t39*t8*t42*t11*t7





			t1 = p1*p2
			  t2 = -teta-1
			  t3 = t1^t2
			  t5 = 1/p1
			  t6 = p1^-teta
			  t7 = p2^-teta
			  t8 = t6+t7-1
			  t9 = 1/teta
			  t10 = -2-t9
			  t11 = t8^t10
			  t12 = t5*t11
			  t16 = t6*teta
			  t17 = 1/t8
			  t18 = t5*t17
			  t21 = -t3*t2
			  t22 = t21*t2
			  t23 = log(t1)
			  t35 = teta^2
			  t37 = log(t8)
			  t39 = log(p1)
			  t41 = log(p2)
			  t44 = t10*(-t6*t39-t7*t41)
			  t46 = 1/t35*t37+t44*t17
			  t62 = t8*t8
der2c.derp1.derrho <- t3*t2*t12-t3*t11*t10*t16*t18-t22*t5*t23*t11-t21*t12+t21*t23*t11*t10*t6*teta*t5*t17+t22*t12*t46-t21*t11*t10*t16*t18*t46+t21*t11*(-t9*t6*t18+t10*(t16*t5*t39-t6*t5)*t17+t44/t62*t16*t5)





			t1 = p1*p2
			  t2 = -teta-1
			  t3 = t1^t2
			  t5 = 1/p2
			  t6 = p2^-teta
			  t7 = p1^-teta
			  t8 = t6+t7-1
			  t9 = 1/teta
			  t10 = -2-t9
			  t11 = t8^t10
			  t12 = t5*t11
			  t16 = t6*teta
			  t17 = 1/t8
			  t18 = t5*t17
			  t21 = -t3*t2
			  t22 = t21*t2
			  t23 = log(t1)
			  t35 = teta^2
			  t37 = log(t8)
			  t39 = log(p2)
			  t41 = log(p1)
			  t44 = t10*(-t6*t39-t7*t41)
			  t46 = 1/t35*t37+t44*t17
			  t62 = t8*t8
der2c.derp2.derrho <- t3*t2*t12-t3*t11*t10*t16*t18-t22*t5*t23*t11-t21*t12+t21*t23*t11*t10*t6*teta*t5*t17+t22*t12*t46-t21*t11*t10*t16*t18*t46+t21*t11*(-t9*t6*t18+t10*(t16*t5*t39-t6*t5)*t17+t44/t62*t16*t5)


}




















}





if(VC$BivD == "F"){

der2h.derp2teta <- -(((1 + p2 * teta) * exp(2 * teta) + (1 + teta * (1 - 
    p2)) * exp(teta * (3 * p1 + p2)) + (1 + teta * (p1 - p2)) * 
    exp(teta * (2 + p1 + p2)) + (1 + teta * (p2 - p1)) * exp(2 * 
    (teta * (1 + p1))) + (2 - teta * (2 * p2 + p1 - 1)) * exp(teta * 
    (1 + p1 + p2)) + (2 + teta * (2 * p2 + p1 - 2)) * exp(teta * 
    (1 + 2 * p1)) + exp(teta * (1 + 3 * p1)) * (teta * (1 - p2) - 
    1) + exp(teta * (2 * p1 + p2)) * (teta * (p1 + p2 - 1) - 
    1) + exp(teta * (2 + p1)) * (teta * (p1 - 2 * p2) - 2) + 
    exp(teta * (2 + p2)) * (p2 * teta - 1) - ((1 + teta * (p1 + 
    p2 - 1)) * exp(teta * (1 + p1)) + (2 + teta * (1 + p1 - 2 * 
    p2)) * exp(teta * (1 + 2 * p1 + p2)))) * exp(teta * (1 + 
    p2))/((exp(p1 * teta) + exp(p2 * teta) - 1) * exp(teta) - 
    exp(teta * (p1 + p2)))^3)
        
        
        
t1 = exp(teta)
t2 = teta*p1
t3 = exp(t2)
t5 = t1*(t3-1)
t6 = teta*p2
t8 = exp(t6+t2)
t10 = exp(t6+teta)
t12 = exp(t2+teta)
t13 = t8-t10-t12+t1
t14 = t13*t13
t20 = (teta*t8-teta*t10)^2
t24 = teta^2

       der2h.derp2p2 <- -2.0*t5/t14/t13*t20+t5/t14*(t24*t8-t24*t10)


t1 = exp(teta);
t2 = teta*p1;
t3 = exp(t2);
t5 = t1*(t3-1.0);
t6 = teta*p2;
t8 = exp(t6+t2);
t10 = exp(t6+teta);
t12 = exp(t2+teta);
t13 = t8-t10-t12+t1;
t14 = 1/t13;
t16 = t1*p1;
t18 = t3*t14;
t20 = t13*t13;
t21 = 1/t20;
t23 = p2+p1;
t25 = p2+1.0;
t26 = p1+1.0;
t28 = t23*t8-t25*t10-t26*t12+t1;
t32 = p1*p1;
t42 = t28*t28;
t44 = t23*t23;
t47 = t25*t25;
t49 = t26*t26;

der2h.derteta.teta.st <- -t5*t14-2.0*t16*t18+2.0*t5*t21*t28-t1*t32*t18+2.0*t16*t3*t21*t28-2.0*t5/t20/t13*t42+t5*t21*(t44*t8-t47*t10-t49*t12+t1);



			t1 = teta^2;
			t2 = exp(teta);
			t3 = t2 - 1;
			t5 = teta*p1;
			t6 = teta*p2;
			t8 = exp(t5+t6+teta);
			t10 = exp(t5+t6);
			t12 = exp(t5+teta);
			t14 = exp(t6+teta);
			t15 = t10-t12-t14+t2;
			t16 = t15*t15;
			
der2h.derp1p2 <- t1*t3*t8/t16-2*teta*t3*t8/t16/t15*(teta*t10-teta*t14);




			t2 = exp(teta);
			t3 = t2-1.0;
			t4 = teta*p2;
			t5 = teta*p1;
			t7 = exp(t4+t5+teta);
			t10 = exp(t4+t5);
			t12 = exp(t4+teta);
			t14 = exp(t5+teta);
			t15 = t10-t12-t14+t2;
			t16 = t15*t15;
			t17 = 1/t16;
			t21 = teta*t3;
der2h.derp1teta <- t3*t7*t17+teta*t2*t7*t17+t21*(p2+p1+1.0)*t7*t17-2.0*t21*t7/t15/t16*((p2+p1)*t10-(p2+1.0)*t12-(p1+1.0)*t14+t2);
		


			t1 = teta^2;
			t2 = exp(teta);
			t3 = t2 - 1;
			t5 = teta*p2;
			t6 = teta*p1;
			t8 = exp(t5+t6+teta);
			t10 = exp(t5+t6);
			t12 = exp(t5+teta);
			t14 = exp(t6+teta);
			t15 = t10-t12-t14+t2;
			t16 = t15*t15;
			
der2h.derp1p1 <-  t1*t3*t8/t16-2*teta*t3*t8/t16/t15*(teta*t10-teta*t14);





if(Cont == TRUE){


			t1 = exp(teta)
			t2 = teta*p2
			t3 = teta*p1
			t5 = exp(t2+t3+teta)
			t8 = exp(t2+t3)
			t10 = exp(t2+teta)
			t12 = exp(t3+teta)
			t13 = t8-t10-t12+t1
			t14 = t13*t13
			t15 = 1/t14
			t18 = t1-1.0
			t19 = p2+p1+1.0
			t21 = t5*t15
			t26 = 1/t14/t13
			t27 = p2 + p1
			t29 = p2+1.0
			t31 = p1+1.0
			t33 = t27*t8-t29*t10-t31*t12+t1
			t37 = teta*t1
			t43 = t5*t26
			t44 = t43*t33
			t47 = teta*t18
			t48 = t19*t19
			t11 = t14*t14
			t9 = t33*t33
			t7 = t27*t27
			t6 = t29*t29
			t4 = t31*t31
der2c.derrho.derrho <-	2.0*t1*t5*t15+2.0*t18*t19*t21-4.0*t18*t5*t26*t33+t37*t21+2.0*t37*t19*t5*t15-4.0*t37*t44+t47*t48*t5*t15-4.0*t47*t19*t44+6.0*t47*t5/t11*t9-2.0*t47*t43*(t7*t8-t6*t10-t4*t12+t1)






			  t1 = teta*teta
			  t3 = exp(teta)
			  t4 = t3-1.0
			  t6 = teta*p2
			  t7 = teta*p1
			  t9 = exp(t6+t7+teta)
			  t11 = exp(t6+t7)
			  t13 = exp(t6+teta)
			  t15 = exp(t7+teta)
			  t16 = t11-t13-t15+t3
			  t17 = t16*t16
			  t24 = t9/t17/t16
			  t27 = teta*t11-teta*t15
			  t31 = teta*t4
			  t32 = t17*t17
			  t35 = t27*t27
der2c.derp1.derp1 <- t1*teta*t4*t9/t17-4.0*t1*t4*t24*t27+6.0*t31*t9/t32*t35-2.0*t31*t24*(t1*t11-t1*t15)


 

 			  t1 = teta*teta
 			  t3 = exp(teta)
 			  t4 = t3-1.0
 			  t6 = teta*p1
 			  t7 = teta*p2
 			  t9 = exp(t6+t7+teta)
 			  t11 = exp(t6+t7)
 			  t13 = exp(t6+teta)
 			  t15 = exp(t7+teta)
 			  t16 = t11-t13-t15+t3
 			  t17 = t16*t16
 			  t24 = t9/t17/t16
 			  t27 = teta*t11-teta*t15
 			  t31 = teta*t4
 			  t32 = t17*t17
 			  t35 = t27*t27
der2c.derp2.derp2 <- t1*teta*t4*t9/t17-4.0*t1*t4*t24*t27+6.0*t31*t9/t32*t35-2.0*t31*t24*(t1*t11-t1*t15)





			t1 = teta*teta;
			t3 = exp(teta);
			t4 = t3-1.0;
			t5 = t1*teta*t4;
			t6 = teta*p2;
			t7 = teta*p1;
			t9 = exp(t6+t7+teta);
			t11 = exp(t6+t7);
			t13 = exp(t6+teta);
			t15 = exp(t7+teta);
			t16 = t11-t13-t15+t3;
			t17 = t16*t16;
			t21 = t1*t4;
			t24 = t9/t17/t16;
			t25 = teta*t11;
			t27 = t25-teta*t13;
			t32 = t25-teta*t15;
			t38 = t17*t17;
der2c.derp1.derp2 <- t5*t9/t17-2.0*t21*t24*t27-2.0*t21*t24*t32+6.0*teta*t4*t9/t38*t32*t27-2.0*t5*t24*t11


 
			  t1 = exp(teta);
			  t2 = t1-1.0;
			  t3 = teta*t2;
			  t4 = teta*p2;
			  t5 = teta*p1;
			  t7 = exp(t4+t5+teta);
			  t9 = exp(t4+t5);
			  t11 = exp(t4+teta);
			  t13 = exp(t5+teta);
			  t14 = t9-t11-t13+t1;
			  t15 = t14*t14;
			  t16 = 1/t15;
			  t17 = t7*t16;
			  t22 = 1/t15/t14;
			  t25 = teta*t9-teta*t13;
			  t29 = teta*teta;
			  t33 = t7*t22;
			  t34 = t33*t25;
			  t37 = t29*t2;
			  t38 = p2+p1+1.0;
			  t46 = p2+p1;
			  t49 = p1+1.0;
			  t51 = t46*t9-(p2+1.0)*t11-t49*t13+t1;
			  t56 = t15*t15;
der2c.derp1.derrho <- 2.0*t3*t17-2.0*t2*t7*t22*t25+t29*t1*t17-2.0*teta*t1*t34+t37*t38*t7*t16-2.0*t3*t38*t34-2.0*t37*t33*t51+6.0*t3*t7/t56*t51*t25-2.0*t3*t33*(t9+t46*teta*t9-t13-t49*teta*t13) 
 
 

 
			  t1 = exp(teta);
			  t2 = t1-1.0;
			  t3 = teta*t2;
			  t4 = teta*p1;
			  t5 = teta*p2;
			  t7 = exp(t4+t5+teta);
			  t9 = exp(t4+t5);
			  t11 = exp(t4+teta);
			  t13 = exp(t5+teta);
			  t14 = t9-t11-t13+t1;
			  t15 = t14*t14;
			  t16 = 1/t15;
			  t17 = t7*t16;
			  t22 = 1/t15/t14;
			  t25 = teta*t9-teta*t13;
			  t29 = teta*teta;
			  t33 = t7*t22;
			  t34 = t33*t25;
			  t37 = t29*t2;
			  t38 = p1+p2+1.0;
			  t46 = p1+p2;
			  t49 = p2+1.0;
			  t51 = t46*t9-(p1+1.0)*t11-t49*t13+t1;
			  t56 = t15*t15;
der2c.derp2.derrho <- 2.0*t3*t17-2.0*t2*t7*t22*t25+t29*t1*t17-2.0*teta*t1*t34+t37*t38*t7*t16-2.0*t3*t38*t34-2.0*t37*t33*t51+6.0*t3*t7/t56*t51*t25-2.0*t3*t33*(t9+t46*teta*t9-t13-t49*teta*t13) 
 


}







}



if(VC$BivD %in% c("G0","G90","G180","G270")){

 der2h.derp2teta <-  -((-log(p2))^(-2 + teta) * ((-log(p1))^teta + (-log(p2))^teta)^(1/teta - 
    3) * (((-log(p1))^teta * (((-log(p1))^teta + (-log(p2))^teta)^(1/teta) - 
    1) * (teta - (1 + log(p2))) + (-log(p2))^teta * (log(p2) - 
    (((-log(p1))^teta + (-log(p2))^teta)^(1/teta) * (log(p2) - 
        2) + ((-log(p1))^teta + (-log(p2))^teta)^(2/teta)))) * 
    ((-log(p1))^teta + (-log(p2))^teta) * log((-log(p1))^teta + 
    (-log(p2))^teta) + teta * ((-log(p1))^teta * (-log(p2))^teta * 
    (teta - (((-log(p1))^teta + (-log(p2))^teta)^(1/teta) * (3 * 
        teta - (1 + log(p2))) + (1 + teta) * log(p2) + 1 + teta * 
        (teta - 2)) * log(-log(p2))) + (-log(p1))^teta * ((-log(p2))^teta * 
    ((-1 + teta) * (log(p2) + teta) + ((-log(p1))^teta + (-log(p2))^teta)^(1/teta) * 
        (2 * teta + log(p2) - 2) + ((-log(p1))^teta + (-log(p2))^teta)^(2/teta)) - 
    (-log(p1))^teta * (((-log(p1))^teta + (-log(p2))^teta)^(1/teta) + 
        teta - 1) * (teta - (1 + log(p2)))) * log(-log(p1)) + 
    (-log(p2))^(2 * teta) * (((-log(p1))^teta + (-log(p2))^teta)^(1/teta) * 
        (log(p2) - 2) + ((-log(p1))^teta + (-log(p2))^teta)^(2/teta) - 
        log(p2)) * log(-log(p2)) + teta * (-log(p1))^(2 * teta) * 
    (1 + log(-log(p2)) * (teta - (1 + log(p2)))))) * exp(-((-log(p1))^teta + 
    (-log(p2))^teta)^(1/teta))/(p2^2 * teta^2))
        
  
  
t1 = log(p2)
t2 = (-t1)^teta 
t3 = log(p1)
t4 = (-t3)^teta  
t5 = t2+t4
t6 = 1/teta
t7 = t5^t6
t8 = exp(-t7)
t9 = t6-1
t10 = t5^t9
t11 = t8*t10
t12 = p2^2
t14 = 1/t12/p2
t15 = t2*t14
t20 = t1*t1
t21 = 1/t20
t26 = 1/t20/t1
t30 = t7*t7
t31 = t2*t2
t32 = t31*t2
t34 = t5*t5
t36 = 1/t34
t37 = t26*t36
t38 = t37*t11
t40 = t11*t2
t41 = teta*t14
t48 = t7*t31
t49 = t48*t14
t50 = 1/t5
t51 = t21*t50
t55 = t26*t50
t59 = t7*t32
t62 = t14*t26
t65 = t10*teta
t69 = t59*t62
t70 = t36*t8
t79 = t11*t9*t31
t86 = t9*t9
t18 = teta^2
t16 = t18*t14
t13 = t16*t37
 

 
der2h.derp2p2 <- -2*t11*t15/t1-3*t11*t15*t21-2*t11*t15*t26-t30*t32*t14*t38+3*t40*t41*t21+3*t40*t41*t26-3*t49*t51*t11-3*t49*t55*t11+t59*t14*t38+3*t48*t62*t50*t8*t65-t69*t70*t65+2*t69*t70*t10*t9*teta+3*t79*t41*t51+3*t79*t41*t55-t11*t86*t32*t13-3*t79*t16*t55+t11*t9*t32*t13-t40*t16*t26 
  
   
  
			  t1 = log(p2);
			  t2 = (-t1)^teta
			  t3 = log(p1);
			  t4 = (-t3)^teta
			  t5 = t2+t4;
			  t6 = 1/teta;
			  t7 = t5^t6
			  t8 = teta^2
			  t9 = 1/t8;
			  t10 = log(t5);
			  t11 = t9*t10;
			  t12 = log(-t1);
			  t13 = t2*t12;
			  t14 = log(-t3);
			  t16 = t13+t4*t14;
			  t18 = 1/t5;
			  t20 = -t11+t6*t16*t18;
			  t21 = t20^2
			  t23 = exp(-t7);
			  t25 = t6-1.0;
			  t26 = t5^t25
			  t28 = 1/p2
			  t29 = 1/t1;
			  t30 = t28*t29;
			  t31 = t26*t2*t30;
			  t36 = 2.0/t8/teta*t10;
			  t39 = 2.0*t9*t16*t18;
			  t40 = t12^2;
			  t42 = t14^2;
			  t44 = t2*t40+t4*t42;
			  t47 = t16^2
			  t49 = t5^2;
			  t50 = 1/t49;
			  t56 = t7*t7;
			  t61 = t23*t26;
			  t62 = t7*t20*t61;
			  t65 = -t11+t25*t16*t18;
			  t70 = t13*t30;
			  t73 = t65*t65;
			  t15 = t2*t28*t29;
			  
der2h.derteta.teta.st <- t7*t21*t23*t31+t7*(t36-t39+t6*t44*t18-t6*t47*t50)*t23*t31-t56*t21*t23*t31+2.0*t62*t65*t2*t30+2.0*t62*t70-t61*t73*t15-t61*(t36-t39+t25*t44*t18-t25*t47*t50)*t15-2.0*t61*t65*t70-t61*t2*t40*t28*t29;  
  
  
  
  
  
  
  			t3 = log(p2);
  			t4 = (-t3)^teta;
  			t5 = log(p1);
  			t6 = (-t5)^teta
  			t7 = t4+t6;
  			t8 = 1/teta;
  			t9 = t7^t8
  			t11 = p2^2;
  			t12 = 1/t11;
  			t13 = 1/t3;
  			t15 = 1/t7;
  			t18 = exp(-t9);
  			t19 = 1/p1;
  			t21 = -1 + t8;
  			t22 = (t7)^(2*t21)
  			t24 = teta - 1
  			t25 = (t3*t5)^t24
  			t27 = t7^-t8
  			t28 = t24*t27;
  			t29 = 1 + t28;
  			t30 = t22*t25*t29;
  			t33 = t18*t12;
  			t36 = t19*t22;
  			
der2h.derp1p2 <-  -t9*t4*t12*t13*t15*t18*t19*t30-t33*t19*t30+2.0*t33*t36*t21*t4*teta*t13*t15*t25*t29+t33*t36*t25*t24*t13*t29-t33*t36*t25*t28*t4*t13*t15;
  
  
  
  
  			t3 = log(p1);
  			  t4 = (-t3)^teta
  			  t5 = log(p2);
  			  t6 = (-t5)^teta
  			  t7 = t4+t6;
  			  t8 = 1/teta;
  			  t9 = t7^t8
  			  t10 = teta^2;
  			  t12 = log(t7);
  			  t13 = 1/t10*t12;
  			  t14 = log(-t3);
  			  t16 = log(-t5);
  			  t18 = t4*t14+t6*t16;
  			  t20 = 1/t7;
  			  t22 = -t13+t8*t18*t20;
  			  t24 = exp(-t9);
  			  t26 = t24/p1;
  			  t28 = 1/p2;
  			  t29 = -1.0+t8;
  			  t30 = t7^(2*t29)
  			  t32 = t3*t5;
  			  t33 = teta-1.0;
  			  t34 = t32^t33
  			  t35 = t7^(-1.0*t8)
  			  t36 = t33*t35;
  			  t17 = 1.0+t36;
  			  t15 = t34*t17;
  			  t11 = t26*t28;
  			  t2 = t30*t34;
  			  t1 = log(t32);
  			  
 der2h.derp1teta <- -t9*t22*t26*t28*t30*t15+t11*t30*(-2.0*t13+2.0*t29*t18*t20)*t15+t11*t2*t1*t17+t11*t2*(t35-t36*t22);
  

  			t3 = log(p1);
  			t4 = (-t3)^teta;
  			t5 = log(p2);
  			t6 = (-t5)^teta
  			t7 = t4+t6;
  			t8 = 1/teta;
  			t9 = t7^t8
  			t11 = p1^2;
  			t12 = 1/t11;
  			t13 = 1/t3;
  			t15 = 1/t7;
  			t18 = exp(-t9);
  			t19 = 1/p2;
  			t21 = -1 + t8;
  			t22 = (t7)^(2*t21)
  			t24 = teta - 1
  			t25 = (t3*t5)^t24
  			t27 = t7^-t8
  			t28 = t24*t27;
  			t29 = 1 + t28;
  			t30 = t22*t25*t29;
  			t33 = t18*t12;
  			t36 = t19*t22;
  			
der2h.derp1p1 <-  -t9*t4*t12*t13*t15*t18*t19*t30-t33*t19*t30+2.0*t33*t36*t21*t4*teta*t13*t15*t25*t29+t33*t36*t25*t24*t13*t29-t33*t36*t25*t28*t4*t13*t15;  



if(Cont == TRUE){



                        t3 = log(p1);
			t4 = (-t3)^teta
			t5 = log(p2)
			t6 = (-t5)^teta
			t7 = t4+t6
			t8 = 1/teta
			t9 = t7^t8
			t10 = teta*teta
			t11 = 1/t10
			t12 = log(t7)
			t13 = t11*t12
			t14 = log(-t3)
			t16 = log(-t5)
			t18 = t4*t14+t6*t16
			t20 = 1/t7
			t22 = -t13+t8*t18*t20
			t23 = t22*t22
			t25 = exp(-t9)
			t27 = t25/p1
			t29 = 1/p2
			t30 = -1.0+t8
			t31 = t7^(2*t30)
			t32 = t29*t31
			t33 = t3*t5
			t34 = teta-1.0
			t35 = t33^t34
			t36 = t7^-t8
			t37 = t34*t36
			t38 = 1.0+t37
			t39 = t35*t38
			t40 = t32*t39
			t44 = 1/t10/teta*t12
			t47 = t11*t18*t20
			t49 = t14*t14
			t51 = t16*t16
			t53 = t4*t49+t6*t51
			t56 = t18*t18
			t58 = t7*t7
			t59 = 1/t58
			t61 = 2.0*t44-2.0*t47+t8*t53*t20-t8*t56*t59
			t65 = t9*t9
			t70 = t9*t22*t27
			t74 = -2.0*t13+2.0*t30*t18*t20
			t75 = t74*t35
			t80 = log(t33)
			t87 = t36-t37*t22
			t88 = t35*t87
			t17 = t27*t29
			t15 = t74*t74
			t2 = t31*t35
			t1 = t80*t80
der2c.derrho.derrho <- -t9*t23*t27*t40-t9*t61*t27*t40+t65*t23*t27*t40-2.0*t70*t32*t75*t38-2.0*t70*t32*t35*t80*t38-2.0*t70*t32*t88+t17*t31*t15*t39+t17*t31*(4.0*t44-4.0*t47+2.0*t30*t53*t20-2.0*t30*t56*t59)*t39+2.0*t27*t32*t75*t80*t38+2.0*t17*t31*t74*t88+t17*t2*t1*t38+2.0*t17*t2*t80*t87+t17*t2*(-2.0*t36*t22+t37*t23-t37*t61)



			  t3 = log(p1)
			  t4 = (-t3)^teta
			  t5 = log(p2)
			  t6 = (-t5)^teta
			  t7 = t4+t6
			  t8 = 1/teta
			  t9 = t7^t8
			  t10 = exp(-t9)
			  t11 = p1*p1
			  t13 = 1/t11/p1
			  t14 = t10*t13
			  t15 = 1/p2
			  t16 = t14*t15
			  t17 = -1.0+t8
			  t18 = t7^(2*t17)
			  t20 = teta-1.0
			  t21 = (t3*t5)^t20
			  t23 = t7^-t8
			  t24 = t20*t23
			  t25 = t24+1.0
			  t26 = t18*t21*t25
			  t29 = t9*t4
			  t30 = 1/t3
			  t32 = 1/t7
			  t35 = t10*t15
			  t36 = t35*t26
			  t39 = t3*t3
			  t40 = 1/t39
			  t41 = t13*t40
			  t43 = t29*t41*t32
			  t45 = t15*t18
			  t46 = t14*t45
			  t47 = t21*t20
			  t52 = t4*t4
			  t53 = t9*t52
			  t54 = t7*t7
			  t55 = 1/t54
			  t56 = t41*t55
			  t57 = t53*t56
			  t66 = t35*t18
			  t68 = t21*t25*teta
			  t71 = t9*t9
			  t79 = 2.0*t45*t17
			  t83 = t47*t25
			  t87 = t47*t23
			  t91 = t14*t79
			  t92 = t4*teta
			  t95 = t32*t21*t25
			  t100 = t14*t45*t21
			  t106 = 2.0*t16*t26+3.0*t29*t13*t30*t32*t36+t43*t36-3.0*t46*t47*t30*t25-t57*t36-t29*teta*t13*t40*t32*t10*t15*t26+t57*t66*t68+t71*t52*t56*t36-2.0*t53*t13*t40*t55*t10*t79*t68-2.0*t43*t66*t83+2.0*t57*t66*t87-3.0*t91*t92*t30*t95+3.0*t100*t24*t4*t30*t32
			  t113 = teta*teta
			  t118 = t52*t113*t40*t55*t21*t25
			  t125 = 2.0*t18*t17
			  t128 = teta*t40
			  t129 = t128*t32
			  t19 = t128*t55
			  t12 = t40*t25
			  t2 = t20*t20
			  t1 = -t91*t92*t40*t95+4.0*t14*t45*t17*t17*t118+t91*t4*t113*t40*t95-t91*t118+2.0*t16*t125*t4*t129*t83-2.0*t16*t125*t52*t19*t87-t46*t47*t12+t46*t21*t2*t12-2.0*t100*t2*t40*t23*t4*t32+t100*t24*t4*t40*t32+t100*t24*t52*t40*t55-t100*t24*t4*t129+t100*t24*t52*t19

der2c.derp1.derp1 <- t106+t1
 

 

			  t3 = log(p2)
			  t4 = (-t3)^teta
			  t5 = log(p1)
			  t6 = (-t5)^teta
			  t7 = t4+t6
			  t8 = 1/teta
			  t9 = t7^t8
			  t10 = exp(-t9)
			  t11 = p2^2
			  t13 = 1/t11/p2
			  t14 = t10*t13
			  t15 = 1/p1
			  t16 = t14*t15
			  t17 = -1.0+t8
			  t18 = t7^(2*t17)
			  t20 = teta-1.0
			  t21 = (t3*t5)^t20
			  t23 = t7^-t8
			  t24 = t20*t23
			  t25 = t24+1.0
			  t26 = t18*t21*t25
			  t29 = t9*t4
			  t30 = 1/t3
			  t32 = 1/t7
			  t35 = t10*t15
			  t36 = t35*t26
			  t39 = t3*t3
			  t40 = 1/t39
			  t41 = t13*t40
			  t43 = t29*t41*t32
			  t45 = t15*t18
			  t46 = t14*t45
			  t47 = t21*t20
			  t52 = t4*t4
			  t53 = t9*t52
			  t54 = t7*t7
			  t55 = 1/t54
			  t56 = t41*t55
			  t57 = t53*t56
			  t66 = t35*t18
			  t68 = t21*t25*teta
			  t71 = t9*t9
			  t79 = 2.0*t45*t17
			  t83 = t47*t25
			  t87 = t47*t23
			  t91 = t14*t79
			  t92 = t4*teta
			  t95 = t32*t21*t25
			  t100 = t14*t45*t21
			  t106 = 2.0*t16*t26+3.0*t29*t13*t30*t32*t36+t43*t36-3.0*t46*t47*t30*t25-t57*t36-t29*teta*t13*t40*t32*t10*t15*t26+t57*t66*t68+t71*t52*t56*t36-2.0*t53*t13*t40*t55*t10*t79*t68-2.0*t43*t66*t83+2.0*t57*t66*t87-3.0*t91*t92*t30*t95+3.0*t100*t24*t4*t30*t32
			  t113 = teta*teta
			  t118 = t52*t113*t40*t55*t21*t25
			  t125 = 2.0*t18*t17
			  t128 = teta*t40
			  t129 = t128*t32
			  t19 = t128*t55
			  t12 = t40*t25
			  t2 = t20*t20
			  t1 = -t91*t92*t40*t95+4.0*t14*t45*t17*t17*t118+t91*t4*t113*t40*t95-t91*t118+2.0*t16*t125*t4*t129*t83-2.0*t16*t125*t52*t19*t87-t46*t47*t12+t46*t21*t2*t12-2.0*t100*t2*t40*t23*t4*t32+t100*t24*t4*t40*t32+t100*t24*t52*t40*t55-t100*t24*t4*t129+t100*t24*t52*t19

der2c.derp2.derp2 <- t106+t1


			  t3 = log(p1);
			  t4 = (-t3)^teta
			  t5 = log(p2);
			  t6 = (-t5)^teta
			  t7 = t4+t6;
			  t8 = 1/teta;
			  t9 = t7^t8
			  t10 = exp(-t9);
			  t11 = p1*p1;
			  t12 = 1/t11;
			  t13 = t10*t12;
			  t14 = p2*p2;
			  t15 = 1/t14;
			  t16 = t13*t15;
			  t17 = -1.0+t8;
			  t18 = t7^(2*t17)
			  t20 = teta-1.0;
			  t21 = (t3*t5)^t20
			  t22 = t18*t21;
			  t23 = t7^-t8
			  t24 = t20*t23;
			  t25 = 1.0+t24;
			  t26 = t22*t25;
			  t28 = t9*t4;
			  t29 = 1/t3;
			  t30 = t12*t29;
			  t31 = 1/t7;
			  t34 = t10*t15;
			  t37 = t9*t6;
			  t38 = t37*t15;
			  t39 = 1/t5;
			  t40 = t7*t7;
			  t41 = 1/t40;
			  t48 = t28*t12;
			  t49 = t29*t41;
			  t51 = t48*t49*t10;
			  t52 = t15*t18;
			  t53 = t52*t21;
			  t55 = teta*t39;
			  t59 = t9*t9;
			  t64 = t15*t39;
			  t70 = 2.0*t18*t17;
			  t71 = t70*t6;
			  t72 = t21*t25;
			  t84 = t6*t39;
			  t85 = t24*t84;
			  t94 = 2.0*t13*t52*t17;
			  t98 = t31*t21*t25;
			  t101 = t13*t52;
			  t102 = t21*t20;
			  t104 = t102*t39*t25;
			  t106 = t13*t53;
			  t107 = t84*t31;
			  t110 = t4*teta;
			  t114 = t16*t26+t28*t30*t31*t34*t26-t38*t39*t41*t4*t30*t10*t26+t51*t53*t25*t6*t55+t59*t4*t12*t49*t6*t64*t10*t26-2.0*t48*t49*t34*t71*t55*t72-t48*t29*t31*t10*t53*t20*t39*t25+2.0*t51*t53*t85+t37*t64*t31*t13*t26-t94*t6*teta*t39*t98-t101*t104+t106*t24*t107-t94*t110*t29*t98;
			  t118 = teta*teta;
			  t128 = t4*t29;
			  t129 = t16*t70*t4;
			  t32 = teta*t29*t31*t104;
			  t27 = t20*t20;
			  t19 = t128*t31;
			  t2 = t16*t22*t20;
			  t1 = 4.0*t16*t18*t17*t17*t6*t118*t39*t41*t128*t72-t129*t118*t29*t41*t72*t84+t129*t32-2.0*t16*t70*t110*t49*t21*t85-t101*t102*t29*t25-t38*t39*t31*t10*t12*t18*t21*t20*t29*t25+t16*t71*t32+t101*t21*t27*t39*t29*t25-t106*t27*t29*t23*t107+t106*t24*t19-t106*t27*t39*t23*t19+t2*t23*t6*t39*t41*t4*t29+t2*t23*t4*t29*t41*t6*t55

der2c.derp1.derp2 <- t114+t1
 



 			  t3 = log(p1);
 			  t4 = (-t3)^teta
 			  t5 = log(p2);
 			  t6 = (-t5)^teta
 			  t7 = t4+t6;
 			  t8 = 1/teta;
 			  t9 = t7^t8
 			  t11 = p1*p1;
 			  t12 = 1/t11;
 			  t13 = 1/t3;
 			  t15 = 1/t7;
 			  t17 = t9*t4*t12*t13*t15;
 			  t18 = teta*teta;
 			  t20 = log(t7);
 			  t21 = 1/t18*t20;
 			  t22 = log(-t3);
 			  t24 = log(-t5);
 			  t26 = t4*t22+t6*t24;
 			  t29 = -t21+t8*t26*t15;
 			  t30 = exp(-t9);
 			  t32 = 1/p2;
 			  t34 = -1.0+t8;
 			  t35 = t7^(2*t34)
 			  t36 = t3*t5;
 			  t37 = teta-1.0;
 			  t38 = t36^t37
 			  t39 = t35*t38;
 			  t40 = t7^-t8
 			  t41 = t37*t40;
 			  t42 = t41+1.0;
 			  t43 = t39*t42;
 			  t47 = 1/p1;
 			  t48 = t47*t13;
 			  t49 = t48*t15;
 			  t50 = t8*t4*t49;
 			  t51 = t4*teta;
 			  t55 = t4*t47*t13;
 			  t56 = t51*t48*t22+t55;
 			  t59 = t7*t7;
 			  t60 = 1/t59;
 			  t63 = -t50+t8*t56*t15-t26*t60*t55;
 			  t65 = t30*t47;
 			  t67 = t32*t35;
 			  t68 = t38*t42;
 			  t69 = t67*t68;
 			  t71 = t9*t9;
 			  t80 = t9*t29;
 			  t81 = t30*t12;
 			  t87 = t80*t30*t12*t32*t35;
 			  t94 = t81*t32;
 			  t97 = t37*t13*t42;
 			  t100 = t38*t37;
 			  t103 = t4*t13*t15;
 			  t104 = t100*t40*t103;
 			  t106 = t30*t32;
 			  t107 = t106*t35;
 			  t109 = 2.0*t34*t26;
 			  t111 = -2.0*t21+t109*t15;
 			  t112 = t111*t38;
 			  t113 = t112*t42;
 			  t121 = 2.0*t94*t35*t34*t4;
 			  t123 = teta*t13*t15;
 			  t126 = t65*t32;
 			  t137 = t81*t67;
 			  t140 = -t17*t29*t30*t32*t43-t9*t63*t65*t69+t71*t29*t4*t12*t13*t15*t30*t32*t43+t80*t81*t69-2.0*t87*t34*t4*teta*t13*t15*t68-t80*t94*t39*t97+t87*t104-t17*t107*t113-t94*t35*t111*t68+t121*t123*t113+t126*t35*(-2.0*t50+2.0*t34*t56*t15-t109*t60*t51*t48)*t68+t137*t112*t97;
 			  t144 = log(t36);
 			  t146 = t38*t144*t42;
 			  t10 = t40-t41*t29;
 			  t2 = t39*t10;
 			  t1 = -t81*t67*t111*t104-t17*t107*t146-t94*t39*t144*t42+t121*t123*t146+t137*t100*t13*t144*t42+t94*t39*t13*t42-t81*t67*t38*t144*t37*t40*t103-t17*t106*t2-t94*t2+2.0*t81*t67*t34*t51*t13*t15*t38*t10+t137*t100*t13*t10+t126*t39*(-t40*t4*t49+t41*t4*t48*t15*t29-t41*t63);

der2c.derp1.derrho <- t140+t1






 			  t3 = log(p2);
 			  t4 = (-t3)^teta
 			  t5 = log(p1);
 			  t6 = (-t5)^teta
 			  t7 = t4+t6;
 			  t8 = 1/teta;
 			  t9 = t7^t8
 			  t11 = p2^2;
 			  t12 = 1/t11;
 			  t13 = 1/t3;
 			  t15 = 1/t7;
 			  t17 = t9*t4*t12*t13*t15;
 			  t18 = teta*teta;
 			  t20 = log(t7);
 			  t21 = 1/t18*t20;
 			  t22 = log(-t3);
 			  t24 = log(-t5);
 			  t26 = t4*t22+t6*t24;
 			  t29 = -t21+t8*t26*t15;
 			  t30 = exp(-t9);
 			  t32 = 1/p1;
 			  t34 = -1.0+t8;
 			  t35 = t7^(2*t34)
 			  t36 = t3*t5;
 			  t37 = teta-1.0;
 			  t38 = t36^t37
 			  t39 = t35*t38;
 			  t40 = t7^-t8
 			  t41 = t37*t40;
 			  t42 = t41+1.0;
 			  t43 = t39*t42;
 			  t47 = 1/p2;
 			  t48 = t47*t13;
 			  t49 = t48*t15;
 			  t50 = t8*t4*t49;
 			  t51 = t4*teta;
 			  t55 = t4*t47*t13;
 			  t56 = t51*t48*t22+t55;
 			  t59 = t7*t7;
 			  t60 = 1/t59;
 			  t63 = -t50+t8*t56*t15-t26*t60*t55;
 			  t65 = t30*t47;
 			  t67 = t32*t35;
 			  t68 = t38*t42;
 			  t69 = t67*t68;
 			  t71 = t9*t9;
 			  t80 = t9*t29;
 			  t81 = t30*t12;
 			  t87 = t80*t30*t12*t32*t35;
 			  t94 = t81*t32;
 			  t97 = t37*t13*t42;
 			  t100 = t38*t37;
 			  t103 = t4*t13*t15;
 			  t104 = t100*t40*t103;
 			  t106 = t30*t32;
 			  t107 = t106*t35;
 			  t109 = 2.0*t34*t26;
 			  t111 = -2.0*t21+t109*t15;
 			  t112 = t111*t38;
 			  t113 = t112*t42;
 			  t121 = 2.0*t94*t35*t34*t4;
 			  t123 = teta*t13*t15;
 			  t126 = t65*t32;
 			  t137 = t81*t67;
 			  t140 = -t17*t29*t30*t32*t43-t9*t63*t65*t69+t71*t29*t4*t12*t13*t15*t30*t32*t43+t80*t81*t69-2.0*t87*t34*t4*teta*t13*t15*t68-t80*t94*t39*t97+t87*t104-t17*t107*t113-t94*t35*t111*t68+t121*t123*t113+t126*t35*(-2.0*t50+2.0*t34*t56*t15-t109*t60*t51*t48)*t68+t137*t112*t97;
 			  t144 = log(t36);
 			  t146 = t38*t144*t42;
 			  t10 = t40-t41*t29;
 			  t2 = t39*t10;
 			  t1 = -t81*t67*t111*t104-t17*t107*t146-t94*t39*t144*t42+t121*t123*t146+t137*t100*t13*t144*t42+t94*t39*t13*t42-t81*t67*t38*t144*t37*t40*t103-t17*t106*t2-t94*t2+2.0*t81*t67*t34*t51*t13*t15*t38*t10+t137*t100*t13*t10+t126*t39*(-t40*t4*t49+t41*t4*t48*t15*t29-t41*t63);




der2c.derp2.derrho <- t140+t1


}




}






if(VC$BivD %in% c("J0","J90","J180","J270")){

 der2h.derp2teta <- (((-1 + teta) * (((1 - p1)^teta - 1) * (1 - p2)^teta - 
    (1 - p1)^teta) * log((1 - p1)^teta - ((1 - p1)^teta - 1) * 
    (1 - p2)^teta) + teta * (((-1 + teta)^2 * log(1 - p2) - teta) * 
    ((1 - p1)^teta - 1) * (1 - p2)^teta + teta * ((-1 + teta) * 
    log(1 - p2) + 1) * (1 - p1)^teta)) * ((1 - p1)^teta - 1) + 
    teta * (-1 + teta) * (((1 - p1)^teta + teta - 1) * (1 - p1)^teta - 
        ((1 - p1)^teta - 1) * ((1 - p1)^teta - teta) * (1 - p2)^teta) * 
        log(1 - p1)) * ((1 - p1)^teta - ((1 - p1)^teta - 1) * 
    (1 - p2)^teta)^(1/teta - 3) * (1 - p1)^teta * (1 - p2)^(-2 + 
    teta)/teta^2
 
			t2 = (1 - p1)^teta
			t3 = 1 - p2
			t4 = t3^teta
			t5 = t2*t4
			t6 = t2+t4-t5
			t8 = 1/teta - 1
			t9 = t6^t8
			t10 = t8^2
			t12 = t4*teta
			t13 = 1/t3
			t14 = -t12*t13+t5*teta*t13
			t18 = t14^2
			t20 = t6^2
			t22 = teta - 1
			t23 = t3^t22
			t24 = 1 - t2
			t26 = 1/t20*t23*t24
			t27 = t9*t8
			t29 = teta^2
			t31 = t3*t3
			t32 = 1/t31
			t41 = 1/t6
			t51 = t9*t23
			t55 = t22*t22
			
der2h.derp2p2 <- t9*t10*t18*t26+t27*(t4*t29*t32-t12*t32-t5*t29*t32+t5*teta*t32)*t41*t23*t24-t27*t18*t26-2.0*t27*t14*t41*t23*t22*t13*t24+t51*t55*t32*t24-t51*t22*t32*t24 
 
 
 


			  t1 = 1.0-p1;
			  t2 = t1^teta
			  t3 = 1.0-p2;
			  t4 = t3^teta;
			  t5 = t2*t4;
			  t6 = t2+t4-t5;
			  t8 = 1/teta-1.0;
			  t9 = t6^t8;
			  t10 = teta^2
			  t11 = 1/t10;
			  t12 = log(t6);
			  t14 = log(t1);
			  t15 = t2*t14;
			  t16 = log(t3);
			  t18 = t4*t16;
			  t20 = t15+t18-t15*t4-t5*t16;
			  t21 = 1/t6;
			  t23 = -t11*t12+t8*t20*t21;
			  t25 = t23^2;
			  t28 = t3^(teta-1)
			  t29 = 1 - t2;
			  t30 = t28*t29;
			  t39 = t14^2;
			  t40 = t2*t39;
			  t42 = t16^2;
			  t50 = t20^2;
			  t56 = t6^2;
			  t13 = t9*t23;
			  t7 = t9*t28;
			  
			  
der2h.derteta.teta.st <- t9*t25*t30+t9*(2.0/t10/teta*t12-2.0*t11*t20*t21+t8*(t40+t4*t42-t40*t4-2.0*t15*t18-t5*t42)*t21-t8*t50/t56)*t30+2.0*t13*t28*t16*t29-2.0*t13*t28*t2*t14+t7*t42*t29-2.0*t7*t16*t2*t14-t7*t40;







t1 = 1.0-p2;
t2 = t1^(1.0*teta);
t3 = 1.0-p1;
t4 = t3^(1.0*teta);
t5 = t2*t4;
t6 = t2+t4-t5;
t8 = 1/teta-2.0;
t9 = t6^(1.0*t8);
t11 = t2*teta;
t12 = 1/t1;
t16 = -t11*t12+t11*t12*t4;
t19 = teta-1.0;
t20 = t1^(1.0*t19);
t22 = t3^(1.0*t19);
t23 = teta-1.0+t2+t4-t5;
t27 = t9*t20;

der2h.derp1p2 <- t9*t8*t16/t6*t20*t22*t23-t27*t19*t12*t22*t23+t27*t22*t16;




t1 = 1.0-p1;
t2 = t1^(1.0*teta);
t3 = 1.0-p2;
t4 = t3^(1.0*teta);
t5 = t2*t4;
t6 = t2+t4-t5;
t8 = 1/teta-2.0;
t9 = t6^(1.0*t8);
t10 = teta^2;
t11 = log(t6);
t12 = log(t1);
t13 = t2*t12;
t14 = log(t3);
t15 = t4*t14;
t16 = t13*t4;
t19 = t5*t14;
t21 = teta-1.0;
t27 = t1^(1.0*t21);
t28 = t3^(1.0*t21);
t30 = teta-1.0+t2+t4-t5;
t33 = t9*t27;

der2h.derp1teta <- t9*(-1/t10*t11+t8*(t13+t15-t16-t19)/t6)*t27*t28*t30+t33*t12*t28*t30+t33*t28*t14*t30+t33*t28*(1.0+t13+t15-t16-t19);
            

t1 = 1.0-p1;
t2 = t1^(1.0*teta);
t3 = 1.0-p2;
t4 = t3^(1.0*teta);
t5 = t2*t4;
t6 = t2+t4-t5;
t8 = 1/teta-2.0;
t9 = t6^(1.0*t8);
t11 = t2*teta;
t12 = 1/t1;
t16 = -t11*t12+t11*t12*t4;
t19 = teta-1.0;
t20 = t1^(1.0*t19);
t22 = t3^(1.0*t19);
t23 = teta-1.0+t2+t4-t5;
t27 = t9*t20;

der2h.derp1p1 <- t9*t8*t16/t6*t20*t22*t23-t27*t19*t12*t22*t23+t27*t22*t16;


if(Cont == TRUE){


                        t1 = 1.0-p1
			t2 = t1^teta
			t3 = 1.0-p2
			t4 = t3^teta
			t5 = t2*t4
			t6 = t2+t4-t5
			t8 = 1/teta-2.0
			t9 = t6^t8
			t10 = teta^2
			t11 = 1/t10
			t12 = log(t6)
			t14 = log(t1)
			t15 = t2*t14
			t16 = log(t3)
			t17 = t4*t16
			t18 = t15*t4
			t19 = t5*t16
			t20 = t15+t17-t18-t19
			t22 = 1/t6
			t24 = -t11*t12+t8*t20*t22
			t25 = t24*t24
			t27 = teta-1.0
			t28 = t1^t27
			t29 = t3^t27
			t30 = t28*t29
			t31 = teta-1.0+t2+t4-t5
			t32 = t30*t31
			t41 = t14*t14
			t42 = t2*t41
			t43 = t16*t16
			t49 = t42+t4*t43-t42*t4-2.0*t15*t17-t5*t43
			t51 = t20*t20
			t53 = t6*t6
			t60 = t9*t24
			t61 = t60*t28
			t62 = t14*t29
			t66 = t29*t16
			t67 = t66*t31
			t70 = 1.0+t15+t17-t18-t19
			t74 = t9*t28
der2c.derrho.derrho <- t9*t25*t32+t9*(2.0/t10/teta*t12-2.0*t11*t20*t22+t8*t49*t22-t8*t51/t53)*t32+2.0*t61*t62*t31+2.0*t61*t67+2.0*t60*t30*t70+t74*t41*t29*t31+2.0*t74*t14*t67+2.0*t74*t62*t70+t74*t29*t43*t31+2.0*t74*t66*t70+t74*t29*t49


 




			  t1 = 1.0-p1;
			  t2 = t1^teta
			  t3 = 1.0-p2;
			  t4 = t3^teta
			  t5 = t2*t4;
			  t6 = t2+t4-t5;
			  t8 = 1/teta-2.0;
			  t9 = t6^t8
			  t10 = t8*t8;
			  t12 = t2*teta;
			  t13 = 1/t1;
			  t17 = -t12*t13+t12*t13*t4;
			  t18 = t17*t17;
			  t20 = t6*t6;
			  t22 = teta-1.0;
			  t23 = t1^t22
			  t25 = t3^t22
			  t26 = teta-1.0+t2+t4-t5;
			  t27 = t25*t26;
			  t28 = 1/t20*t23*t27;
			  t30 = t9*t8;
			  t31 = teta*teta;
			  t32 = t2*t31;
			  t33 = t1*t1;
			  t34 = 1/t33;
			  t37 = t34*t4;
			  t40 = t32*t34-t12*t34-t32*t37+t12*t37;
			  t42 = 1/t6;
			  t43 = t42*t23;
			  t46 = t30*t18;
			  t16 = t13*t25;
			  t15 = t9*t23;
			  t14 = t22*t22;
			  t11 = t34*t25*t26;
			  t7 = t15*t22;
der2c.derp1.derp1 <- t9*t10*t18*t28+t30*t40*t43*t27-t46*t28-2.0*t30*t17*t42*t23*t22*t16*t26+2.0*t46*t43*t25+t15*t14*t11-t7*t11-2.0*t7*t16*t17+t15*t25*t40;




			  t1 = 1.0-p2;
			  t2 = t1^teta
			  t3 = 1.0-p1;
			  t4 = t3^teta
			  t5 = t2*t4;
			  t6 = t2+t4-t5;
			  t8 = 1/teta-2.0;
			  t9 = t6^t8
			  t10 = t8*t8;
			  t12 = t2*teta;
			  t13 = 1/t1;
			  t17 = -t12*t13+t12*t13*t4;
			  t18 = t17*t17;
			  t20 = t6*t6;
			  t22 = teta-1.0;
			  t23 = t1^t22
			  t25 = t3^t22
			  t26 = teta-1.0+t2+t4-t5;
			  t27 = t25*t26;
			  t28 = 1/t20*t23*t27;
			  t30 = t9*t8;
			  t31 = teta*teta;
			  t32 = t2*t31;
			  t33 = t1*t1;
			  t34 = 1/t33;
			  t37 = t34*t4;
			  t40 = t32*t34-t12*t34-t32*t37+t12*t37;
			  t42 = 1/t6;
			  t43 = t42*t23;
			  t46 = t30*t18;
			  t16 = t13*t25;
			  t15 = t9*t23;
			  t14 = t22*t22;
			  t11 = t34*t25*t26;
			  t7 = t15*t22;
der2c.derp2.derp2 <- t9*t10*t18*t28+t30*t40*t43*t27-t46*t28-2.0*t30*t17*t42*t23*t22*t16*t26+2.0*t46*t43*t25+t15*t14*t11-t7*t11-2.0*t7*t16*t17+t15*t25*t40;



 
			t1 = 1.0-p1;
			t2 = t1^teta
			t3 = 1.0-p2;
			t4 = t3^teta
			t5 = t2*t4;
			t6 = t2+t4-t5;
			t8 = 1/teta-2.0;
			t9 = t6^t8
			t10 = t8*t8;
			t13 = 1/t3;
			t17 = -t4*teta*t13+t5*teta*t13;
			t18 = t6*t6;
			t19 = 1/t18;
			t22 = t2*teta;
			t23 = 1/t1;
			t27 = -t22*t23+t22*t23*t4;
			t28 = teta-1.0;
			t29 = t1^t28
			t31 = t3^t28
			t32 = teta-1.0+t2+t4-t5;
			t36 = t9*t8;
			t37 = teta*teta;
			t41 = t4*t13;
			t42 = 1/t6;
			t45 = t29*t31;
			t55 = t28*t13;
			t71 = t23*t31;
			t72 = t9*t29;
			t7 = t28*t28;
der2c.derp1.derp2 <- t9*t10*t17*t19*t27*t29*t31*t32-t36*t2*t37*t23*t41*t42*t45*t32-t36*t27*t19*t45*t32*t17-t36*t27*t42*t45*t55*t32+2.0*t36*t27*t42*t29*t31*t17-t36*t17*t42*t29*t28*t71*t32+t72*t7*t71*t13*t32-t72*t28*t71*t17-t72*t31*t55*t27-t72*t31*t2*t37*t23*t41;



			  t1 = 1.0-p1;
			  t2 = t1^teta
			  t3 = 1.0-p2;
			  t4 = t3^teta
			  t5 = t2*t4;
			  t6 = t2+t4-t5;
			  t8 = 1/teta-2.0;
			  t9 = t6^t8
			  t10 = t9*t8;
			  t11 = t2*teta;
			  t12 = 1/t1;
			  t14 = t12*t4;
			  t16 = -t11*t12+t11*t14;
			  t17 = 1/t6;
			  t19 = t10*t16*t17;
			  t20 = teta*teta;
			  t21 = 1/t20;
			  t22 = log(t6);
			  t24 = log(t1);
			  t25 = t2*t24;
			  t26 = log(t3);
			  t27 = t4*t26;
			  t28 = t25*t4;
			  t29 = t5*t26;
			  t31 = t8*(t25+t27-t28-t29);
			  t33 = -t21*t22+t31*t17;
			  t34 = teta-1.0;
			  t35 = t1^t34
			  t37 = t3^t34
			  t38 = teta-1.0+t2+t4-t5;
			  t39 = t37*t38;
			  t44 = t12*t24;
			  t46 = t2*t12;
			  t52 = -t11*t44-t46+t11*t44*t4+t46*t4+t11*t14*t26;
			  t55 = t6*t6;
			  t61 = t35*t37;
			  t64 = t9*t33;
			  t71 = t9*t35;
			  t80 = t71*t34;
			  t81 = t12*t37;
			  t87 = t26*t38;
			  t94 = 1.0+t25+t27-t28-t29;
der2c.derp1.derrho <- t19*t33*t35*t39+t9*(-t21*t16*t17+t8*t52*t17-t31/t55*t16)*t61*t38-t64*t35*t34*t12*t39+t64*t61*t16+t19*t35*t24*t39-t80*t44*t39-t71*t81*t38+t71*t24*t37*t16+t19*t61*t87-t80*t81*t87+t71*t37*t26*t16+t10*t16*t17*t35*t37*t94-t80*t81*t94+t71*t37*t52;
 

 


			  t1 = 1.0-p2;
			  t2 = t1^teta
			  t3 = 1.0-p1;
			  t4 = t3^teta
			  t5 = t2*t4;
			  t6 = t2+t4-t5;
			  t8 = 1/teta-2.0;
			  t9 = t6^t8
			  t10 = t9*t8;
			  t11 = t2*teta;
			  t12 = 1/t1;
			  t14 = t12*t4;
			  t16 = -t11*t12+t11*t14;
			  t17 = 1/t6;
			  t19 = t10*t16*t17;
			  t20 = teta*teta;
			  t21 = 1/t20;
			  t22 = log(t6);
			  t24 = log(t1);
			  t25 = t2*t24;
			  t26 = log(t3);
			  t27 = t4*t26;
			  t28 = t25*t4;
			  t29 = t5*t26;
			  t31 = t8*(t25+t27-t28-t29);
			  t33 = -t21*t22+t31*t17;
			  t34 = teta-1.0;
			  t35 = t1^t34
			  t37 = t3^t34
			  t38 = teta-1.0+t2+t4-t5;
			  t39 = t37*t38;
			  t44 = t12*t24;
			  t46 = t2*t12;
			  t52 = -t11*t44-t46+t11*t44*t4+t46*t4+t11*t14*t26;
			  t55 = t6*t6;
			  t61 = t35*t37;
			  t64 = t9*t33;
			  t71 = t9*t35;
			  t80 = t71*t34;
			  t81 = t12*t37;
			  t87 = t26*t38;
			  t94 = 1.0+t25+t27-t28-t29;
der2c.derp2.derrho <- t19*t33*t35*t39+t9*(-t21*t16*t17+t8*t52*t17-t31/t55*t16)*t61*t38-t64*t35*t34*t12*t39+t64*t61*t16+t19*t35*t24*t39-t80*t44*t39-t71*t81*t38+t71*t24*t37*t16+t19*t61*t87-t80*t81*t87+t71*t37*t26*t16+t10*t16*t17*t35*t37*t94-t80*t81*t94+t71*t37*t52;
 
 


}



}










if( VC$BivD %in% c("C270","J270","G270") ) {

der2h.derp1p2   <- -der2h.derp1p2
der2h.derp1teta <- -der2h.derp1teta

der2c.derp1.derp2	= -der2c.derp1.derp2
der2c.derp1.derrho	= -der2c.derp1.derrho




}


if( VC$BivD %in% c("C90","J90","G90") ) {

der2h.derp1p1         <- -der2h.derp1p1
der2h.derp2p2         <- -der2h.derp2p2
der2h.derp1teta       <- -der2h.derp1teta
der2h.derteta.teta.st <- -der2h.derteta.teta.st 

der2c.derp1.derp2	= -der2c.derp1.derp2
der2c.derp2.derrho 	= -der2c.derp2.derrho


}

if( VC$BivD %in% c("C180","J180","G180") ) {

der2h.derp2p2              = -der2h.derp2p2
der2h.derteta.teta.st      = -der2h.derteta.teta.st
derteta.derteta.st         = -derteta.derteta.st
der2teta.derteta.stteta.st = -der2teta.derteta.stteta.st 
der2h.derp1p2              = -der2h.derp1p2 
der2h.derp1teta            = -der2h.derp1teta                                   
der2h.derp2teta            = -der2h.derp2teta 
der2h.derp1p1              = -der2h.derp1p1


}





ifef <- function(dv){

epsilon <- 0.0000001 
dv <- ifelse(is.na(dv), epsilon, dv ) 
dv <- ifelse(dv == Inf ,  8.218407e+20, dv )
dv <- ifelse(dv == -Inf ,  -8.218407e+20, dv )
dv 

}


# stuff below is for safety and it generally does nothing

der2h.derp2p2              =  ifef(der2h.derp2p2             )
der2h.derteta.teta.st      =  ifef(der2h.derteta.teta.st     )
derteta.derteta.st         =  ifef(derteta.derteta.st        )
der2teta.derteta.stteta.st =  ifef(der2teta.derteta.stteta.st)
der2h.derp1p2              =  ifef(der2h.derp1p2       )      
der2h.derp1teta            =  ifef(der2h.derp1teta     )      
der2h.derp2teta            =  ifef(der2h.derp2teta     )      
der2h.derp1p1              =  ifef(der2h.derp1p1       )      
der2c.derrho.derrho        =  ifef(der2c.derrho.derrho )      
der2c.derp1.derp1          =  ifef(der2c.derp1.derp1 )        
der2c.derp2.derp2          =  ifef(der2c.derp2.derp2 )        
der2c.derp1.derp2          =  ifef(der2c.derp1.derp2 )        
der2c.derp1.derrho         =  ifef(der2c.derp1.derrho)        
der2c.derp2.derrho         =  ifef(der2c.derp2.derrho)        







list(
der2h.derp2p2              = der2h.derp2p2, 
der2h.derteta.teta.st      = der2h.derteta.teta.st,  
derteta.derteta.st         = derteta.derteta.st, 
der2teta.derteta.stteta.st = der2teta.derteta.stteta.st,  
der2h.derp1p2              = der2h.derp1p2,  
der2h.derp1teta            = der2h.derp1teta,                                     
der2h.derp2teta            = der2h.derp2teta,  
der2h.derp1p1              = der2h.derp1p1,
der2c.derrho.derrho        = der2c.derrho.derrho,
der2c.derp1.derp1          = der2c.derp1.derp1, 
der2c.derp2.derp2          = der2c.derp2.derp2, 
der2c.derp1.derp2          = der2c.derp1.derp2, 
der2c.derp1.derrho         = der2c.derp1.derrho, 
der2c.derp2.derrho         = der2c.derp2.derrho)     




}




     























