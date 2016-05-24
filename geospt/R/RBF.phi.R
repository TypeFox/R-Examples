assign("RBF.phi",
 function(distance, eta, func){
switch(func,                                      
ST = ifelse(distance>0,log(eta*distance/2)+ besselK(distance*eta,0) + 0.5772161, 0),                          
CRS = ifelse(distance>0,log((eta*distance/2)^2) + ifelse(expint_E1((eta*distance/2)^2)=="NaN",0,expint_E1((eta*distance/2)^2)) + 0.5772161, 0),  
M = (sqrt(distance^2+ eta^2))^1,
IM = (sqrt(distance^2+ eta^2))^(-1),
TPS = ifelse(distance>0,((distance*eta)^2)*log(distance*eta),0),
GAU = exp(-eta*(distance^2)),
EXPON = exp(-(distance*eta)),
TRI = sin(distance*eta)
)
}
)