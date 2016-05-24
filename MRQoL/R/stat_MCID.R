stat_MCID <-
function(score_1, score_2, X ) {
     
for (j in 1:length (X)){
    if (is.na(X[j])) (dataghs$A_1[j]= NA) else { 
     if (X[j] ==1) (dataghs$A_1[j]= score_1[j]-score_2[j])
     if (X[j] ==2 | X[j] ==3 | X[j] ==4 | X[j] ==5) (dataghs$A_1[j]=NA) 
     } }
for (j in 1:length (X)){
    if (is.na(X[j])) (dataghs$B_1[j]= NA) else{ 
     if (X[j] ==2) (dataghs$B_1[j]= score_1[j]-score_2[j])
     if (X[j] ==1 | X[j] ==3 | X[j] ==4 | X[j] ==5) (dataghs$B_1[j]=NA) 
     } }
for (j in 1:length (X)) {
    if (is.na(X[j])) (dataghs$C_1[j]= NA) else { 
     if (X[j] ==3) (dataghs$C_1[j]=score_1[j]-score_2[j])
     if (X[j] ==1 | X[j] ==2 | X[j] ==4 | X[j] ==5) (dataghs$C_1[j]=NA) 
     } }
for (j in 1:length (X))  {
    if (is.na(X[j])) (dataghs$D_1[j]= NA) else  { 
     if (X[j] ==4) (dataghs$D_1[j]=score_1[j]-score_2[j])
     if (X[j] ==1 | X[j] ==2 | X[j] ==3 | X[j] ==5) (dataghs$D_1[j]=NA) 
     }}
for (j in 1:length (X)){
    if (is.na(X[j])) (dataghs$E_1[j]= NA) else { 
     if (X[j] ==5) (dataghs$E_1[j]= score_1[j]-score_2[j])
     if (X[j] ==2 | X[j] ==3 | X[j] ==4 | X[j] ==1) (dataghs$E_1[j]=NA) 
      }}

a_1=length(na.omit(dataghs$A_1))
b_1=length(na.omit(dataghs$B_1))
c_1=length(na.omit(dataghs$C_1))
d_1=length(na.omit(dataghs$D_1))
e_1=length(na.omit(dataghs$E_1))


a_3=sd(dataghs$A_1,na.rm=TRUE)
b_3=sd(dataghs$B_1,na.rm=TRUE)
c_3=sd(dataghs$C_1,na.rm=TRUE)
d_3=sd(dataghs$D_1,na.rm=TRUE)
e_3=sd(dataghs$E_1,na.rm=TRUE) 

SD_dim=(a_3*a_1+b_3*b_1+c_3*c_1+d_3*d_1+e_3*e_1)/(a_1+b_1+c_1+d_1+e_1)

SD=c(SD_dim,a_3, b_3, c_3, d_3, e_3)

IC_1=c(t.test(dataghs$A_1, conf.level=0.9)$conf.int, t.test(dataghs$B_1, conf.level=0.9)$conf.int, 
       t.test(dataghs$C_1, conf.level=0.9)$conf.int,t.test(dataghs$D_1, conf.level=0.9)$conf.int,
       t.test(dataghs$E_1, conf.level=0.9)$conf.int)
      
IC=round(IC_1,0)

LCI=c(" ",IC[c(1,3,5,7,9)])
UCI=c(" ", IC[c(2,4,6,8,10)])


data.frame(ID=c("Dimension","MW", "LW", "NC", "LB", "MB"), N=c((a_1+b_1+c_1+d_1+e_1), a_1, b_1, c_1, d_1, e_1), 
SD, LCI, UCI) 
     
         }
