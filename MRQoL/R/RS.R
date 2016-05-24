RS <-
function(score_1, score_2,X ) {
     
for (j in 1:length (X)){
    if (is.na(X[j])) (dataghs$A_3[j]= NA) else { 
     if (X[j] ==1) (dataghs$A_3[j]= score_1[j])
     if (X[j] ==2 | X[j] ==3 | X[j] ==4 | X[j] ==5) (dataghs$A_3[j]=NA) 
     } }
for (j in 1:length (X)){
    if (is.na(X[j])) (dataghs$B_3[j]= NA) else{ 
     if (X[j] ==2) (dataghs$B_3[j]= score_1[j])     
if (X[j] ==1 | X[j] ==3 | X[j] ==4 | X[j] ==5) (dataghs$B_3[j]=NA) 
     } }
for (j in 1:length (X)) {
    if (is.na(X[j])) (dataghs$C_3[j]= NA) else { 
     if (X[j] ==3) (dataghs$C_3[j]=score_1[j])
     if (X[j] ==1 | X[j] ==2 | X[j] ==4 | X[j] ==5) (dataghs$C_3[j]=NA) 
     } }
for (j in 1:length (X))  {
    if (is.na(X[j])) (dataghs$D_3[j]= NA) else  { 
     if (X[j] ==4) (dataghs$D_3[j]=score_1[j])
     if (X[j] ==1 | X[j] ==2 | X[j] ==3 | X[j] ==5) (dataghs$D_3[j]=NA) 
     }}
for (j in 1:length (X)){
    if (is.na(X[j])) (dataghs$E_3[j]= NA) else { 
     if (X[j] ==5) (dataghs$E_3[j]= score_1[j])
     if (X[j] ==2 | X[j] ==3 | X[j] ==4 | X[j] ==1) (dataghs$E_3[j]=NA) 
      }}



for (j in 1:length (X)){
    if (is.na(X[j])) (dataghs$A_4[j]= NA) else { 
     if (X[j] ==1) (dataghs$A_4[j]= score_2[j])
     if (X[j] ==2 | X[j] ==3 | X[j] ==4 | X[j] ==5) (dataghs$A_4[j]=NA) 
     } }
for (j in 1:length (X)){
    if (is.na(X[j])) (dataghs$B_4[j]= NA) else{ 
     if (X[j] ==2) (dataghs$B_4[j]= score_2[j])     
     if (X[j] ==1 | X[j] ==3 | X[j] ==4 | X[j] ==5) (dataghs$B_4[j]=NA) 
     } }
for (j in 1:length (X)) {
    if (is.na(X[j])) (dataghs$C_4[j]= NA) else { 
     if (X[j] ==3) (dataghs$C_4[j]=score_2[j])
     if (X[j] ==1 | X[j] ==2 | X[j] ==4 | X[j] ==5) (dataghs$C_4[j]=NA) 
     } }
for (j in 1:length (X))  {
    if (is.na(X[j])) (dataghs$D_4[j]= NA) else  { 
     if (X[j] ==4) (dataghs$D_4[j]=score_2[j])
     if (X[j] ==1 | X[j] ==2 | X[j] ==3 | X[j] ==5) (dataghs$D_4[j]=NA) 
     }}
for (j in 1:length (X)){
    if (is.na(X[j])) (dataghs$E_4[j]= NA) else { 
     if (X[j] ==5) (dataghs$E_4[j]= score_2[j])
     if (X[j] ==2 | X[j] ==3 | X[j] ==4 | X[j] ==1) (dataghs$E_4[j]=NA) 
      }}


     
dataghs$A_1=dataghs$A_3 - dataghs$A_4
dataghs$B_1=dataghs$B_3 - dataghs$B_4
dataghs$C_1=dataghs$C_3 - dataghs$C_4
dataghs$D_1=dataghs$D_3 - dataghs$D_4
dataghs$E_1=dataghs$E_3 - dataghs$E_4

a_1=length(na.omit(dataghs$A_1))
b_1=length(na.omit(dataghs$B_1))
c_1=length(na.omit(dataghs$C_1))
d_1=length(na.omit(dataghs$D_1))
e_1=length(na.omit(dataghs$E_1)) 

a_3=sd(dataghs$A_3,na.rm=TRUE)
b_3=sd(dataghs$B_3,na.rm=TRUE)
c_3=sd(dataghs$C_3,na.rm=TRUE)
d_3=sd(dataghs$D_3,na.rm=TRUE)
e_3=sd(dataghs$E_3,na.rm=TRUE) 

a_4=mean(dataghs$A_1,na.rm=TRUE)
b_4=mean(dataghs$B_1,na.rm=TRUE)
c_4=mean(dataghs$C_1,na.rm=TRUE)
d_4=mean(dataghs$D_1,na.rm=TRUE)
e_4=mean(dataghs$E_1,na.rm=TRUE) 

N=c((a_1+b_1+c_1+d_1+e_1), a_1, b_1, c_1, d_1, e_1)

RS_bis=(a_4*a_1+b_4*b_1+c_4*c_1+d_4*d_1+e_4*e_1)/(a_1+b_1+c_1+d_1+e_1)
RS=c(RS_bis, a_4, b_4, c_4, d_4, e_4)

ES=c(RS_bis/sd(score_1,na.rm=TRUE), a_4/a_3, b_4/b_3, c_4/c_3, d_4/d_3, e_4/e_3)

IC_1=c(t.test(dataghs$A_1, conf.level=0.9)$conf.int, t.test(dataghs$B_1, conf.level=0.9)$conf.int, 
       t.test(dataghs$C_1, conf.level=0.9)$conf.int,t.test(dataghs$D_1, conf.level=0.9)$conf.int,
       t.test(dataghs$E_1, conf.level=0.9)$conf.int)
      
IC=round(IC_1,0)

LCI=c(" ",IC[c(1,3,5,7,9)])
UCI=c(" ", IC[c(2,4,6,8,10)])

p_value=c(wilcox.test(score_1 , score_2)$p.value, wilcox.test(dataghs$A_3 , dataghs$A_4)$p.value, 
wilcox.test(dataghs$B_3 , dataghs$B_4)$p.value,wilcox.test(dataghs$C_3 , dataghs$C_4)$p.value, 
wilcox.test(dataghs$D_3 , dataghs$D_4)$p.value, wilcox.test(dataghs$E_3 , dataghs$E_4)$p.value )



data.frame(ID=c("Dimension", "MW", "LW", "NC", "LB", "MB"), N, RS, LCI, UCI, p_value, ES)
     
      }
