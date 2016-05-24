MCID <-
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


MCID=data.frame(MW=mean(dataghs$A_1, na.rm=TRUE),LW= mean(dataghs$B_1,na.rm=TRUE), 
NC= mean(dataghs$C_1,na.rm=TRUE ),LB=mean(dataghs$D_1,na.rm=TRUE ),MB=mean(dataghs$E_1,na.rm=TRUE))

MCID
     }
