VarInt.coh <- 
function(XY){
X1=XY[1];  X2=XY[2]; X3=XY[3]; X4=XY[4];
Y1=XY[5];  Y2=XY[6]; Y3=XY[7]; Y4=XY[8];

res <- (X1*Y1)/(X1 + Y1)^3 + (X2*Y2)/(X2 + Y2)^3 + (X3*Y3)/(X3 + Y3)^3 + (X4*Y4)/(X4 + Y4)^3
return(res)
}

