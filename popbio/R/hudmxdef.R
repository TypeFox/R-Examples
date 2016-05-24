hudmxdef<- function(vrs)
{
matrix(c(
vrs[14],      0,             vrs[15],                 vrs[16],                 vrs[17],               vrs[18], 
vrs[19],      0,             vrs[20],                 vrs[21],                 vrs[22],               vrs[23], 
      0, vrs[24], vrs[10]*(1-vrs[1]),      vrs[11]*(1-vrs[2]),      vrs[12]*(1-vrs[5]),                     0, 
      0,       0,     vrs[10]*vrs[1], vrs[11]*(vrs[2]-vrs[3]), vrs[12]*(vrs[5]-vrs[6]),      vrs[13]*(1-vrs[8]), 
      0,       0,                  0, vrs[11]*(vrs[3]-vrs[4]), vrs[12]*(vrs[6]-vrs[7]), vrs[13]*(vrs[8]-vrs[9]), 
      0,       0,                  0,          vrs[11]*vrs[4],          vrs[12]*vrs[7],         vrs[13]*vrs[9]
   ), nrow=6, byrow=TRUE)    
} 
