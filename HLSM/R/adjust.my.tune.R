###Adjust my Tune####

#This function adjusts tuning parameters (1D and higher)#

adjust.my.tune=function(tune, acc, dim){
n=length(acc)
if(dim==1){
for(i in 1:n){
if(acc[i]>0.9){tune[i]=tune[i]*3}
if(acc[i]>0.8){tune[i]=tune[i]*2}
if(acc[i]>0.7){tune[i]=tune[i]*1.5}
if(acc[i]>0.6){tune[i]=tune[i]*1.25}
if(acc[i]>0.5){tune[i]=tune[i]*1.1}
if(acc[i]>0.48){tune[i]=tune[i]*1.05}
if(acc[i]<0.4){tune[i]=tune[i]*.95}
if(acc[i]<0.35){tune[i]=tune[i]*0.9}
if(acc[i]<0.25){tune[i]=tune[i]*0.75}
if(acc[i]<0.1){tune[i]=tune[i]*0.5}
}
}

if(dim>1){
for(i in 1:2){
  if(acc[i]>0.9){tune[i]=tune[i]*3}
if(acc[i]>0.8){tune[i]=tune[i]*2}
if(acc[i]>0.7){tune[i]=tune[i]*1.5}
if(acc[i]>0.6){tune[i]=tune[i]*1.25}
if(acc[i]>0.5){tune[i]=tune[i]*1.1}
if(acc[i]>0.48){tune[i]=tune[i]*1.05}
if(acc[i]<0.4){tune[i]=tune[i]*.95}
if(acc[i]<0.35){tune[i]=tune[i]*0.9}
if(acc[i]<0.25){tune[i]=tune[i]*0.75}
if(acc[i]<0.1){tune[i]=tune[i]*0.5}
}

for(i in 3:n){
if(acc[i]>0.9){tune[i]=tune[i]*3}
if(acc[i]>0.8){tune[i]=tune[i]*2.5}
if(acc[i]>0.7){tune[i]=tune[i]*2}
if(acc[i]>0.55){tune[i]=tune[i]*1.5}
if(acc[i]>0.4){tune[i]=tune[i]*1.25}
if(acc[i]>0.3){tune[i]=tune[i]*1.1}
if(acc[i]>0.26){tune[i]=tune[i]*1.05}
if(acc[i]<0.2){tune[i]=tune[i]*.9}
if(acc[i]<0.18){tune[i]=tune[i]*0.85}
if(acc[i]<0.15){tune[i]=tune[i]*0.8}
if(acc[i]<0.05){tune[i]=tune[i]*0.75}
}
}
return(tune)
}


