"UPMEpiktildefrompik" <-function(pik,eps=1e-6) 
{ 
n=sum(pik) 
n=.as_int(n)
pikt=pik 
arr=1 
while(arr>eps) 
{ 
w=(pikt)/(1-pikt) 
q=UPMEqfromw(w,n)
pikt1=pikt+pik-UPMEpikfromq(q) 
arr=sum(abs(pikt-pikt1)) 
pikt=pikt1 
} 
pikt
}

