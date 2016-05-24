DMCEE <-
function(x,y,b){
n=length(y)
Y=circ.reg(x,y)$fi
MCEE=(abs(MCE(y, Y, n) - MCe(cos(y - Y))))
k=A1inv(sum(cos(y-Y))/n)
if(b==0.90) cc=2*k-1 else cc=2*k
DMCE=DMCE[n,cc]
k = MCEE > DMCE
plot(MCEE, type = "l", ylab= "" , main = paste("Outliers:", 
            paste((1:n)[k], collapse = ",")))
abline(h=DMCE,col=6)
}
