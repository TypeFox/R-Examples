Pvero_ratio <-
function (x,n,rho) {
ene=sum(n);a1=sum(x)
b=ene-(n[2]-x[2])+(ene-(n[1]-x[1]))*rho
raizp=b^2-4*ene*a1*rho
pvero=(b-sqrt(raizp))/(2*ene*rho)
pvero[pvero>1]=1;pvero[pvero<0]=0
pvero}
