power.iv <-
function(n,lambda,gamma,var.z,sigmau,sigmav,rho,alpha=.05){
Lambda=(lambda^2)/((sigmau/sigmav)^2+2*rho*(sigmau/sigmav)*lambda+lambda^2);
ncp=gamma^2*(n*var.z)*Lambda/sigmav^2;
fquantile=qf(1-alpha,1,n-2)
power=1-pf(fquantile,1,n-2,ncp);
list(power=power);
}
