sigmav.func <-
function(prob.d1.given.z1,prob.d1.given.z0,prob.z1){
probvector=c((1-prob.z1)*prob.d1.given.z0,(1-prob.z1)*(1-prob.d1.given.z0),prob.z1*prob.d1.given.z1,prob.z1*(1-prob.d1.given.z1));
value.vector=c(1-prob.d1.given.z0,-prob.d1.given.z0,1-prob.d1.given.z1,-prob.d1.given.z1);
expected.value=sum(probvector*value.vector);
variance=sum(probvector*value.vector^2)-expected.value^2;
list(sigmav=sqrt(variance))
}
