mlica <-
function(prNCP, nruns=10, tol=0.0001, maxit=300, fail.th=5, learn.mu=1){

print("Entering mlica");
print("Performing preliminary run");
# Perform one run to learn ncp and initialise logL of compared
a <- mlicaMAIN(prNCP,tol=0.0001,maxit=10,mu=learn.mu);
ncp <- dim(a$S)[2];
max.logL <- a$LL;
a.best <- a;
print("Finished preliminary run");

# Performs runs 
# only use converged runs and find the one of maximum likelihood
print("Starting runs");
run.n <- 0 ;
fail.count <- 0;
v.logL <- vector();
v.NC <- vector();

while ( run.n < nruns ){

 a <- mlicaMAIN(prNCP,tol=0.0001,maxit=maxit,mu=learn.mu); 
 v.logL <- c(v.logL,a$LL);
 v.NC <- c(v.NC,a$NC);
 
 if ( a$NC == 0 ){ # if converged
  fail.count <- 0;
  run.n <- run.n + 1;

  if ( a$LL > max.logL ){
    a.best <- a;
  }
  
 }

 else { # did not converge
  fail.count <- fail.count + 1;
 }

 if( fail.count >= fail.th ){
   print("Stopping: Five consecutive runs failed to converge!");
   print("Consider either increasing the threshold for pca eigenvalues to perform ICA on a smaller subspace or increasing maxit");
   stop;
 }
 
} # matches run loop
print("End of runs");

return(a.best);

}
