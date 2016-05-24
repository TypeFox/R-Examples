
# Start things up.
# print("Demo for BayesBridge package.");
# readline("Press <Return> to continue...");

# Load the library.
# library("BayesBridge", lib.loc="Test/");

# Load the data.
data(diabetes, package="BayesBridge");
cov.name = colnames(diabetes$x);
y = diabetes$y;
X = diabetes$x;
id = colnames(X);

# Center.
y = y - mean(y);
mX = colMeans(X);
for(i in 1:442){ X[i,] = X[i,] - mX; }

# The least squares solution.
LS = solve(t(X) %*% X, t(X) %*% y);

# Set parameters.
alpha = 0.5;
sig2  = 2500;
tau   = sqrt(sig2)*1e-5;
N = 20000;

# Print parameters
print(paste("alpha, sig2, tau are set to", alpha, sig2, tau));

# Expectation maximization
beta.EM = bridge.EM(y, X, alpha, tau/sqrt(sig2), 10e8, 10e-9, 30, use.cg=FALSE);

print("Bridge Regression when only alpha is known!");

# Gibbs with R routine.
print(paste("Bridge Regression using R routine. ", N, "samples.  This will take a while."));
readline("Press <Return> to continue...");

g.R = bridge.reg.R(y, X, N, alpha, 0.5, 0.5, 500);

# Gibbs with C routine.
print(paste("Bridge Regression using C routine. ", N, "samples.  Same number of samples."));
readline("Press <Return> to continue...");
# gb1 = bridge.reg.know.sig2(y, X, 10000, alpha, sig2, tau, 500);
# gb2 = bridge.reg.know.tau(y, X, 10000,  alpha, tau, 0.0, 0.0, 500);
g.C = bridge.reg(y, X, N, alpha, 0.0, 0.0, 0.5, 0.5, 500);

readline("Histograms: Press <Return> to continue.");

par(mfrow=c(1,2));
par(mai=c(1.02, 0.82, 0.82, 0.42));
for(i in 1:10){
    # Marginal Posterior.
    hist(g.R$beta[,i], breaks=100, prob=TRUE, ylab="dens", xlab=id[i],
         main=paste(id[i], "using R"));
    hist(g.C$beta[i,], breaks=100, prob=TRUE, ylab="dens", xlab=id[i],
         main=paste(id[i], "using C"));
    # Summary statistics.
    print("Mean and SD using R and C.");
    print(paste(mean(g.R$beta[,i]), mean(g.C$beta[i,])));
    print(paste(sd(g.R$beta[,i]), sd(g.C$beta[i,])));
    # Pause.
    readline("Press <Return> to continue.");
}

# Now do this with lots of samples.
N = 500000;
print(paste("Bridge Regression using C routine. ", N, "samples.  (Lots more samples.)"));
readline("Press <Return> to continue...");
g.C = bridge.reg(y, X, N, alpha, 0.0, 0.0, 0.5, 0.5, 500);

readline("Histograms: Press <Return> to continue.");
par(mfrow=c(5,2));
par(mai=c(0.4,0.4,0.1,0.1));
for(i in 1:10){
    # Marginal Posterior.
    hist(g.C$beta[i,], breaks=100, prob=TRUE, ylab="dens", xlab=id[i],
         main=id[i]);
    # Joint Modes.
    abline(v=LS[i], col=2);
    abline(v=beta.EM[i], col=3);
    legend("topright", legend=c("LS", "EM"), col=c(2,3), lty=c(1,1));
}

# How many solves do we save using conjugate gradient method?
# Expectation maximization
readline("Expectation Maximization: Press <Return> to continue.");

beta.EM.f = bridge.EM(y, X, alpha, tau/sqrt(sig2),
                      10e8, 10e-9, 30, use.cg=FALSE, ret.solves=TRUE);
beta.EM.t = bridge.EM(y, X, alpha, tau/sqrt(sig2),
                      10e8, 10e-9, 30, use.cg=TRUE, ret.solves=TRUE);
beta.EM.R = bridge.EM.R(y, X, alpha, tau/sqrt(sig2));

print("Expectation Maximization using CG:");
print(paste("Num solves:", beta.EM.t$num.solves));
print("beta:")
print(beta.EM.t$beta);

print("Expectation Maximization not using CG:");
print(paste("Num solves:", beta.EM.f$num.solves));
print("beta:")
print(beta.EM.f$beta);

print("Simple Expectation Maximization using R code:");
print("beta:")
print(beta.EM.R);

# Go back to a sane parameter.
par(mai=c(1.02, 0.82, 0.82, 0.42));
