ABCEE <-
function(X, Y, U, omega, forX = NA, niter = 5000, nburn = 500, nthin = 10, maxmodelY = NA, OR = 20, family.X = "gaussian")
{
n = length(Y); #Sample size
n_cov = ncol(as.matrix(U)); #Number of covariate, potential confounders
priorX = rep(0, n_cov);
if(is.na(forX[1])) forX = 1:n_cov;
priorX[forX] = 0.5;
if(is.na(maxmodelY)) maxmodelY = min(niter + nburn, 2**n_cov);
model.X = bic.glm(y = X, x = U, glm.family = family.X, OR = OR, prior.param = priorX);
models.X = cbind(model.X$which, model.X$postprob);
alpha_X = model.X$probne0/100;
alpha_Y = as.numeric(bic.glm(y = Y, x = cbind(X, U), glm.family = "gaussian", OR = 1.0000001, 
prior.param = c(1, rep(0.5, n_cov)))$which[1, 2:(n_cov+1)]); #Initial outcome model
if(sum(alpha_Y) == 0)
{
model_Y0 = lm(Y ~ X);
}
else
{
model_Y0 = lm(Y~ X + U[,alpha_Y == 1]);
}
bic_y = BIC(model_Y0);
betas = numeric(ceiling(niter/nthin)); #objects that will contain results
models_Y = matrix(0, nrow = ceiling(niter/nthin), ncol = n_cov); #objects that will contain results
tested_models_Y = matrix(-1, nrow = maxmodelY, ncol = 1); #objects that will contain information for models already tested
pY_tested = matrix(NA, nrow = maxmodelY, ncol = 1); #objects that will contain information for models already tested
pY_tested1 = matrix(NA, nrow = maxmodelY, ncol = n_cov);
pY_tested2 = matrix(NA, nrow = maxmodelY, ncol = n_cov);
j = 1; #Index to be filled next
k_Y = 1; #Index of pY_tested to be filled next
tested_models_Y[k_Y] = paste0(alpha_Y, collapse = ""); #Recording info about candidate model
sy = sd(Y);
su = apply(U, 2, sd);
if(sum(alpha_Y != 0))
{
a = coef(model_Y0)[-(1:2)];
b = sqrt(diag(vcov(model_Y0))[-(1:2)]);
}
else
{
a = NA;
b = NA;
}
pY_tested1[k_Y, 1:sum(alpha_Y)] = a;
pY_tested2[k_Y, 1:sum(alpha_Y)] = b;
pY_tested[k_Y] = bic_y;
k_Y = k_Y + 1; #Changing index
for(i in 1:(nburn+niter))
{
alpha_Y_0 = alpha_Y; #Actual alpha_Y
alpha_Y_1 = alpha_Y; temp = sample(n_cov, 1); alpha_Y_1[temp] = (alpha_Y_1[temp] + 1)%%2; #Candidate alpha_Y
bic_y_0 = bic_y;
a_0 = a;
b_0 = b;

char_alpha_Y_1 = paste0(alpha_Y_1, collapse = "");
tested = which(tested_models_Y == char_alpha_Y_1);
if(length(tested)) #The candidate model was already tested
{
a_1 = pY_tested1[tested,];
b_1 = pY_tested2[tested,];
bic_y_1 = pY_tested[tested];
}
else #The candidate model was never tested
{
if(sum(alpha_Y_1) == 0)
{
model_Y1 = lm(Y ~ X);
a_1 = NA;
b_1 = NA;
}
else
{
model_Y1 = lm(Y ~ X + U[, alpha_Y_1 == 1]);
a_1 = coef(model_Y1)[-(1:2)];
b_1 = sqrt(diag(vcov(model_Y1))[-(1:2)]);
}
pY_tested1[k_Y, 1:sum(alpha_Y_1) ] = a_1;
pY_tested2[k_Y, 1:sum(alpha_Y_1) ] = b_1;
bic_y_1 = BIC(model_Y1);
tested_models_Y[k_Y] = char_alpha_Y_1; #Recording info about candidate model
pY_tested[k_Y] = bic_y_1; #Recording info about candidate model
k_Y = k_Y + 1; #Changing index
}
change_Y = alpha_Y_1 - alpha_Y_0; 
px = alpha_X[change_Y != 0];
add_Y = sum(change_Y); #Is it proposed to add or remove a variable
if(add_Y == 1)
{
delta = omega*(rnorm(1, a_1[change_Y[alpha_Y_1 != 0] != 0], sd = b_1[change_Y[alpha_Y_1 != 0] != 0])*su[change_Y != 0]/sy)**2;
}
else
{
delta = omega*(rnorm(1, a_0[change_Y[alpha_Y_0 != 0] != 0],sd = b_0[change_Y[alpha_Y_0 != 0] != 0])*su[change_Y != 0]/sy)**2;
}

if(add_Y == 1)
{
ratio_P_alpha = (px*(delta/(1 + delta)) + (1 - px)*0.5)/(px*(1/(1 + delta)) + (1 - px)*0.5);
}
if(add_Y == -1)
{
ratio_P_alpha = (px*(1/(1 + delta)) + (1 - px)*0.5)/(px*(delta/(1 + delta)) + (1 - px)*0.5);
}
if(ratio_P_alpha == 0)
{
ratio = 0;
}
else if(ratio_P_alpha == Inf)
{
ratio = 1;
}
else
{
ratio = exp(-bic_y_1/2 + bic_y_0/2)*ratio_P_alpha;
}
if(sample(c(1,0), size = 1, prob = c(min(1, ratio), 1 - min(1, ratio))))
{
alpha_Y = alpha_Y_1;
bic_y = bic_y_1;
a = a_1;
b = b_1;
}
if(i > nburn & (i - nburn)%%nthin == 1)
{
if(sum(alpha_Y) != 0)
{
model_b = lm(Y ~ X + U[,alpha_Y == 1]);
}
else
{
model_b = lm(Y ~ X);
}
beta_alphaY = rnorm(1, mean = model_b$coefficients[2], sd = sqrt(vcov(model_b)[2,2]));
betas[j] = beta_alphaY;
models_Y[j, ] = alpha_Y;
j = j + 1;
}
}
return(list(betas = betas, models.X = models.X, models.Y = models_Y));
}
