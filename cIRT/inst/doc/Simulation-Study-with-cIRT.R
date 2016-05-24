## ----sim_setup, eval = F-------------------------------------------------
#  ### Variables
#  # Y = trial matix
#  # C = KN vector of binary choices
#  # N = #of subjects
#  # J = # of items
#  # K = # of choices
#  # atrue = true item discriminations
#  # btrue = true item locations
#  # thetatrue = true thetas/latent performance
#  # gamma = fixed effects coefficients
#  # Sig = random-effects variance-covariance
#  # subid = id variable for subjects
#  
#  # install mvtnorm is necessary
#  # install.packages("mvtnorm")
#  
#  # Load the Library
#  library(cIRT)
#  
#  # Simulate 2PNO Data
#  N = 1000  # Students
#  J = 20    # Total numbers of possible items per SA
#  
#  # Set a seed for the random generation of a's and b's.
#  set.seed(1337)
#  
#  # Randomly pick a's and b's
#  
#  # Generate as, bs
#  atrue=runif(J)+1
#  btrue=2*runif(J)-1
#  
#  save(atrue, btrue, file="a_b_true.rda")
#  
#  
#  # 2 Level Probit Data
#  K = 30
#  
#  gam_notheta = c(.5,1)
#  gam_theta   = c(3,.25)
#  gamma = c(gam_notheta,gam_theta)
#  
#  Sig = matrix(c(.25,0,0,.125),2,2)
#  
#  # Number of replications
#  B = 100
#  
#  # Begin the Simulation Study
#  for(b in 1:B){
#  
#    # Provide user with state information
#    cat(paste0("On Iteration:",b,"\n"))
#    set.seed(b + 1234)
#  
#    # True theta and etay
#    thetatrue = rnorm(N)
#    etay = outer(rep(1,N),atrue) * thetatrue - outer(rep(1,N),btrue)
#  
#    # Generate Y for 2PNO model
#    p.correct = pnorm(etay)
#    Y = matrix(rbinom(N*J, 1, p.correct),N,J)
#  
#    #################################################
#    # Simulating 2 level probit data
#    #################################################
#  
#    subid = expand.grid(cid = 1:K,sid = 1:N)[,2]
#  
#    pred = rnorm(K*N,0,1) # Pred
#  
#    center_pred = center_matrix(as.matrix(pred))
#  
#    Xnotheta = cbind(1,center_pred)
#  
#    Xtheta = rep(thetatrue,each=K)*Xnotheta
#    X = cbind(Xnotheta,Xtheta)
#  
#  
#    zetas = mvtnorm::rmvnorm(N,mean=c(0,0),sigma=Sig) # mvtnorm environment accessed
#    W_veczeta = apply(Xnotheta*zetas[rep(1:N,each=K),],1,sum)
#  
#    etac = X%*%gamma + W_veczeta
#    Zc = rnorm(N*K,mean=etac,sd=1)
#    C = 1*(Zc>0)
#  
#    # Run the Choice Item Response Model
#    out1 = cIRT(subid,
#                Xnotheta,
#                c(1,2),
#                Xnotheta,
#                Y,
#                C,
#                5000)
#  
#    mname = paste0("model_",b)
#  
#    # Assign the data set name
#    assign(mname, out1)
#  
#    # Save the out object
#    save(list=mname, file=paste0(mname,".rda"))
#  
#    # Clean up Export
#    rm(list = c(mname,"out1"))
#  }

## ----sim_results, eval = F-----------------------------------------------
#  # E[theta] - theta
#  bias = function(theta.star, theta.true){
#    matrix(mean(theta.star) - theta.true,ncol=1)
#  }
#  
#  bias2 = function(theta.star, theta.true){
#    matrix(theta.star - theta.true,ncol=1)
#  }
#  
#  # sqrt ( 1/n * sum( (y_i - y.hat_i)^2 )
#  RMSE = function(y,y.hat){
#    sqrt(  (y-y.hat)^2 )
#  }
#  
#  # Change true values if needed
#  # True Values
#  gam_notheta = c(.5,1)
#  gam_theta   = c(3,.25)
#  gamma = c(gam_notheta,gam_theta)
#  # Loads a and b values
#  load("a_b_true.rda")
#  
#  Sig = as.numeric(matrix(c(.25,0,0,.125),2,2))[c(1,2,4)]
#  B = 100
#  
#  # Storage to hold bootstrap replications
#  a_result = matrix(0, B, 20)
#  
#  b_result = matrix(0, B, 20)
#  
#  gs0_result = matrix(0, B, 2)
#  
#  beta_result = matrix(0, B, 2)
#  
#  sig_result = array(NA, dim=c(2,2,B))
#  
#  
#  for(b in 1:B){
#  
#    mname = paste0("model_",b)
#  
#    load(paste0(mname, ".rda"))
#  
#    d = get(mname)
#  
#    a_result[i,] = apply(d$as, 2, FUN = mean)
#    b_result[i,] = apply(d$bs, 2, FUN = mean)
#  
#    gs0_result[i,] = apply(d$gs0, 2, FUN = mean)
#  
#    beta_result[i,] = apply(d$betas, 2, FUN = mean)
#  
#    sig_result[,,i] = solve(apply(d$Sigma_zeta_inv, c(1,2), FUN = mean))
#  }
#  
#  # Obtain an overall mean for each of the following:
#  
#  m_a_result = apply(a_result, 2, FUN = mean)
#  
#  m_b_result = apply(b_result, 2, FUN = mean)
#  
#  m_gs0_result = apply(gs0_result, 2, FUN = mean)
#  
#  m_beta_result = apply(beta_result, 2, FUN = mean)
#  
#  m_sig_result = as.numeric(apply(sig_result, c(1,2), FUN = mean))[c(1,2,4)]
#  
#  
#  # Perform a bias evaluation given the true values:
#  
#  a_bias = bias2(m_a_result,atrue)
#  
#  b_bias = bias2(m_b_result,btrue)
#  
#  gs0_bias = bias2(m_gs0_result,gamma[1:2])
#  
#  beta_bias = bias2(m_beta_result,gamma[3:4])
#  
#  sig_bias = bias2(m_sig_result, Sig)
#  
#  # Perform the RMSE under the supplied results:
#  
#  a_RMSE = RMSE(m_a_result,atrue)
#  
#  b_RMSE = RMSE(m_b_result,btrue)
#  
#  gs0_RMSE = RMSE(m_gs0_result,gamma[1:2])
#  
#  beta_RMSE = RMSE(m_beta_result,gamma[3:4])
#  
#  sig_RMSE = RMSE(m_sig_result, Sig)

## ----out, eval = F-------------------------------------------------------
#  # Make a results export
#  results_a = cbind(m_a_result, atrue, a_bias, a_RMSE)
#  
#  rownames(results_a) = paste0("Item ",1:length(atrue))
#  
#  colnames(results_a) = c("$a$ Estimate", "$a$ True", "$a$ Bias", "$a$ RMSE")
#  
#  results_b = cbind(m_b_result,  btrue, b_bias, b_RMSE)
#  
#  rownames(results_b) = paste0("Item ",1:length(btrue))
#  colnames(results_b) = c("$b$ Estimate", "$b$ True", "$b$ Bias", "$b$ RMSE")
#  
#  results_gs0 = cbind(m_gs0_result,gamma[1:2],gs0_bias,gs0_RMSE)
#  rownames(results_gs0) = paste0("$\\gamma_",1:2,"$")
#  colnames(results_gs0) = c("$\\beta$ Estimate", "$\\beta$ True", "$\\beta$ Bias", "$\\beta$ RMSE")
#  
#  results_beta = cbind(m_beta_result,gamma[3:4],beta_bias,beta_RMSE)
#  rownames(results_beta) = paste0("$\\beta_",1:2,"$")
#  colnames(results_beta) = c("$\\beta$ Estimate", "$\\beta$ True", "$\\beta$ Bias", "$\\beta$ RMSE")
#  
#  results_sig_zeta = cbind(m_sig_result,Sig,sig_bias,sig_RMSE)
#  rownames(results_sig_zeta) = c("$\\Zeta_{1,1}$","$\\Zeta_{2,1} = $\\Zeta_{1,2}$", "$\\Zeta_{2,2}$")
#  colnames(results_sig_zeta) = c("$\\Sigma_{\\Zeta}$ Estimate", "$\\beta$ True", "$\\beta$ Bias", "$\\beta$ RMSE")
#  
#  save(results_a, results_b, results_gs0, results_beta, results_sig_zeta, file="sim_results.rda")

