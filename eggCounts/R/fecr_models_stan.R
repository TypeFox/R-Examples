# set default values for the priors 
fecr_setPrior <- function(muPrior, kappaPrior, deltaPrior, phiPrior){
  if(missing(muPrior)) muPrior = list(priorDist = "gamma",hyperpars=c(1,0.001))
  if(is.null(muPrior[["priorDist", exact = TRUE]])) muPrior$priorDist = "gamma"
  if(is.null(muPrior[["hyperpars", exact = TRUE]])) muPrior$hyperpars = c(1,0.001)
  
  if(missing(kappaPrior)) kappaPrior = list(priorDist = "gamma",hyperpars=c(1,0.7))
  if(is.null(kappaPrior[["priorDist", exact = TRUE]])) kappaPrior$priorDist = "gamma"
  if(is.null(kappaPrior[["hyperpars", exact = TRUE]])) kappaPrior$hyperpars = c(1,0.7)
  
  if(missing(deltaPrior)) deltaPrior <- list(priorDist="beta", hyperpars=c(1,1))
  if(is.null(deltaPrior[["priorDist", exact = TRUE]])) deltaPrior$priorDist = "beta"
  if(is.null(deltaPrior[["hyperpars", exact = TRUE]])) deltaPrior$hyperpars = c(1,1)
  
  if(missing(phiPrior)) phiPrior = list(priorDist = "beta",hyperpars=c(1,1))
  if(is.null(phiPrior[["priorDist", exact = TRUE]])) phiPrior$priorDist = "beta"
  if(is.null(phiPrior[["hyperpars", exact = TRUE]])) phiPrior$hyperpars = c(1,1)
  
  return(list(mu=muPrior,kappa=kappaPrior,delta=deltaPrior, phi = phiPrior))
}

# Stan model code for paired model without zero inflation 
paired_stan <- function(priors){
  #hyperparameters for pre-treatment mean mu
  a.mu <- priors$mu$hyperpars[1]
  b.mu <- priors$mu$hyperpars[2]
  dist.mu <- priors$mu$priorDist
  #hyperparameters  for overdispersion parameter kappa
  a.kappa <- priors$kappa$hyperpars[1]
  b.kappa <- priors$kappa$hyperpars[2]
  dist.kappa <- priors$kappa$priorDist
  #hyperparameters  for change in mean delta
  a.delta <- priors$delta$hyperpars[1]
  b.delta <- priors$delta$hyperpars[2]
  dist.delta <- priors$delta$priorDist
  paste0('data {
           int J; // number of animals
           int ystararaw[J]; // after treatment McMaster count
           int ystarbraw[J]; // before treatment McMaster count
           int fpre[J];
           int fpost[J];
           }
           parameters {
           real<lower=0> kappa;
           real<lower=0> mu;
           real<lower=0,upper=1> delta;
           real<lower=0> mub[J];
           }
           transformed parameters{
           real lambdaa[J];
           real lambdab[J];
           real kappamu;
           for (i in 1:J){
           lambdab[i] <- mub[i]/fpre[i];
           lambdaa[i] <- delta*mub[i]/fpost[i];
           }
           kappamu <- kappa/mu;
           }
           model {
           mu ~ ',dist.mu,'(',a.mu,',',b.mu,');            // prior
           kappa ~ ',dist.kappa,'(',a.kappa,',',b.kappa,');
           delta ~ ',dist.delta,'(',a.delta,',',b.delta,');
           mub ~ gamma(kappa,kappamu);           // likelihoods
           ystarbraw ~ poisson(lambdab);
           ystararaw ~ poisson(lambdaa);
           }')
}

# Stan model code for unpaired model without zero inflation  --------------

unpaired_stan <- function(priors){
  #hyperparameters for pre-treatment mean mu
  a.mu <- priors$mu$hyperpars[1]
  b.mu <- priors$mu$hyperpars[2]
  dist.mu <- priors$mu$priorDist
  #hyperparameters  for overdispersion parameter kappa
  a.kappa <- priors$kappa$hyperpars[1]
  b.kappa <- priors$kappa$hyperpars[2]
  dist.kappa <- priors$kappa$priorDist
  #hyperparameters  for change in mean delta
  a.delta <- priors$delta$hyperpars[1]
  b.delta <- priors$delta$hyperpars[2]
  dist.delta <- priors$delta$priorDist
  paste0('data {
           int Ja; // number of animals
           int Jb;
           int ystararaw[Ja]; // after treatment McMaster count
           int ystarbraw[Jb]; // before treatment McMaster count
           int fpre[Ja];
           int fpost[Jb];
         }
           parameters {
             real<lower=0> kappa;
             real<lower=0> mu;
             real<lower=0,upper=1> delta;
             real<lower=0> mub[Jb]; # true epg before treatment
             real<lower=0> mua[Ja]; # true epg after treatment
           }
           transformed parameters{
             real lambdaa[Ja];
             real lambdab[Jb];
             real kappamu;
             for (i in 1:Jb){
             lambdab[i] <- mub[i]/fpre[i];
             }
             for (i in 1:Ja){
             lambdaa[i] <- delta*mua[i]/fpost[i];
             }
             kappamu <- kappa/mu;
           }
           model {
             // prior
           mu ~ ',dist.mu,'(',a.mu,',',b.mu,'); 
           kappa ~ ',dist.kappa,'(',a.kappa,',',b.kappa,');
           delta ~ ',dist.delta,'(',a.delta,',',b.delta,');
           increment_log_prob(gamma_log(mub,kappa,kappamu)+gamma_log(mua,kappa,kappamu));             // likelihoods
           ystarbraw ~ poisson(lambdab);
           ystararaw ~ poisson(lambdaa);
           }')
}

# Stan model code for unpaired model with zero inflation  ---------------

ZI_unpaired_stan <- function(priors){
  #hyperparameters for pre-treatment mean mu
  a.mu <- priors$mu$hyperpars[1]
  b.mu <- priors$mu$hyperpars[2]
  dist.mu <- priors$mu$priorDist
  #hyperparameters  for overdispersion parameter kappa
  a.kappa <- priors$kappa$hyperpars[1]
  b.kappa <- priors$kappa$hyperpars[2]
  dist.kappa <- priors$kappa$priorDist
  #hyperparameters  for change in mean delta
  a.delta <- priors$delta$hyperpars[1]
  b.delta <- priors$delta$hyperpars[2]
  dist.delta <- priors$delta$priorDist
  #hyperparameters  for prevalence 
  a.phi <- priors$phi$hyperpars[1]
  b.phi <- priors$phi$hyperpars[2]
  dist.phi <- priors$phi$priorDist
  paste0('data {
           int Ja; // number of animals
           int Jb;
           int ystararaw[Ja]; // after treatment McMaster count
           int ystarbraw[Jb]; // before treatment McMaster count
           int fpost[Ja];
           int fpre[Jb];
        }
        parameters {
           real<lower=0> kappa;
           real<lower=0> mu;
           real<lower=0,upper=1> delta;
           real<lower=0> mua[Ja];
           real<lower=0> mub[Jb];
           real<lower=0,upper=1> phi;
        }
        transformed parameters{
          real lambdaa[Ja];
          real lambdab[Jb];
          real kappamu;
          for (i in 1:Jb){
            lambdab[i] <- mub[i]/fpre[i];
          }
          for (i in 1:Ja){
            lambdaa[i] <- delta*mua[i]/fpost[i];
          }
          kappamu <- kappa/mu;
        }
        model {
          // prior
          mu ~ ',dist.mu,'(',a.mu,',',b.mu,'); 
          kappa ~ ',dist.kappa,'(',a.kappa,',',b.kappa,');
          delta ~ ',dist.delta,'(',a.delta,',',b.delta,');
          phi ~ ',dist.phi,'(',a.phi,',',b.phi,');
          // likelihoods
          ystarbraw ~ poisson(lambdab); 
          ystararaw ~ poisson(lambdaa); 
          for (n in 1:Jb) {
            if (mub[n] == 0)
                 increment_log_prob(bernoulli_log(1,phi));
                 else
                 increment_log_prob(bernoulli_log(0,phi) + gamma_log(mub[n],kappa,kappamu));
           }
          for (n in 1:Ja) {
             if (mua[n] == 0)
                increment_log_prob(bernoulli_log(1,phi));
          else 
                increment_log_prob(bernoulli_log(0,phi) + gamma_log(mua[n],kappa,kappamu));
          }
}')
}

# Stan model code for paired model with zero inflation  -------------------

ZI_paired_stan <- function(priors){
  #hyperparameters for pre-treatment mean mu
  a.mu <- priors$mu$hyperpars[1]
  b.mu <- priors$mu$hyperpars[2]
  dist.mu <- priors$mu$priorDist
  #hyperparameters  for overdispersion parameter kappa
  a.kappa <- priors$kappa$hyperpars[1]
  b.kappa <- priors$kappa$hyperpars[2]
  dist.kappa <- priors$kappa$priorDist
  #hyperparameters  for change in mean delta
  a.delta <- priors$delta$hyperpars[1]
  b.delta <- priors$delta$hyperpars[2]
  dist.delta <- priors$delta$priorDist
  #hyperparameters  for zero-inflation
  a.phi <- priors$phi$hyperpars[1]
  b.phi <- priors$phi$hyperpars[2]
  dist.phi <- priors$phi$priorDist
  paste0('data {
           int J; // number of animals
           int ystararaw[J]; // after treatment McMaster count
           int ystarbraw[J]; // before treatment McMaster count
           int fpre[J];
           int fpost[J];
        }
         parameters {
           real<lower=0> kappa;
           real<lower=0> mu;
           real<lower=0,upper=1> delta;
           real<lower=0> mub[J];
           real<lower=0,upper=1> phi;
         }
        transformed parameters{
          real lambdaa[J];
          real lambdab[J];
          real kappamu;
          for (i in 1:J){
            lambdab[i] <- mub[i]/fpre[i];
            lambdaa[i] <- delta*mub[i]/fpost[i];
          }
          kappamu <- kappa/mu;
        }
        model {
          mu ~ ',dist.mu,'(',a.mu,',',b.mu,');           // prior
          kappa ~ ',dist.kappa,'(',a.kappa,',',b.kappa,');
          delta ~ ',dist.delta,'(',a.delta,',',b.delta,');
          phi ~ ',dist.phi,'(',a.phi,',',b.phi,');
          mub ~ gamma(kappa,kappamu);           // likelihoods
          for (n in 1:J) {
          if (ystarbraw[n] == 0)
            increment_log_prob(log_sum_exp(bernoulli_log(1,phi), bernoulli_log(0,phi)+poisson_log(ystarbraw[n],lambdab[n])));
            else
            increment_log_prob(bernoulli_log(0,phi) + poisson_log(ystarbraw[n],lambdab[n]));
          }
          for (n in 1:J) {
            if (ystararaw[n] == 0)
            increment_log_prob(log_sum_exp(bernoulli_log(1,phi), bernoulli_log(0,phi)+poisson_log(ystararaw[n],lambdaa[n])));
            else 
            increment_log_prob(bernoulli_log(0,phi) + poisson_log(ystararaw[n],lambdaa[n]));
          }
        }
')
}
