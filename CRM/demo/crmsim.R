
prior1 <- c(0.05,0.1,0.2,0.3,0.5,0.7)
true1 <- c(0.1,0.15,0.2,0.4,0.5,0.8)

# simulations using model 1 (hyperbolic tangent model)
crmsim(target=0.2,prior=prior1,true=true1,rate=0.1,cycle=21,cohort=1,nsubject=24,nsim=100,
       model=1,a0=1,b=3,jump=FALSE,start.dose=1,seed=777)

crmsiminc1(target=0.2,prior=prior1,true=true1,rate=0.1,cycle=21,nsubject=24,nsim=100,
           model=1,a0=1,b=3,jump=FALSE,start.dose=1,seed=777)

crmsiminc2(target=0.2,prior=prior1,true=true1,rate=0.1,cycle=21,nsubject=24,nsim=100,
           model=1,a0=1,b=3,jump=FALSE,start.dose=1,seed=777)

# simulations using model 2 (one-parameter logistic model)
crmsim(target=0.2,prior=prior1,true=true1,rate=0.1,cycle=21,cohort=1,nsubject=24,nsim=100,
       model=2,a0=1,b=3,jump=FALSE,start.dose=1,seed=777)

crmsiminc1(target=0.2,prior=prior1,true=true1,rate=0.1,cycle=21,nsubject=24,nsim=100,
           model=2,a0=1,b=3,jump=FALSE,start.dose=1,seed=777)

crmsiminc2(target=0.2,prior=prior1,true=true1,rate=0.1,cycle=21,nsubject=24,nsim=100,
           model=2,a0=1,b=3,jump=FALSE,start.dose=1,seed=777)

