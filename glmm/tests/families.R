library(glmm)
theta<-1
delta<-.01

#check derivatives for bernoulli using finite differences
this<-bernoulli.glmm()$cp(theta)*delta
that<-bernoulli.glmm()$cum(theta+delta)-bernoulli.glmm()$cum(theta)
all.equal(this,that)

this<-bernoulli.glmm()$cpp(theta)*delta
that<-bernoulli.glmm()$cp(theta+delta)-bernoulli.glmm()$cp(theta)
all.equal(this,that)

#check derivatives for poisson using finite differences
this<-poisson.glmm()$cp(theta)*delta
that<-poisson.glmm()$cum(theta+delta)-poisson.glmm()$cum(theta)
all.equal(this,that)

this<-poisson.glmm()$cpp(theta)*delta
that<-poisson.glmm()$cp(theta+delta)-poisson.glmm()$cp(theta)
all.equal(this,that)


