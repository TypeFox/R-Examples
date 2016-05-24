test.normal_default.converge <- function() {
    data(x_gideon)
    ## distlist=NULL defaults to test mix of normals, 
    ## pars.init=NULL defaults to estimated initial conditions
    checkEquals( mix.mle(x_gideon$horn)$convergence, 0, msg='Convergence, default parameters (2 normals, estimated initial conditions')
}

test.normal_method.converge <- function() {
    ## test mix of normals specified
    ## pars.init=NULL defaults to estimated initial conditions
    data(x_gideon)
    checkEquals( mix.mle(x_gideon$horn, method='normal')$convergence, 0, msg='Convergence, 2 normals specified, estimated initial conditions')
}

test.normal_method_mix.prob.converge <- function() {
    ## test mix of normals specified
    ## pars.init=NULL defaults to estimated initial conditions
    data(x_gideon)
    checkEquals( mix.mle(x_gideon$horn, mix.prob=0.4, method='normal')$convergence, 0, msg='Convergence, 2 normals specified, estimated initial conditions')
}

est.facing.gamma <- function() {
    ## define params for facing gamma
    ## values taken from o_taurus
    ## estimate a gamma mix
    .lower <- 0.37
    .upper <- 4.72
    .init <- list(mix.prob=0.55, dist1.par1=1.50, dist1.par2=.4, dist2.par1=3.2, dist2.par2=.5)
    ## simulate synthetic data
    .synthetic <- with(.init, 
        mix.synthetic.facing.gamma( lower=.lower, upper=.upper,  
            shape1=dist1.par1, scale1=dist1.par2, 
            shape2=dist2.par1, scale2=dist2.par2, 
            mix.prob=mix.prob
        ) 
    )
    ## fit and return
    .fit <- with(.init,
        mix.mle(.synthetic,  method='facing.gamma', 
            lower=.lower, upper=.upper, mix.prob=mix.prob, 
            dist1.par1=dist1.par1, dist1.par2=dist1.par2, 
            dist2.par1=dist2.par1, dist2.par2=dist2.par2
        )
    )
}

test.facing_gamma.converge <- function() {
    ## estimate, check for convergence
    test <- est.facing.gamma()
    checkEquals( test$convergence, 0, msg='Convergence, facing gammas, specified initial conditions')
}
