library("RUnit")
library("chngpt")


test.chngptm <- function() {

#tolerance=1e-3
#if (R.Version()$system %in% c("x86_64, mingw32")) {
#    tolerance=1e-6
#} else if (substr(R.Version()$system, 1, 5)=="sparc" | substr(R.Version()$system, 1, 4)=="i686") { 
#    tolerance=1e-1
#}
    
if (R.Version()$system %in% c("x86_64, mingw32")) {

    tolerance=1e-6
    RNGkind("Mersenne-Twister", "Inversion")    
    data=sim.chngpt("sigmoid4", type="step", n=250, seed=1, beta=0, x.distr="norm", e.=3.4, b.=-Inf)

    ########################################
    # main effect
    
    fit = chngptm (y~z, ~x, family="binomial", data, tol=1e-4, maxit=1e3, type="step",      lb.quantile=0.1, ub.quantile=0.9, est.method="sigmoidapprox")
    checkEqualsNumeric(fit$coefficients, c(-0.8522692,0.3352071,0.5255425,3.7278304), tolerance=tolerance)
    fit = chngptm (y~z, ~x, family="binomial", data, tol=1e-4, maxit=1e3, type="step",      lb.quantile=0.1, ub.quantile=0.9, est.method="grid")
    checkEqualsNumeric(fit$coefficients, c( -0.8522692,0.3352071,0.5475595,3.7352043), tolerance=tolerance)
    
    fit = chngptm (y~z, ~x, family="binomial", data, tol=1e-4, maxit=1e3, type="hinge",     lb.quantile=0.1, ub.quantile=0.9, est.method="sigmoidapprox")    
    checkEqualsNumeric(fit$coefficients, c(-0.5201786,0.3036120,0.2908861,5.7181125), tolerance=tolerance)    
    fit = chngptm (y~z, ~x, family="binomial", data, tol=1e-4, maxit=1e3, type="hinge",     lb.quantile=0.1, ub.quantile=0.9, est.method="grid")    
    checkEqualsNumeric(fit$coefficients, c(-0.5198031,0.3037050,0.2920897,5.7275861), tolerance=tolerance)    
    
    # this is an example that test leads to less noise
    fit = chngptm (y~z, ~x, family="binomial", data, tol=1e-4, maxit=1e3, type="segmented", lb.quantile=0.1, ub.quantile=0.9, est.method="sigmoidapprox", verbose=FALSE)
    checkEqualsNumeric(fit$coefficients, c( -0.7608995,0.3125371,0.0616057,0.3969299,6.8184476), tolerance=tolerance)    
    fit = chngptm (y~z, ~x, family="binomial", data, tol=1e-4, maxit=1e3, type="segmented", lb.quantile=0.1, ub.quantile=0.9, est.method="grid", verbose=FALSE)
    checkEqualsNumeric(fit$coefficients, c(-0.76090367,0.31253713,0.06159024,0.39690474,6.81840691), tolerance=tolerance)    
    
#    fit = chngptm (y~z, ~x, family="binomial", data, tol=1e-4, maxit=1e3, type="stegmented", lb.quantile=0.1, ub.quantile=0.9, est.method="sigmoidapprox", verbose=FALSE)
#    checkEqualsNumeric(fit$coefficients, c(), tolerance=tolerance)    
    fit = chngptm (y~z, ~x, family="binomial", data, tol=1e-4, maxit=1e3, type="stegmented", lb.quantile=0.1, ub.quantile=0.9, est.method="grid", verbose=FALSE)
    checkEqualsNumeric(fit$coefficients, c(-0.06705628,0.32633105,-0.28688468,0.76675233,0.32650919,3.73520431), tolerance=tolerance)    
    
    
    
    ########################################
    # interaction effect
    
    fit = chngptm (y~z, ~x*z, family="binomial", data, tol=1e-4, maxit=1e3, type="step",      lb.quantile=0.1, ub.quantile=0.9, est.method="sigmoidapprox")
    checkEqualsNumeric(fit$coefficients, c(-0.5391406,0.2158219,0.4210692,0.5553390,6.0611972), tolerance=tolerance)
    fit = chngptm (y~z, ~x*z, family="binomial", data, tol=1e-4, maxit=1e3, type="step",      lb.quantile=0, ub.quantile=1, est.method="grid", verbose=FALSE)
    checkEqualsNumeric(fit$coefficients, c( -0.4815401,0.3042468,16.0476083,-0.3042468,8.3927654), tolerance=tolerance)
    
    fit = chngptm (y~z, ~x*z, family="binomial", data, tol=1e-4, maxit=1e3, type="hinge",     lb.quantile=0.1, ub.quantile=0.9, est.method="grid")
    checkEqualsNumeric(fit$coefficients, c(-0.5368853,0.2618009,0.2526734,0.1231553,5.3939234), tolerance=tolerance)
# Error in chngpt.test(formula.1, formula.2, data, type = type) : 
#  interaction model for this type not implemented yet: hinge
#    fit = chngptm (y~z, ~x*z, family="binomial", data, tol=1e-4, maxit=1e3, type="hinge",     lb.quantile=0.1, ub.quantile=0.9, est.method="sigmoidapprox")
#    checkEqualsNumeric(fit$coefficients, c(-0.5247333,0.2731954,0.3058931,0.1291960,5.7165258), tolerance=tolerance)
    
    fit = chngptm (y~z, ~x*z, family="binomial", data, tol=1e-4, maxit=1e3, type="segmented", lb.quantile=0.1, ub.quantile=0.9, est.method="grid", verbose=FALSE)
    checkEqualsNumeric(fit$coefficients, c(-1.27487926,2.11715992,0.20647617,-0.07942202,-0.74984522,0.81418233,2.65745247), tolerance=tolerance)
#    fit = chngptm (y~z, ~x*z, family="binomial", data, tol=1e-4, maxit=1e3, type="segmented", lb.quantile=0.1, ub.quantile=0.9, est.method="sigmoidapprox", verbose=FALSE)
#    checkEqualsNumeric(fit$coefficients, c(-0.78767807,0.71956572,0.05555813,0.16831871,-0.11267073,0.30441830,5.31461397), tolerance=tolerance)

    fit = chngptm (y~z, ~x*z, family="binomial", data, tol=1e-4, maxit=1e3, type="stegmented", lb.quantile=0.1, ub.quantile=0.9, est.method="grid", verbose=FALSE)
    checkEqualsNumeric(fit$coefficients, c(-0.77952531,0.89672361,0.05137674,0.41156747,-0.09905192,1.46771433,-0.30845593,-0.17022170,6.00352437), tolerance=tolerance)


    ########################################
    # simulate data from hinge and segmented
    
    data=sim.chngpt("sigmoid2", type="hinge", n=250, seed=1, beta=0, x.distr="norm", e.=3.4, b.=-Inf)
    fit = chngptm (y~z, ~x, family="binomial", data, tol=1e-4, maxit=1e3, type="step",      lb.quantile=0.1, ub.quantile=0.9, est.method="sigmoidapprox")
    checkEqualsNumeric(fit$coefficients, c(-0.6543603,0.3976095,-0.5446341,6.2627189), tolerance=tolerance)
    
    data=sim.chngpt("sigmoid2", type="segmented", n=250, seed=1, beta=0, x.distr="norm", e.=3.4, b.=-Inf)
    fit = chngptm (y~z, ~x, family="binomial", data, tol=1e-4, maxit=1e3, type="step",      lb.quantile=0.1, ub.quantile=0.9, est.method="sigmoidapprox")
    checkEqualsNumeric(fit$coefficients, c( -0.1995570,0.5413582,1.6350509,3.3858460), tolerance=tolerance)
        
    data=sim.chngpt("sigmoid2", type="stegmented", n=250, seed=1, beta=0, x.distr="norm", e.=3.4, b.=-Inf)
    fit = chngptm (y~z, ~x, family="binomial", data, tol=1e-4, maxit=1e3, type="stegmented",      lb.quantile=0.1, ub.quantile=0.9, est.method="grid")
    checkEqualsNumeric(fit$coefficients, c(2.5067137,0.4783413,-0.5892674,1.3710912,-1.6098422,5.8467320), tolerance=tolerance)
        

}
}
