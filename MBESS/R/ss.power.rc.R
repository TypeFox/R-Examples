`ss.power.rc` <-
function(Rho2.Y_X=NULL, Rho2.Y_X.without.k=NULL, K=NULL, 
desired.power=.85, alpha.level=.05, Directional=FALSE, beta.k=NULL, 
sigma.X=NULL, sigma.Y=NULL, Rho2.k_X.without.k=NULL, RHO.XX=NULL, Rho.YX=NULL, 
which.predictor=NULL, Cohen.f2=NULL, Specified.N=NULL, Print.Progress=FALSE)
{
ss.power.reg.coef(Rho2.Y_X=Rho2.Y_X, Rho2.Y_X.without.j=Rho2.Y_X.without.k, p=K, 
desired.power=desired.power, alpha.level=alpha.level, Directional=Directional, beta.j=beta.k, 
sigma.X=sigma.X, sigma.Y=sigma.Y, Rho2.j_X.without.j=Rho2.k_X.without.k, RHO.XX=RHO.XX, Rho.YX=Rho.YX, 
which.predictor=which.predictor, Cohen.f2=Cohen.f2, Specified.N=Specified.N, Print.Progress=Print.Progress)
}

