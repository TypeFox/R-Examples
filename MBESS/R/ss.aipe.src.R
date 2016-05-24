ss.aipe.src <- function(Rho2.Y_X=NULL, Rho2.k_X.without.k=NULL, K=NULL, 
beta.k=NULL, width, which.width="Full", sigma.Y=1, sigma.X.k=1, RHO.XX=NULL, 
Rho.YX=NULL, which.predictor=NULL, alpha.lower=NULL, alpha.upper=NULL, 
conf.level=.95, degree.of.certainty=NULL, assurance=NULL, certainty=NULL,Suppress.Statement=FALSE)
{
if(!is.null(certainty)& is.null(degree.of.certainty)&is.null(assurance)) degree.of.certainty<-certainty
if (is.null(assurance) && !is.null (degree.of.certainty)& is.null(certainty)) assurance <-degree.of.certainty
if (!is.null(assurance) && is.null (degree.of.certainty)& is.null(certainty)) assurance -> degree.of.certainty

if(!is.null(assurance) && !is.null (degree.of.certainty) && assurance!=degree.of.certainty) 
stop("The arguments 'assurance' and 'degree.of.certainty' must have the same value.")

if(!is.null(assurance) && !is.null (certainty) && assurance!=certainty) 
stop("The arguments 'assurance' and 'certainty' must have the same value.")

if(!is.null(degree.of.certainty) && !is.null (certainty) && degree.of.certainty!=certainty) 
stop("The arguments 'degree.of.certainty' and 'certainty' must have the same value.")

return(ss.aipe.reg.coef(Rho2.Y_X=Rho2.Y_X, Rho2.j_X.without.j=Rho2.k_X.without.k, p=K, b.j=beta.k, 
width=width, which.width=which.width, sigma.Y=sigma.Y, sigma.X=sigma.X.k, 
RHO.XX=RHO.XX, Rho.YX=Rho.YX, which.predictor=which.predictor, Noncentral=TRUE, 
alpha.lower=alpha.lower, alpha.upper=alpha.upper, conf.level=conf.level, 
degree.of.certainty=degree.of.certainty, Suppress.Statement=Suppress.Statement))
}
