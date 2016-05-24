ss.aipe.src.sensitivity <-function(True.Var.Y=NULL, True.Cov.YX=NULL, True.Cov.XX=NULL, 
Estimated.Var.Y=NULL, Estimated.Cov.YX=NULL, Estimated.Cov.XX=NULL, Specified.N=NULL, 
which.predictor=1, w=NULL, Noncentral=TRUE, Standardize=TRUE, conf.level=.95, 
degree.of.certainty=NULL, assurance=NULL, certainty=NULL, G=1000, print.iter=TRUE)
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

return(ss.aipe.reg.coef.sensitivity(True.Var.Y=True.Var.Y,True.Cov.YX=True.Cov.YX, True.Cov.XX=True.Cov.XX, 
Estimated.Var.Y=Estimated.Var.Y, Estimated.Cov.YX=Estimated.Cov.YX, Estimated.Cov.XX=Estimated.Cov.XX, Specified.N=Specified.N, 
which.predictor=which.predictor, w=w, Noncentral=Noncentral, Standardize=Standardize, conf.level=conf.level, 
degree.of.certainty=degree.of.certainty, G=G, print.iter=print.iter))
}
