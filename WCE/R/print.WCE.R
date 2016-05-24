print.WCE <- function(x,...){

if (x$analysis == 'Cox'){
if (is.na (sum(x$PL)) == T) {
	for (i in 1:length(x$PL)){ if (is.na(x$PL[i])==T) {cat('Warning : model', i, 'did not converge, and no \npartial log-likelihood was produced. Results \nfor this model should be ignored.\n\n')}}}
for (i in 1:length(x$PL)) {
	if (sum(x$SED[[i]]==0) >0) {cat('Warning : some of the SE for the spline \nvariables in model', i, 'are exactly zero, probably \nbecause the model did not converge. Variable(s)',  names(which(x$SED[[1]]==0)), ' \nhad SE=0. Consider re-parametrizing or increasing \nthe number of iteractions.\n\n')}}
} else {
if (is.na (sum(x$LL)) == T) {
	for (i in 1:length(x$LL)){ if (is.na(x$LL[i])==T) {cat('Warning : model', i, 'did not converge, and no \npartial log-likelihood was produced. Results \nfor this model should be ignored.\n\n')}}}
for (i in 1:length(x$LL)) {
	if (sum(x$SED[[i]]==0) >0) {cat('Warning : some of the SE for the spline \nvariables in model', i, 'are exactly zero, probably \nbecause the model did not converge. Variable(s)',  names(which(x$SED[[1]]==0)), ' \nhad SE=0. Consider re-parametrizing or increasing \nthe number of iteractions.\n\n')}}
}
if (x$constrained == 'Right') {
  cat("\nEstimated right-constrained WCE function(s).\n\n")} 
if (x$constrained == 'Left') {
  cat("\nEstimated left-constrained WCE function(s).\n\n")}
if (x$constrained == FALSE) {
  cat("\nEstimated unconstrained WCE function(s).\n\n")}
print(x$WCEmat)
if (x$analysis == 'Cox') {
	cat("\nNumber of events:\n")
	print(x$ne)
	cat("\nPartial log-Likelihoods:\n")
	print(x$PL)} else {
	cat("\nLog-Likelihoods:\n")
	colnames(x$LL) <- NULL
	print(x$LL)
}
if (x$a == T){cat("\nAIC:\n")} else {cat("\nBIC:\n")}
colnames(x$info.criterion) <- NULL
print(x$info.criterion) 
if (is.null(x$covariates) == F){ 
  cat("\nMatrix of coefficients estimates for the covariates:\n\n")
  print(x$beta.hat.covariates)
  cat("\nMatrix of standard error estimates for the covariates:\n\n")
  print(x$se.covariates)}
cat('\nIf you report these results, please cite Sylvestre MP, Abrahamowicz M. Flexible Modeling of the Effects of Time-Dependent Exposures on the
Hazard. Statistics in Medicine 2009; 28(27):3437-3453.\n')
}