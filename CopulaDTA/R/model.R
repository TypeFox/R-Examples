#' Specify the copula based bivariate beta-binomial distribution to fit to the diagnostic data.
#'
#' @param copula A description of the copula function used to model the correlation between sensitivity and specificty.
#' This is a string naming the copula function. The choices are "fgm", "frank", "gauss", "c90" and "c270".
#' @param modelargs An optional list of control parameter for the prior distributions. The parameters in the list include:
#' \itemize{
#' \item{formula.se:} {An  object of class "formula": A symbolic description of a linear model to be fitted to mean E(x) of sensitivity in the logit scale.
#' the default (when no covariates are included) symbolic description is SID ~ 1 corresponds to the model formula E(x) = mu = exp(a)/(1 + exp(a)) where a is the intercept.
#' When the covariates are categorical and the relative measures are needed it is important to remove the interecept from the model to obtain meaningful parameters. EG for
#' a covariate 'Test' with two levels(A and B) and relative sensitivity of B versus A is needed, then the correct formula is SID ~ Test - 1 or SID ~ Test + 0. See \link[stats]{formula}.
#' For further information on interpretation of parameters in logistic regression see Agresti A(2002) Chapter 5.}
#' \item{formula.sp:} {An object of class "formula": A symbolic description of a linear model to be fitted to specificity data.
#' By default the covariate information for sensitivity is used.}
#' \item{formula.omega:} { An object of class "formula": A symbolic description of a linear model to be fitted to the copula function.
#' By default the covariate information for sensitivity is used.}
#' \item{transform.omega:} { A logical value indicating whether a constrained correlation parameter should be mapped into an non-constrained scale.
#' This applies to all the allowed copula functions except "frank". The default is TRUE.}
#' \item{param:} { Indication of the parameterisation used to map the marginal mean and precision/dispersion to the alpha and beta parameters of the beta distribution.
#' There are two choices: param=1 which uses \deqn{alpha = mu*phi, beta = (1 - mu)*phi} where \deqn{mu = alpha/(alpha + beta), 0<=mu<=1,} and
#' \deqn{phi = alpha + beta, phi >= 0.}
#' param=2 uses \deqn{alpha = ((1 - phi)/phi)*mu}
#' \deqn{beta = ((1 - phi)/phi)*(1 - mu)} where \deqn{mu = alpha/(alpha + beta); 0<=mu<=1,} and
#' \deqn{phi = 1/( 1 + alpha + beta); 0<=phi<=1.}}
#'\item{prior.lse:}{ A description of prior distribution of the marginal mean sensitivity in the logit scale. The default is "normal" distribution.
#'For other distributions see stan documentation at \url{http://mc-stan.org/documentation/}.}
#'\item{par.lse1:}{ A numeric value indicating the location of the prior distribution of the marginal mean sensitivity in the logit scale.
#'The default is 0 which implying a distribution centered around 0.5 in the 0-1 scale.}
#'\item{par.lse2:}{ A numeric value indicating the spread(standard deviation) pf the prior distribution of the marginal mean sensitivity in the logit scale
#'and can be interpreted as the quantity of prior information.
#' vague and non-informative priors are specified by a distribution with large variance. The default is sd=10 implying that the variance is 100.}
#'\item{prior.lsp:}{ A  description of prior distribution of the marginal mean specificity in the logit scale. The default is "normal" distribution.}
#'\item{ par.lsp1:}{ A numeric value indicating the location of the prior distribution of the marginal mean specificity in the logit scale.
#'The default is 0 which implying a distribution centered around 0.5 in the 0-1 scale.}
#'\item{par.lsp2:}{ A numeric value indicating the spread(standard deviation) pf the prior distribution of the marginal mean specificity in the logit scale
#'and can be interpreted as the quantity of prior information.
#' vague and non-informative priors are specified by a distribution with large variance. The default is sd=10 implying that the variance is 100.}
#'\item{prior.omega:}{ A description of prior distribution of the correlation parameter(s). The default is "normal" distribution since "transform.omega=TRUE".
#'When "transform.omega=FALSE" the candidate prior distributions are U[-1, 1] for fgm and gaussian copulas, and half-cauchy(0, 2.5), gamma(0.001, 0.001) for the C90 and C270.}
#'\item{par.omega1:}{ A numeric value indicating the location of the prior distribution of the correlation parameter(s). The default is 0.}
#'\item{par.omega2:}{ A numeric value indicating the scale/spread(standard deviation) of the prior distribution of the correlation parameter(s). The default is sd=10.}
#' }
#'@return An object of cdtamodel class.
#'@examples
#' data(telomerase)
#' model1 <-  cdtamodel(copula = 'fgm')
#'
#' model2 <- cdtamodel(copula = 'fgm',
#'                modelargs=list(param=2,
#'                               prior.lse='normal',
#'                               par.lse1=0,
#'                               par.lse2=5,
#'                               prior.lsp='normal',
#'                               par.lsp1=0,
#'                               par.lsp2=5))
#'
#' model3 <-  cdtamodel(copula = 'fgm',
#'                modelargs = list(formula.se = StudyID ~ Test - 1))
#'
#'@references {Agresti A (2002). Categorical Data Analysis. John Wiley & Sons, Inc.}
#'@references {Clayton DG (1978). A model for Association in Bivariate Life Tables and its Application in
#'Epidemiological Studies of Familial Tendency in Chronic Disease Incidence. Biometrika,65(1), 141-151.}
#'@references {Frank MJ (1979). On The Simultaneous Associativity of F(x, y) and x + y - F(x, y). Aequationes Mathematicae, pp. 194-226.}
#'@references {Farlie DGJ (1960). The Performance of Some Correlation Coefficients for a General Bivariate
#'Distribution. Biometrika, 47, 307-323.}
#'@references {Gumbel EJ (1960). Bivariate Exponential Distributions. Journal of the American Statistical Association, 55, 698-707.}
#'@references {Meyer C (2013). The Bivariate Normal Copula. Communications in Statistics - Theory and Methods, 42(13), 2402-2422.}
#'@references {Morgenstern D (1956). Einfache Beispiele Zweidimensionaler Verteilungen. Mitteilungsblatt furMathematische Statistik, 8, 23 - 235.}
#'@references {Sklar A (1959). Fonctions de Repartition a n Dimensions et Leurs Marges. Publications de l'Institut de Statistique de L'Universite de Paris, 8, 229-231.}
#'@export
#' @author Victoria N Nyaga \email{victoria.nyaga@outlook.com}
#'
#'
cdtamodel <- function(copula,
                  modelargs = list()
                  ) {


    if (is.null(modelargs$formula.se)) modelargs$formula.se = SID ~ 1
    attr(modelargs$formula.se, 'formula')

    if (is.null(modelargs$formula.sp)) modelargs$formula.sp = modelargs$formula.se
    if (is.null(modelargs$formula.omega)) modelargs$formula.omega = modelargs$formula.se
    if (is.null(modelargs$transform.omega)) modelargs$transform.omega=TRUE
    if (is.null(modelargs$param)) modelargs$param=1
    if (is.null(modelargs$prior.lse)) modelargs$prior.lse='normal'
    if (is.null(modelargs$par.lse1)) modelargs$par.lse1=0
    if (is.null(modelargs$par.lse2)) modelargs$par.lse2=10
    if (is.null(modelargs$prior.lsp)) modelargs$prior.lsp='normal'
    if (is.null(modelargs$par.lsp1)) modelargs$par.lsp1=0
    if (is.null(modelargs$par.lsp2)) modelargs$par.lsp2=10
    if (is.null(modelargs$prior.omega)) modelargs$prior.omega='normal'
    if (is.null(modelargs$par.omega1)) modelargs$par.omega1=0
    if (is.null(modelargs$par.omega2)) modelargs$par.omega2=10

#==================================Build model piece by piece =======================
pt1 <- "functions{"

#=================================Choice model ======================================#
GAUSS <- "\n\treal gaussian_log(matrix p_i, matrix alpha, matrix beta, vector omega){
		real f1;
		real f2;
		vector[rows(p_i)] f3;
		vector[rows(p_i)] f4;
		int r;

		r <- rows(p_i);
		f1 <- beta_log(col(p_i, 1), col(alpha,1), col(beta,1));
		f2 <- beta_log(col(p_i, 2), col(alpha,2), col(beta,2));
		for (i in 1:r){
			f3[i] <- (1/sqrt(1 - omega[i]^2))*exp((2*omega[i]*inv_Phi(beta_cdf(p_i[i, 1], alpha[i,1], beta[i,1]))*inv_Phi(beta_cdf(p_i[i, 2], alpha[i,2], beta[i,2])) -
			omega[i]^2*(inv_Phi(beta_cdf(p_i[i, 1], alpha[i,1], beta[i,1]))^2 + inv_Phi(beta_cdf(p_i[i, 2], alpha[i,2], beta[i,2]))^2))/
			(2*(1 - omega[i]^2)));
		}
		return (f1 + f2 + sum(log(f3)));
	}
"

FRANK <- "\n\treal frank_log(matrix p_i, matrix alpha, matrix beta, vector omega){
	real f1;
	real f2;
	int r;
	vector[rows(p_i)] f3;
	vector[rows(p_i)] f4;
	r <- rows(p_i);

	f1 <- beta_log(col(p_i, 1), col(alpha,1), col(beta, 1));
	f2 <- beta_log(col(p_i, 2), col(alpha, 2), col(beta, 2));
	for (i in 1:r){
		f3[i] <- (omega[i]*(1 - exp(-omega[i]))*exp(-omega[i]*(beta_cdf(p_i[i,1], alpha[i,1], beta[i,1]) + beta_cdf(p_i[i,2], alpha[i,2], beta[i,2]))));
		f4[i] <- ((1 - exp(-omega[i])) - (1 - exp(-omega[i]*beta_cdf(p_i[i,1], alpha[i,1], beta[i,1])))*(1 - exp(-omega[i]*beta_cdf(p_i[i,2], alpha[i,2], beta[i,2]))));
	}
	return (f1 + f2 + sum(log((f3)./((f4).*(f4)))));
}
"
FGM <- "\n\treal fgm_log(matrix p_i, matrix alpha, matrix beta, vector omega){
	real f1;
	real f2;
	vector[rows(p_i)] f3;
    int r;

    r <- rows(p_i);

	f1 <- beta_log(col(p_i, 1), col(alpha, 1), col(beta, 1));
	f2 <- beta_log(col(p_i, 2), col(alpha, 2), col(beta, 2));
	for (i in 1:r){
		f3[i] <- log(1 + omega[i]*(1 - 2*beta_cdf(p_i[i,1], alpha[i,1], beta[i,1]))*(1 - 2*beta_cdf(p_i[i,2], alpha[i,2], beta[i,2])));
	}
	return (f1 + f2 + sum(f3));
	}
"
C90 <- "\n\treal clayton90_log(matrix p_i, matrix alpha, matrix beta, vector omega){
	real f1;
	real f2;
	vector[rows(p_i)] f3;
	vector[rows(p_i)] f4;
	int r;
	vector[rows(p_i)] powr;

	r <- rows(p_i);
	f1 <- beta_log(col(p_i, 1), col(alpha,1), col(beta,1));
	f2 <- beta_log(col(p_i, 2), col(alpha,2), col(beta,2));
	powr <- -(2*omega + 1)./(omega);

	for (i in 1:r){
		f3[i] <- (1 + omega[i])*((1 - beta_cdf(p_i[i,1], alpha[i,1], beta[i,1]))^(-(1 + omega[i])))*(beta_cdf(p_i[i,2], alpha[i,2], beta[i,2])^(-(1 + omega[i])))*
		(((1 - beta_cdf(p_i[i,1], alpha[i,1], beta[i,1]))^(-omega[i]) + beta_cdf(p_i[i,2], alpha[i,2], beta[i,2])^(-omega[i]) - 1)^powr[i]);
	}
	return (f1 + f2 + sum(log(f3)));
	}
"
C270 <- "\n\treal clayton270_log(matrix p_i, matrix alpha, matrix beta, vector omega){
	real f1;
	real f2;
	vector[rows(p_i)] f3;
	vector[rows(p_i)] f4;
	int r;
	vector[rows(p_i)] powr;

	r <- rows(p_i);
	f1 <- beta_log(col(p_i, 1), col(alpha,1), col(beta,1));
	f2 <- beta_log(col(p_i, 2), col(alpha,2), col(beta,2));
	powr <- -(2*omega + 1)./(omega);

	for (i in 1:r){
		f3[i] <- (1 + omega[i])*(beta_cdf(p_i[i,1], alpha[i,1], beta[i,1])^(-(1 + omega[i])))*((1 - beta_cdf(p_i[i,2], alpha[i,2], beta[i,2]))^(-(1 + omega[i])))*
		((beta_cdf(p_i[i,1], alpha[i,1], beta[i,1])^(-omega[i]) + (1 - beta_cdf(p_i[i,2], alpha[i,2], beta[i,2]))^(-omega[i]) - 1)^powr[i]);
	}
	return (f1 + f2 + sum(log(f3)));
	}
"

if (copula=="gauss"){
		copkies<-GAUSS
	} else if (copula=="frank") {
		 copkies<-FRANK
	} else if (copula=="fgm") {
		 copkies<-FGM
	} else if (copula=="c90") {
		 copkies<-C90
	} else if (copula=="c270") {
		 copkies<-C270
	} else {
		stop("Invalid copula chosen")}

#========================== Data ============================================#

dat <- "\n}\n data{
		int<lower=0> Ns;
		int<lower=0> tp[Ns];
		int<lower=0> dis[Ns];
		int<lower=0> tn[Ns];
		int<lower=0> nondis[Ns];

		int<lower=0> Npse;
		matrix<lower=0>[Ns,Npse] xse;

		int<lower=0> Npsp;
		matrix<lower=0>[Ns,Npsp] xsp;

		int<lower=0> Npomega;
		matrix<lower=0>[Ns,Npomega] xomega;\n}"

#======================= Parameters ====================================#
params <- "\n parameters{
		vector[Npse] betamuse;
		vector[Npsp] betamusp;
		vector[Npse] betaphise;
		vector[Npsp] betaphisp;
		matrix<lower=0, upper=1>[Ns,2] p_i;
		"

if(copula=="frank" | modelargs$transform.omega==TRUE){
	betaomega <- "vector[Npomega] betaomega;"
} else if (copula=="gauss" |copula=="fgm" ){
	 betaomega <- "vector<lower=-1, upper=1>[Npomega] betaomega;"

} else if (copula=="c90" |copula=="c270" ){
	 betaomega <- "vector<lower=0>[Npomega] betaomega;"}

#======================= Transformed Parameters ====================================#
transf_params.pt1 <- "\n } \n transformed parameters{
		matrix<lower=0, upper=1>[Ns,2] mui;
		vector<lower=0, upper=1>[Ns] musei;
		vector<lower=0, upper=1>[Ns] muspi;
		vector[Npse] MUse;
		vector[Npsp] MUsp;
		vector[Npse] RRse;
		vector[Npsp] RRsp;
		matrix<lower=0>[Ns,2] alpha;
		matrix<lower=0>[Ns,2] beta;"

transf_params.pt2 <- "\n\n\t\tmusei <- exp(xse*betamuse)./(1 + exp(xse*betamuse));
		muspi <- exp(xsp*betamusp)./(1 + exp(xsp*betamusp));
		mui <- append_col(musei, muspi);

		MUse <- exp(betamuse)./(1 + exp(betamuse));
		MUsp <- exp(betamusp)./(1 + exp(betamusp));

		RRse[1] <- 1;
		for (i in 2:Npse)
			RRse[i] <- MUse[i]/MUse[1];

		RRsp[1] <- 1;
		for (i in 2:Npsp)
			RRsp[i] <- MUsp[i]/MUsp[1];"
if(modelargs$param==1){
phi.pt1 <- "\n\t\tvector<lower=0>[Ns] phisei; \n\t\tvector<lower=0>[Ns] phispi; \n\t\tmatrix<lower=0>[Ns,2] phi;\n"
phi.pt2 <- "\n\t\tphisei <- exp(xse*betaphise);
		phispi <- exp(xsp*betaphisp);
		phi <- append_col(phisei, phispi);
		alpha <- (mui).*phi;
		beta <- (1 - mui).*phi;\n"
}else{

phi.pt1 <- "\n\t\tvector<lower=0, upper=1>[Ns] phisei; \n\t\tvector<lower=0, upper=1>[Ns] phispi; \n\t\tmatrix<lower=0, upper=1>[Ns,2] phi;\n"
phi.pt2 <- "\n\t\tphisei <- exp(xse*betaphise)./(1 + exp(xse*betaphise));
		phispi <- exp(xsp*betaphisp)./(1 + exp(xsp*betaphisp));
		phi <- append_col(phisei, phispi);
		alpha <- ((1 - phi)./(phi)).*(mui);
		beta <- ((1 - phi)./(phi)).*(1 - mui);\n"
}

if (copula=="frank"){
	omega.pt1 <- "\t\t vector[Ns] omega;"
	omega.pt2 <- "\n\t\t omega <- xomega*betaomega;"
} else if (copula=="gauss"|copula=="fgm"){
	   if(modelargs$transform.omega==TRUE){
			omega.pt1 <- "\t\tvector[Ns] omegat;\n\t\tvector<lower=-1, upper=1>[Ns] omega;\n\t\tvector[Npomega] betaomegat;"
			omega.pt2 <- "\n\t\tomegat <- xomega*betaomega; \n\t\tfor (s in 1:Ns) \n\t\t\tomega[s] <- tanh(omegat[s]);\n\t\tfor (o in 1:Npomega) \n\t\t\tbetaomegat[o] <- tanh(betaomega[o]);"
		}else{
			omega.pt1 <- "\t\tvector<lower=-1, upper=1>[Ns] omega;"
			omega.pt2 <- "\n\t\tomega <- xomega*betaomega;"
	}
} else {
	if(modelargs$transform.omega==TRUE){
		omega.pt1 <- "\t\tvector[Ns] omegat;\n\t\tvector<lower=0>[Ns] omega;\n\t\tvector[Npomega] betaomegat;"
		omega.pt2 <-"\n\t\tomegat <- xomega*betaomega; \n\t\tfor (s in 1:Ns) \n\t\t\tomega[s] <- exp(omegat[s]);\n\t\tfor (o in 1:Npomega) \n\t\t\tbetaomegat[o] <- exp(betaomega[o]);"
	}else{
		omega.pt1 <- "\t\t vector<lower=0>[Ns]omega;"
		omega.pt2 <- "\n\t\t omega <- xomega*betaomega;"
	}
}

# if (copula=="frank"){
# 	ktau.pt1 <- "\n"
# } else{
	ktau.pt1 <- "\n\t\tvector[Npomega] ktau;"
#}

if (copula=="frank"){
    ktau.pt2 <- "\n\t\tfor (o in 1:Npomega) \n\t\t\tktau[o] <- 1;"
} else if (copula=="gauss"){
	if(modelargs$transform.omega==TRUE){
		ktau.pt2 <- "\n\t\tfor (o in 1:Npomega) \n\t\t\tktau[o] <- (2/pi())*asin(betaomegat[o]);"
	}else{
		ktau.pt2 <- "\n\t\tfor (o in 1:Npomega) \n\t\t\tktau[o] <- (2/pi())*asin(betaomega[o]);"
	}
} else if(copula=="fgm"){
   if(modelargs$transform.omega==TRUE){
		ktau.pt2 <- "\n\t\tfor (o in 1:Npomega) \n\t\t\tktau[o] <- (2*betaomegat[o])/9;"
	}else{
		ktau.pt2 <- "\n\t\tfor (o in 1:Npomega) \n\t\t\tktau[o] <- (2*betaomega[o])/9;"
	}
} else {
	if(modelargs$transform.omega==TRUE){
		ktau.pt2 <- "\n\t\tfor (o in 1:Npomega) \n\t\t\tktau[o] <- -(betaomegat[o])/(betaomegat[o] + 2);"
	}else{
		ktau.pt2 <- "\n\t\tfor (o in 1:Npomega) \n\t\t\tktau[o] <- -(betaomega[o])/(betaomega[o] + 2);"
	}
}

#======================================Priors ===========================================#
pt2 <- "\n}\n model{\n\t #priors \n"

priorse <- paste('\t betamuse ~ ', modelargs$prior.lse,'(', modelargs$par.lse1, ', ', modelargs$par.lse2, ');\n', sep='')
priorsp <- paste('\t betamusp ~ ', modelargs$prior.lsp,'(', modelargs$par.lsp1, ', ', modelargs$par.lsp2, ');\n', sep='')

if(copula=="frank"){if(is.null(modelargs$prior.omega)){prioromega<-paste('\t betaomega ~ ', "normal",'(', 0, ', ', 10, ');\n', sep='')} else{
	prioromega<-paste('\t betaomega ~ ', modelargs$prior.omega,'(', modelargs$par.omega1, ', ', modelargs$par.omega2, ');\n', sep='')

}} else {
	if(modelargs$transform.omega==FALSE){if(is.null(modelargs$prior.omega)){prioromega<-paste('\t betaomega ~ ', "normal",'(', 0, ', ', 10, ');\n', sep='')} else{
		prioromega<-paste('\t betaomega ~ ', modelargs$prior.omega,'(', modelargs$par.omega1, ', ', modelargs$par.omega2, ');\n', sep='')

}} else {
	if(is.null(modelargs$prior.omega)){prioromega<-paste('\t betaomega ~ ', "normal",'(', 0, ', ', 10, ');\n', sep='')} else {
		prioromega=paste('\t betaomega ~ ', modelargs$prior.omega,'(', modelargs$par.omega1, ', ', modelargs$par.omega2, ');\n', sep='')}}}


if (copula=="gauss"){
	priorp <- "\n\t p_i ~ gaussian(alpha, beta, omega);\n"
} else if (copula=="frank"){
	priorp <- "\n\t p_i ~ frank(alpha, beta, omega);\n"
} else if (copula=="fgm"){
	priorp <- "\n\t p_i ~ fgm(alpha, beta, omega);\n"
} else if (copula=="c90"){
	priorp <- "\n\t p_i ~ clayton90(alpha, beta, omega);\n"
} else if(copula=="c270"){
	priorp <- "\n\t p_i ~ clayton270(alpha, beta, omega);\n"
}

#===============================Likelihood=============================================#

lik <- "\n\t tp ~ binomial(dis,col(p_i,1)); \n\t tn ~ binomial(nondis, col(p_i, 2)); \n }"

#===================================Generated Quantities=============================#
GQ <- "\n generated quantities{
	vector[Ns*2] loglik;
	for (i in 1:Ns)
		loglik[i] <- binomial_log(tp[i], dis[i], p_i[i,1]);

	for (i in (Ns+1):(2*Ns))
		loglik[i] <- binomial_log(tn[i-Ns], nondis[i-Ns], p_i[i-Ns,2]);\n}"

modelcode <- paste(pt1,
		  copkies,
		  dat,
		  params,
		  betaomega,
		  transf_params.pt1,
		  phi.pt1,
		  omega.pt1,
		  ktau.pt1,
		  transf_params.pt2,
		  phi.pt2,
		  omega.pt2,
		  ktau.pt2,
		  pt2,
		  priorse,
		  priorsp,
		  prioromega,
		  priorp,
		  lik,
		  GQ,
		  sep='')

out <- new("cdtamodel",
            copula=copula,
            modelcode=modelcode,
            modelargs=modelargs)

return(out)

}


