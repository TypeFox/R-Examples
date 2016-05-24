#'  An R package for the analysis of (genetic) risk prediction studies.
#'
#' Fueled by the substantial gene discoveries from genome-wide association
#' studies, there is increasing interest in investigating the predictive
#' ability of genetic risk models. To assess the performance of genetic risk
#' models, PredictABEL includes functions for the various measures and plots
#' that have been used in empirical studies, including univariate and
#' multivariate odds ratios (ORs) of the predictors, the c-statistic (or AUC),
#' Hosmer-Lemeshow goodness of fit test, reclassification table, net
#' reclassification improvement (NRI) and integrated discrimination
#' improvement (IDI). The plots included are the ROC plot, calibration plot,
#' discrimination box plot, predictiveness curve, and several risk distributions.
#'
#'
#' These functions can be applied to predicted risks that are obtained using
#' logistic regression analysis, to weighted or unweighted risk scores, for
#' which the functions are included in this package. The functions can also be
#' used to assess risks or risk scores that are constructed using other methods, e.g., Cox Proportional
#' Hazards regression analysis, which are not included in the current version.
#' Risks obtained from other methods can be imported into R for assessment
#' of the predictive performance.
#'
#'
#' The functions to construct the risk models using logistic regression analyses
#' are specifically written for models that include genetic variables,
#' eventually in addition to non-genetic factors, but they can also be applied
#' to construct models that are based on non-genetic risk factors only. \cr
#'
#'
#' Before using the functions \code{\link{fitLogRegModel}} for constructing
#' a risk model or \code{\link{riskScore}} for computing risk
#' scores, the following checks on the dataset are advisable to be done:
#'
#' (1) Missing values: The logistic regression analyses and computation of
#' the risk score are done only for subjects that have no missing data. In case
#' of missing values, individuals with missing data can be removed from the
#' dataset or imputation strategies can be used to fill in missing data.
#' Subjects with missing data can be removed with the R function \code{na.omit}
#' (available in \code{stats} package).
#' Example: \code{DataFileNew <- na.omit(DataFile)}
#' will make a new dataset (\code{DataFileNew}) with no missing values;
#'
#'
#' (2) Multicollinearity: When there is strong correlation between the
#' predictor variables, regression coefficients may be estimated imprecisely
#' and risks scores may be biased because the assumption of independent effects
#' is violated. In genetic risk prediction studies, problems with
#' multicollinearity should be expected when single nucleotide polymorphisms
#' (SNPs) located in the same gene are
#' in strong linkage disequilibrium (LD). For SNPs in LD it is common to select
#' the variant with the lowest p-value in the model;
#'
#'
#' (3) Outliers: When the data contain significant outliers, either clinical
#' variables with extreme values of the outcomes or extreme values resulting
#' from errors in the data entry, these may impact the construction of the risk models and
#' computation of the risks scores. Data should be carefully checked and outliers
#' need to be removed or replaced, if justified;
#'
#' (4) Recoding of data: In the computation of unweighted risk scores, it is assumed
#' that the genetic variants are coded \code{0,1,2}
#' representing the number of alleles carried. When variants
#' are coded \code{0,1} representing a dominant or recessive effect of the alleles,
#' the variables need to be recoded before unweighted risk scores can be computed. \cr
#'
#'
#' To import data into R several alternative strategies can be used. Use the
#' \code{Hmisc} package for importing SPSS and SAS data into R.
#' Use "\code{ExampleData <- read.table("DataName.txt", header=T, sep="\t")}" for text
#' files where variable names are included as column headers and data are
#' separated by tabs.
#' Use "\code{ExampleData <- read.table("Name.csv", sep=",", header=T)}"
#' for comma-separated files with variable names as column headers.
#' Use \code{"setwd(dir)"} to set the working directory to "dir". The datafile
#' needs to be present in the working directory. \cr
#'
#'
#' To export datafiles from R tables to a tab-delimited textfile with the first row as
#' the name of the variables,
#' use "\code{write.table(R_Table, file="Name.txt", row.names=FALSE, sep="\t")}"  and
#' when a comma-separated textfile is requested and variable names are provided in the first row,
#' use "\code{write.table(R_Table, file="Name.csv", row.names=FALSE, sep=",")}".
#' When the directory is not specified, the file will be
#' saved in the working directory. For exporting R data into SPSS, SAS and
#' Stata data, use functions in the the \code{foreign} package. \cr
#'
#' Several functions in this package depend on other R packages:
#'
#' (1) \code{Hmisc}, is used to compute NRI and IDI;
#'
#' (2) \code{ROCR}, is used to produce ROC plots;
#'
#' (3) \code{epitools}, is used to compute  univariate odds ratios;
##'
#' (4) \code{PBSmodelling}, is used to produce predictiveness curve.
#'
#' @note The current version of the package includes the basic measures
#' and plots that are used in the assessment of (genetic) risk prediction models and the
#' function to construct a simulated dataset that contains individual genotype
#' data, estimated genetic risk and disease status, used for the evaluation of
#' genetic risk models (see Janssens et al, Genet Med 2006).
#' Planned extensions of the package include functions to construct risk
#' models using Cox Proportional Hazards analysis for prospective data and
#' assess the performance of risk models for time-to-event data.
#'
#'
##' @author  Suman Kundu
##'
##' Yurii S. Aulchenko
##'
##' A. Cecile J.W. Janssens
##'
##' @keywords package
#'
#' @references S Kundu, YS Aulchenko, CM van Duijn, ACJW Janssens. PredictABEL:
#' an R package for the assessment of risk prediction models.
#' Eur J Epidemiol. 2011;26:261-4. \cr
#'
#' ACJW Janssens, JPA Ioannidis, CM van Duijn, J Little, MJ Khoury.
#' Strengthening the Reporting of Genetic Risk Prediction Studies: The GRIPS
#' Statement Proposal. Eur J Epidemiol. 2011;26:255-9. \cr
#'
#' ACJW Janssens, JPA Ioannidis, S Bedrosian, P Boffetta, SM Dolan, N Dowling,
#' I Fortier, AN. Freedman, JM Grimshaw, J Gulcher, M Gwinn, MA Hlatky, H Janes,
#' P Kraft, S Melillo, CJ O'Donnell, MJ Pencina, D Ransohoff, SD Schully,
#' D Seminara, DM Winn, CF Wright, CM van Duijn, J Little, MJ Khoury.
#' Strengthening the reporting of genetic risk prediction studies
#' (GRIPS)-Elaboration and explanation. Eur J Epidemiol. 2011;26:313-37. \cr
#'
#' Aulchenko YS, Ripke S, Isaacs A, van Duijn CM. GenABEL: an R package for genome-wide
#' association analysis. Bioinformatics 2007;23(10):1294-6.
#'
"PredictABEL-package" <- function() {}
#' Function to fit a logistic regression model.
#'
#' The function fits a standard GLM function for the logistic regression model.
#' This function can be used to construct a logistic regression model based on genetic and non-genetic
#' predictors. The function also allows to enter the genetic predictors
#' as a single risk score. For that purpose, the function requires that
#' the dataset additionally includes the risk score.
#' A new dataset can be constructed using
#'  "\code{NewExampleData <- cbind(ExampleData,riskScore)}".
#' The genetic risk scores can be obtained
#' using the function \code{\link{riskScore}} in this package or be
#' imported from other methods.
#'
#' @param data Data frame or matrix that includes the outcome and
#' predictor variables.
#' @param cOutcome Column number of the outcome variable. \code{cOutcome=2}
#'  means that the second column of the dataset is the outcome variable.
#' To fit the logistic regression model, the outcome variable needs to be
#' (re)coded as \code{1} for the presence and \code{0} for the absence of the
#' outcome of interest.
#' @param cNonGenPreds Column numbers of the non-genetic predictors that are
#' included in the model. An example to denote column numbers is
#' \code{c(3,6:8,10)}. Choose \code{c(0)} when no non-genetic predictors
#' are considered.
#' @param cNonGenPredsCat Column numbers of the non-genetic predictors that
#' are entered as categorical variables in the model. When non-genetic
#' predictors are not specified as being categorical they are treated as
#' continuous variables in the model. If no non-genetic predictors are
#' categorical, denote \code{c(0)}.
#' @param cGenPreds Column numbers of the genetic predictors or genetic risk score.
#' Denote \code{c(0)}
#' when the prediction model does not consider
#' genetic predictors or genetic risk score.
#' @param cGenPredsCat Column numbers of the genetic predictors that are
#' entered as categorical variables in the model. When SNPs are considered as
#' categorical, the model
#' will estimate effects per genotype. Otherwise, SNPs are considered as
#' continuous variables for which  the model will estimate an allelic effect.
#' Choose c(0) when no genetic predictors are considered as categorical
#' or when genetic predictors are entered as a risk score into the model.
#'
#' @return No value returned.
#'
#' @keywords models
#'
#'
#' @seealso  \code{\link{predRisk}}, \code{\link{ORmultivariate}}, \code{\link{riskScore}}
#' @examples
#' # specify dataset with outcome and predictor variables
#' data(ExampleData)
#' # specify column number of outcome variable
#'  cOutcome <- 2
#'  # specify column numbers of non-genetic predictors
#'  cNonGenPred <- c(3:10)
#'  # specify column numbers of non-genetic predictors that are categorical
#'  cNonGenPredCat <- c(6:8)
#'  # specify column numbers of genetic predictors
#'  cGenPred <- c(11,13:16)
#'  # specify column numbers of genetic predictors that are categorical
#'  cGenPredCat <- c(0)
#'
#'  # fit logistic regression model
#'  riskmodel <- fitLogRegModel(data=ExampleData, cOutcome=cOutcome,
#'  cNonGenPreds=cNonGenPred, cNonGenPredsCat=cNonGenPredCat,
#'  cGenPreds=cGenPred, cGenPredsCat=cGenPredCat)
#'
#'  # show summary details for the fitted risk model
#'  summary(riskmodel)
#'
"fitLogRegModel" <-
function (data, cOutcome, cNonGenPreds, cNonGenPredsCat,
cGenPreds, cGenPredsCat) {

if (missing(cNonGenPreds)) {cNonGenPreds<- c(0)}
if (missing(cNonGenPredsCat)) {cNonGenPredsCat<- c(0)}
if (missing(cGenPreds)) {cGenPreds<- c(0)}
if (missing(cGenPredsCat)) {cGenPredsCat<- c(0)}
if ((cGenPreds[1]==0) & (cNonGenPreds[1]==0)) {
stop("No predictors have been considered.\n")
}

NonGen_factor <- as.data.frame(data[,cNonGenPredsCat])
Gen_factor<- as.data.frame(data[,cGenPredsCat])
if(dim(Gen_factor)[2] >0 )
{
for(i in 1:dim(Gen_factor)[2])
{
Gen_factor[,i] <- as.factor(Gen_factor[,i])
}
}
if(dim(NonGen_factor)[2] >0 )
{
for(i in 1:dim(NonGen_factor)[2])
{
NonGen_factor[,i] <- as.factor(NonGen_factor[,i])
}
}

p <- cbind(data[,setdiff(cNonGenPreds,cNonGenPredsCat)],NonGen_factor,
data[,setdiff(cGenPreds,cGenPredsCat)],Gen_factor)
colnames(p)<- c(colnames(data)[setdiff(cNonGenPreds,cNonGenPredsCat)],colnames(data)[cNonGenPredsCat],
colnames(data)[setdiff(cGenPreds,cGenPredsCat)],colnames(data)[cGenPredsCat])

outcome<- data[,cOutcome]
genes <- colnames(p)
fmla<- as.formula(paste("outcome~", paste(genes, collapse="+")))
model <- glm(fmla, data=p, family=binomial("logit"))
return(model)

}
#' Function to compute predicted risks for all individuals in the dataset.
#'
#' The function computes predicted risks from a specified logistic regression model.
#' The function \code{\link{fitLogRegModel}} can be used to construct such a model.
#'
#' @param riskModel Name of logistic regression model that can be fitted using
#' the function \code{\link{fitLogRegModel}}.
#' @param data Data frame or matrix that includes the ID number and
#' predictor variables.
#' @param cID Column number of ID variable. The ID number and predicted risks
#' will be saved under \code{filename}. When \code{cID} is not specified, the output is not saved.
#' @param filename Name of the output file in which the ID number and
#' estimated predicted risks will be saved. The file is saved in the working
#' directory as a txt file. Example: filename="name.txt". When no \code{filename}
#' is specified, the output is not saved.
#'
#'
#'
#' @return The function returns a vector of predicted risks.
#'
#'
#' @keywords htest
#'
#'
#' @seealso \code{\link{fitLogRegModel}}, \code{\link{plotCalibration}},
#' \code{\link{plotROC}}, \code{\link{plotPriorPosteriorRisk}}
#'
#' @examples
#' # specify dataset with outcome and predictor variables
#' data(ExampleData)
#' # specify column number of the outcome variable
#'  cOutcome <- 2
#' # specify column number of ID variable
#'  cID <- 1
#'  # specify column numbers of non-genetic predictors
#'  cNonGenPred <- c(3:10)
#'  # specify column numbers of non-genetic predictors that are categorical
#'  cNonGenPredCat <- c(6:8)
#'  # specify column numbers of genetic predictors
#'  cGenPred <- c(11,13:16)
#'  # specify column numbers of genetic predictors that are categorical
#'  cGenPredCat <- c(0)
#'
#'  # fit logistic regression model
#'  riskmodel <- fitLogRegModel(data=ExampleData, cOutcome=cOutcome,
#'  cNonGenPreds=cNonGenPred, cNonGenPredsCat=cNonGenPredCat,
#'  cGenPreds=cGenPred, cGenPredsCat=cGenPredCat)
#'
#'  # obtain predicted risks
#'  predRisk <- predRisk(riskModel=riskmodel)
#'
"predRisk" <-
function(riskModel, data, cID, filename)
{
 if (any(class(riskModel) == "glm"))
  {
   predrisk <- predict(riskModel, newdata=data, type="response")
  }
else
  {
stop("The argument 'riskModel' should be a (GLM)model")
  }

  if (!missing(data)&& !missing(cID)&& !missing(filename))
  {tab <- cbind(ID=data[, cID],PredRisk=predrisk)
  write.table(tab, file=filename, row.names = FALSE,sep = "\t")
  }
   return(predrisk)
  }
#' Function to compute genetic risk scores. The function computes unweighted
#' or weighted genetic risk scores. The relative effects (or weights) of
#' genetic variants can either come from beta coefficients of a risk model
#' or from a vector of beta coefficients imported into R, e.g., when beta cofficients are obtained from meta-analysis.
#'
#'
#' The function calculates unweighted
#' or weighted genetic risk scores. The unweighted genetic risk score is a simple
#' risk allele count assuming that all alleles have the same effect. For this
#' calculation, it is required that the genetic variables are coded as the number of risk
#' alleles. Beta coefficients are used to determine which allele is the risk
#' allele. When the sign of the beta coefficient is negative, the allele coding
#' is reversed. The weighted risk score is a sum of the number of risk alleles
#' multiplied by their beta coefficients.
#'
#' The beta coefficients can come from two different sources, either beta coefficients of a risk model
#' or a vector of beta coefficients imported into R, e.g., when beta cofficients are obtained from meta-analysis.
#' This vector of beta coefficients
#' should be a named vector containing the same names as mentioned in genetic variants.
#' A logistic regression model can be constructed using \code{\link{fitLogRegModel}}
#' from this package.
#'
#' @note When a vector of beta coefficients is imported, it should be checked
#' whether the DNA strands and the coding of the risk alleles are the same
#' as in the study data. The functions are available in the package \code{GenABEL}
#' to accurately compute risk scores when the DNA strands are different or the risk
#' alleles are coded differently in the study data and the data used in meta-analysis.
#'
#' @param weights The vector that includes the weights given to the genetic
#' variants. See details for more informations.
#' @param data Data frame or matrix that includes the outcome
#' and predictors variables.
#' @param cGenPreds Column numbers of the genetic variables on the basis of
#'  which the risk score is computed.
#' @param Type Specification of the type of risk scores that will be computed.
#' Type can be weighted (\code{Type="weighted"}) or
#' unweighted (\code{Type="unweighted"}).
#'
#' @return The function returns a vector of risk scores.
#'
#'
#' @keywords htest
#'
#'
#' @seealso \code{\link{plotRiskDistribution}}, \code{\link{plotRiskscorePredrisk}}
#' @examples
#' # specify dataset with outcome and predictor variables
#'  data(ExampleData)
#'  # specify column numbers of genetic predictors
#'  cGenPred <- c(11:16)
#'
#'  # fit a logistic regression model
#'  # all steps needed to construct a logistic regression model are written in a function
#'  # called 'ExampleModels', which is described on page 4-5
#'  riskmodel <- ExampleModels()$riskModel2
#'
#'  # compute unweighted risk scores
#'  riskScore <- riskScore(weights=riskmodel, data=ExampleData,
#' cGenPreds=cGenPred, Type="unweighted")
#'
"riskScore" <-
function(weights, data, cGenPreds, Type )
{
riskModel <- weights
x <- data[, cGenPreds]
if (any(class(riskModel) == "glm"))
  {
    if(! setequal(intersect(names(riskModel$coef),colnames(x)),colnames(x)))
    {
     stop("The risk model does not contain all the genetic variants as predictors")
    }
     else
    {
     y <- riskModel$coef[intersect(names(riskModel$coef),colnames(x))]
    }
  }

else if(is.vector(riskModel))
 {
   if (length(names(riskModel))!= length(riskModel))
   {
    stop("The argument 'weights' is not a named beta (coefficient) vector")
   }
   if( setequal(intersect(names(riskModel),colnames(x)),colnames(x)))
   {
    y <- riskModel[intersect(colnames(x),names(riskModel))]
   }
   else
   {
 stop("Beta coefficient vector does not contain all the genetic variants")
   }

 }

else
   {
stop("'weights' argument should either be a model or a named beta vector" )
   }

if(Type=="weighted")
{
wrs<- y %*% t(x)
return(as.vector(wrs))
}
else if (Type=="unweighted")
{
Unbetavector <- sign(y)
pro <-   Unbetavector %*% t(x)
num <-sum(Unbetavector==-1)
urs<-pro+(2*num)
return(as.vector(urs))
}
}
#' Function to compute univariate ORs for genetic predictors. The function computes the univariate ORs with 95\% CIs for genetic predictors.
#'
#' The function computes the univariate ORs with 95\% CIs for the specified
#' genetic variants both per allele and per genotype. The ORs are saved with the data from which they are
#' calculated. Genotype frequencies are provided for
#' persons with and without the outcome
#' of interest. The genotype or allele that is coded as \code{'0'} is considered
#' as the reference to computes the ORs.
#'
#' @param data Data frame or matrix that includes the outcome and
#' predictors variables.
#' @param cOutcome Column number of the outcome variable. \code{cOutcome=2}
#' means that the second column of the dataset is the outcome variable.
#' @param cGenPreds Column numbers of genetic variables for which the ORs
#' are calculated.
#' @param filenameGeno Name of the output file in which the univariate ORs
#' and frequencies per genotype will be saved. The file is saved in the working directory as
#' a txt file. When no \code{filenameGeno} is specified, the output is not saved.
#' @param filenameAllele Name of the output file in which the univariate ORs and
#' frequencies per allele will be saved. The file is saved in the working
#' directory as a txt file. When no \code{filenameAllele} is specified, the output is not saved.
#'
#' @return The function returns two different tables. One table contains genotype frequencies
#' and univariate ORs with 95\% CIs and the other contains allele frequencies and
#' univariate ORs with 95\% CIs.
#'
#'
#' @keywords manip
#'
#'
#' @seealso \code{\link{ORmultivariate}}
#' @examples
#' # specify dataset with outcome and predictor variables
#' data(ExampleData)
#' # specify column number of the outcome variable
#'  cOutcome <- 2
#'  # specify column numbers of genetic predictors
#'  cGenPreds <- c(11:13,16)
#'
#'  # compute univariate ORs
#'  ORunivariate(data=ExampleData, cOutcome=cOutcome, cGenPreds=cGenPreds,
#' filenameGeno="GenoOR.txt", filenameAllele="AlleleOR.txt")
#'
"ORunivariate" <- function(data, cOutcome,cGenPreds,filenameGeno, filenameAllele )
  {
  p<- data[,cGenPreds ] # p : A table of Genotype data for all SNP's
  o<- data[,cOutcome] #  o: the binary outcome variable
  m1 <- matrix (nrow=dim(p)[2], ncol=19)
  n1 <- matrix (nrow=dim(p)[2], ncol=12)

  for (i in 1:dim(p)[2])
	 {
s<-table(p[,i],o)
if (dim(s)[1]==1) {s <- rbind(s,c(0,0));s <- rbind(s,c(0,0))}
if (dim(s)[1]==2) {s <- rbind(s,c(0,0))}
a<-oddsratio.wald(s)$measure
b<-oddsratio.wald(s)$data
c <- matrix (nrow=2, ncol=2)
c[1,1] <- 2*b[1,1]+b[2,1]
c[2,1] <- 2*b[3,1]+b[2,1]
c[1,2] <- 2*b[1,2]+b[2,2]
c[2,2] <- 2*b[3,2]+b[2,2]
dimnames(c)[[1]] <- c("Allele-I","Allele-II")
dimnames(c)[[2]] <- c("0","1")
d<-oddsratio.wald(c)$measure
e<-oddsratio.wald(c)$data

	    m1[i,1]<-  colnames(p)[i]
	    m1[i,2]<-  ( b[1,2])
	    m1[i,3]<-  ( round((b[1,2]/b[4,2])*100 ,1))
	    m1[i,4]<-  ( b[2,2])
	    m1[i,5]<-  ( round((b[2,2]/b[4,2])*100 ,1))
	    m1[i,6]<-  ( b[3,2])
	    m1[i,7]<-  ( round((b[3,2]/b[4,2])*100 ,1))
	    m1[i,8]<-  ( b[1,1])
	    m1[i,9]<-  ( round((b[1,1]/b[4,1])*100 ,1))
	    m1[i,10]<- ( b[2,1])
	    m1[i,11]<- ( round((b[2,1]/b[4,1])*100 ,1))
	    m1[i,12]<- ( b[3,1])
	    m1[i,13]<- ( round((b[3,1]/b[4,1])*100 ,1))
	    m1[i,14]<- ( round(a[2,1] ,2))
	    m1[i,15]<- ( round(a[2,2] ,2))
	    m1[i,16]<- ( round(a[2,3] ,2))
	    m1[i,17]<- ( round(a[3,1] ,2))
	    m1[i,18]<- ( round(a[3,2] ,2))
	    m1[i,19]<- ( round(a[3,3] ,2))

	    n1[i,1]<-  colnames(p)[i]
	    n1[i,2]<-  (e[1,2])
	    n1[i,3]<-  round((e[1,2]/e[3,2])*100,1)
	    n1[i,4]<-  e[2,2]
	    n1[i,5]<-  round((e[2,2]/e[3,2])*100,1)
	    n1[i,6]<-   e[1,1]
	    n1[i,7]<-   round((e[1,1]/e[3,1])*100,1)
	    n1[i,8]<-   e[2,1]
	    n1[i,9]<-   round((e[2,1]/e[3,1])*100,1)
	    n1[i,10]<-  round(d[2,1] ,2)
	    n1[i,11]<-  round(d[2,2] ,2)
	    n1[i,12]<-  round(d[2,3] ,2)

	}

       m1 <- as.table(m1)
       dimnames(m1)[[1]] <- c(1:dim(p)[2])
       dimnames(m1)[[2]] <- c("Name", "0 ","0 (%)", "1 ","1 (%)", "2 ","2 (%)",
       "0 ","0 (%)", "1 ","1 (%)", "2 ","2 (%)","OR1","CI-Low","CI-high","OR2",
       "CI-Low","CI-high")
names(dimnames(m1)) <- c("", " Genotype frequencies for Cases & Controls, and OR(95% CI)")
     if (!missing(filenameGeno))
		 write.table(m1,file=filenameGeno, row.names = FALSE,sep = "\t")

       n1 <- as.table(n1)
       dimnames(n1)[[1]] <- c(1:dim(p)[2])
       dimnames(n1)[[2]] <- c("Name", "0","0 (%)","1","1 (%)", "0","0 (%)","1",
       "1 (%)","OR1","CI-Low","CI-high")
names(dimnames(n1)) <- c("", " Allele frequencies for Cases & Controls, and OR(95% CI) ")
       if (!missing(filenameAllele))
       write.table(n1,file=filenameAllele, row.names = FALSE,sep = "\t")

       p<- list(Genotype =m1,  Allelic=n1)
	return(p)
  }
#' Function to obtain multivariate odds ratios from a logistic regression model.
#' The function estimates multivariate (adjusted) odds ratios (ORs) with
#' 95\% confidence intervals (CIs) for all the genetic and non-genetic variables
#' in the risk model.
#'
#' The function requires that first a logistic regression
#' model is fitted either by using \code{GLM} function or the function
#' \code{\link{fitLogRegModel}}. In addition to the multivariate ORs,
#' the function returns summary statistics of model performance, namely the Brier
#' score and the Nagelkerke's \eqn{R^2} value.
#' The Brier score quantifies the accuracy of risk predictions by comparing
#' predicted risks with observed outcomes at individual level (where outcome
#' values are either 0 or 1). The Nagelkerke's \eqn{R^2} value indicates the percentage of variation
#'  of the outcome explained by the predictors in the model.
#'
#' @param riskModel Name of logistic regression model that can be fitted using
#' the function \code{\link{fitLogRegModel}}.
#' @param filename Name of the output file in which the multivariate
#'  ORs will be saved. If no directory is specified, the file is
#' saved in the working directory as a txt file.
#' When \code{filename} is not specified, the output is not saved.
#'
#'
#'  @return
#'   The function returns:
#'   \item{Predictors Summary}{OR with 95\% CI and corresponding  p-values for
#' each predictor in the model}
#'   \item{Brier Score}{Brier score}
#'   \item{Nagelkerke Index}{Nagelkerke's \eqn{R^2} value}
#'
#'
#'
#' @keywords htest
#'
#'
#' @references Brier GW. Verification of forecasts expressed in terms of probability.
#' Monthly weather review 1950;78:1-3.
#'
#'
#' Nagelkerke NJ. A note on a general definition of the coefficient
#' of determination. Biometrika 1991;78:691-692.
#'
#'
#' @seealso \code{\link{fitLogRegModel}}, \code{\link{ORunivariate}}
#'
#' @examples
#' # specify dataset with outcome and predictor variables
#' data(ExampleData)
#' # specify column number of outcome variable
#'  cOutcome <- 2
#'  # specify column numbers of non-genetic predictors
#'  cNonGenPred <- c(3:10)
#'  # specify column numbers of non-genetic predictors that are categorical
#'  cNonGenPredCat <- c(6:8)
#'  # specify column numbers of genetic predictors
#'  cGenPred <- c(11,13:16)
#'  # specify column numbers of genetic predictors that are categorical
#'  cGenPredCat <- c(0)
#'
#'  # fit logistic regression model
#'  riskmodel <- fitLogRegModel(data=ExampleData, cOutcome=cOutcome,
#'  cNonGenPreds=cNonGenPred, cNonGenPredsCat=cNonGenPredCat,
#'  cGenPreds=cGenPred, cGenPredsCat=cGenPredCat)
#'
#'  # obtain multivariate OR(95\% CI) for all predictors of the fitted model
#'   ORmultivariate(riskModel=riskmodel, filename="multiOR.txt")
#'
"ORmultivariate" <- function(riskModel, filename)
{
if (!any(class(riskModel)=="glm"))
stop("riskModel argument should have 'glm' class")
sum.coef<-summary(riskModel)$coef
OR<-exp(sum.coef[,"Estimate"])
Upper.CI<-exp(sum.coef[,"Estimate"]+1.96*sum.coef[,"Std. Error"])
Lower.CI<-exp(sum.coef[,1]-1.96*sum.coef[,2])
tab <-cbind("OR"=round(OR ,4),"Lower.CI"=round(Lower.CI ,4),"Upper.CI"=
round(Upper.CI ,4),"p-value"=round((sum.coef)[, 4],4))

if (!missing(filename))
	write.table(tab,file=filename, row.names=TRUE,sep = "\t")

B <- mean((riskModel$y) * (1-predict(riskModel, type="response"))^2 +
(1-riskModel$y) * (predict(riskModel, type="response"))^2)
# B
Bmax <- mean(riskModel$y) * (1-mean(riskModel$y))^2 +
(1-mean(riskModel$y)) *mean(riskModel$y)^2
Bscaled <- 1 - B/Bmax

LLfull <- (riskModel$deviance)/(-2)
LLnul <- (riskModel$null.deviance)/(-2)
n <- length(riskModel$y)
Mc <- 1-(LLfull/LLnul)   # McFadden's R square
Cox <- 1-exp ((-2)*(LLfull-LLnul)/n)  # Cox-Snell R square
Rmax <- 1-exp((2*LLnul)/n)
Nag <- Cox/Rmax    # Nagelkerke_R2

  p<- list(Predictors_Summary= tab,Brier_Score=round(B ,4),Nagelkerke_R2=round(Nag,4))
 return(p)
}
#' Function to plot predicted risks against risk scores.
#' This function is used to make a plot of predicted risks against risk scores.
#'
#' The function creates a plot of predicted risks against risk scores.
#' Predicted risks can be obtained using the functions
#' \code{\link{fitLogRegModel}} and \code{\link{predRisk}}
#' or be imported from other methods or packages.
#' The function \code{\link{riskScore}} can be
#'  used to compute unweighted or weighted risk scores.
#'
#' @param data Data frame or matrix that includes the outcome and
#' predictors variables.
#' @param riskScore Vector of (weighted or unweighted) genetic risk scores.
#' @param predRisk Vector of predicted risks.
#' @param plottitle  Title of the plot. Specification of \code{plottitle} is optional. Default is "Risk score predicted risk plot".
#' @param xlabel  Label of x-axis. Specification of \code{xlabel} is optional. Default is "Risk score".
#' @param ylabel  Label of y-axis. Specification of \code{ylabel} is optional. Default is "Predicted risk".
#' @param rangexaxis Range of the x axis. Specification of \code{rangexaxis} is optional.
#' @param rangeyaxis Range of the y axis. Specification of \code{rangeyaxis} is optional. Default is \code{c(0,1)}.
#' @param filename Name of the output file in which risk scores and
#' predicted risks for each individual will be saved. If no directory is
#' specified, the file is saved in the working directory as a txt file.
#' When no \code{filename} is specified, the output is not saved.
#' @param fileplot Name of the output file that contains the plot. The file is
#' saved in the working directory in the format specified under \code{plottype}. Example:
#' \code{fileplot="plotname"}. Note that the extension is not specified here.
#' When \code{fileplot} is not specified, the plot is not saved.
#' @param plottype The format in which the plot is saved. Available formats are
#' wmf, emf, png, jpg, jpeg, bmp, tif, tiff, ps,
#' eps or pdf. For example, \code{plottype="eps"} will save the plot in eps format.
#' When \code{plottype} is not specified, the plot will be saved in jpg format.
#'
#' @return
#'  The function creates a plot of predicted risks against risk scores.
#'
#'
#' @keywords hplot
#'
#'
#' @seealso \code{\link{riskScore}}, \code{\link{predRisk}}
#'
#' @examples
#' # specify dataset with outcome and predictor variables
#' data(ExampleData)
#'
#'  # fit a logistic regression  model
#'  # all steps needed to construct a logistic regression model are written in a function
#'  # called 'ExampleModels', which is described on page 4-5
#'  riskmodel <- ExampleModels()$riskModel2
#'
#'  # obtain predicted risks
#'  predRisk <- predRisk(riskmodel)
#'
#'  # specify column numbers of genetic predictors
#'  cGenPred <- c(11:16)
#'
#'  # function to compute unweighted genetic risk scores
#'  riskScore <- riskScore(weights=riskmodel, data=ExampleData,
#'  cGenPreds=cGenPred, Type="unweighted")
#'
#'  # specify range of x-axis
#'  rangexaxis <- c(0,12)
#'  # specify range of y-axis
#'  rangeyaxis <- c(0,1)
#'  # specify label of x-axis
#'  xlabel <- "Risk score"
#'  # specify label of y-axis
#'  ylabel <- "Predicted risk"
#'  # specify title for the plot
#'  plottitle <- "Risk score versus predicted risk"
#'
#'  # produce risk score-predicted risk plot
#' plotRiskscorePredrisk(data=ExampleData, riskScore=riskScore, predRisk=predRisk,
#' plottitle=plottitle, xlabel=xlabel, ylabel=ylabel, rangexaxis=rangexaxis,
#' rangeyaxis=rangeyaxis, filename="RiskscorePredRisk.txt")
#'
"plotRiskscorePredrisk" <-
function(data, riskScore, predRisk, plottitle, xlabel,
ylabel, rangexaxis, rangeyaxis, filename, fileplot, plottype)
  {
  if (missing(plottitle)) {plottitle <- "Risk score predicted risk plot"}
  if (missing(xlabel)) {xlabel<- "Risk score"}
  if (missing(ylabel)) {ylabel<- "Predicted risk"}
  if (missing(rangexaxis)) {rangexaxis<- c(min(riskScore), max(riskScore))}
  if (missing(rangeyaxis)) {rangeyaxis<- c(0,1)}
plot(riskScore,predRisk,main=plottitle,xlab=xlabel,ylab=ylabel,ylim=rangeyaxis,
  xlim = rangexaxis, pch=21, cex.lab= 1.2, cex.axis= 1.1, las=1)
  b <- lm(predRisk~riskScore)
  abline(b,col = 8)

 if (missing(plottype)) {plottype<- "jpg"}
 	if (!missing(fileplot))
      savePlot(filename = fileplot,
	 type =plottype,
	 device = dev.cur())
  tab<- cbind(riskScore,predRisk)
  tab <- as.table(tab)
  dimnames(tab)[[2]] <- c("Risk Score", "Predicted risk ")
  if (!missing(filename))
	write.table(tab,file=filename, row.names=TRUE,sep = "\t")
  }
#' Function to plot posterior risks against prior risks.
#'
#' The function creates a plot of posterior risks (predicted risks using
#' the updated model) against prior risks (predicted risks using the initial
#' model). Predicted risks can be obtained using the functions
#' \code{\link{fitLogRegModel}} and \code{\link{predRisk}} or be
#' imported from other packages or methods.
#'
#' @param data Data frame or matrix that includes the outcome and
#' predictors variables.
#' @param priorrisk Vector of predicted risks based on initial model.
#' @param posteriorrisk  Vector of predicted risks based on updated model.
#' @param cOutcome  Column number of the outcome variable.
#' @param plottitle  Title of the plot. Specification of \code{plottitle} is optional. Default is "PriorPosteriorRisk plot".
#' @param xlabel  Label of x-axis. Specification of \code{xlabel} is optional. Default is "Prior risk".
#' @param ylabel  Label of y-axis. Specification of \code{ylabel} is optional. Default is "Posterior risk".
#' @param rangeaxis  Range of x-axis and y-axis. Specification of \code{rangeaxis} is optional. Default is \code{c(0,1)}.
#' @param  plotAll  \code{plotAll=TRUE} will create one plot for the
#' total population. When  \code{plotAll=FALSE} separate plots will be created
#' for individuals with and without the outcome of interest.
#'  means two separate plots for with and without outcome of interest.
#' @param labels Labels given to the groups of individuals without and with
#' the outcome of interest. Default \code{labels} is
#'  \code{c("without outcome", "with outcome")}. Note that when
#' \code{plotAll=TRUE}, specification of \code{labels} is not necessary.
#' @param  filename Name of the output file in which prior and posterior
#' risks for each individual with the outcome will be saved. If no directory is
#' specified, the file is saved in the working directory as a txt file.
#' When no \code{filename} is specified, the output is not saved.
#' @param fileplot Name of the output file that contains the plot. The file is
#' saved in the working directory in the format specified under \code{plottype}. Example:
#' \code{fileplot="plotname"}. Note that the extension is not specified here.
#' When \code{fileplot} is not specified, the plot is not saved.
#' @param plottype The format in which the plot is saved. Available formats are
#' wmf, emf, png, jpg, jpeg, bmp, tif, tiff, ps,
#' eps or pdf. For example, \code{plottype="eps"} will save the plot in eps format.
#' When \code{plottype} is not specified, the plot will be saved in jpg format.
#'
#' @return
#'  The function creates a plot of posterior risks against prior risks.
#'
#'
#' @keywords hplot
#'
#'
#' @seealso \code{\link{predRisk}}
#' @examples
#' # specify dataset with outcome and predictor variables
#' data(ExampleData)
#' # specify column number of outcome variable
#'  cOutcome <- 2
#'
#'  # fit logistic regression models
#'  # all steps needed to construct a logistic regression model are written in a function
#'  # called 'ExampleModels', which is described on page 4-5
#'  riskmodel1 <- ExampleModels()$riskModel1
#'  riskmodel2 <- ExampleModels()$riskModel2
#'
#'  # obtain predicted risks
#'  predRisk1 <- predRisk(riskmodel1)
#'  predRisk2 <- predRisk(riskmodel2)
#'
#'  # specify label of x-axis
#'  xlabel <- "Prior risk"
#'  # specify label of y-axis
#'  ylabel <- "Posterior risk"
#'  # specify title for the plot
#'  titleplot <- "Prior versus posterior risk"
#'  # specify range of the x-axis and y-axis
#'  rangeaxis <- c(0,1)
#'  # labels given to the groups without and with the outcome of interest
#'  labels<- c("without outcome", "with outcome")
#'
#' # produce prior risks and posterior risks plot
#'  plotPriorPosteriorRisk(data=ExampleData, priorrisk=predRisk1,
#'  posteriorrisk=predRisk2, cOutcome=cOutcome, xlabel=xlabel, ylabel=ylabel,
#'  rangeaxis=rangeaxis, plotAll=TRUE, plottitle=titleplot, labels=labels)
#'
"plotPriorPosteriorRisk" <-
function( data, priorrisk, posteriorrisk, cOutcome, plottitle, xlabel, ylabel,
rangeaxis, plotAll=TRUE, labels, filename, fileplot, plottype)
{
  if (missing(plottitle)) {plottitle <- "PriorPosteriorRisk plot"}
  if (missing(xlabel)) {xlabel<- "Prior risk"}
  if (missing(ylabel)) {ylabel<- "Posterior risk"}
  if (missing(rangeaxis)) {rangeaxis<- c(0,1)}
  if (missing(labels)) {labels<- c("without outcome", "with outcome")}
if((min(priorrisk)<0)|(max(priorrisk)>1))
{stop("Not all predicted risks are lie between 0 to 1")}
if((min(posteriorrisk)<0)|(max(posteriorrisk)>1))
{stop("Not all predicted risks are lie between 0 to 1")}
if(!(length(priorrisk)==length(posteriorrisk)))
{stop("prior and posterior risks have different length")}

  risk1 <- priorrisk
  risk2 <- posteriorrisk

    if( plotAll)
	{
	plot(risk1,risk2,xlab= xlabel, ylab=ylabel,main=plottitle,
	cex.lab=1.2, cex.axis=1.1, las=1,pty='s',xlim=rangeaxis ,
	ylim=rangeaxis,pch=20)
	abline(a=0,b=1, lwd=1,col=8)
	       if (missing(plottype)) {plottype<- "jpg"}
 	       if (!missing(fileplot))
	 savePlot(filename = fileplot,
	 type =plottype,
	 device = dev.cur()
		  )
       }
       else
       {
	op <- par(mfrow=c(1,2),pty="s" )
	plot(risk1,risk2,xlab= xlabel, ylab=ylabel,
	col = (1-(data[,cOutcome]))*1,cex.lab=1.2, cex.axis=1, las=1,pty='s',
	xlim=rangeaxis , ylim=rangeaxis,pch="*")
	abline(a=0,b=1, lwd=1,col=4)
	title(labels[1],cex.main=1)

	plot(risk1,risk2, xlab= xlabel, ylab=ylabel,
	col =(data[,cOutcome])*1, cex.lab=1.2, cex.axis=1, las=1,pty='s',
	xlim=rangeaxis , ylim=rangeaxis,pch="*")
	abline(a=0,b=1, lwd=1,col=4)
	title(labels[2], cex.main=1)
	par(op)
	if (missing(plottype)) {plottype<- "jpg"}
 	      if (!missing(fileplot))
	savePlot(filename = fileplot,
	 type =plottype,
	 device = dev.cur())
       }

	 tab<- cbind(risk1,risk2,data[,cOutcome])
	 tab <- as.table(tab)
       dimnames(tab)[[2]] <- c("Predicted risk 1","Predicted risk 2", "outcome")
	 if (!missing(filename))
	       write.table(tab,file=filename, row.names=TRUE,sep = "\t")
  }
#' Function for calibration plot and Hosmer-Lemeshow goodness of fit test.
#' The function produces a calibration plot and provides Hosmer-Lemeshow
#' goodness of fit test statistics.
#'
#' Hosmer-Lemeshow test statistic is a measure of the fit
#' of the model, comparing observed and predicted risks across subgroups of
#' the population. The default number of groups is 10.
#'
#' The function requires the outcome of interest and predicted risks of
#' all individuals. Predicted risks can be obtained from the
#' functions \code{\link{fitLogRegModel}} and \code{\link{predRisk}} or
#' be imported from other packages or methods.
#'
#' @param data Data frame or numeric matrix that includes the outcome and
#' predictor variables.
#' @param cOutcome  Column number of the outcome variable.
#' @param predRisk Vector of predicted risks of all individuals in the dataset.
#' @param groups  Number of groups considered in
#' Hosmer-Lemeshow test. Specification of \code{groups} is optional (default \code{groups} is 10).
#' @param rangeaxis  Range of x-axis and y-axis. Specification of \code{rangeaxis} is optional. Default is \code{c(0,1)}.
#' @param plottitle  Title of the plot. Specification of \code{plottitle} is optional. Default is "Calibration plot".
#' @param xlabel  Label of x-axis Default. Specification of \code{xlabel} is optional. Default is "Predicted risk".
#' @param ylabel  Label of y-axis. Specification of \code{ylabel} is optional. Default is "Observed risk".
#' @param filename Name of the output file in which the calibration table is saved.
#' The file is saved as a txt file in the working directory. When  no
#' \code{filename} is specified, the output is not saved. Example: filename="calibration.txt"
#' @param fileplot Name of the file that contains the calibation plot.
#'  The file is saved in the working directory in the format specified under \code{plottype}. Example:
#' \code{fileplot="plotname"}. Note that the extension is not specified here.
#' When \code{fileplot} is not specified, the plot is not saved.
#' @param plottype The format in which the plot is saved. Available formats are
#' wmf, emf, png, jpg, jpeg, bmp, tif, tiff, ps,
#' eps or pdf. Foe example, \code{plottype="eps"} will save the plot in eps format.
#' When \code{plottype} is not specified, the plot will be saved in jpg format.
#'
#'  @return
#'   The function creates a calibration plot and returns the following measures:
#'   \item{Chi_square}{Chi square value of  Hosmer-Lemeshow test}
#'   \item{df}{Degrees of freedom, which is \code{(groups-2)} where \code{groups:} number
#'  of groups}
#'   \item{p_value}{p-value of Hosmer-Lemeshow test for goodness of fit}
#'
#'
#' @keywords htest hplot
#'
#'
#' @references Hosmer DW, Hosmer T, Le Cessie S, Lemeshow S. A comparison of
#' goodness-of-fit tests for the logistic regression model. Stat Med 1997;
#' 16:965-980.
#'
#' @seealso \code{\link{predRisk}}
#'
#' @examples
#' # specify dataset with outcome and predictor variables
#' data(ExampleData)
#' # specify column number of the outcome variable
#'  cOutcome <- 2
#'
#'  # fit a logistic regression model
#'  # all steps needed to construct a logistic regression model are written in a function
#'  # called 'ExampleModels', which is described on page 4-5
#'  riskmodel <- ExampleModels()$riskModel2
#'
#'  # obtain predicted risks
#'  predRisk <- predRisk(riskmodel)
#'
#'  # specify range of x-axis and y-axis
#'  rangeaxis <- c(0,1)
#' # specify number of groups for Hosmer-Lemeshow test
#'  groups <- 10
#'
#' # compute calibration measures and produce calibration plot
#'  plotCalibration(data=ExampleData, cOutcome=cOutcome, predRisk=predRisk,
#'  groups=groups, rangeaxis=rangeaxis)
#'
"plotCalibration" <-
function( data, cOutcome, predRisk, groups, rangeaxis, plottitle,
xlabel, ylabel, filename, fileplot, plottype)
{
   if (missing(groups)) {groups<- 10}
  if (missing(plottitle)) {plottitle <- "Calibration plot"}
  if (missing(xlabel)) {xlabel<- "Predicted risk"}
  if (missing(ylabel)) {ylabel<- "Observed risk"}
   if (missing(rangeaxis)) {rangeaxis<- c(0,1)}
 p=predRisk
 y=data[,cOutcome]
 if (length(unique(y))!=2) {
			    stop(" The specified outcome is not a binary variable.\n")
			    }
else{

matres	<-matrix(NA,nrow=groups,ncol=5)
sor	<-order(p)
p	<-p[sor]
y	<-y[sor]
groep	<-cut2(p,g=groups)
total		<-tapply(y,groep,length)
predicted	<-round(tapply(p,groep,sum),2)
observed	<-tapply(y,groep,sum)
meanpred	<-round(tapply(p,groep,mean),3)
meanobs	<-round(tapply(y,groep,mean),3)
matres	<-cbind(total,meanpred,meanobs, predicted, observed)
plot(matres[,2],matres[,3],main=plottitle,xlab=xlabel,ylab=ylabel,
pch=16,ps=2, xlim=rangeaxis , ylim=rangeaxis,cex.lab=1.2,
cex.axis=1.1, las=1)
contr<-((observed-predicted)^2)/(total*meanpred*(1-meanpred))
chisqr<-sum(contr)
df<- (groups-2)
pval<- 1-pchisq(chisqr,df)
lines(x =c(0,1),y=c(0,1))
 if (missing(plottype)) {plottype<- "jpg"}
 	if (!missing(fileplot))
      savePlot(filename = fileplot,
	 type =plottype,
	 device = dev.cur())

if (!missing(filename))
	write.table(matres,file=filename, row.names=TRUE,sep = "\t",dec=",")
  out <- list(Table_HLtest=matres,Chi_square = round(chisqr,3), df=df,
  p_value =round(pval,4))
  return(out)

    }
    }
#' Function for a receiver operating characteristic curve (ROC) plot and area under the ROC curve (AUC) value.
#' The function produces ROC curve and corresponding AUC value with 95\% CI.
#' The function can plot one or  multiple ROC curves in a single plot.
#'
#' The function requirs predicted risks or risk scores and the outcome of
#' interest for all individuals.
#' Predicted risks can be obtained using the functions
#' \code{\link{fitLogRegModel}} and \code{\link{predRisk}}
#' or be imported from other methods or packages.
#'
#' @param data Data frame or matrix that includes the outcome and
#' predictors variables.
#' @param cOutcome  Column number of the outcome variable.
#' @param predrisk  Vector of predicted risk. When multiple curves need to
#' be presented in one plot, specify multiple vectors of predicted
#' risks as \code{predrisk=cbind(predrisk1, predrisk2,...,predriskn)}.
#' @param  labels  Label(s) given to the ROC curve(s). Specification of \code{labels} is optional.
#' When specified, the \code{labels} should be in the same order as specified in \code{predrisk}.
#' @param plottitle  Title of the plot. Specification of \code{plottitle} is optional. Default is "ROC plot".
#' @param xlabel  Label of x-axis. Specification of \code{xlabel} is optional. Default is "1- Specificity".
#' @param ylabel  Label of y-axis. Specification of \code{ylabel} is optional. Default is "Sensitivity".
#' @param fileplot Name of the output file that contains the plot. The file is
#' saved in the working directory in the format specified under \code{plottype}. Example:
#' \code{fileplot="plotname"}. Note that the extension is not specified here.
#' When \code{fileplot} is not specified, the plot is not saved.
#' @param plottype The format in which the plot is saved. Available formats are
#' wmf, emf, png, jpg, jpeg, bmp, tif, tiff, ps,
#' eps or pdf. For example, \code{plottype="eps"} will save the plot in eps format.
#' When \code{plottype} is not specified, the plot will be saved in jpg format.
#'
#'  @return
#'   The function creates ROC plot and returns AUC value with 95\% CI.
#'
#'
#' @keywords htest hplot
#'
#'
#' @references Hanley JA, McNeil BJ. The meaning and use of the area under a
#' receiver operating characteristic (ROC) curve. Radiology 1982;143:29-36.
#'
#'
#' Tobias Sing, Oliver Sander, Niko Beerenwinkel, Thomas Lengauer.
#' ROCR: visualizing classifier performance in R.
#' Bioinformatics 2005;21(20):3940-3941.
#'
#' @seealso \code{\link{predRisk}}, \code{\link{plotRiskDistribution}}
#' @examples
#' # specify the arguments in the function to produce ROC plot
#' # specify dataset with outcome and predictor variables
#' data(ExampleData)
#' # specify column number of the outcome variable
#'  cOutcome <- 2
#'
#'  # fit logistic regression models
#'  # all steps needed to construct a logistic regression model are written in a function
#'  # called 'ExampleModels', which is described on page 4-5
#'  riskmodel1 <- ExampleModels()$riskModel1
#'  riskmodel2 <- ExampleModels()$riskModel2
#'
#'  # obtain predicted risks
#'  predRisk1 <- predRisk(riskmodel1)
#'  predRisk2 <- predRisk(riskmodel2)
#'
#'  # specify label of the ROC curve
#'  labels <- c("without genetic factors", "with genetic factors")
#'
#' # produce ROC curve
#'  plotROC(data=ExampleData, cOutcome=cOutcome,
#'  predrisk=cbind(predRisk1,predRisk2), labels=labels)
#'
"plotROC" <- function( data, cOutcome, predrisk, labels, plottitle,
xlabel, ylabel, fileplot, plottype)
{
   if (missing(plottitle)) {plottitle <- "ROC plot"}
  if (missing(xlabel)) {xlabel<- "1- Specificity"}
  if (missing(ylabel)) {ylabel<- "Sensitivity"}
 if (class(predrisk) == "numeric") {predrisk<- cbind(predrisk)}
  a<-c(1:dim(predrisk)[2])
 	for(i in 1:dim(predrisk)[2])
		{
  rAllele <- rcorr.cens(predrisk[,i], data[,cOutcome], outx=FALSE)
  pred <- prediction(predrisk[,i], data[,cOutcome])
  perf <- performance(pred,"tpr","fpr")

   if (i==1)
   {
plot(perf,xlab=xlabel, ylab=ylabel,col=16+i,lty=i,las=1,lwd =2,
cex.lab=1.2,cex.axis=1.1, main =plottitle )

   }
     else
   {
plot(perf,xlab=xlabel, ylab=ylabel,add=TRUE,col=16+i,lty=i,las=1,
lwd =2 ,cex.lab=1.2,cex.axis=1.1, main = plottitle)

   }
  lines(x=c(0,1), y=c(0,1), lwd=1,col=8)
 cat("AUC [95% CI] for the model",i, ": ", round(rAllele[1],3),
 "[", round(rAllele[1]-1.96/2*rAllele[3],3)," - ",
 round(rAllele[1]+1.96/2*rAllele[3],3), "] \n")
  	}
	  if (!missing(labels)){
  	legend( "bottomright",legend= labels, col=c(17:(16+dim(predrisk)[2])),
    lty=c(1:(dim(predrisk)[2])),lwd =2,cex=1)
    }
 if (missing(plottype)) {plottype<- "jpg"}
 	if (!missing(fileplot))
      savePlot(filename = fileplot,
	 type =plottype,
	 device = dev.cur())

}
#' Function for predictiveness curve. The function creates a plot of cumulative percentage
#' of individuals to the predicted risks.
#'
#' The Predictiveness curve is a plot of cumulative percentage
#' of individuals to the predicted risks. Cumulative percentage indicates
#' the percentage of individual that has a predicted risk equal or lower
#' than the risk value.
#' Predicted risks can be obtained using the functions
#' \code{\link{fitLogRegModel}} and \code{\link{predRisk}}
#' or be imported from other methods or packages.
#'
#'
#' @param predrisk  Vector of predicted risk. When multiple curves need to
#' be presented in one plot, specify multiple vectors of predicted
#' risks as \code{predrisk=cbind(predrisk1, predrisk2,...,predriskn)}.
#' @param  rangeyaxis  Range of the y axis. Default \code{rangeyaxis} is \code{c(0,1)}.
#' @param  labels  Label(s) given to the predictiveness curve(s). Specification of \code{labels} is optional.
#' When specified, the \code{labels} should be in the same order as specified in \code{predrisk}.
#' @param plottitle  Title of the plot. Specification of \code{plottitle} is optional. Default is "Predictiveness curve".
#' @param xlabel  Label of x-axis. Specification of \code{xlabel} is optional. Default is "Cumulative percentage".
#' @param ylabel  Label of y-axis. Specification of \code{ylabel} is optional. Default is "Predicted risks".
#' @param fileplot Name of the output file that contains the plot. The file is
#' saved in the working directory in the format specified under \code{plottype}. Example:
#' \code{fileplot="plotname"}. Note that the extension is not specified here.
#' When \code{fileplot} is not specified, the plot is not saved.
#' @param plottype The format in which the plot is saved. Available formats are
#' wmf, emf, png, jpg, jpeg, bmp, tif, tiff, ps,
#' eps or pdf. For example, \code{plottype="eps"} will save the plot in eps format.
#' When \code{plottype} is not specified, the plot will be saved in jpg format.
#'
#' @return
#'  The function creates a predictiveness curve.
#'
#'
#' @keywords hplot
#'
#'
#' @seealso \code{\link{predRisk}}
#' @references  Pepe MS, Feng Z, Huang Y, et al. Integrating the predictiveness
#' of a marker with its performance as a classifier.
#' Am J Epidemiol 2008;167:362-368.
#'
#' @examples
#' # specify dataset with outcome and predictor variables
#' data(ExampleData)
#'
#'  # fit logistic regression models
#'  # all steps needed to construct a logistic regression model are written in a function
#'  # called 'ExampleModels', which is described on page 4-5
#'  riskmodel1 <- ExampleModels()$riskModel1
#'  riskmodel2 <- ExampleModels()$riskModel2
#'
#'  # obtain predicted risks
#'  predRisk1 <- predRisk(riskmodel1)
#'  predRisk2 <- predRisk(riskmodel2)
#'
#'  # specify range of y-axis
#'  rangeyaxis <- c(0,1)
#'  # specify labels of the predictiveness curves
#'  labels <- c("without genetic factors", "with genetic factors")
#'
#' # produce predictiveness curves
#'  plotPredictivenessCurve(predrisk=cbind(predRisk1,predRisk2),
#'  rangeyaxis=rangeyaxis, labels=labels)
#'
"plotPredictivenessCurve" <-
function(predrisk, rangeyaxis, labels, plottitle, xlabel, ylabel, fileplot, plottype)
 {
   if (missing(plottitle)) {plottitle <- "Predictiveness curve"}
  if (missing(xlabel)) {xlabel<- "Cumulative percentage"}
  if (missing(ylabel)) {ylabel<- "Predicted risks"}
  if (missing(rangeyaxis)) {rangeyaxis<- c(0,1)}
 if (class(predrisk) == "numeric") {predrisk<- cbind(predrisk)}
# if (missing(labels)) {labels <- c(1:dim(predrisk)[2])}
a<-c(1:dim(predrisk)[2])
for(i in 1:dim(predrisk)[2])
		{
  x<-predrisk[,i]
   if (i==1)
   {
    xlim = c(0, 1)
    ylim = c(0, 1)
    xlab = xlabel
    ylab = ylabel
    x <- sort(x)
    n <- length(x)
    y <- (1:n)/n
    z <- y >= ylim[1] & y <= ylim[2]

resetGraph()
evalCall(plot, argu = list(x = y[z], y = x[z], type = "n",
xlab = "", ylab = "", las = 1, mgp = c(0, 0.6, 0),cex.axis =1.1, col=16+i,lty=i,
ylim=rangeyaxis,main =plottitle), checkdef=TRUE, checkpar = TRUE)

    evalCall(lines, argu = list(x = y[z], y = x[z],col=16+i,lty=i, lwd=2),
    checkdef = TRUE, checkpar = TRUE)
    mtext(xlab, side = 1, line = 2.75, cex = 1.2)
    mtext(ylab, side = 2, line = 2.5, cex = 1.2)
    invisible(data.frame(x = x, y = y))
   }
     else
   {
    xlim = c(0, 1)
    ylim = c(0, 1)
    x <- sort(x)
    n <- length(x)
    y <- (1:n)/n
    z <- y >= ylim[1] & y <= ylim[2]

    evalCall(lines, argu = list(x = y[z], y = x[z],col=16+i,lty=i,lwd=2,
    add= TRUE), checkdef = TRUE, checkpar = TRUE)
    mtext(xlab, side = 1, line = 2.75, cex = 1.2)
    mtext(ylab, side = 2, line = 2.5, cex = 1.2)
    invisible(data.frame(x = x, y = y))
   }
  	}
  		  if (!missing(labels)){
legend("bottomright",legend= labels,  col = c(17:(16+dim(predrisk)[2])),
lty=c(1:(dim(predrisk)[2])),lwd = 2,cex=1)
}
 if (missing(plottype)) {plottype<- "jpg"}
 	if (!missing(fileplot))
      savePlot(filename = fileplot,
	 type =plottype,
	 device = dev.cur())


 }
#' Function to plot histogram of risks separated for individuals with and without the outcome of interest.
#'
#'
#' @param data Data frame or numeric matrix that includes the outcome and
#' predictor variables.
#' @param cOutcome  Column number of the outcome variable.
#' @param risks   Risk of each individual. It
#' is specified by either a vector of risk scores or a vector of predicted risks.
#' @param interval  Size of the risk intervals. For example, \code{interval=.1}
#' will construct the following intervals for predicted risks:
#' \code{0-0.1, 0.1-0.2,..., 0.9-1}.
#' @param rangexaxis  Range of the x-axis. Specification of \code{rangexaxis} is optional.
#' @param rangeyaxis  Range of the y-axis.
#' @param plottitle  Title of the plot. Specification of \code{plottitle} is optional. Default is "Histogram of risks".
#' @param xlabel  Label of x-axis. Specification of \code{xlabel} is optional. Default is "Risk score".
#' @param ylabel  Label of y-axis. Specification of \code{ylabel} is optional. Default is "Percentage".
#' @param labels   Labels given to the groups of individuals without and
#'  with the outcome of interest. Specification of \code{labels} is optional. Default is
#'  \code{c("Without outcome", "With outcome")}.
#' @param fileplot Name of the output file that contains the plot. The file is
#' saved in the working directory in the format specified under \code{plottype}. Example:
#' \code{fileplot="plotname"}. Note that the extension is not specified here.
#' When \code{fileplot} is not specified, the plot is not saved.
#' @param plottype The format in which the plot is saved. Available formats are
#' wmf, emf, png, jpg, jpeg, bmp, tif, tiff, ps,
#' eps or pdf. For example, \code{plottype="eps"} will save the plot in eps format.
#' When \code{plottype} is not specified, the plot will be saved in jpg format.
#'
#' @return
#'  The function creates the histogram of risks separated for individuals
#'  with and without the outcome of interest.
#'
#'
#'
#' @keywords hplot
#'
#'
#' @seealso \code{\link{plotROC}}, \code{\link{riskScore}}
#' @examples
#'     # specify dataset with outcome and predictor variables
#'     data(ExampleData)
#'     # specify column number of the outcome variable
#'     cOutcome <- 2
#'
#'  # fit a logistic regression model
#'  # all steps needed to construct a logistic regression model are written in a function
#'  # called 'ExampleModels', which is described on page 4-5
#'  riskmodel <- ExampleModels()$riskModel2
#'
#'  # obtain predicted risks
#'  predRisk <- predRisk(riskmodel)
#'
#'     # specify the size of each interval
#'     interval <- .05
#'     # specify label of x-axis
#'     xlabel <- "Predicted risk"
#'     # specify label of y-axis
#'     ylabel <- "Percentage"
#'     # specify range of x-axis
#'     xrange <- c(0,1)
#'     # specify range of y-axis
#'     yrange <- c(0,40)
#'     # specify title for the plot
#'     maintitle <- "Distribution of predicted risks"
#'     # specify labels
#'     labels <- c("Without outcome", "With outcome")
#'
#'     # produce risk distribution plot
#'     plotRiskDistribution(data=ExampleData, cOutcome=cOutcome,
#'     risks=predRisk, interval=interval, plottitle=maintitle, rangexaxis=xrange,
#'     rangeyaxis=yrange, xlabel=xlabel, ylabel=ylabel, labels=labels)
#'
"plotRiskDistribution" <-
function(data, cOutcome, risks, interval, rangexaxis, rangeyaxis, plottitle,
xlabel, ylabel, labels, fileplot, plottype)
 {
   if (missing(plottitle)) {plottitle <- "Histogram of risks"}
  if (missing(xlabel)) {xlabel<- "Risk score"}
  if (missing(ylabel)) {ylabel<- "Percentage"}
    if (missing(labels)) {labels<- c("Without outcome", "With outcome")}
if (missing(data)){stop("The argument 'data' is missing")}
if (missing(cOutcome)){stop("The argument 'cOutcome' is missing")}
if (missing(risks)){stop("The argument 'risks' is missing")}
if (missing(rangeyaxis)){stop("The argument 'rangeyaxis' is missing")}
if (missing(rangexaxis)){rangexaxis<-c(min(risks),max(risks))}
if((rangexaxis[1]> min(risks))||(rangexaxis[2]< max(risks)))
{ stop("Please specify the range of x-axis that includes the values in 'risks'")}

   m <- table(data[,cOutcome], cut(risks,seq(rangexaxis[1], rangexaxis[2]+interval,interval)))
   p <- apply(m, 1, function(x) (x/sum(x))*100)

   a<-round(seq(rangexaxis[1]+interval,rangexaxis[2]+interval,interval),2)
   b<-round(seq(rangexaxis[1],rangexaxis[2],2*interval),2)

    for(i in 1:length(a))
		{
		if(any(a[i]==b)){a[i]<-a[i] }
		else{a[i] <- ""}
    }

   dimnames(p)[[1]] <- a
   dimnames(p)[[2]] <- labels
   barplot(height=t(p), beside=TRUE,col=c(24,1),ylab=ylabel,
    xlab=xlabel, legend.text=colnames(p), ylim = rangeyaxis,main=plottitle,
    cex.lab = 1.2,cex.axis=1.1, las=1,space = c(0,.10))

  if (missing(plottype)) {plottype<- "jpg"}
 	if (!missing(fileplot))
      savePlot(filename = fileplot,
	 type =plottype,
	 device = dev.cur())

    }
#' Function for reclassification table and statistics.
#' The function creates a reclassification table and provides statistics.
#'
#' The function creates a reclassification table and computes the net
#' reclassification improvement (\code{NRI}) and integrated discrimination
#' improvement (\code{IDI}). A reclassification table indicates the number
#' of individuals who move to another risk category or remain in the same
#' risk category as a result of updating the risk model. NRI equal to \code{x\%}
#' means that compared with individuals without outcome,
#'   individuals with outcome were almost \code{x\%} more likely to move up a category than down.
#' IDI equal to \code{x\%} means that the difference in average
#'   predicted risks between the individuals with and without the outcome
#'   increased by \code{x\%} in the updated model.
#' The function requires predicted risks estimated by using two separate risk
#' models. Predicted risks can be obtained using the functions
#' \code{\link{fitLogRegModel}} and \code{\link{predRisk}}
#' or be imported from other methods or packages.
#'
#' @param data Data frame or matrix that includes the outcome and
#' predictors variables.
#' @param cOutcome  Column number of the outcome variable.
#' @param predrisk1  Vector of predicted risks of all individuals using initial
#' model.
#' @param predrisk2  Vector of predicted risks of all individuals using updated
#' model.
#' @param cutoff  Cutoff values for risk categories.
#' Define the cut-off values as \code{c(0,...,1)}.
#' Multiple values can be defined and always specify \code{0} and \code{1}.
#' Example: \code{c(0,.20,.30,1)}
#'
#'  @return
#'   The function returns the reclassification table, separately
#'  for individuals with and without the outcome of interest and the following measures:
#'  \item{NRI}{Net Reclassification Improvement with 95\% CI and \code{p-value} of the test}
#'  \item{IDI}{Integrated Discrimination Improvement with 95\% CI and \code{p-value}
#'  of the test}
#'
#'
#' @keywords htest
#'
#'
#' @references Cook NR. Use and misuse of the receiver operating characteristic
#' curve in risk prediction. Circulation 2007;115(7):928-935.
#'
#' Pencina MJ, D'Agostino RB Sr, D'Agostino RB Jr, Vasan RS.
#' Evaluating the added predictive ability of a new marker: from
#' area under the ROC curve to reclassification and beyond. Stat
#' Med 2008;27(2):157-172; discussion 207-212.
#'
#'
#' @seealso \code{\link{plotDiscriminationBox}}, \code{\link{predRisk}}
#' @examples
#' # specify dataset with outcome and predictor variables
#' data(ExampleData)
#' # specify column number of the outcome variable
#'  cOutcome <- 2
#'
#'  # fit logistic regression models
#'  # all steps needed to construct a logistic regression model are written in a function
#'  # called 'ExampleModels', which is described on page 4-5
#'  riskmodel1 <- ExampleModels()$riskModel1
#'  riskmodel2 <- ExampleModels()$riskModel2
#'
#'  # obtain predicted risks
#'  predRisk1 <- predRisk(riskmodel1)
#'  predRisk2 <- predRisk(riskmodel2)
#' # specify cutoff values for risk categories
#'  cutoff <- c(0,.10,.30,1)
#'
#' # compute reclassification measures
#'  reclassification(data=ExampleData, cOutcome=cOutcome,
#' predrisk1=predRisk1, predrisk2=predRisk2, cutoff)
#'
"reclassification" <-
function(data,cOutcome,predrisk1,predrisk2, cutoff) {

c1 <- cut(predrisk1,breaks = cutoff ,include.lowest=TRUE,right= FALSE)
c2 <- cut(predrisk2,breaks = cutoff ,include.lowest=TRUE,right= FALSE)
tabReclas <- table("Initial Model"=c1, "Updated Model"=c2)
cat(" _________________________________________\n")
cat(" \n     Reclassification table    \n")
cat(" _________________________________________\n")

 ta<- table(c1, c2, data[,cOutcome])

  cat ("\n Outcome: absent \n  \n" )
  TabAbs <- ta[,,1]
  tab1 <- cbind(TabAbs, " % reclassified"= round((rowSums(TabAbs)-diag(TabAbs))/rowSums(TabAbs),2)*100)
  names(dimnames(tab1)) <- c("Initial Model", "Updated Model")
  print(tab1)

  cat ("\n \n Outcome: present \n  \n" )
  TabPre <- ta[,,2]
  tab2 <- cbind(TabPre, " % reclassified"= round((rowSums(TabPre)-diag(TabPre))/rowSums(TabPre),2)*100)
  names(dimnames(tab2)) <- c("Initial Model", "Updated Model")
  print(tab2)

  cat ("\n \n Combined Data \n  \n" )
  Tab <- tabReclas
  tab <- cbind(Tab, " % reclassified"= round((rowSums(Tab)-diag(Tab))/rowSums(Tab),2)*100)
  names(dimnames(tab)) <- c("Initial Model", "Updated Model")
  print(tab)
cat(" _________________________________________\n")

c11 <-factor(c1, levels = levels(c1), labels = c(1:length(levels(c1))))
c22 <-factor(c2, levels = levels(c2), labels = c(1:length(levels(c2))))

  x<-improveProb(x1=as.numeric(c11)*(1/(length(levels(c11)))),
  x2=as.numeric(c22)*(1/(length(levels(c22)))), y=data[,cOutcome])
  

  y<-improveProb(x1=predrisk1, x2=predrisk2, y=data[,cOutcome])

cat("\n NRI(Categorical) [95% CI]:", round(x$nri,4),"[",round(x$nri-1.96*x$se.nri,4),"-",
 round(x$nri+1.96*x$se.nri,4), "]", "; p-value:", round(2*pnorm(-abs(x$z.nri)),5), "\n" )

 cat(" NRI(Continuous) [95% CI]:", round(y$nri,4),"[",round(y$nri-1.96*y$se.nri,4),"-",
 round(y$nri+1.96*y$se.nri,4), "]", "; p-value:", round(2*pnorm(-abs(y$z.nri)),5), "\n" )

cat(" IDI [95% CI]:", round(y$idi,4),"[",round(y$idi-1.96*y$se.idi,4),"-",
 round(y$idi+1.96*y$se.idi,4), "]","; p-value:", round(2*pnorm(-abs(y$z.idi)),5), "\n")
}
#' Function for box plots of predicted risks separately for individuals with and
#' without the outcome of interest.
#'  The function produces box plots of predicted risks for individuals with
#'  and without the outcome of interest and calculates the discrimination slope.
#'
#'  The discrimination slope is the difference between the mean predicted risks
#'  of individuals with and without the outcome of interest. Predicted risks
#'  can be obtained using the
#' \code{\link{fitLogRegModel}} and \code{\link{predRisk}} or be
#' imported from other programs. The difference between discrimination
#' slopes of two separate risk models is equivalent
#' to (\code{IDI}) which is discussed
#' in details in the \code{\link{reclassification}} function.
#'
#' @param data Data frame or matrix that includes the outcome and
#' predictors variables.
#' @param cOutcome  Column number of the outcome variable.
#' @param predrisk  Vector of predicted risks.
#' @param labels   Labels given to the groups of individuals without and with
#' the outcome of interest. Specification of \code{label} is optional.
#' Default is \code{c("Without disease", "With disease")}.
#' @param plottitle  Title of the plot. Specification of \code{plottitle}
#' is optional. Default is "Box plot".
#' @param ylabel  Label of y-axis. Specification of \code{ylabel}
#' is optional. Default is "Predicted risks".
#' @param fileplot Name of the file that contains the plot. The file is
#' saved in the working directory in the format specified under \code{plottype}.
#' Example: \code{fileplot="name"}. Note that the extension is not specified here.
#' When \code{fileplot} is not specified, the plot is not saved.
#' @param plottype The format in which the plot is saved. Available formats are
#' wmf, emf, png, jpg, jpeg, bmp, tif, tiff, ps,
#' eps or pdf. Foe example, \code{plottype="eps"} will save the plot in eps format.
#' When \code{plottype} is not specified, the plot will be saved in jpg format.
#'
#'  @return
#'  The function creates a box plots of predicted risks for
#'  individuals with and without the outcome of interest and returns the discrimination slope.
#'
#'
#' @keywords hplot
#'
#'
#' @references Yates JF. External correspondence: decomposition of the mean
#' probability score. Organizational Behavior and Human Performance 1982;30:132-156.
#'
#' @seealso \code{\link{reclassification}}, \code{\link{predRisk}}
#' @examples
#' # specify dataset with outcome and predictor variables
#' data(ExampleData)
#' # specify column number of outcome variable
#'  cOutcome <- 2
#'
#'  # fit a logistic regression model
#'  # all steps needed to construct a logistic regression model are written in a function
#'  # called 'ExampleModels', which is described on page 4-5
#'  riskmodel <- ExampleModels()$riskModel2
#'
#'  # obtain predicted risks
#'  predRisk <- predRisk(riskmodel)
#'  # specify labels for the groups without and with the outcome of interest
#'  labels <- c("Without disease", "With disease")
#'
#' # produce discrimination box plot
#'  plotDiscriminationBox(data=ExampleData, cOutcome=cOutcome, predrisk=predRisk,
#'  labels=labels)
#'
plotDiscriminationBox <- function(data, cOutcome, predrisk, labels, plottitle, ylabel, fileplot, plottype)
 {
  if (missing(labels)) {label <- c("Without disease", "With disease")}
  if (missing(plottitle)) {plottitle <- "Box plot"}
  if (missing(ylabel)) {ylabel<- "Predicted risks"}
risk <- predrisk
a<-0;b<-1
if((max(predrisk)>1)|(min(predrisk)<0)){a<-min(predrisk);b<-max(predrisk)}
boxplot(risk~data[,cOutcome],ylab=ylabel,ylim=c(a,b),cex.lab=1.2,las=1 ,
 main= plottitle, cex.axis=1.1, names=labels)
boxplot(c(mean(risk[data[,cOutcome]==0]),mean(risk[data[,cOutcome]==1]))~c(0,1),
add=TRUE, boxlty=0, staplelty=0, medlty=0, medlwd=0, medpch=15,las=1,
cex.axis=1.1, xaxt='n')
p<- list(Discrim_Slope = round(mean(risk[data[,cOutcome]==1]) -
 mean(risk[data[,cOutcome]==0]),3))

 if (missing(plottype)) {plottype<- "jpg"}
 	if (!missing(fileplot))
      savePlot(filename = fileplot,type =plottype,device = dev.cur())
return(p)
}
#' An example code to construct a risk model using logistic regression analysis.
#' \code{ExampleModels} constructs two risk models using logistic regression analysis.
#' Most of the functions in this package require a logistic regression model as an input and
#' estimate predicted risks from this fitted model.
#' To illustrate these functions without repeating the construction of a
#' logistic regression model, this example code has been created.
#' The function returns two different risk models, riskModel1 which is based
#' on non-genetic predictors and riskModel2 which includes genetic and non-genetic predictors.
#'
#' @examples
#' # specify dataset with outcome and predictor variables
#' data(ExampleData)
#' # specify column number of the outcome variable
#'  cOutcome <- 2
#' # specify column numbers of non-genetic predictors
#'  cNonGenPred1 <- c(3:10)
#' cNonGenPred2 <- c(3:10)
#'  # specify column numbers of non-genetic predictors that are categorical
#'  cNonGenPredCat1 <- c(6:8)
#'  cNonGenPredCat2 <- c(6:8)
#'  # specify column numbers of genetic predictors
#'  cGenPred1 <- c(0)
#'  cGenPred2 <- c(11:16)
#'  # specify column numbers of genetic predictors that are categorical
#'  cGenPredsCat1 <- c(0)
#'  cGenPredsCat2 <- c(0)
#'
#'  # fit logistic regression models
#'  riskmodel1 <- fitLogRegModel(data=ExampleData, cOutcome=cOutcome,
#'  cNonGenPreds=cNonGenPred1, cNonGenPredsCat=cNonGenPredCat1,
#'  cGenPreds=cGenPred1, cGenPredsCat=cGenPredsCat1)
#'  riskmodel2 <- fitLogRegModel(data=ExampleData, cOutcome=cOutcome,
#'  cNonGenPreds=cNonGenPred2, cNonGenPredsCat=cNonGenPredCat2,
#'  cGenPreds=cGenPred2, cGenPredsCat=cGenPredsCat2)
#'
#'  # combine output in a list
#'  ExampleModels <- list(riskModel1=riskmodel1, riskModel2=riskmodel2)
#'
"ExampleModels" <- function()
 {
  data(ExampleData, envir = environment())
  cOutcome <- 2
  cNonGenPred1 <- c(3:10)
  cNonGenPredCat1 <- c(6:8)
  cGenPred1 <- c(0)
  cGenPredsCat1 <- c(0)
  cNonGenPred2 <- c(3:10)
  cNonGenPredCat2 <- c(6:8)
  cGenPred2 <- c(11:16)
  cGenPredsCat2 <- c(0)
  riskmodel1 <- fitLogRegModel(data=ExampleData, cOutcome=cOutcome,
  cNonGenPreds=cNonGenPred1, cNonGenPredsCat=cNonGenPredCat1,
  cGenPreds=cGenPred1, cGenPredsCat=cGenPredsCat1)
  riskmodel2 <- fitLogRegModel(data=ExampleData, cOutcome=cOutcome,
  cNonGenPreds=cNonGenPred2, cNonGenPredsCat=cNonGenPredCat2,
  cGenPreds=cGenPred2, cGenPredsCat=cGenPredsCat2)
  out<- list(riskModel1=riskmodel1, riskModel2=riskmodel2)
  return(out)
  }
#' Function to construct a simulated dataset containing individual genotype data, 
#' genetic risks and disease status for a hypothetical population.
#' Construct a dataset that contains individual genotype data, genetic risk,
#' and disease status for a hypothetical population.
#' The dataset is constructed using simulation in such a way that the frequencies
#' and odds ratios (ORs) of the genetic variants and the population disease risk
#' computed from this dataset are the same as specified by the input parameters.
#'
#' The function will execute when the matrix with odds ratios and frequencies,
#' population disease risk and the number of individuals are specified. \cr
#'
#' The simulation method is described in detail in the references. \cr
#'
#'
#' The method assumes that (i) the combined effect of the genetic variants
#' on disease risk follows a multiplicative (log additive) risk model;
#' (ii) genetic variants inherit independently, that is no linkage disequilibrium
#' between the variants; (iii) genetic variants have independent effects on the
#' disease risk, which indicates no interaction among variants; and (iv) all
#' genotypes and allele proportions are in Hardy-Weinberg equilibrium.
#' Assumption (ii) and (iv) are used to generate the genotype data, and assumption
#'(i), (ii) and (iii) are used to calculate disease risk.
#'
#'
#' Simulating the dataset involves three steps: (1) modelling genotype data,
#' (2) modelling disease risks, and (3) modelling disease status. Brief
#' descriptions of these steps are as follows:
#'
#'
#' (1) Modelling genotype data: For each variant the genotype
#' frequencies are either specified or calculated from the allele frequencies
#' using Hardy-Weinberg equilibrium. Then, the genotypes for each genetic
#' variant are randomly distributed without replacement over all individuals.
#'
#'
#' (2) Modelling disease risks: For the calculation of the individual disease
#' risk, Bayes' theorem is used, which states that the posterior odds of disease
#' are obtained by multiplying the prior odds by the likelihood ratio (LR) of
#' the individual genotype data. The prior odds are calculated from the
#' population disease risk or disease prevalence
#' (prior odds= prior risk/ (1- prior risk)) and the posterior odds are converted
#' back into disease risk (disease risk= posterior odds/ (1+ posterior odds)).
#' Under the no linkage disequilibrium (LD) assumption, the LR of a genetic profile
#' is obtained by multiplying the LRs of the single genotypes that are included in
#' the risk model. The LR of a single genotype is calculated using frequencies
#' and ORs of genetic variants and population disease risk. See references
#' for more details.
#'
#'
#' (3) Modelling disease status: To model disease status, we used a procedure
#' that compares the estimated disease risk of each subject to a randomly drawn value
#' between 0 and 1 from a uniform distribution. A subject is assigned to the
#' group who will develop the disease when the disease risk is higher than the
#' random value and to the group who will not develop the disease when the risk
#' is lower than the random value.
#'
#' This procedure ensures that for each genomic profile, the percentage of
#' people who will develop the disease equals the population disease risk
#' associated with that profile, when the subgroup of individuals with that
#' profile is sufficiently large.
#'
#'
#' @param ORfreq Matrix with ORs and frequencies of the genetic variants.
#' The matrix contains four columns in which the first two describe ORs and the
#' last two describe the corresponding frequencies. The number of rows in this
#' matrix is same as the number of genetic variants included. Genetic variants
#' can be specified as per genotype, per allele, or as dominant/ recessive
#' effect of the risk allele. When per genotype data are used, OR of the
#' heterozygous and homozygous risk genotypes are mentioned in the first two
#' columns and the corresponding genotype frequencies are mentioned in the last
#' two columns. When per allele data are used, the OR and frequency of the risk
#' allele are specified in the first and third column and the remaining two cells
#' are coded as '1'.  Similarly, when dominant/ recessive effects of the risk
#' alleles are used, the
#' OR and frequency of the dominant/ recessive variant are specified in the first
#' and third column, and the remaining two cells are coded as '0'.
#' @param poprisk Population disease risk (expressed in proportion).
#' @param popsize  Total number of individuals included in the dataset.
#' @param filename Name of the file in which the dataset will be saved.
#' The file is saved in the working directory as a txt file. When no filename
#' is specified, the output is not saved.
#'
#'  @return
#'   The function returns:
#'   \item{Dataset}{A data frame or matrix that includes genotype data,
#' genetic risk and disease status for a hypothetical population.
#' The dataset contains (4 + number of genetic variants included) columns,
#' in which the first column is the un-weighted risk score, which is the sum
#' of the number of risk alleles for each individual, the third column is the
#' estimated genetic risk, the forth column is the individual disease status expressed
#' as '0' or '1', indicating without or with the outcome of interest, and the fifth until
#' the end column are genotype data for the variants expressed as '0', '1' or '2',
#' which indicate the number of risk alleles present in each individual for the genetic variants.}
#'
#'
#' @keywords models
#'
#'
#' @references Janssens AC, Aulchenko YS, Elefante S, Borsboom GJ, Steyerberg EW,
#' van Duijn CM. Predictive testing for complex diseases using multiple genes:
#' fact or fiction? Genet Med. 2006;8:395-400.
#'
#' Kundu S, Karssen LC, Janssens AC: Analytical and simulation methods for 
#' estimating the potential predictive ability of genetic profiling: a comparison 
#' of methods and results. Eur J Hum Genet. 2012 May 30.
#' 
#' van Zitteren M, van der Net JB, Kundu S, Freedman AN, van Duijn CM,
#' Janssens AC. Genome-based prediction of breast cancer risk in the general
#' population: a modeling study based on meta-analyses of genetic associations.
#' Cancer Epidemiol Biomarkers Prev. 2011;20:9-22.
#'
#' van der Net JB, Janssens AC, Sijbrands EJ, Steyerberg EW. Value of genetic
#' profiling for the prediction of coronary heart disease.
#' Am Heart J. 2009;158:105-10.
#'
#' Janssens AC, Moonesinghe R, Yang Q, Steyerberg EW, van Duijn CM, Khoury MJ.
#' The impact of genotype frequencies on the clinical validity of genomic
#' profiling for predicting common chronic diseases. Genet Med. 2007;9:528-35.
#'


#'
#'

#'
#'
#' @examples
#' # specify the matrix containing the ORs and frequencies of genetic variants
#' # In this example we used per allele effects of the risk variants
#' ORfreq<-cbind(c(1.35,1.20,1.24,1.16), rep(1,4), c(.41,.29,.28,.51),rep(1,4))
#'
#' # specify the population disease risk
#'  popRisk <- 0.3
#' # specify size of hypothetical population
#'  popSize <- 10000
#'
#' # Obtain the simulated dataset
#' Data <- simulatedDataset(ORfreq=ORfreq, poprisk=popRisk, popsize=popSize)
#'
#' # Obtain the AUC and produce ROC curve
#' plotROC(data=Data, cOutcome=4, predrisk=Data[,3])
#'
"simulatedDataset" <- function(ORfreq, poprisk, popsize, filename) 
{
if (missing(poprisk)) {stop("Population disease risk is not specified")}
if (missing(popsize)) {stop("Total number of individuals is not mentioned")}

g <- nrow(ORfreq)
reconstruct.2x2table <- function(p,d,OR,s)
{
a <- 0
b <- 0
c <- (OR*p*s*(1-d)*d*s)/((1-p)*s*(1-d)+OR*p*s*(1-d))
dd <- p*s-c
e <- d*s-c
f <- (1-p)*s-e
tabel <- cbind(a,b,c,dd,e,f,g,OR)
tabel
}

reconstruct.2x3table <- function(OR1,OR2,p1,p2,d,s){
	a	<- 1
	eOR	<- 0
	while (eOR<=OR2){
		b	<- p2*s*(1-d)
		snew <- s-a-b
		p1new <-p1/(1-p2)
		dnew <- (d-(a/s))/((d-(a/s))+ ((1-d)-b/s))
		c	<-	(OR1*p1new*snew*(1-dnew)*dnew*snew)/((1-p1new)*snew*(1-dnew)+OR1*p1new*snew*(1-dnew))
		dd	<-	p1new*((1-d)-b/s)*s
		e	<-	(d-(a/s))*s-c
		f	<- ((1-d)-b/s)*s-dd
		eOR	<- (a*f)/(b*e)
		tabel <- cbind(a,b,c,dd,e,f,g,OR1,OR2)
		a	<- a+1
		tabel
		}
	tabel
}

reconstruct.2x3tableHWE <- function(OR,p,d,s){
  OR1 <- OR
  OR2 <- OR^2
	p1 <-  2*p*(1-p)
	p2 <-  p*p

	a	<- 1
	eOR	<- 0
	while (eOR<=OR2){
		b	<- p2*s*(1-d)
		snew <- s-a-b
		p1new <-p1/(1-p2)
		dnew <- (d-(a/s))/((d-(a/s))+ ((1-d)-b/s))
		c	<-	(OR1*p1new*snew*(1-dnew)*dnew*snew)/((1-p1new)*snew*(1-dnew)+OR1*p1new*snew*(1-dnew))
		dd	<-	p1new*((1-d)-b/s)*s
		e	<-	(d-(a/s))*s-c
		f	<- ((1-d)-b/s)*s-dd
		eOR	<- (a*f)/(b*e)
		tabel <- cbind(a,b,c,dd,e,f,g,OR1,OR2)
		a	<- a+1
		tabel
		}
	tabel
}


adjust.postp <- function (pd, LR){		
	odds.diff <- 0
	prior.odds <- pd/(1-pd)	
	for (i in (1:100000)) {
	Postp <- (prior.odds*LR)/(1+(prior.odds*LR))
	odds.diff <- (pd-mean(Postp))/ (1-(pd-mean(Postp)))
	prior.odds	<- prior.odds+odds.diff
	if (odds.diff < .0001) break
	}
	Postp
}

func.data <- function(p,d,OR,s,g){
  Data <- matrix (NA,s,4+g)
  Data[,1] <- rep(0,s)                   
	Data[,2] <- rep(1,s)									
	Data[,3] <- rep(0,s)
	i <- 0
	while (i < g){
    i <- i+1
    cells2x3 <- rep(NA,9)
    cells2x3 <- if(p[i,2]==0) {reconstruct.2x2table(p=p[i,1],d,OR=OR[i,1],s)} else {if(p[i,2]==1) {reconstruct.2x3tableHWE(OR=OR[i,1],p=p[i,1],d,s)}
  else {reconstruct.2x3table(OR1=OR[i,1],OR2=OR[i,2],p1=p[i,1],p2=p[i,2],d,s)}}   
      LREE 	  <- ((cells2x3[1]/d*s)/(cells2x3[2]/(1-d)*s))			
      LREe	  <- ((cells2x3[3]/d*s)/(cells2x3[4]/(1-d)*s))
      LRee	  <- ((cells2x3[5]/d*s)/(cells2x3[6]/(1-d)*s))
 
 Gene <- if(p[i,2]==0){c(rep(0,((1-p[i,1]-p[i,2])*s)),rep(1,p[i,1]*s),rep(2,p[i,2]*s))}
 else {if(p[i,2]==1) {c(rep(0,(((1-p[i,1])^2)*s)),rep(1,2*p[i,1]*(1-p[i,1])*s),rep(2,p[i,1]*p[i,1]*s))}
  else {c(rep(0,((1-p[i,1]-p[i,2])*s)),rep(1,p[i,1]*s),rep(2,p[i,2]*s))}}  
		Filler <- s-length(Gene)                               
		Gene <- sample(c(Gene,rep(0,Filler)),s,replace=FALSE)
    Data[,4+i] <- Gene
    GeneLR <- ifelse(Gene==0,LRee,ifelse(Gene==1,LREe,LREE))
   
    Data[,1] <- Data[,1]+Gene
    Data[,2] <- Data[,2]*GeneLR	
			
#	 cat(i,"")
		}
		
		Data[,3] <- adjust.postp(pd=d, LR=Data[,2])				
		Data[,4]  <- ifelse(runif(s)<=(Data[,3]), 1, 0)  					          	
    Data <- as.data.frame(Data)
    Data
    }
  
 simulatedData <- func.data  (p=ORfreq[,c(3,4)],d=poprisk,OR=ORfreq[,c(1,2)],s=popsize,g=nrow(ORfreq))   


if (!missing(filename)) 
	{write.table( simulatedData,file=filename, row.names=TRUE,sep = "\t")  }

 return(simulatedData)
}
