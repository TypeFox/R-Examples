#' Bootstrapping t-statistics
#'
#' This function is used to create non-parametrically bootstrapped samples of marginal effects
#' calculated from a model with interaction. The function takes data, model, an index of observations, 
#' and identities of the interacted variables and returns a vector of t-statistics
#' corresponding to a single bootstrap deviate for marginal effects in the direction of 
#' x.name and z.name.
#'
#' @param dat The data set on which the original model was run
#' @param K Index of observations in dat that will be selected for a particular bootstrap run
#' @param form The formula for the GLM model
#' @param fam The family for the GLM model
#' @param x.name The identity of the first interacted variable
#' @param z.name the identity of the second interacted variable
#' @author Justin Esarey and Jane Lawrence Sumner
#'
#' @return The bootstrap-t statistics from the function for dy/dx and dy/dz, respectively.
#' @note Form should include x.name*z.name, in that order.
#' @examples 
#' \dontrun{
#'  data(legfig)                # Clark and Golder 2006 replication data
#' 
#'  set.seed(1231124)
#'  
#   # limit to established democracies from the 1990s
#'  dat<-subset(legfig, subset=(nineties==1 & old==1))
#'
#'  # create bootstrap samples of marginal effects of eneg and logmag on enep1
#'  library(boot)
#'  boot.t.dist <- boot(data = dat, statistic = bootFun, R = 1000, 
#'           form=enep1 ~ eneg * logmag + uppertier_eneg + uppertier + proximity1 + 
#'           proximity1_enpres + enpres, fam="gaussian", x.name="eneg", 
#'           z.name="logmag")$t
#'  boot.t.x.dist<-boot.t.dist[,1:10]
#'  
#'  # calculate critical t-statistic that sets familywise error rate to 10%
#'  # for statistical significance of marginal effect of of eneg at any value of logmag
#'  findMultiLims(boot.t.x.dist, type="any", err=0.1)$minimum         # answer: 2.593086
#'  }
#' @references Clark, William R., and Matt Golder. 2006. "Rehabilitating Duverger's Theory." \emph{Comparative Political Studies} 39(6): 679-708.
#' @references Esarey, Justin, and Jane Lawrence Sumner. 2015. "Marginal Effects in Interaction Models: Determining and Controlling the False Positive Rate." URL: http://jee3.web.rice.edu/interaction-overconfidence.pdf.
#' @importFrom stats na.omit glm vcov
#' @export
 
bootFun<-function(dat, K, form, fam="gaussian", x.name, z.name){
    
  # remove all extraneous variables from the data, delete missing observations
  dat<-subset(dat, select=all.vars(form))
  dat<-na.omit(dat)
  
  # anticipate what the name of the interaction term will be
  int.name<-paste(x.name, ":", z.name, sep="")
  
  # fit the model to the full data set
  int.mod<-glm(formula=form, family=fam, data=dat)
  
  # determine what MEs to calculate in the x dimension
  # calculate at every unique value of z, OR ten evenly spaced values of z, depending on the nature of z
  if(length(unique(dat[[z.name]]))<=10)
  {z.test<-sort(unique(dat[[z.name]]))}else{z.test<-seq(from=min(dat[[z.name]]), to=max(dat[[z.name]]), length.out=10)}
  # calculate dy/dx | z for the values of z.test
  coef.x.o<-summary(int.mod)$coefficients[x.name,1]+summary(int.mod)$coefficients[int.name,1]*z.test
  
  # determine what MEs to calculate in the z dimension
  # calculate at every unique value of x, OR ten evenly spaced values of x, depending on the nature of x
  if(length(unique(dat[[x.name]]))<=10)
  {x.test<-sort(unique(dat[[x.name]]))}else{x.test<-seq(from=min(dat[[x.name]]), to=max(dat[[x.name]]), length.out=10)}
  # calculate dy/dz | x for the values of z.test
  coef.z.o<-summary(int.mod)$coefficients[z.name,1]+summary(int.mod)$coefficients[int.name,1]*x.test
  
  # create the bootstrapped data set using the indices K
  dat.boot<-dat[K,]
  
  # run the model on the bootstrapped data
  int.mod.store<-glm(formula=form, family=fam, data=dat.boot)
  
  # calculate dy/dx | z for the values of z.test
  coef.x<-summary(int.mod.store)$coefficients[x.name,1]+summary(int.mod.store)$coefficients[int.name,1]*z.test
  se.x<-sqrt(vcov(int.mod.store)[x.name,x.name]+(z.test^2)*vcov(int.mod.store)[int.name,int.name]+2*z.test*vcov(int.mod.store)[x.name,int.name])
  
  # calculate the bootstrap-t values for dy/dx | z
  t.stat.x<-(coef.x-coef.x.o)/se.x
  
  # calculate dy/dz | x for the values of z.test
  coef.z<-summary(int.mod.store)$coefficients[z.name,1]+summary(int.mod.store)$coefficients[int.name,1]*x.test
  se.z<-sqrt(vcov(int.mod.store)[z.name,z.name]+(x.test^2)*vcov(int.mod.store)[int.name,int.name]+2*x.test*vcov(int.mod.store)[z.name,int.name])
  
  # calculate the bootstrap-t values for dy/dz | x
  t.stat.z<-(coef.z-coef.z.o)/se.z
  
  # return the bootstrap-t statistics from the function
  return(c(t.stat.x, t.stat.z))
  
}

