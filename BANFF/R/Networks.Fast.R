#'Fast Algorithm for Bayesian Network Discovery
#'@param pvalue a vector of p-values obtained from large scale statistical hypothesis testing
#'@param net a n by n network configuration, n is the length of pvalue
#'@param iter number of iterations. The default is 5000
#'@param nburn number of burn-in. The default is 2000
#'@param initials the way you choose to get the inital posteritor samples of parameters produced by DPM fitting. Generating initals by function Mclust() when you set "mclust" or by DPdensity() when you set "DPdensity". It is recommended to choose "Mclust" when the dimension of data is large, and to choose "DPdensity" when the dimenstion is small.
#'@param v number of iterations set for DPM fitting by "DPdensity". v is only valid when you choose initals as "DPdensity"
#'@param piall a vector of selections of pi0. The default vector is 0.75, 0.8, 0.85, 0.9. The selections of pi0 should be placed in sequence, from smaller to larger.
#'@param rhoall a vector of selections of rho0 and rho1. The default vector is 1, 2, 5, 10, 15. The selections of rho0 and rho1 should be placed in sequence, from smaller to larger.
#'@details This generic function fits a Bayesian Nonparametric Mixture Model for gene selection incorporating network information (Zhao et al., 2014):
#' \itemize{
#' \item  r_i| g_i, \strong{theta} ~ N(mu_{g_i}, sigma_{g_i}),
#' \item  g_i | z_i=k, \strong{q}_k ~ Discrete(\strong{a}_k, \strong{q}_k),
#' \item	\strong{theta} ~ G_{0k}, for g in \strong{a}_k,
#' \item	\strong{q}_k ~ Dirichlet(tau_k \strong{1}_{L_k}/L_k),
#' \item	\strong{theta}={\strong{theta}_g}_{g in \strong{a}_0 and \strong{a}_1} 
#' \item	\strong{theta}_g=(mu_g, sigma_g)
#' }
#' where we define 
#' \describe{
#' \item{Index}{\strong{a}_0=(-L_0+1,-L_0+2,...,0) , \strong{a}_1=(1,2,...,L_1) and the correspondent probability q_0=(q_{-L_0+1}, q_{-L_0+2}, ...,q_0), q_1=(q_1, q_2, ..., q_{L_1}), according to the defination of Discrete(\strong{a}_k, \strong{b}_k), for example, Pr(g_i={L_0+2})=q_{-L_0+2}. }
#' \item{Assumption}{We have an assumption that "selected" gene or image pixel should have larger statiscs comparing to "unselected" ones without the loss of generality. In this regard, we set the restriction mu_g<mu_{g+1} for g=-L_0+1, -L_0+2,...,L_1.}
#' }
#' For this function, The NET-DPM-3, considered as Fast function is applied , and more details about the algorithm can be referred from Appnendix B.3 of Zhao et al., 2014
#'@return a list cantaining Bayesian Statistics information of the distribution of zi
#'\describe{
#'\item{trace}{a n by (iter-nburns) matrix (n is the length of elements and iter-nburns is the length of the saved chain), showing the evolution of the chain for each element}
#'\item{mean}{the mean of the distribution for each element}
#'\item{median}{the median of the distribution for each element}
#'\item{var}{the variance of the distribution for each element}
#'\item{quantile}{the quantiles of the distribution for each element}
#'}
#'@examples ####Example1. For a 10X10 image with 5X5 signal for example
#' ##Creating the network of 10X10 image
#' library(igraph)
#' library(BayesNetDiscovery)
#' g <- graph.lattice(length=10,dim=2)
#' net=as(get.adjacency(g,attr=NULL),"matrix")##this is the input of argument \code{net}
#' ##Assign the signal elements with signal intenstion as normal distribution N(1,0.2). While noise is set as N(0,0.2) 
#' newz=rep(0,100)
#' for (i in 3:7)
#' {
#'  newz[(i*10+3):(i*10+7)]=1
#' }
#' testcov<-0
#' for(i in 1:100){
#'  if(newz[i]==0){
#'    testcov[i]<-rnorm(1,mean=0,sd=0.2)
#'  }else{
#'   testcov[i]<-rnorm(1,mean=1,sd=0.2)
#'  }
#' }
#' ##The profile of the impage
#' image(matrix(testcov,10,10),col=gray(seq(0,1,length=255)))
#' ##Transform the signals into pvalue form and begin identification
#' pvalue=pnorm(-testcov)
#' total2=Networks.Fast(pvalue,net,iter=5000,initials="mclust",piall=c(0.8, 0.85, 0.9, 0.95),rhoall=c(0.5,1,5,10,15))
#' 
#' 
#' 
#' 
#' ####Example.2 Gene Network discovery
#' ##Generating Scale free Gene Network
#' library(igraph)
#' library(BayesNetDiscovery)
#' g <- barabasi.game(50, power=1, zero.appeal=1.5,directed = F)
#' tkplot(g, layout=layout.kamada.kawai)
#' net=as(get.adjacency(g,attr=NULL),"matrix")
#' ##Random assign selected genes and make the signal intension as gaussian mixture
#' newz=rep(c(1,0,0,1,0),10)
#' Simnorm=function(n){
#' weight = c(0.4, 0.6)
#'mu = c(5,4)
#'sigma = c(1,0.5)
#'z = sample(c(1,2),size=n, prob=weight,replace=TRUE)
#'r = rnorm(n,mean=mu[z],sd=sigma[z])
#'return(r)
#'}
#'testcov<-0
#'for(i in 1:50){
#'  if(newz[i]==0){
#'    testcov[i]<-rnorm(1,mean=0,sd=1)
#'  }else{
#'   testcov[i]<-Simnorm(1)
#'  }
#'}
#'pvalue=pnorm(-testcov)
#'total1=Networks.Fast(pvalue,net,iter=5000,v=20,initials="DPdensity",piall=c(0.8, 0.85, 0.9, 0.95),rhoall=c(1, 2, 5, 10, 15))
#'@export
Networks.Fast=function(pvalue,net,iter=5000,nburns=2000,algorithms=c("EM","DPM"),v=20,DPM.mcmc=list(nburn=2000,nsave=1,nskip=0,ndisplay=10),DPM.prior=list(a0=2,b0=1,m2=rep(0,1),s2=diag(100000,1),
                                                                                                                                                           psiinv2=solve(diag(0.5,1)),
                                                                                                                                                           nu1=4,nu2=4,tau1=1,tau2=100),DPparallel=FALSE,n.cores=1,piall=c(0.8, 0.85, 0.9, 0.95),rhoall=c(1, 2, 5, 10, 15),show.steps=10,showlikelihood=FALSE, likelihood.frequency=100)
{
  ####Preparing Stage:Defining Parameters, 
  print("NOW_Transferring p-values to testing statistics")
  rstat=Transfer(pvalue)###Step1 : Transferring 
  print("NOW_Getting initials by Kmeans")
  wholeindex=Kmeans(rstat)###Step2: Using Kmeans for giving the initial values
  
  #initial_mu_var=Initial_mu_var(rstat,wholeindex) ###Step3: Genering initial mu and var based on intialindex
  
  pirhopair=Pi_rho_gen(piall,rhoall)####Step4.1 Generating all the possible selections of pi and rho
  print("NOW_Generating_Zi for likelihood comparison")
  znew=Generating_Zi(net,pirhopair,wholeindex)####Step4.2 Generating Zi by different setting of pi and rho.
  print("NOW_Comparing the likelihood to select the best set")
  choice=Pi_rho_selecting(znew,rstat,piall,rhoall)####Step 4.3 Using likelihood to select best set of pi and rho
  #parameter=Parameter_Define(initial_mu_var,pirhopair,choice)####Defining parameters
  
  if(length(choice)!=1){choice=sample(choice,1,prob=c(rep(1/length(choice),length(choice))))}###in case there are multiple choices
  
  print("Iteration Begins~")
  if(algorithms=="EM")
  {
    mclust=Mclust(rstat)
    if (length(mclust$parameter$mean)==1) {warning("warning: the input is not appropriate for mclust since only one cluster was detected by the function Mclust,Please try other methods" )}
    if (length(mclust$parameter$mean)==1) break
    hodcmclust=HODCMclust(mclust,rstat)        
    total=Iteration3_Mclust(iter,wholeindex,hodcmclust,net,pirhopair,choice,rstat,show.steps,showlikelihood,likelihood.frequency)
  }else if (algorithms=="DPM"){
    dpdensitycluster=DPdensitycluster(v,rstat,DPM.mcmc,DPM.prior)
    if (DPparallel==FALSE) {total=Iteration3_DPdensity(iter,wholeindex,dpdensitycluster,net,pirhopair,choice,rstat,v,show.steps)
    }else{total=Iteration3_DPdensity_Par(iter,wholeindex,dpdensitycluster,net,pirhopair,choice,rstat,v,show.steps,n.cores)}
  }
  mcmcdata=sapply(1:iter, function(kk) return(t(matrix(total[kk,]))%*%matrix(rstat)))
  mcmcfile=mcmc(data=mcmcdata)
  convergence=heidel.diag(mcmcfile, eps=0.1, pvalue=0.05)
  graph <- graph.adjacency(net,mode="undirected")
  networks.fast=list()
  networks.fast$trace=total[(nburns+1):iter,]
  networks.fast$HyperParameter$pi0=pirhopair$pi0[choice]
  networks.fast$HyperParameter$rho0=pirhopair$rho0[choice]
  networks.fast$HyperParameter$rho1=pirhopair$rho1[choice]
  networks.fast$convergence=convergence
  networks.fast$graph=graph
  networks.fast$statistics$mean=sapply(1:length(pvalue), function(kk) return(mean(total[,kk])))
  networks.fast$statistics$median=sapply(1:length(pvalue), function(kk) return(median(total[,kk])))
  networks.fast$statistics$var=sapply(1:length(pvalue), function(kk) return(networks.fast$parameter$mean[kk]*(1-networks.fast$parameter$mean[kk])))
  networks.fast$statistics$quantile=sapply(1:length(pvalue), function(kk) return(quantile(total[,kk])))
  
  #  if (trace==FALSE){
  #    show(networks.fast$parameter)       
  #  }else if (trace==TRUE){
  #    if (show.steps=="All") {
  #      show(networks.fast)
  #    }else{steps=c()
  #          for (kk in 1:nrow(networks.fast$trace))
  #          {
  #            if(kk%%show.steps==0){steps=rbind(steps,total[kk,])}
  #          }
  #          showdata=list()
  #          showdata$steps=steps
  #          show(showdata)
  #          show(networks.fast$parameter)
  #          
  #    }
  #  }
  
  
  structure(networks.fast,class="Networks.Fast")
  
  
  
  
}
