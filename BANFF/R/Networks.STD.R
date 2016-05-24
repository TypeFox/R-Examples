#'Standard Algorithm for Bayesian Network Discovery
#'@param pvalue a vector of p-values obtained from large scale statistical hypothesis testing
#'@param net a n by n network configuration, n is the length of pvalue
#'@param iter number of iterations. The default is 5000
#'@param nburns number of burn-in. The default is 2000
#'@param piall a vector of selections of pi0. The default vector is 0.75, 0.8, 0.85, 0.9. The selections of pi0 should be placed in sequence, from smaller to larger.
#'@param rhoall a vector of selections of rho0 and rho1. The default vector is 0.5, 1, 5, 10, 15. The selections of rho0 and rho1 should be placed in sequence, from smaller to larger.
#'@details This generic function fits a Bayesian Nonparametric Mixture Model for gene selection incorporating network information (Zhao et al., 2014):
#' \itemize{
#' \item  r_i| g_i, \strong{theta} ~ N(mu_{g_i}, sigma_{g_i}),
#' \item	g_i | z_i=k, \strong{q}_k ~ Discrete(\strong{a}_k, \strong{q}_k),
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
#' For this function, The NET-DPM-1, considered as standard function is applied , and more details about the algorithm can be referred from Appnendix B.1 of Zhao et al., 2014
#'@return The trace of gi showing the evolution of the Monte Carlo Markov Chain
#'@examples #' ##Creating the network of 10X10 image
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
#'  
#'  }else{
#'   testcov[i]<-rnorm(1,mean=1,sd=0.2)
#'    
#'  }
#' }
#' ##The profile of the impage
#' image(matrix(testcov,10,10),col=gray(seq(0,1,length=255)))
#' ##Transform the signals into pvalue form and begin identification
#' pvalue=pnorm(-testcov)
#' total=Networks.STD(pvalue,net,iter=5000,piall=c(0.8, 0.85, 0.9, 0.95),rhoall=c(0.5,1,5,10,15))
#'@references Zhao, Y.*, Kang, J., Yu, T. A Bayesian nonparametric mixture model for gene and gene-sub network selection
#' Annals of Applied Statistics, In press: 2014. 
#'@export
Networks.STD=function(pvalue,net,iter=5000, nburns=2000, 
                      piall=c(0.75, 0.8, 0.85, 0.9),rhoall=c(0.5, 1, 5, 10, 15),
                      status=FALSE,fit,show.steps=1,showlikelihood=FALSE,likelihood.frequency=100){
  if(status==FALSE){
    ####Preparing Stage:Defining Parameters, 
    rstat=Transfer(pvalue)###Step1 : Transferring 
    wholeindex=Kmeans(rstat)###Step2: Using Kmeans for giving the initial values
    initial_mu_var=Initial_mu_var(rstat,wholeindex) ###Step3: Genering initial mu and var based on intialindex
    pirhopair=Pi_rho_gen(piall,rhoall)####Step4.1 Generating all the possible selections of pi and rho
    znew=Generating_Zi(net,pirhopair,wholeindex)####Step4.2 Generating Zi by different setting of pi and rho.
    choice=Pi_rho_selecting(znew,rstat,piall,rhoall)####Step 4.3 Using likelihood to select best set of pi and rho
    parameter=Parameter_Define(initial_mu_var,pirhopair,choice)####Defining parameters
    #parameter$pi0=0.85
    #parameter$pi1=0.15
    #parameter$rho0=0.5
    #parameter$rho1=1
    #return(parameter)
    ####Iteration
    model=Model(wholeindex,net,rstat,parameter,initial_mu_var)#####Defining the values which would be used in the later iterations
  }else{model=fit}
  ####iteration begin~
  sample.save = matrix(NA,ncol=length(model$wholeindex),nrow=iter-nburns)
  for (jj in 1:iter)
  {
    
    if(jj%%show.steps==0){
      cat("iteration:",jj,"\n")
      flush.console()}
    for (num1 in 1:length(model$wholeindex)){
      model=Step1_1_Update_gi_zi(model,num1)
      model=Step1_2_Check(model=model,num1)
      model=Step_2_Update_mu_var(model,num1)
      model=Step_3_switch_the_label(model,num1)
    }
    if(jj>nburns){sample.save[jj-nburns,]=model$wholeindex}
    
    if(showlikelihood==TRUE){
      if(jj%%likelihood.frequency==0){
        mylog<-0
        mu0=mean(rstat[which(model$wholeindex<=0)])
        mu1=mean(rstat[which(model$wholeindex>=1)])
        var0=var(rstat[which(model$wholeindex<=0)])
        var1=var(rstat[which(model$wholeindex>=1)])
        for(num1 in 1:length(rstat)){
          
          if(model$wholeindex[num1]==0){
            mylog[num1]<--(rstat[num1]-mu0)^2/(2*var0)-log(sqrt(var0))
          }else{
            mylog[num1]<--(rstat[num1]-mu1)^2/(2*var1)-log(sqrt(var1))
          }
          
        }
        cat("Now for the step:" ,jj, "the log-likelihood value is" ,sum(mylog) , "\n")
        flush.console()
      }
    }
    
  }
  
  graph <- graph.adjacency(net,mode="undirected")
  mcmcdata=sapply(1:nrow(sample.save), function(kk) return(t(matrix(sample.save[kk,]))%*%matrix(rstat)))
  mcmcfile=mcmc(data=mcmcdata)
  convergence=heidel.diag(mcmcfile, eps=0.1, pvalue=0.05)
  results=list()
  results$trace=sample.save
  results$convergence=convergence
  results$model=model
  results$graph=graph
  
  
  
  structure(results,class="Networks.STD")
  
}


