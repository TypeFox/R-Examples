loglikPC<-function(lambda,X,nobj,partsList,ENV)
{
      # print dot at each iteration (max 50/line)
      if(ENV$pr.it) {
          if (ENV$iter>50){
            cat("\n")
            ENV$iter<-0
          }
          ENV$iter<-ENV$iter+1
          cat(".")
          flush.console()
      }

      # calculates normalisation if NI==TRUE
      if(ENV$NI) ENV$nrm.all<-normalize(lambda,X,nobj,ENV)

      # calculates contribution to loglik for each Cov group
      ccontrib<-lapply(partsList, NIgroupscontrib, lambda, X,nobj,ENV)

      ll<-sum(unlist(lapply(ccontrib,function(x)lapply(x,function(y) y$ll))))
      ENV$fl<-sum(unlist(lapply(ccontrib,function(x)lapply(x,function(y) y$fl))))
      ENV$ll<-ll
     -ll
}
