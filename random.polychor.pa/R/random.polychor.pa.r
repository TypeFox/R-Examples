random.polychor.pa <- function (nvar = "NULL", n.ss = "NULL", nrep, nstep = "NULL", data.matrix, q.eigen, r.seed = "NULL", diff.fact=FALSE, distr="NULL", comparison = "random", fit.pa=FALSE, print.all=FALSE) 
{
  ### PRELIMINAR CHECK: START
  start.t <- Sys.time()
  cat("***** RESULTS FOR PARALLEL ANALYSIS *****", "\n")
  cat("*** computation starts at:", format(start.t, "%X"), "\n")
  flush.console()
  #    cat("\n")
  
  if (distr!="uniform" & distr!="multinomial" & distr!="NULL" ) { 
    stop("Bad parameter specification. Please provide a valid option for comparison parameter")
  }
  else{}
  if (distr=="uniform" & diff.fact==TRUE) { 
    stop("Bad parameter specification. Uniform and difficulty factor not implemented")
  }
  else{}
  if (distr=="multinomial" & diff.fact==TRUE) { 
    stop("Bad parameter specification. Multinomial and difficulty factor not implemented")
  }
  else{}
  if (comparison=="bootstrap" & diff.fact==TRUE) { 
    stop("Bad parameter specification. Bootstrap and difficulty factor not implemented")
  }
  else{}
  
  if (comparison!="random" & comparison!="bootstrap" & comparison!="random-mg" & comparison!="bootstrap-mg") { 
    stop("Bad parameter specification. Please provide a valid option for comparison parameter")
  }
  else{}
  

  data.matrix.0 <- na.omit(data.matrix)
  n.ss <- nrow(data.matrix)
  cat("*** number of units (rows) in data.matrix:", n.ss, "\n")
  if (nrow(data.matrix.0) < nrow(data.matrix)) {
    cat("*** LISTWISE deletion needed. The new sample size is:", nrow(data.matrix.0), "\n")
    data.matrix <- data.matrix.0
    n.ss <- nrow(data.matrix.0)
  }
  else{      
    cat("*** No missing values found", "\n")
  }
  
  if (comparison=="random-mg" | comparison=="bootstrap-mg") { # if multisample is TRUE then the first variable of the dataset is the factor variable that will be used to distinguish among sub-sample
    multi.sample<-data.matrix[,1]
    cat("*** Number of samples and size of samples to be compared: ", "\n") # prints the number of samples to be used
    print(table(as.factor(multi.sample)))
    cat(" ", "\n") # prints the number of samples to be used
    data.matrix <- data.matrix[,-1]
  }
  else{}
  
  if (comparison=="random" | comparison=="bootstrap") { # if random is TRUE then a constant will be added to the dataset for compatibility with the multi-sample algorithm
    multi.sample<-cbind(rep(1, nrow(data.matrix)))
    cat("*** SINGLE sample Parallel Analysis", "\n") # prints the number of samples to be used
  }
  else{}
  if (comparison=="bootstrap" | comparison=="bootstrap-mg") {
    cat("*** simulation method: BOOTSTRAP (permutation)", "\n")
    cat("*** difficulty factor: FALSE","\n")
  }
  else{
    cat("*** simulation method: RANDOM", "\n")      
    if(distr=="uniform" | distr=="NULL") {
      cat("*** distribution: UNIFORM","\n")
      if(diff.fact==TRUE) {
        cat("*** difficulty factor: TRUE","\n")
      }
      else{
        cat("*** difficulty factor: FALSE","\n")
      }
    }
    else{
      cat("*** distribution: MULTINOMIAL","\n")
      cat("*** difficulty factor: FALSE","\n")
    }
  }  
  
  for (z in 1:ncol(data.matrix)) {
    if (is.numeric(data.matrix[, z]) == FALSE) {
      data.matrix[, z] <- as.numeric(data.matrix[, z])
    }
    else{}
  }
  
  data.matrix <- as.matrix(data.matrix)
  if (!is.null(dimnames(data.matrix))) {
    dimnames(data.matrix) <- list(NULL, NULL)
    data.matrix
  }
  else{}
  if (is.null(q.eigen)) {
    stop("the quantile to be extracted is not declared. Please provide a valid quantile")
  }
  else{}
  
  if (q.eigen < 0 | q.eigen >= 1) {
    stop("the quantile to be extracted is out of bound. Please provide a valid quantile")
  }
  else{}
  
  if (r.seed == "NULL") {
    set.seed(1335031435)
  }
  else {
    set.seed(r.seed)
  }
  
  nvar <- ncol(data.matrix)
  cat("*** number of variables (cols) in data.matrix:", nvar, "\n")
  flush.console()
  item.tab.ex <- matrix(, ncol(data.matrix), 4) 
  for (h in 1:ncol(data.matrix)) {
    item.tab.ex[h, 2] <- (max(data.matrix[, h])-min(data.matrix[ , h])+1)
    item.tab.ex[h, 3] <- min(data.matrix[, h])
    item.tab.ex[h, 4] <- max(data.matrix[, h])
    item.tab.ex[h, 1] <- h
  }
  flag<-matrix(0, nrow(item.tab.ex),2) 
  flag[,1]<-1:nrow(item.tab.ex)
  item.tab<-matrix(0,1,4) 
  item.tab.row<-matrix(0,1,4) 
  i<-1
  item.tab[i,2:4] <- item.tab.ex[i,2:4] 
  while(sum(flag[,2]) < nrow(item.tab.ex)) {
    if(flag[i,2] == 0) {  
      for (j in 1:nrow(item.tab.ex)) {
        if((flag[j,2] == 0) & (item.tab.ex[j,2]==item.tab[nrow(item.tab),2]) & (item.tab.ex[j,3]==item.tab[nrow(item.tab),3]) & (item.tab.ex[j,4]==item.tab[nrow(item.tab),4])) { 
          item.tab[nrow(item.tab),1]<-(item.tab[nrow(item.tab),1])+1 
          flag[j,2] <- 1   
        }
        else{}
      }
    }
    else{}
    if(i+1 == nrow(item.tab.ex)+1) break
    else { 
      if(flag[i+1,2] == 0) {  
        i<-i+1  
        item.tab<-rbind(item.tab,item.tab.row) 
        item.tab[nrow(item.tab),2:4] <- item.tab.ex[i,2:4] 
      }
      else i<-i+1
    }
  }
  colnames(item.tab) <- c("Items", "Categories", "Min.Cat", "Max.Cat")
  row.labels <- paste(c(1:nrow(item.tab)), "GROUP")
  rownames(item.tab) <- c(row.labels)
  cat("*** Groups of items with diffent number of categories found in your data.matrix:", "\n")
  print(item.tab)
  flush.console()
  cat("\n")
  
  item.tab <- as.matrix(item.tab)
  if (!is.null(dimnames(item.tab))) {
    dimnames(item.tab) <- list(NULL, NULL)
    item.tab
  }
  else{}
  
  ### PRELIMINAR CHECK: END
  
  table(multi.sample)
  
  sub.sample<-matrix(, , sum(dim(table(multi.sample)))) # matrix that will contain the names of the sub-groups
  print.perm.eigen.polyc.fa<-list()
  print.perm.eigen.polyc.pca<-list()
  print.perm.eigen.pear.fa<-list()
  print.perm.eigen.pear.pca<-list()
  print.diff.weight<-list() # print all differential weights used
  print.nr.fact <- list() # number of factors selected by PA
  print.matrix.multisample.pol<-list() # matrices for polychoric FA
  print.matrix.multisample.pea<-list() # matrices for pearson FA
  print.matrix.pca.multisample.pol<-list() # matrices for polychoric PCA
  print.matrix.pca.multisample.pea<-list() # matrices for pearson PCA
  map.result.multisample<-list() # matrices for MAP
  table.fit.res<-matrix(,sum(dim(table(multi.sample)))*2, 10) # matrices for fit indexes
  index2<-0 # index for the FIT INDEXES
  table.pa.res<-matrix(,8,sum(dim(table(multi.sample)))) # matrices for MAP and PA solutions
  
  ### BEGINNING the loop that is repeated for each sub-sample: START
  for (w in 1:sum(dim(table(multi.sample)))) {  
    
    categ<-levels(as.factor(multi.sample))[w]  # selecting each value of the multisample factor 
    
    ### selecting the sub-sample
    data.matrix.sub<-subset(data.matrix, multi.sample==categ)
    ### computing the number of subjects and variables
    nvar.sub <- ncol(data.matrix.sub)
    n.ss.sub <- nrow(data.matrix.sub)
    
    ######### RANDOM MULTISAMPLE POLYCHORIC.PA: START        
    if (comparison=="random" | comparison=="random-mg") {
      
      ### COMPUTING DIFFERENTIAL WEIGHTS: START
      diff.weigth<-function(dataset) {
        ### just for printing
        for (g in 1:nvar.sub) {
          ### START: loop for each variable
          for (q in 1:length(table(dataset[,g]))){
            p.val<-table(dataset[,g])[q]/n.ss.sub
            categ<-q-1
            item<-g
            terna<-cbind(p.val,categ,item)
            if(g==1 & q==1){
              pre.peso<-terna
            }
            else {
              pre.peso<-rbind(pre.peso, terna)
            }
          }
        }
        return(pre.peso)
      }
      peso<-diff.weigth(dataset = data.matrix.sub)
      colnames(peso)<-c("p-value", "cateogry", "variable")
      ### COMPUTING DIFFERENTIAL WEIGHTS: END     
      
      ### generating random samples: START
      sim.random.matrix<-function() {
        pre.matrix.3 <- matrix(0, n.ss.sub, nvar.sub)
        u <- 1
        t <- 1
        if(diff.fact==TRUE & distr=="NULL"){
          ### calculating p-value for each item: POLYTHOMOUS
          for (q in 1:nvar.sub) {
            ### START: loop for each variable
            peso<-matrix(,length(table(data.matrix.sub[,q])),3)
            for (p in 1:length(table(data.matrix.sub[,q]))){
              peso[p,1]<-table(data.matrix.sub[,q])[p]/n.ss.sub
              peso[p,2]<-p-1
            }
            if (length(table(data.matrix.sub[,q])) == 2) {
              pre.matrix.3[,q]<-sample(x=peso[,2], size = n.ss.sub, replace = TRUE, prob=peso[,1])
            }
            else{
              pre.matrix.3[,q]<-sample(x=peso[,2]+1, size = n.ss.sub, replace = TRUE, prob=peso[,1])
            }
            ### END: for each variable
          }
        }
        else {
          if((distr=="uniform" | distr=="NULL") & diff.fact==FALSE) {
            #cat("UNIFORM distributions will be simulated","\n")
            while (u <= nrow(item.tab)) {
              for (z in 1:(item.tab[u, 1])) {
                pre.matrix.3[, t] <- round(runif(n.ss.sub, min = 1, max = item.tab[u, 2]))
                if (t < sum(item.tab[, 1])) t <- t + 1
              }
              u <- u + 1
            }
          }
          else{}
          if(distr=="multinomial") {
            #cat("MULTINOMIAL distributions will be simulated","\n")
            for (q in 1:nvar.sub) {
              ### START: loop for each variable
              peso<-matrix(,length(table(data.matrix.sub[,q])),3)
              for (p in 1:length(table(data.matrix.sub[,q]))){
                peso[p,1]<-table(data.matrix.sub[,q])[p]/n.ss.sub
                peso[p,2]<-p-1
              }
              draws <- rmultinom(n=1, size=nrow(data.matrix.sub), prob=peso[,1])
              pre.matrix.3[,q]<-sample(x=peso[,2]+1, size = n.ss.sub, replace = TRUE, prob=(draws/nrow(data.matrix.sub)))              
              ### END: for each variable
            }
          }
        }
        return(pre.matrix.3)
      }      
      ### generating random samples: START
        
      
      ### COMPUTING FA eigenvalue distributons from RANDOM polychoric and pearson correlations: START
      random.eigen.fa.distr <- function(dataset) {
        eigen.data <- matrix(0, nvar.sub, nrep)   
        eigen.data1 <- matrix(0, nvar.sub, nrep)
        f1.poly.cor <- matrix(0, nvar.sub, )
        f1.cor <- matrix(0, nvar.sub, )
        pre.st.matrix <- matrix(0, nvar.sub, 9)
        
        for (j in 1:nrep) {
          matrix.3<-sim.random.matrix()
          f1.poly.cor <- suppressMessages(polychoric(matrix.3, global=FALSE)$rho)   # simulated polychoric corr
          f1.cor <- cor(matrix.3)   # simulated pearson corr
          eigen.data[, j] <- eigen(corFA(f1.poly.cor))$values
          eigen.data1[, j] <- eigen(corFA(f1.cor))$values
          if (j == 1 & w == 1) {
            end.pt <- Sys.time()
            estimated.t <- difftime(end.pt, start.t, units="auto")
            estimated.total <- estimated.t * nrep * sum(dim(table(multi.sample)))
            estimated.t <- as.numeric(estimated.t, units = "secs")
            estimated.total <- as.numeric(estimated.total, units = "secs")
            cat("The first simulation for FA took:", round(estimated.t,3), "secs.", "\n")
#            cat("The whole simulation will take no less than:", round(estimated.total/60), "min.", "to terminate", "\n")
            flush.console()
          }
          else{}
        }
        for (col in 1:nvar.sub) {
          eigen.data.t <- t(eigen.data)
          pre.st.matrix[col, 1] <- mean(eigen.data.t[, col])
          pre.st.matrix[col, 2] <- sd(eigen.data.t[, col])
          pre.st.matrix[col, 3] <- quantile(eigen.data.t[, col], 0.95)
          pre.st.matrix[col, 4] <- quantile(eigen.data.t[, col], q.eigen)
          pre.st.matrix[col, 5] <- (col)
          eigen.data1.t <- t(eigen.data1)
          pre.st.matrix[col, 6] <- mean(eigen.data1.t[, col])
          pre.st.matrix[col, 7] <- sd(eigen.data1.t[, col])
          pre.st.matrix[col, 8] <- quantile(eigen.data1.t[, col], 0.95)
          pre.st.matrix[col, 9] <- quantile(eigen.data1.t[, col], q.eigen)
        }
        return(pre.st.matrix)
      }
      st.matrix <- random.eigen.fa.distr(dataset = data.matrix.sub)
      ### COMPUTING FA eigenvalue distributons from RANDOM polychoric and pearson correlations: END
      
      ### COMPUTING PCA eigenvalue distributons from RANDOM polychoric and pearson correlations: START    
      random.eigen.pca.distr <- function(dataset) {
        eigen.data.pca <- matrix(0, nvar.sub, nrep)
        eigen.data1.pca <- matrix(0, nvar.sub, nrep)
        f1.poly.cor <- matrix(0, nvar.sub, )
        f1.cor <- matrix(0, nvar.sub, )
        pre.st.matrix.pca <- matrix(0, nvar.sub, 9)
        
        for (j in 1:nrep) {
          matrix.3<-sim.random.matrix()
          f1.poly.cor <- suppressMessages(polychoric(matrix.3, global=FALSE)$rho)   # simulated polychoric corr
          f1.cor <- cor(matrix.3)   # simulated pearson corr
          eigen.data.pca[, j] <- eigen(f1.poly.cor)$values
          eigen.data1.pca[, j] <- eigen(f1.cor)$values
          if (j == 1 & w == 1) {
            end.pt <- Sys.time()
            estimated.t <- difftime(end.pt, start.t, units="auto")
            estimated.total <- estimated.t * nrep * sum(dim(table(multi.sample)))
            estimated.t <- as.numeric(estimated.t, units = "secs")
            estimated.total <- as.numeric(estimated.total, units = "secs")
            cat("The first simulation for PCA took:", round(estimated.t,3), "secs.", "\n")
#            cat("The whole simulation will take no less than:", round(estimated.total/60), "min.", "to terminate", "\n")
            flush.console()
          }
          else{}
        }
        for (col in 1:nvar.sub) {
          eigen.data.pca.t <- t(eigen.data.pca)
          pre.st.matrix.pca[col, 1] <- mean(eigen.data.pca.t[, col])
          pre.st.matrix.pca[col, 2] <- sd(eigen.data.pca.t[, col])
          pre.st.matrix.pca[col, 3] <- quantile(eigen.data.pca.t[, col], 0.95)
          pre.st.matrix.pca[col, 4] <- quantile(eigen.data.pca.t[, col], q.eigen)
          pre.st.matrix.pca[col, 5] <- (col)
          eigen.data1.pca.t <- t(eigen.data1.pca)
          pre.st.matrix.pca[col, 6] <- mean(eigen.data1.pca.t[, col])
          pre.st.matrix.pca[col, 7] <- sd(eigen.data1.pca.t[, col])
          pre.st.matrix.pca[col, 8] <- quantile(eigen.data1.pca.t[, col], 0.95)
          pre.st.matrix.pca[col, 9] <- quantile(eigen.data1.pca.t[, col], q.eigen)
        }
        return(pre.st.matrix.pca)
      }
      st.matrix.pca <- random.eigen.pca.distr(dataset = data.matrix.sub)
    }
    ### COMPUTING PCA eigenvalue distributons from RANDOM polychoric and pearson correlations: START
    
    
    ### random PERMUTATIONS of cases: START
    ######### RANDOM MULTISAMPLE POLYCHORIC.PA: START        
    if (comparison=="bootstrap" | comparison=="bootstrap-mg") {
      
      ### this function is called within the boot function        
      bootstrap.eigen.polyc.fa <- function(d, i) {
        boot.matrix[i,1]<-d[,1][i]
        for (j in 2:ncol(d)) {
          column<-d[,j]
          boot.matrix[,j]<-column[i]
        }
        # cat("\n", "BOOTSTRAP sampled FA eigenvalues from POLYCHORIC correlations", "\n")
        boot.polyc<-suppressMessages(polychoric(boot.matrix, global=FALSE)$rho) ### polychoric correlation
        eigen.polyc.fa<-eigen(corFA(boot.polyc))$values       ### FA
        return(eigen.polyc.fa)
      }
      
      bootstrap.eigen.polyc.pca <- function(d, i) {
        boot.matrix[i,1]<-d[,1][i]
        for (j in 2:ncol(d)) {
          column<-d[,j]
          boot.matrix[,j]<-column[i]
        }
        # cat("\n", "BOOTSTRAP sampled PCA eigenvalues from POLYCHORIC correlations", "\n")
        boot.polyc<-suppressMessages(polychoric(boot.matrix, global=FALSE)$rho) ### polychoric correlation
        eigen.polyc.pca<-eigen(boot.polyc)$values             ### PCA
        return(eigen.polyc.pca)
      }
      
      bootstrap.eigen.pear.fa <- function(d, i) {
        boot.matrix[i,1]<-d[,1][i]
        for (j in 2:ncol(d)) {
          column<-d[,j]
          boot.matrix[,j]<-column[i]
        }
        # cat("\n", "BOOTSTRAP sampled FA eigenvalues from PEARSON correlations", "\n")
        boot.pear<-cor(boot.matrix)                   ### pearson correlation
        eigen.pear.fa<-eigen(corFA(boot.pear))$values ### FA
        return(eigen.pear.fa)
      }
      
      bootstrap.eigen.pear.pca <- function(d, i) {
        boot.matrix[i,1]<-d[,1][i]
        for (j in 2:ncol(d)) {
          column<-d[,j]
          boot.matrix[,j]<-column[i]
        }
        # cat("\n", "BOOTSTRAP sampled PCA eigenvalues from PEARSON correlations", "\n")
        boot.pear<-cor(boot.matrix)             ### pearson correlation
        eigen.pear.pca<-eigen(boot.pear)$values ### PCA
        return(eigen.pear.pca)
      }
      ### RANDOM PERMUTATIONS of cases: END
      
      ### BOOTSTRAP: START
      boot.matrix<-matrix(0, nrow(data.matrix.sub), ncol(data.matrix.sub))    
      b.perm.eigen.polyc.fa <- boot(data.matrix.sub, bootstrap.eigen.polyc.fa, R=nrep, sim="permutation", stype="i")
      boot.matrix<-matrix(0, nrow(data.matrix.sub), ncol(data.matrix.sub))
      b.perm.eigen.polyc.pca <- boot(data.matrix.sub, bootstrap.eigen.polyc.pca, R=nrep, sim="permutation", stype="i")
      boot.matrix<-matrix(0, nrow(data.matrix.sub), ncol(data.matrix.sub))
      b.perm.eigen.pear.fa <- boot(data.matrix.sub, bootstrap.eigen.pear.fa, R=nrep, sim="permutation", stype="i")
      boot.matrix<-matrix(0, nrow(data.matrix.sub), ncol(data.matrix.sub))
      b.perm.eigen.pear.pca <- boot(data.matrix.sub, bootstrap.eigen.pear.pca, R=nrep, sim="permutation", stype="i")
      ### BOOTSTRAP: END
      
      sub.sample[,w]<-paste('sample',categ, sep=".")      # list of names for sub-groups
      print.perm.eigen.polyc.fa[sub.sample[,w]]<-list(capture.output(print(b.perm.eigen.polyc.fa)))
      print.perm.eigen.polyc.pca[sub.sample[,w]]<-list(capture.output(print(b.perm.eigen.polyc.pca)))
      print.perm.eigen.pear.fa[sub.sample[,w]]<-list(capture.output(print(b.perm.eigen.pear.fa)))
      print.perm.eigen.pear.pca[sub.sample[,w]]<-list(capture.output(print(b.perm.eigen.pear.pca)))
      
      ############ computing the percentile from boostraped eigenvalues: START      
      perm.eigen.polyc.fa<-matrix(0,ncol(b.perm.eigen.polyc.fa$t),1)
      perm.eigen.pear.fa<-matrix(0,ncol(b.perm.eigen.pear.fa$t),1)
      perm.eigen.polyc.pca<-matrix(0,ncol(b.perm.eigen.polyc.pca$t),1)
      perm.eigen.pear.pca<-matrix(0,ncol(b.perm.eigen.pear.pca$t),1)
      
      for(n in 1:ncol(data.matrix.sub)) {
        perm.eigen.polyc.fa[n,1] <- quantile(b.perm.eigen.polyc.fa$t[,n], q.eigen)
        perm.eigen.pear.fa[n,1] <- quantile(b.perm.eigen.pear.fa$t[,n], q.eigen)
        perm.eigen.polyc.pca[n,1] <- quantile(b.perm.eigen.polyc.pca$t[,n], q.eigen)
        perm.eigen.pear.pca[n,1] <- quantile(b.perm.eigen.pear.pca$t[,n], q.eigen)
      }
      ############ computing the 95th percentile from boostraped eigenvalues: END          
    }
    ### permutations: END
    
    
    ###### PUT THE FOLLOWING RESULTS IN A NEW NAMED COMMON MATRIX (COMMON TO RANDOM AND BOOSTRAP)
    
    ### COMPUTING EMPIRICAL CORRELATION MATRICES: START 
    matrix.cor1 <- suppressMessages(polychoric(data.matrix.sub, global=FALSE)$rho)
    matrix.cor2 <- cor(data.matrix.sub)
    ### COMPUTING EMPIRICAL CORRELATION MATRICES: END
    ### COMPUTING EMPIRICAL EIGENVALUES FA/PCA RESULTS: START
    raw.eigen.poly.fa<- eigen(corFA(matrix.cor1))$values  # Empirical FA with polychoric corr 
    raw.eigen.pear.fa<- eigen(corFA(matrix.cor2))$values  # Empirical FA with Pearson corr
    raw.eigen.poly.pca<- eigen(matrix.cor1)$values        # Empirical PCA with polychoric corr
    raw.eigen.pear.pca<- eigen(matrix.cor2)$values        # Empirical PCA with Pearson corr
    
    if(comparison=="random" | comparison=="random-mg") {
      st.matrix <- cbind(st.matrix, raw.eigen.poly.fa, raw.eigen.pear.fa) # Empirical+FA+polychoric/pearson
      colnames(st.matrix) <- c("P.SimMeanEigen", "P.SimSDEigen", "P.Sim95perc", "P.SimQuant", "Factor", 
                               "C.SimMeanEigen", "C.SimSDEigen", "C.Sim95perc", "C.SimQuant",
                               "Emp.Polyc.Eigen", "Emp.Pears.Eigen")  #  st.matrix
      st.matrix.pca <- cbind(st.matrix.pca, raw.eigen.poly.pca, raw.eigen.pear.pca) # Empirical+PCA+polychoric/pearson
      colnames(st.matrix.pca) <- c("P.SimMeanEig.PCA", "P.SimSDEig.PCA", "P.Sim95perc.PCA", "P.SimQuant.PCA", "Comp", 
                                   "C.SimMeanEig.PCA", "C.SimSDEig.PCA", "C.Sim95perc.PCA", "C.SimQuant.PCA",
                                   "Emp.Pol.Eig.PCA", "Emp.Pear.Eig.PCA")  #  st.matrix.pca
    }
    else{
      st.matrix <- cbind(perm.eigen.polyc.fa, perm.eigen.pear.fa, raw.eigen.poly.fa, raw.eigen.pear.fa) # Empirical+FA+polychoric/pearson
      colnames(st.matrix) <- c("P.SimQuant", "C.SimQuant", "Emp.Polyc.Eigen", "Emp.Pears.Eigen")  #  st.matrix
      st.matrix.pca <- cbind(perm.eigen.polyc.pca, perm.eigen.pear.pca, raw.eigen.poly.pca, raw.eigen.pear.pca) # Empirical+PCA+polychoric/pearson
      colnames(st.matrix.pca) <- c("P.SimQuant.PCA", "C.SimQuant.PCA", "Emp.Pol.Eig.PCA", "Emp.Pear.Eig.PCA")  #  st.matrix.pca
    }
    
    ### COMPUTING EMPIRICAL EIGENVALUES FA/PCA RESULTS: END
    
    
    ### COMPARE EIGENVALUES AND COMPUTE THE NUMBER OF FACTORS TO RETAIN WITH PA: START
    pa.polychor <- function(poly.matrix, corr.matrix, fa.matrix, pca.matrix) {
      # poly.matrix: empirical polychoric correlation
      # corr.matrix: empirical pearson correlation
      # fa.matrix:   sim FA eigenvalues
      # pca.matrix:  sim PCA eigenvalues
      
      ### COMPUTING NUMBER OF FACTORS TO RETAIN WITH PA: START
      diff.eigen<-matrix(, nvar.sub, 4) # differenze tra gli eigenvalue emp-sim
      nr.fact <- matrix(0,4,1) # numero di fattori scelti con PA
      
      ### differences between empirical and simulated eigenvalues: START
      diff.eigen[, 1]<-((fa.matrix[, 3] - fa.matrix[, 1])>0)*1   # FA+polychoric
      diff.eigen[, 2]<-((fa.matrix[, 4] - fa.matrix[, 2])>0)*1   # FA+Pearson
      diff.eigen[, 3]<-((pca.matrix[, 3] - pca.matrix[, 1])>0)*1 # PCA+polychoric
      diff.eigen[, 4]<-((pca.matrix[, 4] - pca.matrix[, 2])>0)*1 # PCA+Pearson
      ### differences between empirical and simulated eigenvalues: END
      ### counting the number of factors to be retained: START
      nr.fact <- matrix(0,4,1)
      for(count in 1:4){
        coin <- 1
        righetta <- 1
        while (coin == 1) {
          if(righetta<nvar.sub){
            nr.fact[count, 1] <- righetta
            coin <- diff.eigen[righetta, 1]>0
            righetta <- righetta + 1
          }
          else{
            coin<-0
          }
        }
      }
      ### counting the number of factors to be retained: END
      
      rownames(nr.fact) <- c("nr.FA.polyc", "nr.FA.pear", "nr.PCA.polyc", "nr.PCA.pear")
      colnames(nr.fact) <- c("nr.factors")
      return(nr.fact-1)
    }
    if(comparison=="random" | comparison=="random-mg") {
      pa.result <- pa.polychor(poly.matrix = matrix.cor1, 
                               corr.matrix = matrix.cor2,
                               fa.matrix = cbind(st.matrix[, 4],st.matrix[, 9:11]),
                               pca.matrix = cbind(st.matrix.pca[, 4],st.matrix.pca[, 9:11]))
    }
    else{
      pa.result <- pa.polychor(poly.matrix = matrix.cor1, 
                               corr.matrix = matrix.cor2,
                               fa.matrix = cbind(perm.eigen.polyc.fa, perm.eigen.pear.fa, raw.eigen.poly.fa, raw.eigen.pear.fa),
                               pca.matrix = cbind(perm.eigen.polyc.pca, perm.eigen.pear.pca, raw.eigen.poly.pca, raw.eigen.pear.pca))          
    }
    ### COMPARE EIGENVALUES AND COMPUTE THE NUMBER OF FACTORS TO RETAIN WITH PA: END
    
    
    ### COMPUTING map: START
    map.polychor <- function(poly.matrix, corr.matrix) {
      loadings.map <- matrix(0, nvar.sub, nvar.sub)
      loadings.map1 <- matrix(0, nvar.sub, nvar.sub)
      fm <- matrix(0, nvar.sub, 5)
      eigen.map <- eigen(poly.matrix)
      eigen.map1 <- eigen(corr.matrix)
      loadings.map <- (eigen.map$vectors %*% (sqrt(diag(eigen.map$values))))
      loadings.map1 <- (eigen.map1$vectors %*% (sqrt(diag(eigen.map1$values))))
      fm[1, 2] <- (sum(poly.matrix^2) - nvar.sub)/(nvar.sub * (nvar.sub - 1))
      fm[1, 3] <- (sum(poly.matrix^4) - nvar.sub)/(nvar.sub * (nvar.sub - 1))
      fm[1, 4] <- (sum(corr.matrix^2) - nvar.sub)/(nvar.sub * (nvar.sub - 1))
      fm[1, 5] <- (sum(corr.matrix^4) - nvar.sub)/(nvar.sub * (nvar.sub - 1))
      for (m in 1:(nvar.sub - 1)) {
        A <- loadings.map[, 1:m]
        partcov <- poly.matrix - (A %*% t(A))
        d <- diag(1/sqrt(diag(partcov)))
        pr <- d %*% partcov %*% d
        fm[m + 1, 2] <- (sum(pr^2) - nvar.sub)/(nvar.sub * (nvar.sub - 1))
        fm[m + 1, 3] <- (sum(pr^4) - nvar.sub)/(nvar.sub * (nvar.sub - 1))
        A1 <- loadings.map1[, 1:m]
        partcov1 <- corr.matrix - (A1 %*% t(A1))
        d1 <- diag(1/sqrt(diag(partcov1)))
        pr1 <- d1 %*% partcov1 %*% d1
        fm[m + 1, 4] <- (sum(pr1^2) - nvar.sub)/(nvar.sub * (nvar.sub - 1))
        fm[m + 1, 5] <- (sum(pr1^4) - nvar.sub)/(nvar.sub * (nvar.sub - 1))
      }
      minfm.map <- fm[1, 2]
      minfm4.map <- fm[1, 3]
      minfm.map1 <- fm[1, 4]
      minfm4.map1 <- fm[1, 5]
      nfactors.map <- 0
      nfactors4.map <- 0
      nfactors.map1 <- 0
      nfactors4.map1 <- 0
      for (s in 1:nrow(fm)) {
        fm[s, 1] <- s - 1
        if (fm[s, 2] < minfm.map) {
          minfm.map = fm[s, 2]
          nfactors.map = s - 1
        }
        else{}
      }
      for (s in 1:nrow(fm)) {
        fm[s, 1] <- s - 1
        if (fm[s, 3] < minfm4.map) {
          minfm4.map = fm[s, 3]
          nfactors4.map = s - 1
        }
        else{}
      }
      for (s in 1:nrow(fm)) {
        fm[s, 1] <- s - 1
        if (fm[s, 4] < minfm.map1) {
          minfm.map1 = fm[s, 4]
          nfactors.map1 = s - 1
        }
        else{}
      }
      for (s in 1:nrow(fm)) {
        fm[s, 1] <- s - 1
        if (fm[s, 5] < minfm4.map1) {
          minfm4.map1 = fm[s, 5]
          nfactors4.map1 = s - 1
        }
        else{}
      }
      colnames(fm) <- c("Factor", "POLY.MAP.squared", "POLY.MAP.4th", "CORR.MAP.squared", "CORR.MAP.4th")
      max.nr <- rbind(nfactors.map, nfactors4.map, nfactors.map1, nfactors4.map1)
      return(fm[(1:(max(max.nr) + 2)), ])
    }
    map.result <- map.polychor(poly.matrix = matrix.cor1, corr.matrix = matrix.cor2)
    ### COMPUTING map: START
    
    
    
    ### STORING FIT STATISTICS FOR DIFFERENT FACTOR SOLUTIONS: START
    if(fit.pa==TRUE){
      q<-max(1,(pa.result)) # number of factors for PA
      if(q>1)  {table.fit.res<-rbind(table.fit.res, matrix(,(q-1)*2, 10))} # adds 4 rows each extracted factor 
      type.cor<-matrix.cor1 
      type.label<-rbind("emp-polyc-corr","emp-pears-corr")
      ### fit indexs are computed for only the empirical data
      index1<-1
      for (t in 1:2) {     # 1) emp-polyc-corr; 2) sim-polyc-corr; 3) emp-pears-corr; 4) sim-pears-corr; 
        if (t==2) {type.cor<-matrix.cor2} # 3) emp-pears-corr; 
        for (n in 1:q) {   # q is the number factor, when zero then the fit is computed only on the first extracted factor
          fit.fa<-fa(type.cor, nfactors=n, n.obs=n.ss.sub, alpha=0.01)  # Empirical+FA+polychoric
          ind<-index2+index1
          table.fit.res[ind,1]<-paste('sample',categ, sep=".") # list of name for the sub-group
          table.fit.res[ind,2]<-type.label[t,]
          table.fit.res[ind,3]<-as.numeric(n)
          table.fit.res[ind,4]<-round(fit.fa$STATISTIC, 3)
          table.fit.res[ind,5]<-round(fit.fa$dof, 3)
          table.fit.res[ind,6]<-round(fit.fa$PVAL, 3)
          table.fit.res[ind,7]<-round(fit.fa$TLI[1], 3)
          table.fit.res[ind,8]<-round(fit.fa$RMSEA[1], 3)
          table.fit.res[ind,9]<-round(fit.fa$rms, 3)
          table.fit.res[ind,10]<-round(fit.fa$BIC, 3)
          #          table.fit.res[ind,11]<-round(fit.fa$R2, 3)
          index1<-index1+1
        }
      }      
      index2<-index2+q*2
    }
    
    ### STORING FIT STATISTICS FOR DIFFERENT FACTOR SOLUTIONS: END
    
    
    ### storing results for each sub-sample: START 
    sub.sample[,w]<-paste('sample',categ, sep=".")      # list of names for sub-groups
    if(comparison=="random" | comparison=="random-mg"){
      print.diff.weight[sub.sample[,w]]<-list(peso)       # differential weights
      print.matrix.multisample.pol[sub.sample[,w]]<-list(cbind(st.matrix[ , c(5, 10, 1, 2, 4)]))
      print.matrix.multisample.pea[sub.sample[,w]]<-list(cbind(st.matrix[ , c(5, 11, 6, 7, 9)]))
      print.matrix.pca.multisample.pol[sub.sample[,w]]<-list(cbind(st.matrix.pca[ , c(5, 10, 1, 2, 4)]))
      print.matrix.pca.multisample.pea[sub.sample[,w]]<-list(cbind(st.matrix.pca[ , c(5, 11, 6, 7, 9)]))
    }
    else{
      print.matrix.multisample.pol[sub.sample[,w]]<-list(cbind(st.matrix[ , c(1, 3)]))
      print.matrix.multisample.pea[sub.sample[,w]]<-list(cbind(st.matrix[ , c(2, 4)]))
      print.matrix.pca.multisample.pol[sub.sample[,w]]<-list(cbind(st.matrix.pca[ , c(1, 3)]))
      print.matrix.pca.multisample.pea[sub.sample[,w]]<-list(cbind(st.matrix.pca[ , c(2, 4)]))
    }
    print.nr.fact[sub.sample[,w]]<-list(pa.result)      # for each sub-group
    map.result.multisample[sub.sample[,w]]<-list(map.result)
    ### storing results for each sub-sample: END         
    
    ### SAVING RESULTS matrix: START
    end.t <- Sys.time()
    elapsed.t <-as.numeric(difftime(end.t, start.t), units = "secs")
    cat("computation ended at:", format(end.t, "%X"), "\n")
    cat("Elapsed Time:", round(elapsed.t/60), "min", "\n")
    
    table.pa.res[1,w]<-which.min(map.result[,2])-1
    table.pa.res[2,w]<-which.min(map.result[,3])-1
    table.pa.res[3,w]<-which.min(map.result[,4])-1
    table.pa.res[4,w]<-which.min(map.result[,5])-1
    table.pa.res[5,w]<-(pa.result[3,1])
    table.pa.res[6,w]<-(pa.result[4,1])
    table.pa.res[7,w]<-(pa.result[1,1])
    table.pa.res[8,w]<-(pa.result[2,1])
    
    ### SAVING RESULTS matrix: END
    
  }
  ### ENDING the loop that is repeated for each sub-sample: END
  
  ### START PRINTING DEAFULT OUTPUT: START
  if(comparison=="random"){
    cat("\n", "Comparison between RANDOM eigenvalues and EMPIRICAL eigenvalues", "\n")
  }
  if(comparison=="bootstrap"){
    cat("\n", "Comparison between BOOTSTRAP eigenvalues and EMPIRICAL eigenvalues", "\n")
  }
  if(comparison=="random-mg"){
    cat("\n", "Comparison between MULTI-GROUP RANDOM eigenvalues and EMPIRICAL eigenvalues", "\n")
  }
  if(comparison=="bootstrap-mg"){
    cat("\n", "Comparison between MULTI-GROUP BOOTSTRAP eigenvalues and EMPIRICAL eigenvalues", "\n")
  }
  
  rownames(table.pa.res)<-c("# of factors (PCA) for Velicer MAP criterium (Pearson corr)...: ", 
                            "# of factors (PCA) for Velicer MAP(4th power)(Polychoric corr): ", 
                            "# of factors (PCA) for Velicer MAP criterium (Polychoric corr): ", 
                            "# of factors (PCA) for Velicer MAP(4th power)(Pearson corr)...: ", 
                            "# of factors (PCA) for PA method (Polychoric Corr.)...........: ", 
                            "# of factors (PCA) for PA method (Pearson Corr.)..............: ", 
                            "# of factors for PA method (Polychoric Corr.).................: ",
                            "# of factors for PA method (Pearson Corr.)....................: ")      
  colnames(table.pa.res)<-sub.sample
  cat("\n******* RESULTS for PARALLEL ANALYSIS: ", "\n")
  print(table.pa.res)
  if(fit.pa==TRUE){  
    colnames(table.fit.res)<-c("sample", "Emp/Sim", "#.factors", "Chi-sqrd", "dof", "p-val", "TLI", "RMSEA", "RMS", "BIC")
    cat("\n******* FIT INDEXES from FACTOR SOLUTIONS INDICATED BY PARALLEL ANALYSIS: ", "\n")
    print(as.data.frame(table.fit.res))
  }
  ### START PRINTING DEAFULT OUTPUT: END
  
  ### PLOTTING RESULTS: START
  if(comparison=="random-mg" | comparison=="bootstrap-mg")  {
    ## for a complete list of graphs check project sheet
    
    # PLOT 1a: FA=T; Polychoric=T; Empirical=T; Multi-Sample=T    
    categ.first<-levels(as.factor(multi.sample))[1]
    bullet.sample <-matrix(,sum(dim(table(multi.sample))),) 
    bullet.sample[1,] <- paste("sample", min(levels(as.factor(multi.sample))), sep = " = ")
    op <- par(no.readonly = TRUE)
    par(mfrow=c(2,2), mar=c(4, 4, 5, 2), oma=c(0,0,0,5), xpd=F) # configure multiple plots per page: in this case 2 rows and 1 column
    
    ### START: finding the maximum value of first eigenvalues extracted from all samples
    first.eigenv.max<-max(print.matrix.multisample.pol[[1]][, "Emp.Polyc.Eigen"])
    max.value<-first.eigenv.max
    for (w in 2:sum(dim(table(multi.sample)))) {  # this loop is repeated for each sub-sample
      next.eigenv.max<-max(print.matrix.multisample.pol[[w]][, "Emp.Polyc.Eigen"])    
      if(next.eigenv.max > max.value)  {
        max.value<-next.eigenv.max
      }
      else{max.value<-max.value}
    }
    ### END: finding the maximum value of first eigenvalues extracted from all samples
    
    plot(1:nvar, print.matrix.multisample.pol[[1]][,"Emp.Polyc.Eigen"], type = "b", 
         xlim = c(1, max(1:nvar)), 
         ylim = c((min(print.matrix.multisample.pol[[1]])), max.value), 
         xlab = "# factors", 
         ylab = "eigenvalues", 
         main = "FA - Empirical Data - Polychoric corr.", cex=1)
    mtext("***** Parallel Analysis - Multisample *****", side=3, col=c("blue"), line=3, cex=1.5)    
    points(1:nvar, print.matrix.multisample.pol[[1]][, "P.SimQuant"], type = "b", pch = 1, col=1+sum(dim(table(multi.sample))))
    #      mtext("Parallel Analysis - Multisample", side=3, line=3, cex=1.5, at=5.5)    
    for (w in 2:sum(dim(table(multi.sample)))) {  # this loop is repeated for each sub-sample
      categ<-levels(as.factor(multi.sample))[w]  # selecting each value of the multisample factor
      bullet.sample[w,] <- paste('sample', categ, sep = " = ")
      points(1:nvar, print.matrix.multisample.pol[[w]][, "Emp.Polyc.Eigen"], type = "b", pch = w, col=w)
      points(1:nvar, print.matrix.multisample.pol[[w]][, "P.SimQuant"], type = "b", pch = w, col=w+sum(dim(table(multi.sample))))
    }
    abline(h = 1)
    abline(h = 0)
    ### PLOTTING RESULTS: END
    
    # PLOT 2a: PCA=T; Polychoric=T; Empirical=T; Multi-Sample=T        
    ### START: finding the maximum value of first eigenvalues extracted from all samples
    first.eigenv.max<-max(print.matrix.pca.multisample.pol[[1]][, "Emp.Pol.Eig.PCA"])
    max.value<-first.eigenv.max
    for (w in 2:sum(dim(table(multi.sample)))) {  # this loop is repeated for each sub-sample
      next.eigenv.max<-max(print.matrix.pca.multisample.pol[[w]][, "Emp.Pol.Eig.PCA"])    
      if(next.eigenv.max > max.value)  {
        max.value<-next.eigenv.max
      }
      else{max.value<-max.value}
    }
    ### END: finding the maximum value of first eigenvalues extracted from all samples
    
    plot(1:nvar, print.matrix.pca.multisample.pol[[1]][,"Emp.Pol.Eig.PCA"], type = "b", 
         xlim = c(1, max(1:nvar)), 
         ylim = c((min(print.matrix.pca.multisample.pol[[1]])), max.value), 
         xlab = "# factors", 
         ylab = "eigenvalues", 
         main = "PCA - Empirical Data - Polychoric corr.", cex=1) 
    points(1:nvar, print.matrix.pca.multisample.pol[[1]][, "P.SimQuant.PCA"], type = "b", pch = 1, col=1+sum(dim(table(multi.sample))))
    for (w in 2:sum(dim(table(multi.sample)))) {  # this loop is repeated for each sub-sample
      categ<-levels(as.factor(multi.sample))[w]  # selecting each value of the multisample factor
      points(1:nvar, print.matrix.pca.multisample.pol[[w]][, "Emp.Pol.Eig.PCA"], type = "b", pch = w, col=w)
      points(1:nvar, print.matrix.pca.multisample.pol[[w]][, "P.SimQuant.PCA"], type = "b", pch = w, col=w+sum(dim(table(multi.sample))))
    }
    abline(h = 1)
    abline(h = 0)
    # PLOT: END
    
    
    # PLOT 1b: FA=T; Pearson=T; Empirical=T; Multi-Sample=T        
    ### START: finding the maximum value of first eigenvalues extracted from all samples
    first.eigenv.max<-max(print.matrix.multisample.pea[[1]][, "Emp.Pears.Eigen"])
    max.value<-first.eigenv.max
    for (w in 2:sum(dim(table(multi.sample)))) {  # this loop is repeated for each sub-sample
      next.eigenv.max<-max(print.matrix.multisample.pea[[w]][, "Emp.Pears.Eigen"])    
      if(next.eigenv.max > max.value)  {
        max.value<-next.eigenv.max
      }
      else{max.value<-max.value}
    }
    ### END: finding the maximum value of first eigenvalues extracted from all samples
    
    plot(1:nvar, print.matrix.multisample.pea[[1]][,"Emp.Pears.Eigen"], type = "b", 
         xlim = c(1, max(1:nvar)), 
         ylim = c((min(print.matrix.multisample.pea[[1]])), max.value), 
         xlab = "# factors", 
         ylab = "eigenvalues", 
         main = "FA - Empirical Data - Pearson corr.", cex=1)
    points(1:nvar, print.matrix.multisample.pea[[1]][, "C.SimQuant"], type = "b", pch = 1, col=1+sum(dim(table(multi.sample))))
    for (w in 2:sum(dim(table(multi.sample)))) {  # this loop is repeated for each sub-sample
      categ<-levels(as.factor(multi.sample))[w]  # selecting each value of the multisample factor
      points(1:nvar, print.matrix.multisample.pea[[w]][, "Emp.Pears.Eigen"], type = "b", pch = w, col=w)
      points(1:nvar, print.matrix.multisample.pea[[w]][, "C.SimQuant"], type = "b", pch = w, col=w+sum(dim(table(multi.sample))))
    }
    abline(h = 1)
    abline(h = 0)
    
    # PLOT 2b: PCA=T; Pearson=T; Empirical=T; Multi-Sample=T        
    first.eigenv.max<-max(print.matrix.pca.multisample.pea[[1]][, "Emp.Pear.Eig.PCA"])
    max.value<-first.eigenv.max
    for (w in 2:sum(dim(table(multi.sample)))) {  # this loop is repeated for each sub-sample
      next.eigenv.max<-max(print.matrix.pca.multisample.pea[[w]][, "Emp.Pear.Eig.PCA"])    
      if(next.eigenv.max > max.value)  {
        max.value<-next.eigenv.max
      }
      else{max.value<-max.value}
    }
    ### END: finding the maximum value of first eigenvalues extracted from all samples
    
    plot(1:nvar, print.matrix.pca.multisample.pea[[1]][,"Emp.Pear.Eig.PCA"], type = "b", 
         xlim = c(1, max(1:nvar)), 
         ylim = c((min(print.matrix.pca.multisample.pea[[1]])), max.value), 
         xlab = "# factors", 
         ylab = "eigenvalues", 
         main = "PCA - Empirical Data - Pearson corr.", cex=1) 
    points(1:nvar, print.matrix.pca.multisample.pea[[1]][, "C.SimQuant.PCA"], type = "b", pch = 1, col=1+sum(dim(table(multi.sample))))
    for (w in 2:sum(dim(table(multi.sample)))) {  # this loop is repeated for each sub-sample
      categ<-levels(as.factor(multi.sample))[w]  # selecting each value of the multisample factor
      points(1:nvar, print.matrix.pca.multisample.pea[[w]][, "Emp.Pear.Eig.PCA"], type = "b", pch = w, col=w)
      points(1:nvar, print.matrix.pca.multisample.pea[[w]][, "C.SimQuant.PCA"], type = "b", pch = w, col=w+sum(dim(table(multi.sample))))
    }
    abline(h = 1)
    abline(h = 0)
    # PLOT: END
    
    par(xpd=NA)
    legend.bullet.sample<-matrix(,sum(dim(table(multi.sample)))*2,1)
    for (w in 1:sum(dim(table(multi.sample)))) {  # this loop is repeated for each sub-sample
      categ<-levels(as.factor(multi.sample))[w]  # selecting each value of the multisample factor
      legend.bullet.sample[w,1]<-paste(categ,"emp", sep="=")
      legend.bullet.sample[w+sum(dim(table(multi.sample))),1]<-paste(categ,"sim", sep="=")      
    }
    legend("topright", inset=c(-0.5,-0.8), legend=legend.bullet.sample, bty="n", x.intersp=0.2, y.intersp=0.2, col = 1:length(legend.bullet.sample), pch = 1:(sum(dim(table(multi.sample)))))
    
    par(op)
    par(mfrow=c(1,1)) # re-configure just one plot per page      
  }
  else {
    #### PLOT FOR SINGLE SAMPLE: START
    ### PLOT RANDOM VS EMPIRICAL
    if(comparison=="random"){
      if(diff.fact==TRUE) {
        text.title<-c("Parallel Analysis: diff.fact")
      }
      else {
        text.title<-c("Parallel Analysis")
      }
      plot(st.matrix[, 5], st.matrix[, 10], type = "b", xlim = c(1, max(st.matrix[, 5])), ylim = c((min(st.matrix)), (max(st.matrix[, c(4, 9, 10, 11)]))), xlab = "# factors", ylab = "eigenvalues", main = text.title)
      points(st.matrix[, 5], st.matrix[, 4], type = "b", pch = 8)
      points(st.matrix[, 5], st.matrix[, 11], type = "b", pch = 20, col = "red")
      points(st.matrix[, 5], st.matrix[, 9], type = "b", pch = 2, col = "red")
      perc <- paste(q.eigen * 100, "* perc. Polychoric corr. Sim. FA", sep = "")
      perc1 <- paste(q.eigen * 100, "* perc. Pearson corr. Sim. FA", sep = "")
      res <- (pa.result[1,])
      res1 <- (pa.result[2,])
      res.pca <- (pa.result[3,])
      res1.pca <- (pa.result[4,])
      ris <- paste("# factors with Polyc.PA: ", res, sep = "")
      ris1 <- paste("# factors with Pear.PA: ", res1, sep = "")
      
      legend(x = "topright", c("Polychoric corr. Empirical FA", "Pearson corr. Empirical FA", perc, perc1, ris, ris1), 
             col = c(1, 2, 1, 2, 1, 2), pch = c(1, 20, 8, 2), y.intersp=0.4, cex=0.95)
      
      abline(h = 1)
      abline(h = 0)
    }
    else{
      ### PLOT BOOSTRAP VS EMPIRICAL
      ### Scree-plot of FA Empirical vs Bootstrap for polychoric and pearson corr: START
      par(pch=20)
      plot(1:ncol(data.matrix.sub), raw.eigen.poly.fa, type="b", xlim=c(1, max(ncol(data.matrix.sub))), ylim=c(-1, 5), xlab="# factors", ylab="eigenvalues", main="FA: Empirical vs Bootstrap")
      points(1:ncol(data.matrix.sub), raw.eigen.pear.fa, type="b", pch=8)
      points(1:ncol(data.matrix.sub), perm.eigen.polyc.fa, type="b", pch=20, col="red")
      points(1:ncol(data.matrix.sub), perm.eigen.pear.fa, type="b", pch=8, col="red")
      legend.boot.poly.fa<-paste("Bootstrap Polychoric corr. FA (", q.eigen*100, "\u00B0)")
      legend.boot.pear.fa<-paste("Bootstrap Pearson corr. FA (", q.eigen*100, "\u00B0)")
      legend(x = "topright", c("Empirical Polychoric corr. FA", "Empirical Pearson corr. FA",  legend.boot.poly.fa, legend.boot.pear.fa), 
             col = c("black", "black", "red", "red"), pch = c(20, 8, 20, 8),
             cex = 0.8, y.intersp=0.4)
      abline(h = 1)
      abline(h = 0)
      ### plotta lo scree-test: END
      
      ### Scree-plot of PCA Empirical vs Bootstrap for polychoric and pearson corr: START
      par(pch=20)
      plot(1:ncol(data.matrix.sub), raw.eigen.poly.pca, type="b", xlim=c(1, max(ncol(data.matrix.sub))), ylim=c(-1, 5), xlab="# factors", ylab="eigenvalues", main="PCA: Empirical vs Bootstrap")
      points(1:ncol(data.matrix.sub), raw.eigen.pear.pca, type="b", pch=8)
      points(1:ncol(data.matrix.sub), perm.eigen.polyc.pca, type="b", pch=20, col="red")
      points(1:ncol(data.matrix.sub), perm.eigen.pear.pca, type="b", pch=8, col="red")
      legend.boot.poly.pca<-paste("Bootstrap Polychoric corr. PCA (", q.eigen*100, "\u00B0)")
      legend.boot.pear.pca<-paste("Bootstrap Pearson corr. PCA (", q.eigen*100, "\u00B0)")
      legend(x = "topright", c("Empirical Polychoric corr. PCA", "Empirical Pearson corr. PCA",  legend.boot.poly.pca, legend.boot.pear.pca), 
             col = c("black", "black", "red", "red"), pch = c(20, 8, 20, 8),
             cex=0.8, y.intersp=0.4)
      abline(h = 1)
      abline(h = 0)
      ### plotta lo scree-test: END
    }
  }
  ### PLOTTING RESULTS: END
  
  
  ### PRINTING RESULTS matrix: START
  if(print.all==TRUE) {
    if(diff.fact==TRUE & (comparison=="random" | comparison=="random-mg")){
      cat("\n******* DIFFERENTIAL WEIGHTS used for simulating random datasets: ")      
      cat("\nThe diff.fact paramether is set to TRUE, so random dataset will be\n")
      cat("simulated by taking into account the weights of each category for each variable\n")
      cat("within the provided dataset. Weights are listed in the following table:\n\n")
      print(print.diff.weight)
      flush.console()
      cat("\n\n")
    }
    else{}
    if(comparison=="bootstrap" | comparison=="bootstrap-mg") {
      cat("\n******* PA RESULTS from polychoric FA: ", "\n")
      print(print.matrix.multisample.pol) # print results from FA-multisample
      #        print(as.data.frame(print.perm.eigen.polyc.fa),justify = c("left"))
      print(noquote(print.perm.eigen.polyc.fa))
      cat("\n******* PA RESULTS from polychoric PCA: ", "\n")
      print(print.matrix.pca.multisample.pol) # print results from PCA-multisample
      print(noquote(print.perm.eigen.polyc.pca))
      cat("\n******* PA RESULTS from pearson FA: ", "\n")
      print(print.matrix.multisample.pea) # print results from FA-multisample
      print(noquote(print.perm.eigen.pear.fa),justify = c("left"))
      cat("\n******* PA RESULTS from pearson PCA: ", "\n")
      print(print.matrix.pca.multisample.pea) # print results from PCA-multisample
      print(noquote(print.perm.eigen.pear.pca),justify = c("left"))
      cat("\n******* MAP RESULTS from PCA: ", "\n")
      print(map.result.multisample)  # print results from MAP-multisample
    }
    else{
      cat("\n******* PA RESULTS from polychoric FA: ", "\n")
      print(print.matrix.multisample.pol) # print results from FA-multisample
      cat("\n******* PA RESULTS from polychoric PCA: ", "\n")
      print(print.matrix.pca.multisample.pol) # print results from PCA-multisample
      cat("\n******* PA RESULTS from pearson FA: ", "\n")
      print(print.matrix.multisample.pea) # print results from FA-multisample
      cat("\n******* PA RESULTS from pearson PCA: ", "\n")
      print(print.matrix.pca.multisample.pea) # print results from PCA-multisample
      cat("\n******* MAP RESULTS from PCA: ", "\n")
      print(map.result.multisample)  # print results from MAP-multisample
    }
  }
  ### PRINTING RESULTS matrix: END
}