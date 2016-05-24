stochasticProfilingML <-
function() {

   # definition of variables (necessary for CMD R check)
   # (these variables will be initialized later, but they are not visible as global functions/data)
   get.range <- NULL
   rm(get.range)


   cat("This function performs maximum likelihood estimation for the stochastic\nprofiling model. In the following, you are asked to enter your data and\nspecify some settings. By pressing 'enter', you choose the default option.\n\n")

   model.default <- 1
   TY.default <- 2
   n.default <- 10
   m.default <- 1

   options(warn=-2)

   # enter data
   cat("---------\n")
   continue <- F
   this.text <- "How would you like to input your data?\n 1: enter manually\n 2: read from file\n 3: enter the name of a variable\n(default: 1).\n"
   while (!continue) {
      manually <- readline(this.text)
      if (manually=="") {
         manually <- 1
      }
      if (manually %in% 1:3) {
         continue <- T
      }
      else {
         this.text <- "Invalid choice. Please try again.\n"
      }
   }

   # enter data manually
   if (manually==1) {
      # choose m
      cat("---------\n")
      continue <- F
      this.text <- paste("Please enter the number of genes your dataset contains measurements for\n(default: 1).\n",sep="")
      while (!continue) {
         m <- readline(this.text)
         if (m=="") { m <- m.default }
         if (is.na(as.numeric(m))) {
            this.text <- "Invalid choice. Please enter a finite natural number.\n"
         }
         else {
            m <- as.numeric(m)
            if ((round(m)==m) && ((m>0) && (m<Inf))) {
               continue <- T
            }
            else {
               this.text <- "Invalid choice. Please enter a finite natural number.\n"
            }
         }
      }

      # enter data
      cat("---------\n")
      for (g in 1:m) {
         continue <- F
         add.text <- ""
         if (m>1) {
            add.text <- paste("\nFor gene number ",g,":\n",sep="")
         }
         if (g==1) {
            this.text <- paste("Please enter your measurements, separated by either commas or spaces,\ne.g. in case of five samples\n100.3, 150.7, 110.0, 80.2, 201.9\n or\n100.3 150.7 110.0 80.2 201.9.\n",add.text,sep="")
         }
         else {
            this.text <- add.text
         }
         while (!continue) {
            this.data <- readline(this.text)

            if (this.data=="") {
               this.text <- "Please enter some numbers as there is no default value.\n"
            }
            else {
               # try comma separation
               this.data.tmp <- as.numeric(unlist(strsplit(this.data, ",")))
               # try space separation
               if (is.na(this.data.tmp)) {
                  this.data <- as.numeric(unlist(strsplit(this.data, " ")))
               }
               else {
                  this.data <- this.data.tmp
               }
               if (any(is.na(this.data))) {
                  this.text <- "This is not a numeric vector. Please try again.\n"
               }
               else {
                  if (any(abs(this.data)==Inf)) {
                     this.text <- "There are infinite values. Please choose finite ones.\n"
                  }
                  else if (any(this.data<=0)) {
                     this.text <- "There are non-positive values. Please choose positive ones.\n"
                  }
                  else if (any(is.na(this.data))) {
                     this.text <- "There are missing values. Please choose real ones.\n"
                  }
                  else {
                     if (g==1) {
                        dataset <- matrix(NA,ncol=m,nrow=length(this.data))
                        dataset[,1] <- this.data
                        continue <- T
                     }
                     else {
                        if (length(this.data)!=nrow(dataset)) {
                           this.text <- "The number of samples for this gene does not agree with the number of\nsamples for the previous genes. Please try again.\n"
                        }
                        else {
                           dataset[,g] <- this.data
                           continue <- T
                        }
                     }
                  }
               }
            }
         }
      }
      cnames <- "no"
      rnames <- "no"
      dims <- 1
   }
   else if (manually==2) {
      # choose filename
      cat("---------\n")
      continue <- F
      cat("The file should contain a data matrix with one dimension standing for genes\nand the other one for samples. Fields have to be separated by tabs or white\nspaces, but not by commas. If necessary, please delete the commas in the\ntext file using the \'replace all\' function of your text editor.\n")
      cat("\nPlease enter a valid path and filename, either a full path, e.g.\n",getwd(),"/mydata.txt\nor just a file name, e.g.\nmydata.txt.\nThe current directory is\n",getwd(),".\n",sep="")
      this.text <- ""

      while (!continue) {
         path <- readline(this.text)

         if (!file.exists(path)) {
            this.text <- "This file does not exist. Please try again.\n"
         }
         else {
            continue <- T
         }
      }

      continue <- F
      this.text <- "Does the file contain column names? Please enter 'yes' or 'no'.\n"
      while (!continue) {
         cnames <- readline(this.text)
         if ((cnames=="") || (!(cnames %in% c("yes","no")))) {
            this.text <- "Please enter 'yes' or 'no'.\n"
         }
         else {
            continue <- T
         }
      }

      continue <- F
      this.text <- "Does the file contain row names? Please enter 'yes' or 'no'.\n"
      while (!continue) {
         rnames <- readline(this.text)
         if ((rnames=="") || (!(rnames %in% c("yes","no")))) {
            this.text <- "Please enter 'yes' or 'no'.\n"
         }
         else {
            continue <- T
         }
      }

      continue <- F
      this.text <- "Do the columns stand for different genes or different samples?\n 1: genes\n 2: samples.\n"
      while (!continue) {
         dims <- readline(this.text)
         if ((dims=="") || (!(dims %in% 1:2))) {
            this.text <- "Invalid choice. Please try again.\n"
         }
         else {
            continue <- T
         }
      }

      # read file

      if (rnames=="yes") {
         dataset <- read.table(file=path,header=(cnames=="yes"),row.names=1)
      }
      else {
         dataset <- read.table(file=path,header=(cnames=="yes"))
      }
      if (dims==2) {
         dataset <- t(dataset)
      }
   }
   if (manually==3) {
      cat("---------\n")
      continue <- F
      this.text <- "The variable should be a matrix with one dimension standing for genes\nand the other one for samples. Please enter the name of the variable.\n"
      while (!continue) {
         variablename <- readline(this.text)
         variable <- try(eval(parse(text=variablename)),silent=T)

         if (class(variable)=="try-error") {
            this.text <- "This variable does not exist. Please try again.\n"
         }
         else {
            continue <- T
            if (is.data.frame(variable)) {
               cat("Data frame has been converted to matrix.\n")
               variable <- as.matrix(variable)
            }
            if (!is.matrix(variable)) {
               cat(paste("This is not a matrix. The variable is converted to a matrix with 1 column\nand ",length(variable)," rows.\n\n",sep=""))
               variable <- matrix(variable,ncol=1)
            }
         }
      }

      continue <- F
      this.text <- "Do the columns stand for different genes or different samples?\n 1: genes\n 2: samples.\n"
      while (!continue) {
         dims <- readline(this.text)
         if ((dims=="") || (!(dims %in% 1:2))) {
            this.text <- "Invalid choice. Please try again.\n"
         }
         else {
            continue <- T
            if (dims==1) {
               dataset <- variable
            }
            else {
               dataset<- t(variable)
            }
         }
      }
      if (is.null(colnames(dataset))) {
         cnames <- "no"
      }
      else {
         cnames <- "yes"
      }
      rnames <- "no"
   }

   dataset <- as.matrix(dataset)
   m <- ncol(dataset)

   cat("\nThis is the head of the dataset (columns contain different genes):\n")
   print(head(dataset))

   # check whether  data matrix is correct
   if (any(is.na(as.numeric(dataset)))) {
      stop("\nError: The dataset contains non-numeric elements.")
   }
   if (any(is.na(dataset))) {
      stop("\nError: The dataset contains missing values.")
   }
   if (any(is.null(dataset))) {
      stop("\nError: The dataset contains missing values.")
   }
   if (any(dataset<=0)) {
      stop("\nError: The dataset contains non-positive elements.")
   }

   cat("\nIf the matrix does not look correct to you, there must have been an error\nin the answers above. In this case, please quit by pressing 'escape' and call\nstochasticProfilingML() again.\n\n")

   # names of genes
   continue <- F
   if (((cnames=="yes") && (dims==1)) || ((rnames=="yes") && (dims==2))) {
      genenames <- colnames(dataset)
      cat("The file contained the following gene names:\n")
      cat(genenames)
      this.text <- "\nWould you like to enter different names? Please enter 'yes' or 'no'.\n"
   }
   else {
      genenames <- NULL
      this.text <- "Would you like to enter names for the genes? Otherwise, genes will\nsimply be numbered. Please enter 'yes' or 'no'.\n"
   }
   continue <- F
   while (!continue) {
      enter.genenames <- readline(this.text)
      if ((enter.genenames=="") || (!(enter.genenames %in% c("yes","no")))) {
         this.text <- "Please enter 'yes' or 'no'.\n"
      }
      else {
         continue <- T
      }
   }

   if (enter.genenames=="yes") {
      continue <- F
      this.text <- "Please enter the names of the genes, e.g.\n POT1, CLDN7, GCS1.\nThe names have to be separated by commas. Spaces within names are allowed.\n"
      while (!continue) {
         genenames <- readline(this.text)

         if (genenames=="") {
            this.text <- "Please enter the names of the genes.\n"
         }
         else {
            # comma separation
            genenames <- unlist(strsplit(genenames, ","))
            if (length(genenames)!=m) {
               this.text <- "The number of names does not coincide with the number of genes.\nPlease try again.\n"
            }
            else {
               continue <- T
            }
         }
      }
   }
   else if (is.null(genenames)) {
      genenames <- paste("Gene",1:m)
   }
   colnames(dataset) <- genenames

   # choose model
   cat("---------\n")
   continue <- F
   this.text <- paste("Please choose the model you would like to estimate:\n 1: LN-LN\n 2: rLN-LN\n 3: EXP-LN\n(default: ",model.default,")\n",sep="")
   while (!continue) {
      model <- readline(this.text)
      if (model=="") { model <- model.default }
      if (model=="LN-LN") { model <- 1 }
      else if (model=="rLN-LN") { model <- 2 }
      else if (model=="EXP-LN") { model <- 3 }
      if (model %in% c(1,2,3)) {
         continue <- T
      }
      else {
         this.text <- "Invalid choice. Please choose 1, 2 or 3 oder simply press 'enter' for the\ndefault option.\n"
      }
   }

   # choose TY
   cat("---------\n")
   continue <- F
   this.text <- paste("Please enter the number of different populations you would like to estimate:\n(default: ",TY.default,")\n",sep="")
   while (!continue) {
      TY <- readline(this.text)
      if (TY=="") { TY <- TY.default }
      if (is.na(as.numeric(TY))) {
         this.text <- "Invalid choice. Please enter a finite natural number.\n"
      }
      else {
         TY <- as.numeric(TY)
         if (TY %in% 1:50) {
            continue <- T
         }
         else if (TY>50) {
            this.text <- "More than 50 populations are theoretically possible but not meaningful.\nPlease choose a smaller number.\n"
         }
         else {
            this.text <- "Invalid choice. Please enter a natural number.\n"
         }
      }
   }

   # choose n
   cat("---------\n")
   continue <- F
   this.text <- paste("Please enter the number of cells that entered each sample:\n(default: ",n.default,")\n",sep="")
   while (!continue) {
      n <- readline(this.text)
      if (n=="") { n <- n.default }
      if (is.na(as.numeric(n))) {
         this.text <- "Invalid choice. Please enter a finite natural number.\n"
      }
      else {
         n <- as.numeric(n)
         if ((round(n)==n) && ((n>0) && (n<Inf))) {
            continue <- T
         }
         else {
            this.text <- "Invalid choice. Please enter a finite natural number.\n"
         }
      }
   }


   ############
   # settings #
   ############

   genes <- genenames

   # length of vector mu
   if (model==1) {
      mu.length <- TY*m
      model <- "LN-LN"

   }
   else if (model==2) {
      mu.length <- TY*m
      model <- "rLN-LN"
   }
   else if (model==3) {
      model <- "EXP-LN"
      if (TY>1) {
         mu.length <- (TY-1)*m
      }
      else {
         mu.length <- 0
      }
   }
   set.model.functions(model)
   use.constraints <- F
   show.plots <- T


   # ask for preanalysis and 'main analysis'
   if ((m==1) || (mu.length==0)) {
      preanalyze <- F
      mainanalyze <- F
   }
   else {
      # ask for preanalysis
      cat("---------\n")
      continue <- F
      this.text <- "Would you like to do a preanalysis in which estimation is initially\ncarried out for each gene individually? This may give a first rough\nidea about the single genes and speed up the estimation for the entire\ncluster. Please enter 'yes' or 'no'.\n"
      while (!continue) {
         preanalyze.answer <- readline(this.text)
         if ((preanalyze.answer=="") || (!(preanalyze.answer %in% c("yes","no")))) {
            this.text <- "Please enter 'yes' or 'no'.\n"
         }
         else {
            preanalyze <- (preanalyze.answer=="yes")
            continue <- T
         }
      }

      if (m>4) {
         # ask for analysis of subclusters
         continue <- F
         cat("\nThe dataset contains more than four genes. It is hence recommended to\ncarry out estimation for subgroups of size four first. This will speed up\nthe final estimation and yield more reliable results. Would you like\nto consider subgroups first? Please enter 'yes' or 'no'.\n")
         this.text <- ""
         while (!continue) {
            mainanalyze.answer <- readline(this.text)
            if ((mainanalyze.answer=="") || (!(mainanalyze.answer %in% c("yes","no")))) {
               this.text <- "Please enter 'yes' or 'no'.\n"
            }
            else {
               mainanalyze <- (mainanalyze.answer=="yes")
               continue <- T
            }
         }
      }
      else {
         mainanalyze <- F
      }
   }

   options(warn=0)

   ##############################################
   # Prompting finished!! Start estimation now. #
   ##############################################

   cat("\n***** Estimation started! *****\n\n")

   #################
   # optional step #
   #################
   # Consider all genes separately first in order to get a rough idea about
   # the location of the parameters. This might speed up the analysis of the larger clusters.

   if ((preanalyze) && (mu.length>0)) {
      cat("###################################################\n")
      cat("## Preanalysis: Consider all genes individually. ##\n")
      cat("###################################################\n")
      # number of parameters
      if (model=="LN-LN") {
         npar <- TY*2
         mu.indices <- TY:(npar-1)
         mu.names <- paste("mu_",1:TY,sep="")
      }
      else if (model=="rLN-LN") {
         npar <- TY*3-1
         mu.indices <- TY:(2*TY-1)
         mu.names <- paste("mu_",1:TY,sep="")
      }
      else if (model=="EXP-LN") {
         if (TY>1) {
            npar <- TY*2
            mu.indices <- TY:(2*TY-2)
            mu.names <- paste("mu_",1:(TY-1),sep="")
         }
         else {
            npar <- 1
            mu.indices <- NULL
            mu.names <- NULL
         }
      }

      # store the results for the single genes in this matrix:
      # lower and upper bounds of confidence intervals
      single.gene.lower <- matrix(nrow=npar,ncol=length(genes))
      single.gene.upper <- matrix(nrow=npar,ncol=length(genes))

      # single gene estimation
      for (i in 1:length(genes)) {
         this.gene <- genes[i]

         this.text <- paste("|",this.gene,"|")
         this.length <- nchar(this.text)
         this.line <- NULL
         for (l in 1:this.length) {
            this.line <- paste(this.line,"-",sep="")
         }
         cat(paste(this.line,"\n",sep=""))
         cat(paste(this.text,"\n",sep=""))
         cat(paste(this.line,"\n",sep=""))

         single.res <- stochprof.loop(model=model,dataset=dataset[,i,drop=F],n=n,TY=TY,genenames=this.gene,fix.mu=F,loops=3,until.convergence=T,print.output=F,show.plots=show.plots,plot.title=paste("Gene",this.gene),use.constraints=use.constraints)

         if (!(is.null(single.res$ci))) {
            single.gene.lower[,i] <- (single.res$ci)[,1]
            single.gene.upper[,i] <- (single.res$ci)[,2]
         }
         else {
            other.interval <- get.range(method="quant",prev.result=single.res$pargrid,TY=TY,fix.mu=F)
            single.gene.lower[,i] <- other.interval[,1]
            single.gene.upper[,i] <- other.interval[,2]
         }
      }

      colnames(single.gene.lower) <- paste("gene",genes)
      cn <- colnames(single.res$pargrid)
      rownames(single.gene.lower) <- cn[-length(cn)]
      if (length(mu.indices)>0) {
         rownames(single.gene.lower)[mu.indices] <- mu.names
      }
      if (model=="EXP-LN") {
         rownames(single.gene.lower)[nrow(single.gene.lower)] <- "lambda"
      }
      dimnames(single.gene.upper) <- dimnames(single.gene.lower)

      cat("--------------------------\n")
      cat("| Result of preanalysis: |\n")
      cat("--------------------------\n")
      cat("lower bounds:\n")
      print(single.gene.lower)
      cat("\n")
      cat("upper bounds:\n")
      print(single.gene.upper)
   }

   #############
   # main step #
   #############

   if ((mu.length>0) && (mainanalyze)) {
      cat("#################################################\n")
      cat("## Main analysis: Consider subgroups of genes. ##\n")
      cat("#################################################\n")

      set.set <- list()
      subgroup.size <- 4
      no.of.subgroups <- ceiling(m/subgroup.size)
      if (no.of.subgroups==1) {
         set.set <- list(genes)
      }
      else {
         for (i in 1:no.of.subgroups) {
            if (i < no.of.subgroups) {
               set.set[[i]] <- genes[(i-1)*subgroup.size + 1:subgroup.size]
            }
            else {
               set.set[[i]] <- genes[m-(3:0)]
            }
         }
      }

      # according sets of indices with respect to vector "genes"
      index.list <- list(length=length(set.set))
      for (j in 1:length(set.set)) {
         this.set <- set.set[[j]]
         k.vector <- NULL
         for (i in 1:length(this.set)) {
            this.gene <- this.set[i]
            k <- which(genes==this.gene)
            k.vector <- c(k.vector,k)
         }
         index.list[[j]] <- k.vector
      }

      # for each of the subgroups, estimate all parameters
      result.set <- list(length=length(set.set))

      for (set.index in 1:length(set.set)) {
         # consider the subsets of genes
         this.set <- set.set[[set.index]]
         # size of this subset
         this.m <- length(this.set)


         # strings for printing the results
         label <- "| "
         for (v in 1:length(this.set)) {
            label <- paste(label,this.set[v])
         }
         label <- paste(label,"|")
         label2 <- NULL
         for (v in 1:nchar(label)) {
            label2 <- paste(label2,"-",sep="")
         }
         cat(label2,"\n")
         cat(label,"\n")
         cat(label2,"\n")

         if (preanalyze) {
            # range from which parameters should be drawn
            if (model=="LN-LN") {
               this.npar <- TY*(this.m+1)
               sigma.indices <- this.npar
               single.sigma.indices <- nrow(single.gene.lower)
               single.mu.indices <- TY:(nrow(single.gene.lower)-1)
            }
            else if (model=="rLN-LN") {
               this.npar <- TY*(this.m+2)-1
               sigma.indices <- ((this.m+1)*TY):((this.m+2)*TY-1)
               single.sigma.indices <- (2*TY):(3*TY-1)
               single.mu.indices <- TY:(2*TY-1)
            }
            else if (model=="EXP-LN") {
               if (TY>1) {
                  this.npar <- TY*(this.m+1)
                  sigma.indices <- (this.m+1)*(TY-1)+1
                  single.sigma.indices <- 2*(TY-1)+1
                  single.mu.indices <- TY:(2*(TY-1))
               }
               else {
                  this.npar <- this.m
                  sigma.indices <- NULL
                  single.sigma.indices <- NULL
                  single.mu.indices <- NULL
               }
            }
            par.range <- matrix(nrow=this.npar,ncol=2)
            # initialize bounds for p
            if (TY>1) {
               par.range[1:(TY-1),1] <- 1
               par.range[1:(TY-1),2] <- 0
            }
            # initialize bounds for sigma
            if (length(sigma.indices)>0) {
               par.range[sigma.indices,1] <- 10^7
               par.range[sigma.indices,2] <- 0
            }
            if (model=="EXP-LN") {
               # initialize bounds for lambda
               if (TY>1) {
                  par.range[(this.m+1)*(TY-1)+1+(1:this.m),1] <- 10^7
                  par.range[(this.m+1)*(TY-1)+1+(1:this.m),2] <- 0
               }
               else {
                  par.range[,1] <- 10^7
                  par.range[,2] <- 0
               }
            }
            # derive values from preanalysis
            for (i in 1:length(this.set)) {
               # consider each gene in the subset
               this.gene <- this.set[i]
               # position of this gene in the overall dataset
               k <- which(genes==this.gene)
               # p: union of all confidence intervals for p in this subset
               if (TY>1) {
                  par.range[1:(TY-1),1] <- pmin(par.range[1:(TY-1),1],single.gene.lower[1:(TY-1),k])
                  par.range[1:(TY-1),2] <- pmax(par.range[1:(TY-1),2],single.gene.upper[1:(TY-1),k])
               }
               # sigma: union of all confidence intervals for sigma in this subset
               if (length(sigma.indices)>0) {
                  par.range[sigma.indices,1] <- min(par.range[sigma.indices,1],single.gene.lower[single.sigma.indices,k])
                  par.range[sigma.indices,2] <- max(par.range[sigma.indices,2],single.gene.upper[single.sigma.indices,k])
               }
               # mu: equal to confidence interval for mu for the respective gene
               if (model %in% c("LN-LN","rLN-LN")) {
                  mu.indices <- TY-1+(0:(TY-1))*this.m+i
               }
               else if (model=="EXP-LN") {
                  if (TY>1) {
                     mu.indices <- TY-1+(0:(TY-2))*this.m+i
                  }
                  else {
                     mu.indices <- NULL
                  }
               }
               if (length(mu.indices)>0) {
                  par.range[mu.indices,1] <- single.gene.lower[single.mu.indices,k]
                  par.range[mu.indices,2] <- single.gene.upper[single.mu.indices,k]
               }

               if (model=="EXP-LN") {
                  # lambda
                  if (TY>1) {
                     par.range[(this.m+1)*(TY-1)+1+i,1] <- pmin(par.range[(this.m+1)*(TY-1)+1+i,1],single.gene.lower[nrow(single.gene.lower),k])
                     par.range[(this.m+1)*(TY-1)+1+i,2] <- pmax(par.range[(this.m+1)*(TY-1)+1+i,2],single.gene.upper[nrow(single.gene.lower),k])
                  }
                  else {
                     par.range[i,1] <- pmin(par.range[i,1],single.gene.lower[nrow(single.gene.lower),k])
                     par.range[i,2] <- pmax(par.range[i,2],single.gene.upper[nrow(single.gene.lower),k])
                  }
               }
            }
         }
         else {
            par.range <- NULL
         }

         # set of indices
         indices <- index.list[[set.index]]

         # estimation for entire subgroup
         this.result <- NULL
         while (is.null(this.result)) {
            this.result <- stochprof.loop(model=model,dataset=dataset[,indices,drop=F],n=n,TY=TY,par.range=par.range,genenames=this.set,fix.mu=F,loops=10,until.convergence=T,print.output=F,show.plots=show.plots,plot.title=paste("Subgroups of genes: group number",set.index),use.constraints=use.constraints)
         }
         result.set[[set.index]] <- this.result
      }
   }


   ##############
   # final step #
   ##############

   # Fix mu to the values determined in the subclusters.
   # There remain the following parameters to estimate:
   # p and sigma in the LN-LN and rLN-LN model, and
   # p, sigma and lambda in the EXP-LN model.
   if ((mu.length>0) && (mainanalyze)) {
      cat("################################################################\n")
      cat("## Final analysis: Fix mu, estimate p and sigma (and lambda). ##\n")
      cat("################################################################\n")

      fixed.mu <- vector(length=mu.length)
      for (i in 1:length(result.set)) {
         # result of this subgroup
         this.result <- result.set[[i]]
         this.mle <- this.result$mle
         # indices of genes in this subgroup
         indices <- index.list[[i]]
         # according indices in fixed.mu of entire cluster
         if (model=="LN-LN") {
            this.mu.indices <- TY:(length(this.mle)-1)
         }
         else if (model=="rLN-LN") {
            this.mu.indices <- TY:(length(this.mle)-TY)
         }
         else if (model=="EXP-LN") {
            if (TY>1) {
               this.mu.indices <- TY:(length(this.mle)-1-this.m)
            }
            else {
               this.mu.indices <- NULL
            }
         }
         fm.indices <- NULL
         for (j in 1:(mu.length/m)) {
            fm.indices <- c(fm.indices,(j-1)*m+indices)
         }
         fixed.mu[fm.indices] <- this.mle[this.mu.indices]
      }
      # estimation of p, sigma and lambda for mu fixed to the values estimated in the subgroups
      result <- stochprof.loop(model=model,dataset=dataset[,genes,drop=F],n=n,TY=TY,genenames=genes,fix.mu=T,fixed.mu=fixed.mu,loops=10,until.convergence=T,print.output=F,show.plots=show.plots,plot.title="Entire Cluster",use.constraints=use.constraints,subgroups=set.set)
   }
   else {
      if (preanalyze) {
         cat("###################################################\n")
         cat("## Final analysis: Estimate model for all genes. ##\n")
         cat("###################################################\n")
      }

      # estimation of p and lambda (there is no mu)
      result <- stochprof.loop(model=model,dataset=dataset[,genes,drop=F],n=n,TY=TY,genenames=genes,fix.mu=F,loops=10,until.convergence=F,print.output=F,show.plots=F,plot.title="Entire Cluster",use.constraints=use.constraints)
   }

   return(invisible(result))
}
