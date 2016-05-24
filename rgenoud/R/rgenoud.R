#
#  RGENOUD
#
#  Walter R. Mebane, Jr.
#  University of Michigan
#  http://www-personal.umich.edu/~wmebane
#  <wmebane@umich.edu>
#
#  Jasjeet Singh Sekhon 
#  UC Berkeley
#  http://sekhon.berkeley.edu
#  <sekhon@berkeley.edu>
#


genoud <- function(fn, nvars, max=FALSE, pop.size=1000, max.generations=100, wait.generations=10,
                   hard.generation.limit=TRUE, starting.values=NULL, MemoryMatrix=TRUE, 
                   Domains=NULL, default.domains=10, solution.tolerance=0.001,
                   gr=NULL, boundary.enforcement=0, lexical=FALSE, gradient.check=TRUE, BFGS=TRUE, 
                   data.type.int=FALSE, hessian=FALSE, unif.seed=812821, int.seed=53058,
                   print.level=2, share.type=0, instance.number=0,
                   output.path="stdout", output.append=FALSE, project.path=NULL, 
                   P1=50, P2=50, P3=50, P4=50, P5=50, P6=50, P7=50, P8=50, P9=0,
                   P9mix=NULL, BFGSburnin=0, BFGSfn=NULL, BFGShelp = NULL,
                   control = list(), optim.method=ifelse(boundary.enforcement < 2, "BFGS", "L-BFGS-B"),
                   transform=FALSE, debug=FALSE, cluster=FALSE, balance=FALSE,  ...)
{

  if (transform!=FALSE & transform!=TRUE)
    {
      warning("'transform' must be TRUE or FALSE.  Setting to FALSE")
      transform <- FALSE
    }
  if(transform) {
    return(do.call(what  = genoud_transform,
                   args  = as.list(match.call())[-1],
                   quote = FALSE,
                   envir = parent.frame()))
  }

  if(!is.null(BFGShelp) && !is.function(BFGShelp)) stop("'BFGShelp' must be NULL or a function")

  #These used to be option. From >5.7-3, they are no longer available
  #because R CMD check complains about calls to 'stdout'
  if (output.path!="stdout")
    {
      output.path <- "stdout"              
      warning("'output.path' can no longer be changed. Please use 'sink'. Option is only provided for backward compatibility of the API.")
    }
  if(output.append!=FALSE)
    {
      output.append <- FALSE
      warning("'output.append' can no longer be changed. Please use 'sink'. Option is only provided for backward compatibility of the API.")      
    }
  
  if(!is.null(P9mix) && !is.double(P9mix))  {
    stop("'P9mix' must be NULL or a number between 0 and 1")
  } else {
    if(is.null(P9mix)) {
      P9mix <- -1
    } else {
      if(! ( (1 >= P9mix) && (P9mix > 0) ))
        stop("'P9mix' must be NULL or a number between 0 and 1 (it may be equal to 1)")
    }
  }

  if( (BFGSburnin < 0) & !gradient.check )
    {
      warning("If 'BFGSburnin' is negative, gradient.check must be TRUE for 'BFGSburnin' to have any effect.")
    }

  if(!is.list(control))
    stop("'control' must be a list, see ?optim")
    g.scale <- control$fnscale
    if(!is.null(g.scale)) {
      if(g.scale > 0 & max) {
        stop("positive control$fnscale is inconsistent with maximization")
      }
      else if(g.scale < 0 & !max) {
        stop("negative control$fnscale is inconsistent with minimization")
      }
      else if(g.scale == 0) {
        stop("optim divides the function value by control$fnscale ",
             "setting control$fnscale to zero is therefore impossible")
      }
      FiniteBadFitValue <- ifelse(max, -.Machine$double.xmax, .Machine$double.xmax)
    }
    else { # NULL g.scale
      if (max == FALSE) {
        g.scale <- 1
        FiniteBadFitValue <- .Machine$double.xmax
      }
      else {
        g.scale <- -1
        FiniteBadFitValue <- -.Machine$double.xmax
      }
    }
    control$fnscale <- g.scale

  if(!lexical & !is.null(BFGSfn))
    {
      stop("'BFGSfn' can only be provided with lexical optimization or when 'transform=TRUE'")
    }
  if (!is.null(BFGSfn) & BFGS==FALSE)
    {
      if (!is.function(BFGSfn))
        stop("IF 'BFGSfn' is not a function, it must be NULL")
      warning("setting BFGS==TRUE because 'BFGSfn' is not null")
      BFGS <- TRUE
    }

  fn1 <- function(par) {
    fit <- fn(par, ...)

    if(is.null(fit))
      fit <- FiniteBadFitValue

    if(length(fit)==1)
      if(!is.finite(fit))
        fit <- FiniteBadFitValue

    return(fit)
  }#end of fn1

  if(!is.null(BFGShelp)) {
      if (!is.null(gr)) {
        gr1 <- function(par, helper = do.call(BFGShelp, 
					args = list(initial = par, done = TRUE))) {
				gr(par, helper, ...)
			}
      } else gr1 <- NULL
  }
  else {  
    if (!is.null(gr)) {
           gr1 <- function(par) gr(par, ...)
    } else gr1 <- NULL
  }

  #setpath to tempdir
  if(is.null(project.path))
    {
      project.path <- file.path(tempdir(), "genoud.pro")
    }
  
  #do we have stating values?
  if (is.null(starting.values)) {
    nStartingValues <- 0;
  }
  else if(is.matrix(starting.values)) {
    if(any(dim(starting.values) == nvars)) {
      if(nrow(starting.values) == nvars & ncol(starting.values) !=nvars) starting.values <- t(starting.values)
       nStartingValues <- nrow(starting.values)
       if(nStartingValues > pop.size) {
         warning("increasing 'pop.size' because too many starting.values were provided")
         pop.size <- nStartingValues
       }
    }
    else {
      warning("ignoring 'starting.values' because the wrong number of parameters was provided")
      nStartingValues <- 0
    }
  }
  else if(is.numeric(starting.values) | is.logical(starting.values)) {
    nStartingValues <- 1;

    if(length(starting.values)!=nvars)
      {
        nStartingValues <- 0
        warning("Ignoring 'starting.values' because length(staring.values)!=nvars")
      }
    else starting.values <- matrix(starting.values, nrow = 1)
  }
  else stop("starting.values must be NULL, a vector, or a matrix")

  #set output.type
  if (output.path=="stdout")
    {
      output.type <- 0;
    }
  else
    {
      if (output.append)
        {
          output.type <- 2;
        }
      else
        {
          output.type <- 1;
        }
    }

  # let's create the Domains if none have been passed.
  if (!(is.matrix(Domains)))
    {
      Domains <- matrix(nrow=nvars, ncol=2);
      for (i in 1:nvars)
        {
          Domains[i,1] <- -1*default.domains;
          Domains[i,2] <- default.domains;
        } # end of for loop
    } 
  else if(nrow(Domains) != nvars) {
    stop("number of rows in Domains must match 'nvars'")
  }
  else if(ncol(Domains) != 2) {
    stop("number of cols in Domains must be 2")
  }

  if(!all(Domains[,1] <= Domains[,2])) {
    stop("Domains[,1] must be less than or equal to Domains[,2]")
  }
  if(any(Domains[,1] == Domains[,2])) {
    warning("some Domains[,1]==Domains[,2]")
  }  

  # BG: now check all starting values are sane
  if(nStartingValues > 0 && any(is.na(starting.values))) {
	stop("Some starting values are NA")
  }
  if(nStartingValues > 0 && boundary.enforcement != 0 && 
     !all(apply(starting.values, 1, FUN = function(x) 
               Domains[,1] <= x & x <= Domains[,2])) )
        warning("'starting.values' which are outside of the bounds have been provided.
           Continuing, but unexpected behavior can occur with 'boundary.enforcement!=0'")

  # has the user provided any seeds?
  if (unif.seed==812821 && int.seed==53058)
    provide.seeds <- FALSE
  else
    provide.seeds <- TRUE;

  #if lexical==TRUE we need to know how many items will be returned
  if (lexical < 0)
    {
      warning("lexical < 0.  Resetting to FALSE\n")
      lexical <- 0
    }
  if (lexical>=1) 
    {
      #creating visible binding for variable "indx"; although it is
      #actually defined in fnLexicalSort() via an eval and paste
      if (!exists("indx"))
        {
          indx <- NULL; rm(indx)
        }
      if(share.type > 0) {
        warning("'share.type' being set to 0 because of lexical optimization")
        share.type <- 0
      }
      
      if(nStartingValues)
        {
          foo <- fn1(starting.values[1,])          
        } else {
          rfoo <- stats::runif(nrow(Domains), Domains[,1], Domains[,2])
          if(data.type.int)
            rfoo <- as.integer(round(rfoo))
          foo <- fn1(rfoo)          
        }
	foo.length <- length(as.vector(foo))
	if(lexical > 1 && foo.length != lexical) {
	  warning(paste("Function returns a vector of length", foo.length, 
                        "\nbut you specified lexical =", lexical))
	}
        if(foo.length == 1) {
          lexical <- 0
          warning("you specified lexical = TRUE but the function returns a scalar")
	}
        else lexical <- foo.length
    }
  else foo.length <- 1

  if (lexical > 0)
    {
      if(is.null(BFGSfn))
         {
           #All derivative stuff is turned off if we are going to do lexical if BFGSfn is not provided
           BFGS=FALSE
           gradient.check=FALSE
           if(hessian) {
             warning("'hessian' being set to false because of lexical optimization.  See 'BFGSfn' for workaround")
             hessian=FALSE             
           }

           P9 = 0
         } else {
             if(!is.null(BFGShelp)) {
               fn1.bfgs <-  function(par, helper = do.call(BFGShelp, 
                                     args = list(initial = par, done = TRUE), 
                                     envir = environment(fn))) {
                 fit <- BFGSfn(par, helper, ...) 
             
                 if(is.null(fit)) fit <- FiniteBadFitValue
             
                 if(length(fit)==1) if(!is.finite(fit)) fit <- FiniteBadFitValue
             
                 return(fit)
               }#end of fn1.bfgs
             } else {
               fn1.bfgs <-  function(par) {
                 fit <- BFGSfn(par, ...) 
             
                 if(is.null(fit)) fit <- FiniteBadFitValue
             
                 if(length(fit)==1) if(!is.finite(fit)) fit <- FiniteBadFitValue
             
                 return(fit)
               }#end of fn1.bfgs
             } # end else

           if(is.null(gr)) {
             gr <- function(par, helper = NA, ...)
               {
                  gr.fn1.bfgs <- function(par, helper, FBFV) {
                    fit <- if(is.null(BFGShelp)) BFGSfn(par, ...) else BFGSfn(par, helper, ...) 
                    
                    if(is.null(fit))
                      fit <- FBFV
                    
                    if(length(fit)==1)
                      if(!is.finite(fit))
                        fit <- FBFV
                    
                    return(fit)
                  }  # end of gr.fn1.bfgs               
                 genoud.wrapper101.env <- new.env()
                 assign("x", par, envir = genoud.wrapper101.env)
                 assign("helper", helper, envir = genoud.wrapper101.env)
                 assign("FiniteBadFitValue", FiniteBadFitValue, envir = genoud.wrapper101.env)
                 foo <- as.double(attr(stats::numericDeriv(quote(gr.fn1.bfgs(x, helper, FiniteBadFitValue)), theta=c("x"), genoud.wrapper101.env), "gradient"))
                 return(foo)
               } #end of gr
             gr1 <- if(is.null(BFGShelp)) function(par, ...) gr(par) else
                    function(par, helper = do.call(BFGShelp, args = list(initial = par,
					    done = TRUE), envir = environment(fn))) {
				gr(par, helper, ...)
			} # end of gr1
	   } # end of if(!is.null(gr))
           gr1func <- gr1
         }# end of else
    }#if lexical > 0
  if (lexical==0)
    lexical <- 1

  #optim st
  if(is.null(BFGSfn))
     {
       if(optim.method != "L-BFGS-B") {
         genoud.optim.wrapper101 <- function(foo.vals)
           {
             ret <- stats::optim(foo.vals, fn=fn1, gr=gr1, method=optim.method,
                          control=control);
             return(c(ret$value,ret$par));
           } # end of genoud.optim.wrapper101
       }
       else {
         genoud.optim.wrapper101 <- function(foo.vals)
           {
             ret <- stats::optim(foo.vals, fn=fn1, gr=gr1, method="L-BFGS-B",
                          lower = Domains[,1], upper = Domains[,2],
                          control=control);
             return(c(ret$value,ret$par));
           } # end of genoud.optim.wrapper101
       }
     } else {
       if(optim.method != "L-BFGS-B") {
         genoud.optim.wrapper101 <- function(foo.vals)
           {
             if(print.level > 2) {
               fit <- fn1(foo.vals)
               cat("\nPre-BFGS Complete Lexical Fit:\n")
               print(fit)
             }
             if(is.null(BFGShelp)) {
               ret <- stats::optim(foo.vals, fn=fn1.bfgs, gr=gr1, method=optim.method,
                            control=control);
             }
             else {
               ret <- stats::optim(foo.vals, fn=fn1.bfgs, gr=gr1, method=optim.method,
                            control=control,
                            helper = do.call(BFGShelp, args = list(initial = foo.vals), envir = environment(fn)) );
             }
             
             if(print.level > 2)
               {
                 cat("BFGS Fit:",ret$value,"\n")
                 
                 fit <- fn1(ret$par)
                 cat("Post-BFGS Complete Lexical Fit:\n")
                 print(fit)
               }           
             
             foo <- c(ret$value,ret$par)
             return(foo);
           } # end of genoud.optim.wrapper101
       } else { # "L-BFGS-B"
         genoud.optim.wrapper101 <- function(foo.vals)
           {
             if(print.level > 2) {
               fit <- fn1(foo.vals)
               cat("\nPre-BFGS Complete Lexical Fit:\n")
               print(fit)
             }
             if(is.null(BFGShelp)) {
               ret <- stats::optim(foo.vals, fn=fn1.bfgs, gr=gr1, method="L-BFGS-B",
                            lower = Domains[,1],
                            upper = Domains[,2],
                            control=control);
             }
             else {
               ret <- stats::optim(foo.vals, fn=fn1.bfgs, gr=gr1, method="L-BFGS-B",
                            lower = Domains[,1],
                            upper = Domains[,2],
                            control=control,
                            helper = do.call(BFGShelp, args = list(initial = foo.vals), envir = environment(fn)) );
             }
             
             if(print.level > 2)
               {
                 cat("BFGS Fit:",ret$value,"\n")
                 
                 fit <- fn1(ret$par)
                 cat("Post-BFGS Complete Lexical Fit:\n")
                 print(fit)
               }           
             
             foo <- c(ret$value,ret$par)
             return(foo);
           } # end of genoud.optim.wrapper101
       } 
     }
  
  # create the P vector
  P <- vector(length=9, mode="numeric");
  P[1] <- P1; P[2] <- P2; P[3] <- P3; P[4] <- P4;
  P[5] <- P5; P[6] <- P6; P[7] <- P7; P[8] <- P8;
  P[9] <- P9;

  clustertrigger=1
  if (is.logical(cluster))
    {
      if (cluster==FALSE)  {
        clustertrigger=0
      } else {
        stop("cluster option must be either FALSE, an object of the 'cluster' class (from the 'parallel' package) or a list of machines so 'genoud' can create such an object")
      }
    }

  if(clustertrigger) {
    parallel.exists = requireNamespace("parallel")
    if (!parallel.exists) {
      stop("The 'cluster' feature cannot be used unless the package 'parallel' can be loaded.")
    }
  }

  if(clustertrigger)
    {
      if(!MemoryMatrix)
        {
          MemoryMatrix = TRUE
          warning("MemoryMatrix has been set to 'TRUE' because the cluster option has been requested.")
        }

      GENclusterExport <- function (cl, list) 
        {
          gets <- function(n, v) {
            assign(n, v)
            NULL
          }
          for (name in list) {
            parallel::clusterCall(cl, gets, name, get(name))
          }
        }

      ExpandDots  <- function(...)
        {
          return(match.call())
        }      

      dots  <- ExpandDots(...)
      if ( (length(dots) > 1) & balance==TRUE)
        {
          balance  <- FALSE
          warning("load balancing has been turned off because the function to be optimized requires extra arguments")
        }
      
      if (class(cluster)[1]=="SOCKcluster" | class(cluster)[1]=="PVMcluster" | class(cluster)[1]=="spawnedMPIcluster" | class(cluster)[1]=="MPIcluster") {
        clustertrigger=1
        cl <- cluster
        GENclusterExport(cl, "fn")
        GENclusterExport(cl, "fn1")
      } else {
        clustertrigger=2
        cluster <- as.vector(cluster)
        cat("Initializing Cluster\n")
        cl <- parallel::makePSOCKcluster(cluster)
        GENclusterExport(cl, "fn")        
        GENclusterExport(cl, "fn1")
      }

      if (length(cl) < 2 )
        {
          stop("You only have one node. You probably should not be using the cluster option")
        }
    }

  fnLexicalSort <- function(mat, parms)
    {
      # parms = (1) MinMax, (2) nvars, (3) lexical_end, (4) type (0, vars, 1, obj)

      decreasing=FALSE
      if(parms[1]==1)
        decreasing=TRUE

      if(parms[4]==0)
         {
           #on nvars not on obj function
           foo = "indx <- order("
           for(i in (2:(parms[2]+1)) )
             {
               foo = paste(foo,"mat[,",i,"], ",sep="")
             }
           foo = paste(foo,"mat[,",parms[2]+2,"], ",sep="")
           foo = paste(foo,"decreasing=FALSE)",sep="")
           eval(parse(text=foo))
           mat = mat[indx,]           
         } else {
           #lexical on obj function
           foo = "indx <- order(mat[,1], "
           for(i in (parms[2]+3):parms[3] )
             {
               foo = paste(foo,"mat[,",i,"], ",sep="")
             }
           foo = paste(foo,"decreasing=",decreasing,")",sep="")
           eval(parse(text=foo))
           mat = mat[indx,]
         }
      return(mat)
    }

  fnMemoryMatrixEvaluate <- function(Memory, population, parms)
    {
      EVALUATE = -93813381

      MinMax = parms[1]
      nvars = parms[2]
      lexical = parms[3]
      
      lexical.end = ncol(population)
      pop.size = nrow(population)

      vars.indx <- 2:(nvars+1)
      lexical.indx <- c(1,(nvars+3):lexical.end)

      FIRSTTIME = TRUE
      if (nrow(Memory) > 1)
        FIRSTTIME = FALSE

      nevaluate = 0

      mfunc <- function(pop,memory)
        {
          f <- function(...) paste(..., sep=":")
          a2 <- do.call("f", as.data.frame(pop))
          b2 <- do.call("f", as.data.frame(memory))          
          
          return(match(a2,b2))
        }

      population.mat = matrix(population[,vars.indx], ncol=nvars)
      if (!FIRSTTIME)
        {
          Memory.mat = matrix(Memory[,vars.indx], ncol=nvars)
          match.matrix.indx <- mfunc(population.mat,Memory.mat)

          for (i in 1:pop.size)
            {
              match.indx = match.matrix.indx[i]
              found.match <- FALSE
              if ( is.finite(match.indx) )
                {
                  #this is needed because mfunc truncates floats (by the call to paste) to 15 digits
#                  if (all(population.mat[i,]==Memory.mat[match.indx,]))
#                    {
                      found.match <- TRUE
                      population[i,] <- Memory[match.indx,]
#                    }
                }
              if(!found.match)
                {
                  if (population[i,nvars+2] != 0) {
                    population[i,nvars+2] = EVALUATE
                    nevaluate = nevaluate+1
                  }
                }
            }
        } else {
          for (i in 1:pop.size)
            {
              population[i,nvars+2] = EVALUATE
              nevaluate = nevaluate+1
            }
        }

      #evaluation loop
      if (nevaluate > 0)
        {
          eval.indx <- population[,nvars+2]==EVALUATE
          ret=0
          if (clustertrigger==0)
            {
              in.mat = matrix(population.mat[eval.indx,], ncol=nvars)
              ret <- matrix(t(apply(in.mat, 1, fn1)), ncol=lexical)
            } else {
              if (balance==TRUE) {
                in.mat = t(matrix(population.mat[eval.indx,], ncol=nvars))
                cl.in <- as.list(as.data.frame(in.mat))
                cl.out <- parallel::clusterApplyLB(cl, cl.in, fn1)
                try(ret <- matrix(t(data.frame(cl.out)), ncol=lexical), TRUE)
                if (!is.matrix(ret)) {
                  if (!debug) {
                    stop("Cluster returned an object which could not be turned into a data.frame.  Cluster may be in bad state.  Please consider restarting it. To see what the cluster returned please turn on the 'debug' option.")
                  } else {
                    cat("Cluster returned an object which could not be turned into a data.frame:\n")
                    print(cl.out)
                    stop("Cluster may be in bad state.  Please consider restarting it.")
                  }
                }
              } else {
                in.mat = matrix(population.mat[eval.indx,], ncol=nvars)
                if(length(cl) > 1 )
                  {
                    cl.out = parallel::parRapply(cl, in.mat, fn1)
                  } else {
                    stop("You only have one node. Cluster option should not be used")
                  } 

                try(ret <- matrix(cl.out, byrow=TRUE, ncol=lexical), TRUE)
                if (!is.matrix(ret)) {
                  if (!debug) {
                    stop("Cluster returned an object which could not be turned into a data.frame.  Cluster may be in bad state.  Please consider restarting it. To see what the cluster returned please turn on the 'debug' option.")
                  } else {
                    cat("Cluster returned an object which could not be turned into a data.frame:\n")
                    print(cl.out)
                    stop("Cluster may be in bad state.  Please consider restarting it.")
                  }
                } 
              }
            } # else clustertrigger==0

          if (lexical < 2)
            {
              population[eval.indx,1] = ret[,1]
            } else {
              population[eval.indx,lexical.indx] = ret
            }
          population[eval.indx,nvars+2] = 0
          
          if(!FIRSTTIME)
            {
              Memory = rbind(Memory,population[eval.indx,])
            } else {
              Memory = matrix(population[eval.indx,], ncol=lexical.end)
              FIRSTTIME = FALSE
            }
        }#end of nevaluate
      
      if (lexical > 1)
        {
          population <- fnLexicalSort(population, c(MinMax,nvars,lexical.end,1))
        } else {
          if (MinMax==0)
            {
              population <- population[order(population[,1]),]
            } else {
              population <- population[order(population[,1], decreasing=TRUE),]
            }
        }

      return(as.vector(c(nrow(Memory), Memory, population)))
    } #end of fnMemoryMatrixEvaluate

  if (!is.null(gr))
    {
      UserGradient = TRUE
      gr1func <- gr1
    } else {
      UserGradient = FALSE      
      gr1func <- function() {}
    }

  if(data.type.int)
    {
      BFGS = FALSE
      gradient.check=FALSE
    }

  if(is.matrix(starting.values))
    starting.values <- t(starting.values)

  #C++ code checks if at least Generation 0 has been run and if
  #print.level>0 and project.path!="/dev/null".  Otherwise, 'interrupted'
  #remains FALSE
  interrupted <- FALSE
  interrupt.message <- paste("genoud interrupted:\none may recover the best solution found so far by executing")
  interrupt.expression <- paste("pop <- read.table('",project.path, 
				"', comment.char = 'G')", sep = "")
  interrupt.expression2 <- "best <- pop[pop$V1 == 1,, drop = FALSE]"
  interrupt.expression3 <- paste("very.best <- as.matrix(best[nrow(best), ", 
				foo.length + 2, ":ncol(best)])", sep = "")
  on.exit(if(interrupted) cat(interrupt.message, "\n", interrupt.expression, "\n",
				interrupt.expression2, "\n", interrupt.expression3, "\n"))
  
  gout <- .Call("rgenoud", as.function(fn1), new.env(),
                as.integer(nvars), as.integer(pop.size), as.integer(max.generations),
                as.integer(wait.generations),
                as.integer(nStartingValues), as.double(starting.values),
                as.vector(P), as.matrix(Domains),
                as.integer(max), as.integer(gradient.check), as.integer(boundary.enforcement),
                as.double(solution.tolerance), as.integer(BFGS), as.integer(data.type.int),
                as.integer(provide.seeds), as.integer(unif.seed), as.integer(int.seed),
                as.integer(print.level), as.integer(share.type), as.integer(instance.number),
                as.integer(MemoryMatrix), as.integer(debug),
                as.character(output.path), as.integer(output.type), as.character(project.path),
                as.integer(hard.generation.limit),
                as.function(genoud.optim.wrapper101), 
                as.integer(lexical), as.function(fnLexicalSort), as.function(fnMemoryMatrixEvaluate),
                as.integer(UserGradient), as.function(gr1func), as.double(P9mix),
                as.integer(BFGSburnin), as.integer(transform),
                PACKAGE="rgenoud");

  indx1 <- 4;
  indx2 <- (indx1+lexical-1);
  value=gout[indx1:indx2];

  indx1 <- indx2+1
  indx2 <- indx1+nvars-1
  par = gout[indx1:indx2]

  indx1 <- indx2+1
  indx2 <- indx1+nvars-1
  if(!gradient.check & !BFGS )
    {
      gradients= gout[indx1:indx2]
      gradients = rep(NA, length(gradients))
    } else {
      gradients = gout[indx1:indx2]
    }

  indx1 <- indx2+1
  indx2 <- indx1+8
  operators=gout[indx1:indx2]

  if (hessian==TRUE)
    {
      con <- list(trace = 0, fnscale = g.scale,
                  parscale = rep.int(1, length(par)), ndeps = rep.int(0.001, length(par)), maxit = 100,
                  abstol = -Inf, reltol = sqrt(.Machine$double.eps), alpha = 1,
                  beta = 0.5, gamma = 2, REPORT = 10, type = 1, lmm = 5,
                  factr = 1e+07, pgtol = 0, tmax = 10, temp = 10)
      
      nm <- names(gout[5:(nvars+4)])
      
      if(lexical == 1) {
        hess.fn <- function(par) fn1(par, ...)
        hess <- stats::optimHess(par=gout[5:(nvars+4)], fn=fn1, gr=gr1, control=con)
      }
      else {
        help.stuff <- do.call(BFGShelp, args = list(initial = gout[5:(nvars+4)], done = TRUE), 
                              envir = environment(fn))
        hess.fn <- function(par, helper = help.stuff) fn1.bfgs(par, helper, ...)
        hess <- stats::optimHess(par=gout[5:(nvars+4)], fn=hess.fn, gr=NULL, control=con)
      }
      
      hes <- 0.5 * (hess + t(hess))
      if (!is.null(nm)) dimnames(hes) <- list(nm, nm)
      
      
      ret <- list(value=value, par=par, gradients=gradients,
                  generations=gout[1], peakgeneration=gout[2], popsize=gout[3],
                  operators=operators,
                  hessian=hes)
    }
  else
    {
      ret <- list(value=value, par=par, gradients=gradients,
                  generations=gout[1], peakgeneration=gout[2], popsize=gout[3],
                  operators=operators)
    }
  
  if (clustertrigger==2)
    parallel::stopCluster(cl)

  interrupted <- FALSE
  
  return(ret)  
} #end of genoud()

genoud_transform <- function(fn, nvars, max=FALSE, pop.size=1000, max.generations=100, wait.generations=10,
                             hard.generation.limit=TRUE, starting.values=NULL, MemoryMatrix=TRUE,
                             Domains=NULL, default.domains=10, solution.tolerance=0.001,
                             gr=NULL, boundary.enforcement=0, lexical=FALSE, gradient.check=TRUE, BFGS=TRUE,
                             data.type.int=FALSE, hessian=FALSE, unif.seed=812821, int.seed=53058,
                             print.level=2, share.type=0, instance.number=0,
                             output.path="stdout", output.append=FALSE, project.path=NULL,
                             P1=50, P2=50, P3=50, P4=50, P5=50, P6=50, P7=50, P8=50, P9=0,
                             P9mix=NULL, BFGSburnin=0, BFGSfn=NULL, BFGShelp = NULL,
                             control = list(), optim.method=ifelse(boundary.enforcement < 2, "BFGS", "L-BFGS-B"),
                             transform=FALSE, debug=FALSE, cluster=FALSE, balance=FALSE,  ...)
{
  #wrapper for the transform option. beta changes should go here as well
  if(data.type.int)
    stop("transform option cannot be used when 'data.type.int=TRUE'")

  if(MemoryMatrix)
    warning("setting MemoryMatrix to FALSE")
  
  MemoryMatrix <- FALSE
  if(!identical(cluster, FALSE)) warning("setting cluster to FALSE")
  cluster <- FALSE

  if(!is.null(BFGShelp) && !is.function(BFGShelp)) stop("'BFGShelp' must be NULL or a function")

  #creating visible binding for variable "indx"; although it is
  #actually defined in fnLexicalSort() via an eval and paste
  if (!exists("indx"))
    {
      indx <- NULL; rm(indx)
    }

  
  if(!is.null(P9mix) && !is.double(P9mix))  {
    stop("'P9mix' must be NULL or a number between 0 and 1")
  } else {
    if(is.null(P9mix)) {
      P9mix <- -1
    } else {
      if(! ( (1 >= P9mix) && (P9mix > 0) ))
        stop("'P9mix' must be NULL or a number between 0 and 1 (it may be equal to 1)")
    }
  }

  if( (BFGSburnin < 0) & !gradient.check )
    {
      warning("If 'BFGSburnin' is negative, gradient.check must be TRUE for 'BFGSburnin' to have any effect.")
    }

  if(!is.list(control))
    stop("'control' must be a list, see ?optim")
    g.scale <- control$fnscale
    if(!is.null(g.scale)) {
      if(g.scale > 0 & max) {
        stop("positive control$fnscale is inconsistent with maximization")
      }
      else if(g.scale < 0 & !max) {
        stop("negative control$fnscale is inconsistent with minimization")
      }
      else if(g.scale == 0) {
        stop("optim divides the function value by control$fnscale ",
             "setting control$fnscale to zero is therefore impossible")
      }
      FiniteBadFitValue <- ifelse(max, -.Machine$double.xmax, .Machine$double.xmax)
    }
    else { # NULL g.scale
      if (max == FALSE) {
        g.scale <- 1
        FiniteBadFitValue <- .Machine$double.xmax
      }
      else {
        g.scale <- -1
        FiniteBadFitValue <- -.Machine$double.xmax
      }
    }
    control$fnscale <- g.scale

  ## BFGSfn can (and probably should) be provided if transform = TRUE
#   if(!lexical & !is.null(BFGSfn))
#     {
#       stop("'BFGSfn' can only be provided with lexical optimization or when 'transform=TRUE'")
#     }
  if (!is.null(BFGSfn) & BFGS==FALSE)
    {
      if (!is.function(BFGSfn))
        stop("IF 'BFGSfn' is not a function, it must be NULL")
      warning("setting BFGS==TRUE because 'BFGSfn' is not null")
      BFGS <- TRUE
    }

  fn1 <- function(par) {
    fit <- fn(par, ...)
    if(!is.finite(fit[1])) fit[1] <- FiniteBadFitValue ## only catch [1]
    return(fit)
  }

  if(!is.null(BFGShelp)) {
      if (!is.null(gr)) {
        gr1 <- function(par, helper = do.call(BFGShelp,
					args = list(initial = par, done = TRUE))) {
				gr(par, helper, ...)
			}
      } else gr1 <- NULL
  }
  else {
    if (!is.null(gr)) {
           gr1 <- function(par) gr(par, ...)
    } else gr1 <- NULL
  }

  #setpath to tempdir
  if(is.null(project.path))
    {
      project.path <- file.path(tempdir(), "genoud.pro")
    }

  #do we have stating values?
  if (is.null(starting.values)) {
    nStartingValues <- 0;
  }
  else if(is.matrix(starting.values)) {
    if(any(dim(starting.values) == nvars)) {
      if(nrow(starting.values) == nvars & ncol(starting.values) !=nvars) starting.values <- t(starting.values)
      nStartingValues <- nrow(starting.values)
      if(nStartingValues > pop.size) {
        warning("increasing 'pop.size' because too many starting.values were provided")
        pop.size <- nStartingValues
      }
    }
    else {
      warning("ignoring 'starting.values' because the wrong number of parameters was provided")
      nStartingValues <- 0
    }
  }
  else if(is.numeric(starting.values) | is.logical(starting.values)) {
    nStartingValues <- 1;

    if(length(starting.values)!=nvars)
      {
        nStartingValues <- 0
        warning("Ignoring 'starting.values' because length(staring.values)!=nvars")
      }
    else starting.values <- matrix(starting.values, nrow = 1)
  }
  else stop("starting.values must be NULL, a vector, or a matrix")

  #set output.type
  if (output.path=="stdout")
    {
      output.type <- 0;
    }
  else
    {
      if (output.append)
        {
          output.type <- 2;
        }
      else
        {
          output.type <- 1;
        }
    }

  # let's create the Domains if none have been passed.
  if (!(is.matrix(Domains)))
    {
      Domains <- matrix(nrow=nvars, ncol=2);
      for (i in 1:nvars)
        {
          Domains[i,1] <- -1*default.domains;
          Domains[i,2] <- default.domains;
        } # end of for loop
    }
  else if(nrow(Domains) != nvars) {
    stop("number of rows in Domains must match 'nvars'")
  }
  else if(ncol(Domains) != 2) {
    stop("number of cols in Domains must be 2")
  }

  if(!all(Domains[,1] <= Domains[,2])) {
    stop("Domains[,1] must be less than or equal to Domains[,2]")
  }
  if(any(Domains[,1] == Domains[,2])) {
    warning("some Domains[,1]==Domains[,2]")
  }

  # BG: now check all starting values are sane
  if(nStartingValues > 0 && any(is.na(starting.values))) {
    stop("Some starting values are NA")
  }
  if(nStartingValues > 0 && boundary.enforcement != 0 &&
     !all(apply(starting.values, 1, FUN = function(x)
               Domains[,1] <= x & x <= Domains[,2])) )
        warning("'starting.values' which are outside of the bounds have been provided.
           Continuing, but unexpected behavior can occur with 'boundary.enforcement!=0'")

  # has the user provided any seeds?
  if (unif.seed==812821 && int.seed==53058)
    provide.seeds <- FALSE
  else
    provide.seeds <- TRUE;

  # we need to know how many items will be returned bf fn()
  if(nStartingValues)
        {
          foo <- fn1(starting.values[1,])
        } else {
          rfoo <- stats::runif(nrow(Domains), Domains[,1], Domains[,2])
          foo <- fn1(rfoo)
  }
  foo.length <- length(c(foo))
  if(foo.length <= nvars) stop("'fn' must return criteria and all transformed parameters")
  lexical <- foo.length - nvars ## it is possible that lexical == 1 or lexical > 1

  if (lexical >= 1) ## always TRUE
    {
      if(share.type > 0) {
        warning("'share.type' being set to 0 because of 'transform' option")
        share.type <- 0
      }
      if(is.null(BFGSfn))
         {
           #All derivative stuff is turned off if BFGSfn is not provided
           BFGS <- FALSE
           gradient.check <- FALSE
           if(hessian) {
             warning("'hessian' being set to false because of lexical optimization.  See 'BFGSfn' for workaround")
             hessian <- FALSE
           }

           P9 = 0
         } else {
             if(!is.null(BFGShelp)) {
               fn1.bfgs <-  function(par, helper = do.call(BFGShelp,
                                     args = list(initial = par, done = TRUE),
                                     envir = environment(fn))) {
                 fit <- BFGSfn(par, helper, ...)

                 if(is.null(fit)) fit <- FiniteBadFitValue

                 if(length(fit)==1) if(!is.finite(fit)) fit <- FiniteBadFitValue

                 return(fit)
               }#end of fn1.bfgs
             } else {
               fn1.bfgs <-  function(par) {
                 fit <- BFGSfn(par, ...)

                 if(is.null(fit)) fit <- FiniteBadFitValue

                 if(length(fit)==1) if(!is.finite(fit)) fit <- FiniteBadFitValue

                 return(fit)
               }#end of fn1.bfgs
             } # end else

           if(is.null(gr)) { ## should we do numerical gradients when transform = TRUE?
             gr <- function(par, helper = NA, ...)
               {
                  gr.fn1.bfgs <- function(par, helper, FBFV) {
                    fit <- if(is.null(BFGShelp)) BFGSfn(par, ...) else BFGSfn(par, helper, ...)

                    if(is.null(fit))
                      fit <- FBFV

                    if(length(fit)==1)
                      if(!is.finite(fit))
                        fit <- FBFV

                    return(fit)
                  }  # end of gr.fn1.bfgs
                 genoud.wrapper101.env <- new.env()
                 assign("x", par, envir = genoud.wrapper101.env)
                 assign("helper", helper, envir = genoud.wrapper101.env)
                 assign("FiniteBadFitValue", FiniteBadFitValue, envir = genoud.wrapper101.env)
                 foo <- as.double(attr(stats::numericDeriv(quote(gr.fn1.bfgs(x, helper, FiniteBadFitValue)), theta=c("x"), genoud.wrapper101.env), "gradient"))
                 return(foo)
               } #end of gr
             gr1 <- if(is.null(BFGShelp)) function(par, ...) gr(par) else
                    function(par, helper = do.call(BFGShelp, args = list(initial = par,
                            done = TRUE), envir = environment(fn))) {
                                gr(par, helper, ...)
             } # end of gr1
           } # end of if(!is.null(gr))
           gr1func <- gr1
      }# end of else
    }#if lexical > 1

  #optim st
  if(is.null(BFGSfn))
     {
       if(optim.method != "L-BFGS-B") {
         genoud.optim.wrapper101 <- function(foo.vals)
           {
             ret <- stats::optim(foo.vals, fn=fn1, gr=gr1, method=optim.method,
                          control=control);
             return(c(ret$value,ret$par));
           } # end of genoud.optim.wrapper101
       }
       else {
         genoud.optim.wrapper101 <- function(foo.vals)
           {
             ret <- stats::optim(foo.vals, fn=fn1, gr=gr1, method="L-BFGS-B",
                          lower = Domains[,1], upper = Domains[,2],
                          control=control);
             return(c(ret$value,ret$par));
           } # end of genoud.optim.wrapper101
       }
     } else {
       if(optim.method != "L-BFGS-B") {
         genoud.optim.wrapper101 <- function(foo.vals)
           {
             if(print.level > 2) {
               fit <- fn1(foo.vals)
               cat("\nPre-BFGS Complete Lexical Fit:\n")
               print(fit)
             }
             if(is.null(BFGShelp)) {
               ret <- stats::optim(foo.vals, fn=fn1.bfgs, gr=gr1, method=optim.method,
                            control=control);
             }
             else {
               ret <- stats::optim(foo.vals, fn=fn1.bfgs, gr=gr1, method=optim.method,
                            control=control,
                            helper = do.call(BFGShelp, args = list(initial = foo.vals), envir = environment(fn)) );
             }

             if(print.level > 2)
               {
                 cat("BFGS Fit:",ret$value,"\n")

                 fit <- fn1(ret$par)
                 cat("Post-BFGS Complete Lexical Fit:\n")
                 print(fit)
               }

             foo <- c(ret$value,ret$par)
             return(foo);
           } # end of genoud.optim.wrapper101
       } else { # "L-BFGS-B"
         genoud.optim.wrapper101 <- function(foo.vals)
           {
             if(print.level > 2) {
               fit <- fn1(foo.vals)
               cat("\nPre-BFGS Complete Lexical Fit:\n")
               print(fit)
             }
             if(is.null(BFGShelp)) {
               ret <- stats::optim(foo.vals, fn=fn1.bfgs, gr=gr1, method="L-BFGS-B",
                            lower = Domains[,1],
                            upper = Domains[,2],
                            control=control);
             }
             else {
               ret <- stats::optim(foo.vals, fn=fn1.bfgs, gr=gr1, method="L-BFGS-B",
                            lower = Domains[,1],
                            upper = Domains[,2],
                            control=control,
                            helper = do.call(BFGShelp, args = list(initial = foo.vals), envir = environment(fn)) );
             }

             if(print.level > 2)
               {
                 cat("BFGS Fit:",ret$value,"\n")

                 fit <- fn1(ret$par)
                 cat("Post-BFGS Complete Lexical Fit:\n")
                 print(fit)
               }

             foo <- c(ret$value,ret$par)
             return(foo);
           } # end of genoud.optim.wrapper101
       }
     }

  # create the P vector
  P <- vector(length=9, mode="numeric");
  P[1] <- P1; P[2] <- P2; P[3] <- P3; P[4] <- P4;
  P[5] <- P5; P[6] <- P6; P[7] <- P7; P[8] <- P8;
  P[9] <- P9;

  ## cluster is coerced to FALSE at the top so none of this matters
  clustertrigger=1
  if (is.logical(cluster))
    {
      if (cluster==FALSE)  {
        clustertrigger=0
      } else {
        stop("cluster option must be either FALSE, an object of the 'cluster' class (from the 'parallel' package) or a list of machines so 'genoud' can create such an object")
      }
    }

  if(clustertrigger) {
    parallel.exists = requireNamespace("parallel")
    if (!parallel.exists) {
      stop("The 'cluster' feature cannot be used unless the package 'parallel' can be loaded.")
    }
  }

  if(clustertrigger)
    {
      if(!MemoryMatrix)
        {
          MemoryMatrix = TRUE
          warning("MemoryMatrix has been set to 'TRUE' because the cluster option has been requested.")
        }

      GENclusterExport <- function (cl, list)
        {
          gets <- function(n, v) {
            assign(n, v)
            NULL
          }
          for (name in list) {
            parallel::clusterCall(cl, gets, name, get(name))
          }
        }

      ExpandDots  <- function(...)
        {
          return(match.call())
        }

      dots  <- ExpandDots(...)
      if ( (length(dots) > 1) & balance==TRUE)
        {
          balance  <- FALSE
          warning("load balancing has been turned off because the function to be optimized requires extra arguments")
        }

      if (class(cluster)[1]=="SOCKcluster" | class(cluster)[1]=="PVMcluster" | class(cluster)[1]=="spawnedMPIcluster" | class(cluster)[1]=="MPIcluster") {
        clustertrigger=1
        cl <- cluster
        GENclusterExport(cl, "fn")
        GENclusterExport(cl, "fn1")
      } else {
        clustertrigger=2
        cluster <- as.vector(cluster)
        cat("Initializing Cluster\n")
        cl <- parallel::makePSOCKcluster(cluster)
        GENclusterExport(cl, "fn")
        GENclusterExport(cl, "fn1")
      }

      if (length(cl) < 2 )
        {
          stop("You only have one node. You probably should not be using the cluster option")
        }
    }

  fnLexicalSort <- function(mat, parms)
    {
      # parms = (1) MinMax, (2) nvars, (3) lexical_end, (4) type (0, vars, 1, obj)

      decreasing=FALSE
      if(parms[1]==1)
        decreasing=TRUE

      if(parms[4]==0)
         {
           #on nvars not on obj function
           foo = "indx <- order("
           for(i in (2:(parms[2]+1)) )
             {
               foo = paste(foo,"mat[,",i,"], ",sep="")
             }
           foo = paste(foo,"mat[,",parms[2]+2,"], ",sep="")
           foo = paste(foo,"decreasing=FALSE)",sep="")
           eval(parse(text=foo))
           mat = mat[indx,]
         } else {
           #lexical on obj function
           foo = "indx <- order(mat[,1], "
           for(i in (parms[2]+3):parms[3] )
             {
               foo = paste(foo,"mat[,",i,"], ",sep="")
             }
           foo = paste(foo,"decreasing=",decreasing,")",sep="")
           eval(parse(text=foo))
           mat = mat[indx,]
         }
      return(mat)
    }

  fnMemoryMatrixEvaluate <- function(Memory, population, parms)
    {
      stop("this function should not be called when transform = TRUE")
      EVALUATE = -93813381

      MinMax = parms[1]
      nvars = parms[2]
      lexical = parms[3]

      lexical.end = ncol(population)
      pop.size = nrow(population)

      vars.indx <- 2:(nvars+1)
      lexical.indx <- c(1,(nvars+3):lexical.end)

      FIRSTTIME = TRUE
      if (nrow(Memory) > 1)
        FIRSTTIME = FALSE

      nevaluate = 0

      mfunc <- function(pop,memory)
        {
          f <- function(...) paste(..., sep=":")
          a2 <- do.call("f", as.data.frame(pop))
          b2 <- do.call("f", as.data.frame(memory))

          return(match(a2,b2))
        }

      population.mat = matrix(population[,vars.indx], ncol=nvars)
      if (!FIRSTTIME)
        {
          Memory.mat = matrix(Memory[,vars.indx], ncol=nvars)
          match.matrix.indx <- mfunc(population.mat,Memory.mat)

          for (i in 1:pop.size)
            {
              match.indx = match.matrix.indx[i]
              found.match <- FALSE
              if ( is.finite(match.indx) )
                {
                  #this is needed because mfunc truncates floats (by the call to paste) to 15 digits
#                  if (all(population.mat[i,]==Memory.mat[match.indx,]))
#                    {
                      found.match <- TRUE
                      population[i,] <- Memory[match.indx,]
#                    }
                }
              if(!found.match)
                {
                  if (population[i,nvars+2] != 0) {
                    population[i,nvars+2] = EVALUATE
                    nevaluate = nevaluate+1
                  }
                }
            }
        } else {
          for (i in 1:pop.size)
            {
              population[i,nvars+2] = EVALUATE
              nevaluate = nevaluate+1
            }
        }

      #evaluation loop
      if (nevaluate > 0)
        {
          eval.indx <- population[,nvars+2]==EVALUATE
          ret=0
          if (clustertrigger==0)
            {
              in.mat = matrix(population.mat[eval.indx,], ncol=nvars)
              ret <- matrix(t(apply(in.mat, 1, fn1)), ncol=lexical)
            } else {
              if (balance==TRUE) {
                in.mat = t(matrix(population.mat[eval.indx,], ncol=nvars))
                cl.in <- as.list(as.data.frame(in.mat))
                cl.out <- parallel::clusterApplyLB(cl, cl.in, fn1)
                try(ret <- matrix(t(data.frame(cl.out)), ncol=lexical), TRUE)
                if (!is.matrix(ret)) {
                  if (!debug) {
                    stop("Cluster returned an object which could not be turned into a data.frame.  Cluster may be in bad state.  Please consider restarting it. To see what the cluster returned please turn on the 'debug' option.")
                  } else {
                    cat("Cluster returned an object which could not be turned into a data.frame:\n")
                    print(cl.out)
                    stop("Cluster may be in bad state.  Please consider restarting it.")
                  }
                }
              } else {
                in.mat = matrix(population.mat[eval.indx,], ncol=nvars)
                if(length(cl) > 1 )
                  {
                    cl.out = parallel::parRapply(cl, in.mat, fn1)
                  } else {
                    stop("You only have one node. Cluster option should not be used")
                  }

                try(ret <- matrix(cl.out, byrow=TRUE, ncol=lexical), TRUE)
                if (!is.matrix(ret)) {
                  if (!debug) {
                    stop("Cluster returned an object which could not be turned into a data.frame.  Cluster may be in bad state.  Please consider restarting it. To see what the cluster returned please turn on the 'debug' option.")
                  } else {
                    cat("Cluster returned an object which could not be turned into a data.frame:\n")
                    print(cl.out)
                    stop("Cluster may be in bad state.  Please consider restarting it.")
                  }
                }
              }
            } # else clustertrigger==0

          if (lexical < 2)
            {
              population[eval.indx,1] = ret[,1]
            } else {
              population[eval.indx,lexical.indx] = ret
            }
          population[eval.indx,nvars+2] = 0

          if(!FIRSTTIME)
            {
              Memory = rbind(Memory,population[eval.indx,])
            } else {
              Memory = matrix(population[eval.indx,], ncol=lexical.end)
              FIRSTTIME = FALSE
            }
        }#end of nevaluate

      if (lexical > 1)
        {
          population <- fnLexicalSort(population, c(MinMax,nvars,lexical.end,1))
        } else {
          if (MinMax==0)
            {
              population <- population[order(population[,1]),]
            } else {
              population <- population[order(population[,1], decreasing=TRUE),]
            }
        }

      return(as.vector(c(nrow(Memory), Memory, population)))
    } #end of fnMemoryMatrixEvaluate

  if (!is.null(gr))
    {
      UserGradient = TRUE
      gr1func <- gr1
    } else {
      UserGradient = FALSE
      gr1func <- function() {}
    }

  ## impossible
  if(data.type.int)
    {
      BFGS = FALSE
      gradient.check=FALSE
    }

  if(is.matrix(starting.values))
    starting.values <- t(starting.values)

  #C++ code checks if at least Generation 0 has been run and if
  #print.level>0 and project.path!="/dev/null".  Otherwise, 'interrupted'
  #remains FALSE
  interrupted <- FALSE
  interrupt.message <- paste("genoud interrupted:\none may recover the best solution found so far by executing")
  interrupt.expression <- paste("pop <- read.table('",project.path,
				"', comment.char = 'G')", sep = "")
  interrupt.expression2 <- "best <- pop[pop$V1 == 1,, drop = FALSE]"
  interrupt.expression3 <- paste("very.best <- as.matrix(best[nrow(best), ",
				foo.length + 2, ":ncol(best)])", sep = "")
  on.exit(if(interrupted) cat(interrupt.message, "\n", interrupt.expression, "\n",
				interrupt.expression2, "\n", interrupt.expression3, "\n"))

  gout <- .Call("rgenoud", as.function(fn1), new.env(),
                as.integer(nvars), as.integer(pop.size), as.integer(max.generations),
                as.integer(wait.generations),
                as.integer(nStartingValues), as.double(starting.values),
                as.vector(P), as.matrix(Domains),
                as.integer(max), as.integer(gradient.check), as.integer(boundary.enforcement),
                as.double(solution.tolerance), as.integer(BFGS), as.integer(data.type.int),
                as.integer(provide.seeds), as.integer(unif.seed), as.integer(int.seed),
                as.integer(print.level), as.integer(share.type), as.integer(instance.number),
                as.integer(MemoryMatrix), as.integer(debug),
                as.character(output.path), as.integer(output.type), as.character(project.path),
                as.integer(hard.generation.limit),
                as.function(genoud.optim.wrapper101),
                as.integer(lexical), as.function(fnLexicalSort), as.function(fnMemoryMatrixEvaluate),
                as.integer(UserGradient), as.function(gr1func), as.double(P9mix),
                as.integer(BFGSburnin), as.integer(transform),
                PACKAGE="rgenoud");

  indx1 <- 4;
  indx2 <- (indx1+lexical-1);
  value=gout[indx1:indx2];

  indx1 <- indx2+1
  indx2 <- indx1+nvars-1
  par = gout[indx1:indx2]

  indx1 <- indx2+1
  indx2 <- indx1+nvars-1
  if(!gradient.check & !BFGS )
    {
      gradients= gout[indx1:indx2]
      gradients = rep(NA, length(gradients))
    } else {
      gradients = gout[indx1:indx2]
    }

  indx1 <- indx2+1
  indx2 <- indx1+8
  operators=gout[indx1:indx2]

  if (hessian==TRUE) ## only TRUE if BFGSfn is provided
    {
      con <- list(trace = 0, fnscale = g.scale,
                  parscale = rep.int(1, length(par)), ndeps = rep.int(0.001, length(par)), maxit = 100,
                  abstol = -Inf, reltol = sqrt(.Machine$double.eps), alpha = 1,
                  beta = 0.5, gamma = 2, REPORT = 10, type = 1, lmm = 5,
                  factr = 1e+07, pgtol = 0, tmax = 10, temp = 10)

      nm <- names(gout[5:(nvars+4)])

#       if(lexical == 1) {
#         hess.fn <- function(par) fn1(par, ...)
#         hess <- stats::optimHess(par=gout[5:(nvars+4)], fn=fn1, gr=gr1, control=con)
#       }
#       else {
        help.stuff <- do.call(BFGShelp, args = list(initial = gout[5:(nvars+4)], done = TRUE),
                              envir = environment(fn))
        hess.fn <- function(par, helper = help.stuff) fn1.bfgs(par, helper, ...)
        hess <- stats::optimHess(par=gout[5:(nvars+4)], fn=hess.fn, gr=NULL, control=con)
#       }

      hes <- 0.5 * (hess + t(hess))
      if (!is.null(nm)) dimnames(hes) <- list(nm, nm)


      ret <- list(value=value, par=par, gradients=gradients,
                  generations=gout[1], peakgeneration=gout[2], popsize=gout[3],
                  operators=operators,
                  hessian=hes)
    }
  else
    {
      ret <- list(value=value, par=par, gradients=gradients,
                  generations=gout[1], peakgeneration=gout[2], popsize=gout[3],
                  operators=operators)
    }

  if (clustertrigger==2)
    parallel::stopCluster(cl)

  interrupted <- FALSE

  return(ret)
}
