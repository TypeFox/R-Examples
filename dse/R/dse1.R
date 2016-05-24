###
###  There is still some function description here and in EvalEst that
###    should be moved to help files.
###

.onAttach  <- function(libname, pkgname) {
   .DSEflags(list(COMPILED=TRUE))
   invisible(TRUE)
   }

.DSEflags <- local({
   .DSE.flags <- character(0)
   function(new) {
       if(!missing(new))
           .DSE.flags <<- new
       else
           .DSE.flags
       }
})

# With NAMESPACEs use .onAttach instead of .First.Lib 
#.First.lib <- function(library, section) {
#   assign(".DSEflags()$COMPILED", TRUE, env = .GlobalEnv)
#   assign(".DSEDUP",      TRUE, env = .GlobalEnv)
#   invisible(library.dynam("dse"))
#   }
  

DSEversion <- function() 
  {if (!is.R()) return("version cannot be determined.") else
   z <- c("setRNG", "tframe", "dse","EvalEst",  
          "curve", "CDNmoney", "tsfa", "TSdbi")
   z <- z[ z %in% library()$results[,1] ]
   r <- NULL
   for (pac in z ) 
     r <- c(r,  packageDescription(pac, fields="Version"))
   names(r) <- z
   r
  }


##############################################################################

# The code was divided roughly into 
#   the groups listed below, but the organization has changed a little bit.
#   The code grouping can be seen by grep ing on the string '<<<<<<' eg:
#    grep "<<<<<<" dse.R

#    Functions which work on a model (i.e. if a model with data is allowed as
#             an arguement then the data is ignored):
#        -model summary, description, print, comparison and 
#              calculation of properties
#        -model conversion 
             
#    Functions which work on data ( and a model ):
#        -likelihood and residual calculation
#        -statistical tests
#        -parameter estimation 
#        -model reduction (this could work just on models but
#                    comparisons require data)     

#    Utility Functions:
#        -utilities for generating theoretical statistics 
#              (i.e.- model generated, not from data)
#        -data transformations
#        -model and data scaling
#        -utilities for polynomial arithmetic
#        -internal utilities used for updating objects 
#        -data interface functions

#   Of special note in the internal utilities are two programs,
#   setArrays and setTSmodelParameters which take a model list and
#   update the representation and parameter information respectively. 
#   i.e. setArrays uses the parameter information and ensures that the 
#    array representation is consistent while setTSmodelParameters uses the 
#    arrays and ensures that the parameter vector is consistent.

#############################################################################
#############################################################################
residuals.TSestModel <- function(object, ...) object$estimates$pred-outputData(object)

# I don't think this exists really works for namespaces
 if (!exists("acf.default", mode="function")){
   acf.default  <- stats::acf
   acf <- function(x, ...)UseMethod("acf")
   }
 
 acf.TSdata <- function(x, ...)acf(outputData(x))
 acf.TSestModel <- function(x, ...) acf(x$estimates$pred - outputData(x), ...)

#############################################################################

#functions which work on models   <<<<<<<<<<


############################################################

#     functions for model summary, description, print and
#      comparison and functions for calculation of model properties

############################################################


print.TSestModel <- function(x, ...) 
{ cat("neg. log likelihood=",x$estimates$like[1],"\n")
  if(!is.null(x$converged)) {if(!x$converged) cat(" (not converged)") }
  print(x$model,...) 
  invisible(x)
}

print.SS <- function(x, digits=options()$digits, latex=FALSE, ...) 
    {#  (... further arguments, currently disregarded)
     printM <- function(st, M, digits, latex)
	{arrayc <- function(n)
	   {cat("\\left[ \\begin{array}{"); for(j in 1:n) cat("c"); cat("}\n") }
	 if (latex) cat("\\begin{equation}\n")
	 cat("\n", st, "=\n")
	 col <- if (is.matrix(M)) dim(M)[2] else length(M)
	 if (latex) 
	   {arrayc(col)
	    if (!is.matrix(M)) M <- matrix(M, 1, length(M))
	    M <- signif(M, digits=digits)
            for (i in 1:dim(M)[1])
	      {for (j in 1:dim(M)[2]) 
	         {cat(M[i,j])
		  if (j != dim(M)[2]) cat(" & ")
		 }
	       if (i != dim(M)[1]) cat("\\\\ \n")
	      }
	    cat("\n\\end{array} \\right] \n\\end{equation}\n")
	   }
	 else print(M, digits=digits)
	 invisible()
	}
	
     printM("F", x$F, digits, latex)
     if(!is.null(x$G)) printM("G", x$G, digits, latex)
     printM("H", x$H, digits, latex)
     if ( is.nonInnov.SS(x)) 
       {printM("Q", x$Q, digits, latex)
        printM("R", x$R, digits, latex)
       }
     else if (is.innov.SS(x))  printM("K", x$K, digits, latex)
     if(!is.null(x$z0)) printM("initial state", x$z0, digits, latex)
     if(!is.null(x$rootP0))
        printM("square root of initial state tracking error", x$rootP0, digits, latex)
     else if(!is.null(x$P0))
        printM("initial state tracking error", x$P0, digits, latex)
     invisible(x)
    }

print.ARMA <- function(x, digits=options()$digits, latex=FALSE, L=TRUE, fuzz=1e-10, ...) 
#  (... further arguments, currently disregarded)
   {arrayc <- function(n)
	{cat("\\left[ \\begin{array}{"); for(j in 1:n) cat("c"); cat("}\n") }
        
     A <- x$A
     B <- x$B
     C <- x$C
     if(!is.null(x$TREND))
        cat("TREND= ", format(signif(x$TREND,digits=digits)))
     if (latex | L)
       {if (latex) cat("\\begin{equation}")
        cat("\nA(L) =\n")
        if (latex) arrayc(dim(A)[3])
        for(i in 1:dim(A)[2]) 
          {for(j in 1:dim(A)[3]) 
             {cat(format(signif(A[1,i,j],digits=digits)))
              if(dim(A)[1] > 1) for(l in 2:dim(A)[1]) 
                 if (abs(A[l,i,j]) > fuzz)
                  {if(1==sign(A[l,i,j])) cat("+")
                   cat(format(signif(A[l,i,j],digits=digits)))
                   if (latex) cat("L^{",l-1,"}", sep="") else cat("L",l-1,sep="")
                  }
               if (latex & j != dim(A)[3]) cat(" & ") else cat("    ")
             }
           if (latex) cat("\\\\ \n") else cat("\n")
         }
        if (latex) cat("\\end{array} \\right] \n\\end{equation}\n\\begin{equation}")
        cat("\nB(L) =\n")
        if (latex) arrayc(dim(B)[3]) 
        for(i in 1:dim(B)[2]) 
          {for(j in 1:dim(B)[3]) 
             {cat(signif(B[1,i,j],digits=digits))
              if (2 <= dim(B)[1]) for(l in 2:dim(B)[1]) 
                 if (abs(B[l,i,j]) > fuzz)
                  {if(1==sign(B[l,i,j])) cat("+")
                   cat(signif(B[l,i,j],digits=digits))
                   if (latex) cat("L^{",l-1,"}", sep="") else cat("L",l-1,sep="")
                  }
              if (latex & j != dim(B)[3]) cat(" & ") else cat("    ")
             }
           if (latex) cat("\\\\ \n") else cat("\n")
         }
         if (latex) cat("\\end{array} \\right] \n\\end{equation}\n")
         if(!is.null(x$C)) 
          {if (latex) cat("\\begin{equation}")
           cat("\nC(L) =\n")
           if (latex) arrayc(dim(C)[3])
           for(i in 1:dim(C)[2]) 
           {for(j in 1:dim(C)[3]) 
             {cat(signif(C[1,i,j],digits=digits))
              if (2<=dim(C)[1]) for(l in 2:dim(C)[1]) 
                 if (abs(C[l,i,j]) > fuzz)
                  {if(1==sign(C[l,i,j])) cat("+")
                   cat(signif(C[l,i,j],digits=digits))
                   if (latex) cat("L^{",l-1,"}", sep="") else cat("L",l-1,sep="")
                  }
              if (latex & j != dim(C)[3]) cat(" & ") else cat("    ")
             }
             if (latex) cat("\\\\ \n") else cat("\n")
           } 
           if (latex) cat("\\end{array} \\right] \n\\end{equation}\n")
	  }
       }
     else
       {for(l in 1:dim(A)[1]) {cat("\nA(L=",l-1,")\n");print(A[l,,],digits=digits)}
        for(l in 1:dim(B)[1]) {cat("\nB(L=",l-1,")\n");print(B[l,,],digits=digits)}
        if(!is.null(x$C))
          for(l in 1:dim(C)[1]) {cat("\nC(L=",l-1,")\n");print(C[l,,],digits=digits)}
       }
     invisible(x)
} 
 
summary.TSestModel <- function(object, ...)
  {#  (... further arguments, currently disregarded)
   residual <- residuals(object)
   sampleT <- nrow(residual)
   p <- ncol(residual)	
   #Om <- t(residual) %*% residual/sampleT
   Om <- crossprod(residual) /sampleT
   rmse <- matrix( diag(Om)^.5 ,1,p)
   dimnames(rmse) <- list(c("RMSE"), seriesNamesOutput(object))

   classed(list(  # summary.TSestModel constructor
     estimates=list(
        l=object$estimates$like[1],
        rmse=rmse,
        sampleT=sampleT,
        converged= object$converged,
        nlmin.results= (!is.null(object$nlmin.results)),
        conv.type= object$nlmin.results$conv.type,
        filter= (!is.null(object$filter)),
        smooth= (!is.null(object$smooth)) ),
     model=summary(object$model)), "summary.TSestModel")
  }


print.summary.TSestModel <- function(x, digits=options()$digits, ...)
  {#  (... further arguments, currently disregarded)
   cat("neg. log likelihood =",x$estimates$l)
   cat("    sample length ="     ,x$estimates$sampleT, "\n")
   print(x$estimates$rmse)
   if (!is.null(x$estimates$converged))
                           cat("convergence: ",x$estimates$converged,"\n")
   if (x$estimates$nlmin.results)
                           cat("convergence type: ", x$estimates$conv.type,"\n")
   if (x$estimates$filter) cat("Includes  filter  estimates.\n")
   if (x$estimates$smooth) cat("Includes smoother estimates.\n")
   print(x$model, digits=digits)
   invisible(x)
  }



summary.SS <- function(object, ...)
  {#  (... further arguments, currently disregarded)
   m <- nseriesInput(object)
   p <- nseriesOutput(object)
   n <- nstates(object)
   classed(list(  # summary.SS constructor
         description=object$description,
         input.series=seriesNamesInput(object),
         output.series=seriesNamesOutput(object),
         innov=is.innov.SS(object),
         m=m,
         p=p,
         n=n,
         P=n * (m+2*p),  #assumes full rank noise
         P.actual = length(coef(object)),
         P.IC = sum(object$location %in% c("z", "P", "r P")),
	 constants=length(object$const),
         ICs=(!is.null(object$z0)),
         init.track=(!is.null(object$P0)) | (!is.null(object$rootP0))
	 ), "summary.SS")
  }

print.summary.SS <- function(x, digits=options()$digits, ...)
    {#  (... further arguments, currently disregarded)
     if (x$innov) cat("innovations form ")
     cat("state space: ")
     cat(x$description,"\n")
     cat("inputs : ", x$input.series, "\n")
     cat("outputs: ", x$output.series, "\n")
     cat("   input  dimension = ", x$m)
     cat("   state  dimension = ",x$n)
     cat("   output dimension = ",x$p,"\n")
     cat("   theoretical parameter space dimension = ",x$P,"\n")
     cat("  ",x$P.actual, " actual parameters")
     if (0 != x$P.IC) cat(" (of which ",x$P.IC, " are ICs)")
     cat("   ",x$constants," non-zero constants\n")
     if (x$ICs)        cat("   Initial values specified.\n")
     else              cat("   Initial values not specified.\n")
     if (!x$innov)
       {if (x$init.track) cat("   Initial tracking error specified.\n")
        else              cat("   Initial tracking error not specified.\n")
       }
     invisible(x)
    }

summary.ARMA <- function(object, ...)
  {#  (... further arguments, currently disregarded)
   m <- nseriesInput(object)
   p <- nseriesOutput(object)
   classed(list(  # summary.ARMA constructor
         description=object$description,
         input.series=seriesNamesInput(object),
         output.series=seriesNamesOutput(object),

         a=dim(object$A)[1]-1,
         b=dim(object$B)[1]-1,
         c=dim(object$C)[1]-1, 
         m=m,
         p=p,
         P.actual=length(coef(object)),
         constants=length(object$const),
         trend=(!is.null(object$TREND)) ), "summary.ARMA")
}
 
print.summary.ARMA <- function(x, digits=options()$digits, ...)
    {#  (... further arguments, currently disregarded)
     cat("ARMA: ")
     cat(x$description,"\n")
     cat("inputs : ", x$input.series, "\n")
     cat("outputs: ", x$output.series, "\n")
     cat("     input  dimension = ", x$m)
     cat("     output dimension = ", x$p,"\n")
     cat("     order A = ", x$a)
     cat("     order B = ", x$b)
     cat("     order C = ", x$c,"\n") 
     cat("     ",x$P.actual, " actual parameters")
     cat("     ",x$constants," non-zero constants\n")
     if(x$trend) cat("     trend estimated.\n")
     else        cat("     trend not estimated.\n")
     invisible(x)
    }
 
tfplot.TSestModel <- function(x, ..., 
    tf=NULL, start=tfstart(tf), end=tfend(tf), 
    select.inputs=NULL, select.outputs=NULL,
    Title=NULL, xlab=NULL, ylab=NULL, 
    graphs.per.page=5, mar=par()$mar, reset.screen=TRUE) {

  # plot one-step ahead estimates and actual data.
  # ... is a list of models of class TSestModel.
  model <- x
  if (is.null(Title))
     Title <- "One step ahead predictions (dotted) and actual data (solid)"
  p<- nseriesOutput(model)
  if (is.null(select.outputs)) select.outputs <-1:p
  if (all(0==select.outputs)) select.outputs <- NULL
  Ngraphs <- length(select.outputs)
  if (!is.null(select.inputs)) if (all(0==select.inputs)) select.inputs <- NULL
  m<- nseriesInput(model)
  if (is.null(m)) m <-0
  else Ngraphs <- Ngraphs+length(select.inputs)  # NULL is zero
  Ngraphs <- min(Ngraphs, graphs.per.page)
  if(reset.screen) {
     old.par <- par(mfcol = c(Ngraphs, 1), mar = mar, no.readonly=TRUE)
     on.exit(par(old.par)) }
  names <-seriesNamesOutput(model)
  if (m!=0) names <-c(seriesNamesInput(model), names)
  if (is.null(names)) names <- rep(" ", m+p)
  if (1 == length(xlab)) xlab <- rep(xlab, m+p)

  if (m!=0) 
    {for (i in select.inputs) 
      {z <-NULL 
       for (model in append(list(x),list(...)) )
           z<-tbind(z,inputData(model, series=i))
       tframe(z) <-tframe(inputData(model))
       tfOnePlot(z,start=start,end=end, xlab=xlab[i], ylab=names[i]) # tsplot
       if(!is.null(Title) && (i==1) && (is.null(options()$PlotTitles)
                || options()$PlotTitles)) title(main = Title)
    } }
  for (i in select.outputs ) 
    {#z <-c(outputData(model, series=i),
     #      rep(NA,Tobs(model$estimates$pred) - TobsOutput(model)))
     z <- outputData(model, series=i)
     for (model in append(list(x),list(...))){
         if (! is.TSestModel(model)) 
	    stop("Argument in ... to be plotted is not a TSestModel object.")
         #z <- cbind(z,model$estimates$pred[,i,drop=FALSE])}
         z <- tbind(z,selectSeries(model$estimates$pred, series=i))}
     #tframe(z) <- tframe(outputData(model))
     tfOnePlot(z,start=start,end=end, xlab=xlab[m+i], ylab=names[m+i]) # tsplot
     if(!is.null(Title) && (i==1) && (is.null(options()$PlotTitles)
                || options()$PlotTitles)) title(main = Title)
    }
  invisible()
  }
    


testEqual.TSestModel <- function(obj1, obj2, fuzz=0) # this could be better
  { testEqual.TSmodel( obj1$model, obj2$model, fuzz=fuzz) &
        testEqual.TSdata(obj1$data, obj2$data, fuzz=fuzz)
  }

testEqual.TSmodel <- function(obj1, obj2, fuzz=0)
{# return T if models are identical (excluding description)
  r       <- all(class(obj1) == class(obj2))
  if (r) r <-length(coef(obj1)) == length(coef(obj2))
  if (r) r <-all(fuzz >= abs(coef(obj1)   -  coef(obj2)))
  if (r) r <-length(obj1$location) == length(obj2$location)
  if (r) r <-all(obj1$location  ==     obj2$location)
  if (r) r <-length(obj1$i) == length(obj2$i)
  if (r) r <-all(obj1$i  ==     obj2$i)
  if (r) r <-length(obj1$j) == length(obj2$j)
  if (r) r <-all(obj1$j  ==     obj2$j)
  if (r) r <-length(obj1$const) == length(obj2$const)
  if (r) r <-all(obj1$const  ==     obj2$const)
  if (r) r <-length(obj1$const.location) == length(obj2$const.location)
  if (r) r <-all(obj1$const.location  ==     obj2$const.location)
  if (r) r <-length(obj1$const.i) == length(obj2$const.i)
  if (r) r <-all(obj1$const.i  ==     obj2$const.i)
  if (r) r <-length(obj1$const.j) == length(obj2$const.j)
  if (r) r <-all(obj1$const.j  ==     obj2$const.j)
  if (r)
    {if (is.ARMA(obj1))    r <-testEqual.ARMA(obj1,obj2, fuzz=fuzz)
     else if (is.SS(obj1)) r <-testEqual.SS(obj1,obj2, fuzz=fuzz)
    }
  r
}

testEqual.ARMA <- function(obj1, obj2, fuzz=0)
{    r <-length(obj1$l) == length(obj2$l)
     if (r) r <-all(obj1$l  ==     obj2$l)
     if (r) r <-length(obj1$const.l) == length(obj2$const.l)
     if (r) r <-all(obj1$const.l  == obj2$const.l)
     if (r) r <-length(obj1$A) == length(obj2$A)
     if (r) r <-all(fuzz >= abs(obj1$A   -     obj2$A))
     if (r) r <-length(obj1$B) == length(obj2$B)
     if (r) r <-all(fuzz >= abs(obj1$B   -     obj2$B))
     if (r) 
           {if (is.null(obj1$C)) r <- is.null(obj2$C)
             else
               {if (r) r <-length(obj1$C) == length(obj2$C)
                if (r) r <-all(fuzz >= abs(obj1$C   -     obj2$C))
           }  }
  r
}

testEqual.SS <- function(obj1,obj2, fuzz=0)
{    r <-length(obj1$F) == length(obj2$F)
     if (r) r <-all(fuzz >= abs(obj1$F   -     obj2$F))
     if (r) 
         {if(is.null(obj1$G)) r <- is.null(obj2$G)
          else
            {if (r) r <-length(obj1$G) == length(obj2$G)
             if (r) r <-all(fuzz >= abs(obj1$G   -     obj2$G))
         }  }
     if (r) r <-length(obj1$H) == length(obj2$H)
     if (r) r <-all(fuzz >= abs(obj1$H   -     obj2$H))
     if (is.innov.SS(obj1))
       {if (r) r <-length(obj1$K) == length(obj2$K)
        if (r) r <-all(fuzz >= abs(obj1$K   -     obj2$K))
       }
     else
       {if (r) r <-length(obj1$Q) == length(obj2$Q)
        if (r & (0 != length(obj2$Q)) ) r <-all(fuzz >= abs(obj1$Q - obj2$Q))
        if (r) r <-length(obj1$R) == length(obj2$R)
        if (r & (0 != length(obj2$R)) ) r <-all(fuzz >= abs(obj1$R - obj2$R))
       }
     if (r) if(is.null(obj1$z0)) r <- is.null(obj2$z0)
     else
       {if (r) r <-length(obj1$z0) == length(obj2$z0)
        if (r) r <-all(fuzz >= abs(obj1$z0   -     obj2$z0))
       }
     if (r) if(is.null(obj1$P0)) r <- is.null(obj2$P0)
     else
       {if (r) r <-length(obj1$P0) == length(obj2$P0)
        if (r) r <-all(fuzz >= abs(obj1$P0   -     obj2$P0))
       }
     if (r) if(is.null(obj1$rootP0)) r <- is.null(obj2$rootP0)
     else
       {if (r) r <-length(obj1$rootP0) == length(obj2$rootP0)
        if (r) r <-all(fuzz >= abs(obj1$rootP0   -     obj2$rootP0))
       }
  r
}


McMillanDegree <- function(model,  ...) UseMethod("McMillanDegree")

McMillanDegree.TSestModel <- function(model, ...)
 {McMillanDegree(TSmodel(model),...) }

McMillanDegree.ARMA <- function(model, fuzz=1e-4, verbose=TRUE, warn=TRUE, ...)
 {#  (... further arguments, currently disregarded)
    z  <- roots(model, warn=warn)
    gross <- length(z)
    zz <- outer(z,z,FUN="-")
    distinct <-sum(!apply((outer(1:length(z),1:length(z),FUN="<") & (Mod(zz) <fuzz)),2,any))
    deg <- list(gross=gross,distinct=distinct)
    if (verbose)
      {cat("Assuming the model is left prime:\n")
       cat("   Without distinguishing distinct roots the degree det(A(L)) =",deg$gross,"\n")
       cat("   Distinguishing  distinct  roots       the degree det(A(L)) =",deg$distinct,"\n")
       if(!is.null(model$TREND))
          {cat("The trend adds unit roots which are added to the degree. Multiple")
           cat("unit roots are not considered distinct (but probably should be).\n")
          }
      }
    invisible(deg)
 }

McMillanDegree.SS <- function(model, fuzz=1e-4, ...)
 {#  (... further arguments, currently disregarded)
   cat("state dimension = ", nstates(model),"/n")
   invisible()
 }

stability <- function(obj, fuzz=1e-4, eps=1e-15, digits=8, verbose=TRUE) UseMethod("stability")

stability.TSestModel <- function(obj, fuzz=1e-4, eps=1e-15, digits=8, verbose=TRUE)
   {stability(TSmodel(obj), fuzz=fuzz, eps=eps, digits=digits, verbose=verbose)}

stability.TSmodel <- function(obj, fuzz=1e-4, eps=1e-15, digits=8, verbose=TRUE)
  {stability(roots(obj, fuzz=fuzz, randomize=FALSE),
             fuzz=fuzz, eps=eps, digits=digits, verbose=verbose)}

stability.roots <- function(obj, fuzz=1e-4, eps=1e-15, digits=8, verbose=TRUE) 
   {#obj <- roots(model, fuzz=fuzz, randomize=FALSE)
    m <- Mod(obj)
    s <- if (all(m < (1.0 - eps))) TRUE else  FALSE
    if (verbose)
      {if (s) cat("The system is stable.\n")
       else   cat("The system is NOT stable.\n")
      }
    r <- cbind(obj, m)
    dimnames(r) <- list(NULL, c("Eigenvalues of F","moduli"))
    attr(s, "roots") <- r
    s
   }

stability.ARMA <- function(obj, fuzz=1e-4, eps=1e-15, digits=8, verbose=TRUE) 
   {z <- roots(obj, fuzz=fuzz, randomize=FALSE)
    if (is.null(z))
      {warning("The model has only a zero root.")
       return(TRUE)}
    m <- Mod(z)
    s <- if (all(m < (1.0 - eps))) TRUE else  FALSE
    if (verbose) {
        cat("Distinct roots of det(A(L)) and moduli are:\n")
        print(cbind(1/z, Mod(1/z)), digits = digits)
        if (!is.null(obj$TREND)) 
            cat("Trend not taken into account: ")
        if (s) 
            cat("The system is stable.\n")
        else cat("The system is NOT stable.\n")
      }
    r <- cbind(z, m)
    dimnames(r) <- list(NULL, c("Inverse of distinct roots of det(A(L))","moduli"))
    attr(s, "roots") <- r
    s
   }


roots <- function(obj, ...) UseMethod("roots")

roots.TSestModel <- function(obj, ...){roots(TSmodel(obj), ...)}

roots.SS <- function(obj, fuzz=0, randomize=FALSE, ...) 
{#  (... further arguments, currently disregarded)
    z <- eigen(obj$F, symmetric=FALSE, only.values=TRUE)$values
    if (randomize) if (sample(c(TRUE,FALSE),1)) z <- Conj(z)
      #this prevents + - ordering of complex roots (for Monte Carlo evaluations)
    classed(z,"roots")  # constructor (roots.SS)
}


roots.ARMA <- function(obj, fuzz=0, randomize=FALSE, warn=TRUE, by.poly=FALSE, ...) 
{#  (... further arguments, currently disregarded)
    if(1 == dim(obj$A)[1]) {
        warning("An ARMA model with no lags in the AR term has no roots.")
	return(NULL)
	}
    if(by.poly) z <- 1/polyrootDet(obj$A)
    else        z <- roots(toSS(obj))
    if (fuzz!=0)
      {zz <- outer(z,z,FUN="-")  # find distinct roots within fuzz
       z <- z[ !apply((outer(1:length(z),1:length(z),FUN="<")
                 & (Mod(zz) <fuzz)),2,any)]
      }
    # add unit roots for TREND elements.
    if (!is.null(obj$TREND))
      {#z <- c(rep(1,sum(0!=obj$TREND)), z)this depends on nature of TREND
       if (warn)
         warning("Unit roots may need to be added for non-zero trend elements.")
      }
    if (randomize) if (sample(c(TRUE,FALSE),1)) z <- Conj(z)
      #this prevents + - ordering of complex roots (for Monte Carlo evaluations)
    classed(z,"roots")  # constructor (roots.ARMA)
}


plot.roots <- function(x, pch='*', fuzz=0, ...)
{#  (... further arguments, currently disregarded)
 if (is.TSmodel(x))    x <- roots(x, fuzz=fuzz)
 if (is.TSestModel(x)) x <- roots(x, fuzz=fuzz)
 i <- 2*pi*(1:1000)/1000
 if (max(Mod(x)) <1 )
   {plot(sin(i),cos(i),pch='.', xlab='Im', ylab='Re')
    points(Re(x),Im(x), pch=pch)
   }
 else
   {#plot.default(Re(x),Im(x), pch=pch, xlab='Im', ylab='Re')
    plot(Re(x),Im(x), pch=pch, xlab='Im', ylab='Re')
    points(sin(i),cos(i),pch='.')
   }
 points(0,0,pch='+')
 invisible(x)
}


addPlotRoots <- function(v, pch='*', fuzz=0)
{if (is.TSmodel(v))    v <- roots(v, fuzz=fuzz)
 if (is.TSestModel(v)) v <- roots(v, fuzz=fuzz)
 points(Re(v),Im(v), pch=pch)
 invisible(v)
}


observability <- function(model)  UseMethod("observability")

observability.TSestModel <- function(model){observability(TSmodel(model)) }

observability.SS <- function(model){ 
  FF<-    model$F
  O <-    model$H
  HFn <- O
  for (n in 1:dim(FF)[1])  {
    HFn <- HFn %*% FF
    O <- rbind(O,HFn)
    }
  La.svd(O)$d
  }

observability.ARMA <- function(model){
  cat("not applicable to ARMA models\n")
  invisible()
  }


reachability <- function(model) UseMethod("reachability") 

reachability.TSestModel <- function(model){reachability(TSmodel(model))}

reachability.SS <- function(model){
 FF <-    model$F
 C  <-    model$G
 if (!is.null(C))
   {FnG <- C
    for (n in 1:dim(FF)[1])  
      {FnG <- FF %*% FnG
       C <- cbind(C,FnG)
      }
    cat("Singular values of reachability matrix for input: ",La.svd(C)$d)
   }
 if (is.innov.SS(model)) C <- model$K
 else      
  {C <- model$R
   if(dim(C)[1]==1)
     {if (any(C==0))
         {cat("State noise matrix is singular. All states are NOT excited!\n")
          return(C)
         }
      C <- 1/C
     }
   else
     {v<-La.svd(C)
      if (any(v$d==0))
         {cat("State noise matrix is singular. All states are NOT excited!\n")
          return(v$d)
         }
  #   C <-v$v %*% diag(1/v$d) %*% t(v$u) following is equivalent
  #   C <-v$v %*% (t(v$u) * 1/v$d) with svd (vs La.svd)
      C <-(v$u %*% v$vt ) * 1/v$d 
     }
   C <- model$Q %*% C
  }
 FnK <- C
 for (n in 1:dim(FF)[1])  
   {FnK <- FF %*% FnK
    C <- cbind(C,FnK)
   }
 d <- La.svd(C)$d
 cat("Singular values of reachability matrix for noise: ",d,"\n")
 invisible(d)
 }

reachability.ARMA <- function(model){ 
  cat("not applicable to ARMA models\n")
  invisible()
  }


checkBalance <- function(model) UseMethod("checkBalance") 

checkBalance.TSestModel <- function(model){checkBalance(TSmodel(model))}

checkBalance.SS <- function(model){ 
  FF<-    model$F
  O <-    model$H
  HFn <- O
  for (n in 1:dim(FF)[1])  {
    HFn <- HFn %*% FF
    O <- rbind(O,HFn)
    }
  #O <- t(O) %*% O 
  O <- crossprod(O) # observability gramian
  C <-    cbind(model$G,model$K)
  FnG <- C
  for (n in 1:dim(FF)[1])  {
    FnG <- FF %*% FnG
    C <- cbind(C,FnG)
    }
  C <- C %*% t(C) # controllability gramian
  difference <- O-C
  #cat("observability gramian minus controllability gramian:\n")
  #print(difference)
  cat("maximum absolute difference (O-C): ", max(abs(difference)),"\n")
  cat("maximum off-diagonal element of C: ", max(abs(C-diag(diag(C)))),"\n")
  cat("maximum off-diagonal element of O: ", max(abs(O-diag(diag(O)))),"\n")
  invisible()
  }

checkBalance.ARMA <- function(model){ 
  cat("not applicable to ARMA models\n")
  invisible()
  }


checkBalanceMittnik <- function(model) UseMethod("checkBalanceMittnik")

checkBalanceMittnik.TSestModel <- function(model)
   checkBalanceMittnik(TSmodel(model))

checkBalanceMittnik.SS <- function(model){ 
  FF  <-    model$F - model$K %*% model$H
  O   <-    model$H
  HFn <- O
  for (n in 1:dim(FF)[1])  {
    HFn <- HFn %*% FF
    O <- rbind(O,HFn)
    }
  #O <- t(O) %*% O 
  O <- crossprod(O)# observability gramian
  C <-    cbind(model$G,model$K)
  FnG <- C
  for (n in 1:dim(FF)[1])  {
    FnG <- FF %*% FnG
    C <- cbind(C,FnG)
    }
  C <- C %*% t(C) # controllability gramian
  difference <- O-C
  #cat("observability gramian minus controllability gramian:\n")
  #print(difference)
  cat("maximum absolute difference (O-C): ", max(abs(difference)),"\n")
  cat("maximum off-diagonal element of C: ", max(abs(C-diag(diag(C)))),"\n")
  cat("maximum off-diagonal element of O: ", max(abs(O-diag(diag(O)))),"\n")
  invisible()
  }

checkBalanceMittnik.ARMA <- function(model){ 
  cat("not applicable to ARMA models\n")
  invisible()
  }


############################################################

#     functions for model conversion   <<<<<<<<<<

############################################################

toSS <- function(model, ...) UseMethod("toSS")

toSS.TSestModel <- function(model, ...) 
	l(toSS(TSmodel(model), ...),TSdata(model))

toSS.SS <- function(model, ...) {model}
 #  (... further arguments, currently disregarded)
 
toSS.ARMA <- function(model,...){
    # convert an ARMA (or VAR) to a SS (innovations) representation
    if (is.null(model$A)) a<-0
    else a <- dim(model$A)[1] - 1  #order of polynomial arrays
    if (is.null(model$B)) b<-0
    else b <- dim(model$B)[1] - 1
    if (is.null(model$C)) cc<-0
    else cc<- dim(model$C)[1] - 1
    if ((b<=a) & (cc<=(a-1))) model <- toSSaugment(model)
    else  model <-toSSnested(model,...) #  (otherwise best working method) 
                  # A better approach would be an algorithm like Guidorzi's. 
    model
    }


toSSnested <- function(model, ...) UseMethod("toSSnested")

toSSnested.TSestModel <- function(model, ...) toSSnested(TSmodel(model), ...)

toSSnested.SS <- function(model, n=NULL, Aoki=FALSE, ...){
  #  (... further arguments, currently disregarded)
  # convert to a nested-balanced state space model by svd  a la Mittnik (or Aoki)
  if (is.null(n)) n <-ncol(model$F)  
  if (Aoki) stop("Aoki balancing no longer supported.") #return(Aoki.balance(model, n=n))
  else      return(balanceMittnik(model, n=n)) 
  }

toSSnested.ARMA <- function(model, n=NULL, Aoki=FALSE, ...){
  #  (... further arguments, currently disregarded)
  # convert to a nested-balanced state space model by svd  a la Mittnik (or Aoki)
  if (is.null(n)) n <- McMillanDegree(model)$distinct
  if (Aoki) stop("Aoki balancing no longer supported.") #return(Aoki.balance(model, n=n))
  else      return(balanceMittnik(model, n=n)) 
  }


toSSaugment <- function(model, ...) UseMethod("toSSaugment")

toSSaugment.TSestModel <- function(model, ...)
   l(toSSaugment(TSmodel(model), ...), TSdata(model))


toSSaugment.ARMA <- function(model, fuzz=1e-14, ...) {
  #  (... further arguments, currently disregarded)
  # convert by augmentation - state dimension may not be minimal
  # First sets A[1,,] = B[1,,] = I if that is not already the case.
   A <- model$A
   B <- model$B
   C <- model$C 
   if (fuzz  < max(abs(A[1,,]-diag(1,dim(A)[2]) )) )
      {A0.inv <- solve(A[1,,])
       A <-  polyprod(A0.inv,A)
       B <-  polyprod(A0.inv, B)
       if (!is.null(C)) C <- polyprod(A0.inv, C)
#       if (!is.null(TREND)) TREND <- t(A0.inv %*% t(TREND))
       }
   if (fuzz  < max(abs(B[1,,]-diag(1,dim(B)[2]) )) )
          B<- polyprod(solve(B[1,,]), B)
   if (is.null(model$A)) a<-0
   else a <- dim(model$A)[1] - 1  #order of polynomial arrays
   if (is.null(model$B)) b<-0
   else b <- dim(model$B)[1] - 1
   if (is.null(model$C)) cc<-0
   else cc<- dim(model$C)[1] - 1
   p <- nseriesOutput(model)          # Dim of endoenous Variables.
   m <-  nseriesInput(model)          # Dim of exogenous Variables.
   if (b>a) stop("The MA order cannot exceed the AR order to convert with state augmentation.")
   if (cc>(a-1)) stop(
      " The order of the input polynomial cannot exceed the AR order -1 to convert with state augmentation.")   
  #make three parameters A,B and C have convenient order by adding 0's.
   k <- 1 + a  
#  if (b != 0)
    {BB <- array(0,c(k,dim(B)[2:3]))
     BB[1:(b+1),,] <- B
    }
   if (m!=0) 
     {CC <- array(0,c(k,dim(C)[2:3]))  
      CC[1:(cc+1),,] <- C
     }
   FF <- matrix(NA,a*p,p)
   for (i in 1:a) FF[(1+p*(i-1)):(p*i),] <- -A[a-i+2,,]
   if(a>1) FF<-cbind(rbind(matrix(0,p,(a-1)*p),diag(1,(a-1)*p)),FF)
   if (m == 0) G <-NULL
   else
     {G <- matrix(NA,a*p,m)
      for (i in 1:a) G[(1+p*(i-1)):(p*i),] <- CC[a-i+1,,] 
     }
   H <- diag(1,p)
   if(a>1) H <- cbind( matrix(0,p,(p*(a-1))),H)
   K <- matrix(NA,a*p,p)
   for (i in 1:a) K[(1+p*(i-1)):(p*i),] <- -A[a-i+2,,]+BB[a-i+2,,]
   z0 <-NULL
   if(!is.null(model$TREND))   #add a constant state which feeds into the states
     {FF<-rbind(cbind(FF,0),0) # identified with outputs (through H).
      n <-dim(FF)[1]
      FF[n,n] <-1 
      if (p != length(model$TREND)) stop("This fails for matrix TREND.")
      FF[n-p:1,n] <- model$TREND
      z0 <- rep(0,n)
      z0[n] <-1
      H<-cbind(H,0)
      if (m!=0) G <- rbind(G,0)
      K<- rbind(K,0)
     }                     
   descr<-c(model$description,
            " Converted to state space by state augmentation.")
   SS(F.=FF,G=G,H=H,K=K,z0=z0,description=descr,
         input.names= seriesNamesInput(model),
        output.names=seriesNamesOutput(model))        
 }


gmap <- function(g, model) 
{# convert to an equivalent representation using a given matrix
 if (is.TSestModel(model)) model <- TSmodel(model)
 if  (!is.TSmodel(model)) stop("gmap expecting a TSmodel.")
 if ( is.SS(model))# transform State space model by g in GL(n)
  {n <- dim(model$F)[1]
   if (!is.matrix(g)) g<-diag(g,n) # if g is not a matrix make it into a diagonal matrix.
   if ((n !=dim(g)[1]) | (n !=dim(g)[2]) )
      stop("g must be a square matrix of dimensions equal the model state (or a scalar).")
   inv.g <- solve(g)
   model$F <-inv.g%*%model$F%*% g
   if (!is.null(model$G)) model$G <-inv.g %*%model$G
   model$H <-model$H %*% g
   if (!is.null(model$z0)) model$z0 <-c(inv.g %*%model$z0)
   if (is.innov.SS(model)) model$K <-inv.g %*% model$K
   if (is.nonInnov.SS(model)) 
      {model$Q <-inv.g %*% model$Q
       model$R <-model$R
      }       
 }
 if ( is.ARMA(model))
       {# if g is not a matrix make it into a diagonal matrix.
        if (! is.matrix(g)) g<- diag(g,dim(model$A)[2]) 
	for(l in 1:dim(model$A)[1]) model$A[l,  ,  ] <- g %*% model$A[l, ,]	
	for(l in 1:dim(model$B)[1]) model$B[l,  ,  ] <- g %*% model$B[l, ,]
	for(l in 1:dim(model$C)[1]) model$C[l,  ,  ] <- g %*% model$C[l, ,]
	if(!is.null(model$TREND))  model$TREND <- t(g %*% t(model$TREND))
       }
 setTSmodelParameters(model)
}


#findg <- function(model1,model2, minf=nlmin){ 
#  if (is.TSestModel(model1)[1]) model1 <- TSmodel(model1)
#  if   (!is.TSmodel(model1)) stop("findg expecting a TSmodel.")
#  if (is.TSestModel(model2)[1]) model2 <- TSmodel(model2)
#  if  (!is.TSmodel(model2 )) stop("findg expecting a TSmodel.")
#
#  if ( !is.SS(model1)| !is.SS(model2)) 
#      stop("findg only works for state space models")
#   n <- dim(model1$F)[1]
#   if ((n!= dim(model2$F)[1])
#     |(dim(model1$G)[2] != dim(model2$G)[2])
#     |(dim(model1$H)[1] != dim(model2$H)[1]))
#      stop("models must have the same dimensions for findg.")
#   para<- c(diag(1,n))
#   zzz.model1 <<- model1	  # This could be done with assign(frame=1  ??)
#   zzz.model2 <<- model2
#   zzz.n <<-n
#   func <- function(para){
#      gmodel1<- gmap(matrix(para,zzz.n,zzz.n),zzz.model1)
#      error <- 	gmodel1$F-zzz.model2$F
#      error <-c(error,(gmodel1$G-zzz.model2$G))
#      error <-c(error,(gmodel1$H-zzz.model2$H))
#      error <-c(error,(gmodel1$K-zzz.model2$K))
#      error <-c(error,(gmodel1$R-zzz.model2$R))
#      sum(error^2)
#   }
#   para <-minf(func,para)
#   rm(zzz.model1,zzz.model2,zzz.n)
#   matrix(para[[1]],n,n)
#   }


fixConstants <- function(model, fuzz=1e-5, constants=NULL){ 
  if (is.TSestModel(model)) model <- TSmodel(model)
  if  (!is.TSmodel(model)) stop("fixConstants expecting a TSmodel.")
  if (is.null(constants))
    {p <-abs(coef(model) - 1.0) < fuzz
     model$const <- c(model$const,rep(1.0,sum(p)))
     model$const.location <- c(model$const.location,model$location[p])
     model$const.i <- c(model$const.i,model$i[p])
     model$const.j <- c(model$const.j,model$j[p])
     if(is.ARMA(model)) model$const.l <- c(model$const.l,model$l[p])
     p <- (!p) & (abs(coef(model)) > fuzz) 
     model$coefficients <- coef(model)[p]
     model$location <- model$location[p]
     model$i <- model$i[p]
     model$j <- model$j[p]
     if(is.ARMA(model)) model$l <- model$l[p]
     return(setArrays(model))
    }
  else return(setTSmodelParameters(model,constants=constants))
  }


toSSinnov <- function(model, ...){
 #  (... further arguments, currently disregarded)
 data <- NULL
 if (is.TSestModel(model)) 
    {data <- TSdata(model)
     model <- TSmodel(model)
    }
 if  (!is.TSmodel(model)) stop("toSSinnov expecting a TSmodel.")
 if (!is.SS(model))  model <- toSS(model)
 if ( is.nonInnov.SS(model)) 
   {warning("this is not exactly correct.")
    PH  <-  model$Q %*% t(model$Q) %*% t(model$H)# use QQ' as state tracking error ??
    ft    <- (model$H %*% PH) + (model$R %*% t(model$R))        
    ft    <-  (ft + t(ft))/2   # force ft to be symmetric 
    model$K <-  t(solve(ft,t(model$F %*% PH))) 
    model$R <- NULL
    model$Q <- NULL
   }  
 model <- setTSmodelParameters(classed(model, c("innov","SS","TSmodel"))) # bypass constructor
 if (is.null(data)) model else l(model, data)
 }

toSSOform <- function(model) UseMethod("toSSOform")

toSSOform.TSestModel <- function(model) 
   l(toSSOform(TSmodel(model)), TSdata(model))

toSSOform.TSmodel <- function(model){
 if (!is.SS(model))       model <- toSS(model)
 if (!is.innov.SS(model)) model <- toSSinnov(model)
 n <- dim(model$H)[2]
 p <- nseriesOutput(model)
 if (p >= n) 
   {ginv <-model$H[1:n,]
    g <-solve(ginv)
   }
 else
   {sv   <- La.svd(model$H)
# svd    g    <- sv$v %*% diag(1/sv$d, ncol = length(sv$d))  %*%  t(sv$u) #right inv
# svd    ginv <- sv$u %*% diag( sv$d,  ncol = length(sv$d))  %*%  t(sv$v) #left  inv
    g    <- t(sv$vt) %*% diag(1/sv$d, ncol = length(sv$d))  %*%  t(sv$u) #right inv
    ginv <-   sv$u   %*% diag( sv$d,  ncol = length(sv$d))  %*%    sv$vt #left  inv
    g    <- cbind(g,   matrix(0,n, n-p))  # no good. these need to be full rank
    ginv <- rbind(ginv,matrix(0,n-p, n))  # and still convert H to [ I | 0 ]
    stop("This procedure is not working properly yet.") 
    # have fixed only nxp not yet nxn elements
   }
 model$F <- ginv %*% model$F %*% g
 if(!is.null(model$G)) model$G <- ginv %*% model$G
 model$H <- model$H %*% g
 model$K <- ginv %*% model$K  
 fixConstants(setTSmodelParameters(model))
}


fixF <- function(model){
  if (is.TSestModel(model)) model <- TSmodel(model)
  if  (!is.TSmodel(model)) stop("fixF expecting a TSmodel.")
  if (!is.SS(model))         model <- toSS(model)
  if ( is.nonInnov.SS(model))  model <- toSSinnov(model)
  p <-model$location == "f"
  model$const <- c(model$const, coef(model)[p])
  model$const.location <- c(model$const.location,model$location[p])
  model$const.i <- c(model$const.i,model$i[p])
  model$const.j <- c(model$const.j,model$j[p])
  p <- !p
  model$coefficients <- coef(model)[p]
  model$location <- model$location[p]
  model$i <- model$i[p]
  model$j <- model$j[p]
  cat("Remaining parameters: ",sum(p),"\n")
  if (is.null(model$G)) m <-0
  else  m <- dim(model$G)[2]
  n <- dim(model$F)[1]
  p <- dim(model$H)[1]
  cat("Theoretical parameter space dimension: ",n*(m+2*p),"\n")
  setArrays(model)
}


toSSChol <- function(model, ...) UseMethod("toSSChol")

toSSChol.TSestModel <- function(model, Om=NULL, ...) {
   #  (... further arguments, currently disregarded)
   if(is.null(Om)) Om <-model$estimates$cov
   l(toSSChol(TSmodel(model), Om=Om), TSdata(model))
   }

toSSChol.TSmodel <- function(model, Om=diag(1,nseriesOutput(model)), ...){
 #  (... further arguments, currently disregarded)
 if (!is.SS(model))  model <- toSS(model)
 if (is.innov.SS(model)) 
   {model$R <-t(chol(Om) )  # Om = RR'
    model$Q <- model$K %*% model$R
    model$K <- NULL
   }  
 classed(model, c( "nonInnov","SS","TSmodel" ) )  # bypass constructor
 }

toARMA <- function(model, ...) UseMethod("toARMA")

toARMA.TSestModel <- function(model, ...) 
	l(toARMA(TSmodel(model), ...), TSdata(model))

toARMA.ARMA <- function(model, ...) model

toARMA.SS <- function(model, fuzz=1e-10, ...){
    #  (... further arguments, currently disregarded)
    if (is.nonInnov.SS(model)) model <- toSSinnov(model)
    FF<-model$F
    G <-model$G
    H <-model$H
    K <-model$K
    m <-dim(G)[2]
    if (is.null(m)) m <-0
    n <-dim(FF)[1]
    p <-dim(H)[1]
    poly <-  - characteristicPoly(FF) # N.B. sign change  in Aoki vs Kailath
    if ( n != length(poly))
      stop("There is some problem. The characteristic polynomial length should = state dimension.")
    A <- array(NA,c(1+n,p,p))
    A[1,,] <-diag(1,p)
    for (i in 1:n) A[i+1,,] <- diag(-poly[i],p)
    if(any(is.na(A))) stop("error in calculation of A in toARMA.")

#                                       i
#            Fn [i,,] corresponds to   F    
    if (n > 1)
      {Fn <- array(NA,c(n-1,n,n))
       Fn[1,,] <- FF
      }
    if (n > 2) for (i in 2:(n-1))   Fn[i,,] <- FF %*% Fn[i-1,,]

    
    HFnK <- array(NA, c(n+1, p, p))
    HFnK[1, , ] <- diag(1, p)
    HFnK[2, , ] <- H %*% K
    if (n > 1) for (i in 3:(n+1)) HFnK[i, , ] <- H %*% Fn[i-2, , ] %*% K
    B <- array(NA, c(1+n, p, p))
    B[1, , ] <- diag(1, p)
    for (i in 1:n) 
       {B[i+1, , ] <- HFnK[i+1, , ]
        for (j in 1:i) B[i+1, , ] <- B[i+1, , ] - poly[j] *  HFnK[i+1-j, , ]
       }
    if(any(is.na(B))) stop("error in calculation of B in toARMA.")

    if (m == 0) C <- NULL
    else
      {C <- array(NA,c(n,p,m))   
       HFnG <- array(NA,c(n,p,m)) 
       HFnG[1,,] <- H %*% G
       if (n > 1) for (i in 2:n) HFnG[i,,] <- H %*% Fn[i-1,,] %*% G
       C[1,,] <- HFnG[1,,]
       for (i in 2:n)
         {C[i,,] <- HFnG[i,,]
          for (j in 1:(i-1)) C[i,,] <- C[i,,]-poly[j]*HFnG[i-j,,]
         }
       if(any(is.na(C))) stop("error in calculation of C in toARMA.")
      }
 fixConstants(ARMA(A=A,B=B,C=C, input.names =  seriesNamesInput(model),
                  output.names = seriesNamesOutput(model)), fuzz=fuzz)
}


#######################################################################

#                  Utility functions  <<<<<<<<<<

############################################################

#             functions for generating  statistics  
#        from data and for generating theoretical statistics
#           (i.e.- calculated from model not data)  

############################################################



Riccati <- function(A, B, fuzz=1e-10, iterative=FALSE)
{warning("This procedure has not been tested!")
 if (!iterative) 
  {n <- dim(A)[1]
   Atinv <-solve(t(A))    # symplectic matrix Vaughan (10)(12), R=0
   S   <- rbind(cbind(Atinv     , diag(0,n)),
               cbind(B %*% Atinv, A)) 
   Q <- eigen(S, symmetric=FALSE, only.values=FALSE)
   Q <- Q$vectors[,rev(order(Mod(Q$values)))]   # This may have imaginary parts.
   P <- Re( Q[(n+1):(n+n),1:n] %*% solve(Q[1:n,1:n]) ) #This should not have any significant im parts.
  }
 else
  {P<- diag(0,dim(A)[1])
   i <-0
   repeat    # yuk
     {P <- A%*% P %*% t(A) + B
      i <- i+1
      if (i>1000) break
      if (fuzz > max(abs(P-A %*% P %*% t(A) - B))) break
  }  }
 if (fuzz < max(abs(P-A %*% P %*% t(A) - B)))
      warning("Riccati failed! Result does not solve the Riccati equation!")
 P
}


markovParms <- function(model, blocks=NULL) 
{if (is.TSestModel(model)) model <- TSmodel(model)
 if  (!is.TSmodel(model)) stop("markovParms expecting a TSmodel.")
 if ( is.ARMA(model))
  #  M ={ Mi }={ Ci-1|Bi| -Ai}, i=2,...,k. k=max(a,b,cc). Assumes I=A[1,,]=B[1,,]
  {A <- model$A
   B <- model$B
   C <- model$C
   a <- dim(A)[1] - 1      # order of polynomial arrays
   if (is.na(a)) a <- 0
   b <- dim(B)[1] - 1
   if (is.na(b)) b <- 0
   cc <- dim(C)[1] - 1 
   if (is.na(cc)) cc <- 0
   m <- dim(C)[3]          # Dim of exogenous Variables.
   if (is.null(m))  m <- 0                         
   #make three parameters A,B and C have convenient order by adding 0's.    
   if (is.null(blocks))
      blocks <- 1 + max(2,a,cc,b) 
    # if blocks is not at least 3 the Hankel shift does not work
   if (a != 0) 
    {AA <- array(0,c(blocks,dim(A)[2:3]))
     AA[1:(a+1),,] <- A
    }
   if (b != 0)
    {BB <- array(0,c(blocks,dim(B)[2:3]))
     BB[1:(b+1),,] <- B
    }
   if (m != 0)
    {CC <- array(0,c(blocks,dim(C)[2:3]))  
     CC[1:(cc+1),,] <- C
    }
   M <- NULL
   if (b != 0) cat (" Warning: This has only been developed for SS and VARX models.\n")
   if(!is.null(model$TREND)) cat(" Warning: Markov parameter generation does not account for trends.")
   for(i in 2:blocks)  
       {if (m != 0) M <- cbind(M,CC[(i-1),,]) 
        if (b != 0) M <- cbind(M, BB[i,,])      # constant term ignored (=I)
        if (a != 0) M <- cbind(M,-AA[i,,])      # constant term ignored (=I)
       }
  }   
else if ( is.SS(model))
   {FF<-model$F
    G <-model$G
    H <-model$H
    if (is.innov.SS(model)) K <-model$K
    else               K <- model$Q %*% solve(model$R)
    if (is.null(blocks)) blocks <- 1+dim(FF)[1]
    FF <- FF - K %*% H  # model transformed a la Mittnik so lagged outputs are inputs
    FnGK <- cbind(G,K)
    M <- H %*% FnGK
    i<-0 # M should have at least 2 blocks or Hankel shift does not work.
    stopp <- FALSE
    while((i<=blocks) & !stopp)   #  no. of blocks affect Hankel size
       {FnGK <- FF %*% FnGK
        M <- cbind(M,H %*% FnGK)
        i <-i+1   # count should not be necessary, but insures an end.
        stopp <- (i>3) & ( max(abs(FnGK)) < 1e-15)
       }
  }
else stop("markovParms requires an ARMA or SS model.")
M       
}



############################################################

#     polynomial utility functions   <<<<<<<<<<

############################################################


polyprod <- function(a,b)
{ # product of two polynomials.
# The convention used is by polyvalue and polyroot is constant first,
#  highest order coef. last. The reverse convention could also be used for multiplication.
# This function handles scalar (ie. non-matrix) and matrix polynomials.
# Scalar polynomials are vectors of length 1+the polynomial order.
# Polynomial matrices are defined as 3 dimensional arrays with the last 2
# dimensions as the matrix dimension and the first equal 1+the
# polynomial order.

    pprod <- function(a,b)  # local function, product of non-matrix polys.
       {n <- length(a) +length(b) -1
        if (is.null(a))    return(NA)
        if (is.null(b))    return(NA)
        if (0 ==length(a)) return(NA)
        if (0 ==length(b)) return(NA)
        if (any(is.na(a))) return(NA)
        if (any(is.na(b))) return(NA)
	r <- rep(NA, n)
        z <- outer(a, b) 
        zi <- 1 + outer(1:length(a)-1,1:length(b)-1,"+")
	for(i in 1:n) r[i]<- sum(z[i==zi])
	r
       }
   psum <- function(a,b)  # local function, sum of non-matrix polys.
    {if (length(a) < length(b)) return(c(a,rep(0,length(b)-length(a))) + b)
     else                       return(c(b,rep(0,length(a)-length(b))) + a)
    }
   if (is.vector(b)  && (is.array(a) | is.matrix(a)))
     {z <- b; b <- a; a <- z } # scalar multiplication commutes (even for scalar polynomials)
   if (is.null(a))      r <- NULL  
   else if (is.null(b)) r <- NULL   
   else if (is.vector(a))  
          {if      (is.vector(b)) r <-pprod(a,b)
           else if (is.matrix(b))
              {r <- array(NA,c(length(a),dim(b)))
               for (i in 1:(dim(b)[1])) 
                  for (j in 1:(dim(b)[2]))
                     r[,i,j] <- pprod(a,b[i,j])
              }
           else if (is.array(b))
              {r <- array(NA,c(length(a)+dim(b)[1]-1,dim(b)[2:3]))
               for (i in 1:(dim(b)[2])) 
                  for (j in 1:(dim(b)[3]))
                     r[,i,j] <- pprod(b[,i,j],a)
              }
          }
   else if (is.matrix(a))
          {if (is.matrix(b)) r <- a %*% b
           else if (is.array(b))
              {if (dim(a)[2] != dim(b)[2]) 
                  stop("Matrix polynomial dimensions do not conform.")
               r <- array(0,c(dim(b)[1],dim(a)[1], dim(b)[3]))
               #for (i in 1:(dim(a)[1])) for (j in 1:(dim(b)[3]))
                 #for (k in 1:(dim(a)[2])) r[,i,j] <- psum(r[,i,j], pprod(a[i,k],b[,k,j]))
               for (k in 1:(dim(b)[1]))
                 r[k,,] <- a %*% array(b[k,,],dim(b)[2:3])
                 #array above insures b a matrix, drop=FALSE cannot be used
              }
          }
   else if (is.array(a))
        if (is.matrix(b))
          {if (dim(a)[3] != dim(b)[1])
              stop("Matrix polynomial dimensions do not conform.")
           r <- array(0,c(dim(a)[1],dim(a)[2], dim(b)[2]))
           #for (i in 1:(dim(a)[2])) for (j in 1:(dim(b)[2]))
           #  for (k in 1:(dim(a)[3]))  r[,i,j] <- psum(r[,i,j], pprod(a[,i,k],b[k,j]))
           for (k in 1:(dim(a)[1]))
                 r[k,,] <- array(a[k,,],dim(a)[2:3]) %*% b
                 #array above insures b a matrix, drop=FALSE cannot be used
          }
        else if (is.array(b))
          {if (dim(a)[3] != dim(b)[2]) 
              stop("Matrix polynomial dimensions do not conform.")
           r <- array(0,c(dim(b)[1]+dim(a)[1]-1,dim(a)[2], dim(b)[3]))
           for (i in 1:(dim(a)[2])) 
              for (j in 1:(dim(a)[3]))
                 for (k in 1:(dim(a)[3]))
                   r[,i,j] <- psum(r[,i,j], pprod(a[,i,k],b[,k,j]))
          }
   else stop("polynomial product not defined for these objects")
r
}



polysum <- function(a,b){
   # sum of two polynomials (including polynomial arrays)
    psum <- function(a,b)  # local function for non-matrix polys.
 {if (length(a) < length(b))  {r <- b;  r[1:length(a)] <-r[1:length(a)] + a }
  else  {r <- a; r[1:length(b)] <-r[1:length(b)] + b }
  r 
 }
   if (is.null(a) )     r <- b
   else if (is.null(b)) r <- NULL
   else if (is.vector(a) && is.vector(b)) r <-psum(a,b)
   else if (is.array(a) )
        if (is.array(b))
           if ( all(dim(a)[2:3] == dim(b)[2:3]))
             {if (dim(a)[1] < dim(b)[1])
                    {r <- b
                     r[1:(dim(a)[1]),,] <-r[1:(dim(a)[1]),,] + a
                    }
                 else
                    {r <-  a
                     r[1:(dim(b)[1]),,] <-r[1:(dim(b)[1]),,] + b
                    }
             }
           else stop("polynomial matrix dimensions must agree")
        else if (is.vector(b))
          {r <- array(NA,c(max(length(b),dim(a)[1]),dim(a)[2:3]))
           for (i in 1:(dim(a)[2])) 
              for (j in 1:(dim(a)[3]))
                 r[,i,j] <- psum(a[,i,j],b)
          }
   else if (is.vector(a) && is.array(b))
          {r <- array(NA,c(max(length(a),dim(b)[1]),dim(b)[2:3]))
           for (i in 1:(dim(b)[2])) 
              for (j in 1:(dim(b)[3]))
                 r[,i,j] <- psum(b[,i,j],a)
          }
   else stop("polynomial sum not defined for these objects")
r             
}

polyrootDet <- function(a)
{#roots of the determinant of a.  Note: polydet is slow. There is room for improvement here!
z <- polydet(a)
if (length(z)==1) 
  stop(paste("root cannot be calculated for degree 0 determinant = ", as.character(z)))
polyroot(z)
}

polydet <- function(a)
{# recursive pessimist: Life is just one damned thing before another.
 # attributed to Richard Bird in The Mathematical Intelligencer, Winter 1994.
 #Recursively form the determinant of a polynomial matrix a, where the first
 #  dim stores the coefs. of the polynomial for each a[,i,j].
	n <- dim(a)[2]  
        if (n != dim(a)[3])
           stop( "The determinant is only defined for square polynomial matrices")
	if(1 == n) r <- c(a)
	else 
          {r<- 0   # previously NULL
           for (i in 1:n) 
             {if(!all(0==a[,i,1]))#not nec.but faster for sparse arrays
                   r<- polysum(r,(-1)^(i+1)*
          polyprod(a[,i,1],Recall(a[,(1:n)[i!=(1:n)],2:n,drop=FALSE])))
              r[is.na(r)] <- 0
              if (any(0==r)) 
                   {if (all(r==0)) r <- 0
                    else  r <-r[0 == rev(cumprod(rev(r==0)))] #remove trailing zeros
                    if(0==length(r)) r <- 0
             }     } 
          }
r
}

polyvalue <- function(coef, z)
{# evaluate a polynomial given by coef (constant first) at z
 # could be extended for matrix coef.
#           n-1           n-2
#  coef[n]*z   + coef[n-1]*z   + ... + coef[1]  
  coef %*% z^(0:(length(coef)-1))
}


characteristicPoly <- function(a){
	# coefficients of the characteristic polynomial of a matrix
	# return a vector of the coefficients of the characteristic polynomial of a.
	# ref. Kailath "Linear Systems" p657

   tr <- function(a){ # calc the trace of a matrix
     sum(diag(a))
   }
   n <- dim(a)[1]
   if (n != dim(a)[2]) stop(" arguement must be a square matrix.")
   s <-array(0,c(n,n,n))
   if (n==1) a2 <- -a
   else
     {a2 <- rep(0,n)  
      s[1,,] <- diag(1,n)
      for (i in 1:(n-1))
        {a2[i] <- -tr(s[i,,]%*%a)/i
         s[i+1,,] <- (s[i,,]%*%a) + diag(a2[i],n)
        }
      a2[n] <- -tr(s[n,,]%*%a)/n
     }
   a2
}
companionMatrix <- function(a)
{# return the (top) companion matrix for a 3 dim array a (polynomial matrix), 
#  where the 1st dim corresponds to coefs. of polynomial powers (lags).
# ref. Kailath "Linear Systems" p659
  p <- dim(a)[2]
  if (p!= dim(a)[3]) 
    stop("companion matrix can only be computed for square matrix polynomials")
  l <- dim(a)[1]  # 1+ order of a
  C <- rbind(matrix(0,p,l*p), cbind(diag(1,(l-1)*p,(l-1)*p), matrix(0,(l-1)*p,p)))
  for (i in 1:l) C[1:p,((l-1)*p+1):(l*p)] <- -a[l,,]
  C
}

############################################################

#     internal utility functions    <<<<<<<<<<

############################################################

nstates <- function(x) UseMethod( "nstates")
nstates.SS <- function(x)  nrow(x$F) 
nstates.ARMA <- function(x) stop("ARMA models do not have a state vector.")
nstates.TSestModel <- function(x) nstates(TSmodel(x))

nseriesInput <- function(x) UseMethod( "nseriesInput")
nseriesInput.default <- function(x)
   if (is.null(x$input)) 0 else nseries(x$input)
nseriesInput.SS <- function(x)   {if (is.null(x$G)) 0 else  dim(x$G)[2] }
nseriesInput.ARMA <- function(x) {if (is.null(x$C)) 0 else dim(x$C)[3] }
nseriesInput.TSestModel <- function(x){nseriesInput(x$data)}

nseriesOutput <- function(x)UseMethod("nseriesOutput")
nseriesOutput.default <- function(x)
   if (is.null(x$output)) 0 else nseries(x$output)
nseriesOutput.SS <- function(x){dim(x$H)[1] }
nseriesOutput.ARMA <- function(x){dim(x$A)[2] }
nseriesOutput.TSestModel <- function(x){nseriesOutput(x$data)}


checkConsistentDimensions <- function(obj1, obj2=NULL)
   # check data & model dimensions
   UseMethod("checkConsistentDimensions")
   
checkConsistentDimensions.TSdata <- function(obj1, obj2=NULL)
  {if(is.null(obj2))
     stop("A TSmodel must be supplied as the second argument to this method.")
   checkConsistentDimensions(obj2, obj1) #(data,model) -> (model, data)
   }
checkConsistentDimensions.TSestModel <- function(obj1, obj2=NULL)
   {if(is.null(obj2)) obj2 <- obj1$data
    checkConsistentDimensions(obj1$model, obj2)
   }

checkConsistentDimensions.SS <- function(obj1, obj2=NULL) # (model, data)
 {m <- ncol(obj1$G)
  n <- nrow(obj1$F)
  p <- nrow(obj1$H)
  if (n!= ncol(obj1$F)) stop("Model F matrix must be square.")
  if (n!= ncol(obj1$H))
      stop("Model H matrix have second dimension consistent with matrix F.")
  if (!is.null(obj1$G)) if(n!= nrow(obj1$G))
      stop("Model G matrix have first dimension consistent with matrix F.")
  if (!is.null(obj1$K)) if(n!= nrow(obj1$K))
      stop("Model K matrix have first dimension consistent with matrix F.")
  if (!is.null(obj1$K)) if(p!= ncol(obj1$K))
      stop("Model K matrix have second dimension consistent with matrix H.")
  if (!is.null(obj1$Q)) if(n!= nrow(obj1$Q))
      stop("Model Q matrix must have first dimension consistent with matrix F.")
  if (!is.null(obj1$R)) if(p!= nrow(obj1$R))
      stop("Model R matrix must have first dimension consistent with matrix H.")
  if (!is.null(obj1$R)) if(p!= ncol(obj1$R))
      stop("Model R matrix must be square.")

  if (!is.null(obj2))
   {if(dim(obj1$H)[1] != nseriesOutput(obj2))
       stop("Model and data output dimensions do not correspond.\n")
    if(is.null(obj1$G))
      {if(0 != nseriesInput(obj2))
        stop("Model and data input dimensions do not correspond.\n")
      }
    else
      {if(dim(obj1$G)[2] != nseriesInput(obj2))
         stop("Model and data input dimensions do not correspond.\n")
      }
   }
  return(TRUE)
 }

checkConsistentDimensions.ARMA <- function(obj1, obj2=NULL) #(model, data)
 {p <-dim(obj1$A)[2]
  if (p!= dim(obj1$A)[3]) stop("Model A array dim 2 and 3 should be equal.")
  if (p!= dim(obj1$B)[2]) stop("Model B array dim inconsistent with array A.")
  if (p!= dim(obj1$B)[3]) stop("Model B array dim inconsistent with array A.")
  if (!is.null(obj1$C))
      if (p!= dim(obj1$C)[2]) stop("Model C array dim inconsistent with array A.")

  if (!is.null(obj2))
   {if(dim(obj1$A)[2] != nseriesOutput(obj2))
       stop("Model and data output dimensions do not correspond.\n")
    if(is.null(obj1$C))
      {if(0 != nseriesInput(obj2))
        stop("Model and data input dimensions do not correspond.\n")
      }
    else
      {if(dim(obj1$C)[3] != nseriesInput(obj2))
         stop("Model and data input dimensions do not correspond.\n")
      }
   }
  return(TRUE)
 }

checkConsistentDimensions.default <- function(obj1, obj2=NULL)
   {stop(paste("No method for obj1ect of class ", class(obj1), "\n"))}



TSestModel <- function(obj) UseMethod("TSestModel")
TSestModel.TSestModel <- function(obj) {obj} # return everything


TSmodel <- function(obj, ...) UseMethod("TSmodel")
TSmodel.TSmodel <- function(obj, ...){obj}
#  (... further arguments, currently disregarded)

TSmodel.TSestModel <- function(obj, ...)
  {#  (... further arguments, currently disregarded)
   # Return a TSmodel object but also retains convergence info (if not null).
   obj$model$converged <- obj$converged
   obj$model
  }




ARMA <- function(A=NULL, B=NULL, C=NULL, TREND=NULL, 
          constants=NULL,
	  description=NULL, names=NULL, input.names=NULL, output.names=NULL) 
  {if  (is.null(A)) stop("specified structure is not correct for ARMA model.")
   # and fix some simple potential problems
   if(is.null(dim(A)))   A <- array(A, c(length(A),1,1))
   if(is.null(B)) stop("B array must be specified for ARMA class models.")
   if(is.null(dim(B)))   B <- array(B, c(length(B),1,1))
   if(2==length(dim(B))) B <- array(B, c(1, dim(B)))
   if(!is.null(C) && is.null(dim(C))) C <- array(C, c(length(C),1,1))
   model <- list(A=A, B=B, C=C, TREND=TREND, 
                 constants=constants, description=description)
   class(model) <- c("ARMA","TSmodel") # constructor
   if(!is.null(names)) seriesNames(model) <- names
   else
     {if(!is.null( input.names))  seriesNamesInput(model) <- input.names
      if(!is.null(output.names)) seriesNamesOutput(model) <- output.names
     }
   checkConsistentDimensions(model)
   setTSmodelParameters(model)
  }



SS <- function(F.=NULL, G=NULL, H=NULL, K=NULL, Q=NULL, R=NULL,
          z0=NULL, P0=NULL, rootP0=NULL,
          constants=NULL,
          description=NULL, names=NULL, input.names=NULL, output.names=NULL)   
  {if (is.null(F.) | is.null(H))
       stop("specified stucture is not correct for SS model.")
   # and fix some simple potential problems
   if(1==length(F.))  
     {if(!is.matrix(F.))                 F. <- matrix(F.,1,1)
      if(!is.matrix(H))                  H  <- matrix(H,length(H),1)
      if(!is.null(G) && !is.matrix(G))   G  <- matrix(G,1,length(G))
      if(!is.null(K) && !is.matrix(K))   K  <- matrix(K,1,length(K))
      if(!is.null(Q) && !is.matrix(Q))   Q  <- matrix(Q,1,length(Q))
     }
   model <- list(F=F., G=G, H=H, K=K, Q=Q, R=R, z0=z0, P0=P0, rootP0=rootP0,
                 constants=constants, description=description)
   if (!is.null(model$K)) class(model) <- c("innov","SS","TSmodel" ) else
   if (!is.null(model$Q)) class(model) <- c( "nonInnov","SS","TSmodel") # constructor
   else stop("specified stucture is not correct for SS model.")
   if(!is.null(names)) seriesNames(model) <- names
   else
     {if(!is.null( input.names))  seriesNamesInput(model) <-  input.names
      if(!is.null(output.names)) seriesNamesOutput(model) <- output.names
     }
   checkConsistentDimensions(model)
   setTSmodelParameters(model)
  }



"coef<-" <- function(object, value) UseMethod("coef<-")
"coef<-.default" <- function(object, value) {
   object$coefficients <- value
   object
   }
   
coef.TSmodel <- function(object, ...) object$coefficients
#  (... further arguments, currently disregarded)

coef.TSestModel <- function(object, ...) coef(TSmodel(object))
#  (... further arguments, currently disregarded)

is.TSmodel <- function(obj){inherits(obj,"TSmodel")}
is.TSestModel <- function(obj){inherits(obj,"TSestModel")}
is.SS <- function(obj){inherits(obj,"SS")}
is.innov.SS <- function(obj){inherits(obj,"SS")& inherits(obj,"innov")}
is.nonInnov.SS <- function(obj){inherits(obj,"SS")&inherits(obj,"nonInnov")}
is.ARMA <- function(obj){inherits(obj,"ARMA")}

# complete parameter info. based on representation info. 

setTSmodelParameters <- function(model, constants=model$constants)  
   UseMethod("setTSmodelParameters")

setTSmodelParameters.TSestModel <- function(model, constants=TSmodel(model)$constants)  
  setTSmodelParameters(TSmodel(model), constants=constants)


setTSmodelParameters.SS <- function(model, constants=model$constants) { 

 locateSS <- function(A,Ac,a,I,J,plist)# local function for locating parameters
  {indicate <-  (A==1.0)                # constants
   if (!is.null(Ac)) indicate <- indicate | Ac  # Ac is T for fixed entries
   plist$const <- c(plist$const,A[indicate])
   plist$const.location <- c(plist$const.location,rep(a,sum(indicate)))
   if (is.vector(A))  # z is not a matrix, (a=="z") also should work
     {plist$const.i <- c(plist$const.i,(1:length(A))[indicate]) 
      plist$const.j <- c(plist$const.j,(1:length(A))[indicate]) #dup.but ok
     }
   else if(is.matrix(A))
     {plist$const.i <- c(plist$const.i,row(A)[indicate]) 
      plist$const.j <- c(plist$const.j,col(A)[indicate])
     }
   else stop("The dimension of something in the SS model structure is bad.")
   indicate <- (A!=0.0) & (A!=1.0) 
   if (!is.null(Ac)) indicate <- indicate & (!Ac)
   coef(plist) <- c(coef(plist), A[indicate])          # parameters
   plist$location <- c(plist$location,rep(a,sum(indicate)))
   if (is.vector(A))  # z is not a matrix, (a=="z") also should work
     {plist$i <- c(plist$i,(1:length(A))[indicate]) 
      plist$j <- c(plist$j,(1:length(A))[indicate]) #dup.but ok
     }
   else if(is.matrix(A))
     {plist$i <- c(plist$i,row(A)[indicate])
      plist$j <- c(plist$j,col(A)[indicate])
     }
   else stop("The dimension of something in the SS model structure is bad.")
   plist
 }# end locateSS    

    m <-dim(model$G)[2]
    n <-dim(model$F)[1]
    if (n!= dim(model$F)[2]) stop("Model F matrix must be square.")
    if (n!= dim(model$H)[2])
      stop("Model H matrix have second dimension consistent with matrix F.")
    if (!is.null(model$G)) if(n!= dim(model$G)[1])
      stop("Model G matrix have first dimension consistent with matrix F.")
    p <-dim(model$H)[1]
    plist <- list(coefficients=NULL,location=NULL,i=NULL,j=NULL,
                 const=NULL,const.location=NULL,const.i=NULL,const.j=NULL)
    plist <- locateSS(model$F, constants$F,"f",n,n,plist)
    if(!is.null(m)) plist <- locateSS(model$G,constants$G,"G",n,m,plist)
    plist <- locateSS(model$H,constants$H,"H",p,n,plist)
    if(!is.null(model$z0)) plist <- locateSS(model$z0,constants$z0,"z",p,n,plist)
    if(!is.null(model$P0)) plist <- locateSS(model$P0,constants$P0,"P",p,n,plist)
    if(!is.null(model$rootP0)) plist <- locateSS(model$rootP0,constants$rootP0,"rP",p,n,plist)
    if (is.innov.SS(model)) 
      {plist <- locateSS(model$K,constants$K,"K",n,p,plist)
       # note constants$H, etc are logical arrays (if not NULL) to indicate
       #      parameters which are to remain fixed, so that setTSmodelParameters
       #       knows to put them in const. This feature has not been used much.
      }
    else
      {plist <- locateSS(model$Q,constants$Q,"Q",n,n,plist)
       plist <- locateSS(model$R,constants$R,"R",p,p,plist)
      }

#    list.add(model, names(plist) ) <- plist
    model[ names(plist) ] <- plist
    model
}

setTSmodelParameters.ARMA <- function  (model, constants=model$constants) { 

 locateARMA <- function(A,Ac, a,I,J,L,plist){ # local function for locating parameters
  ind <- function(x, i) # equivalent of row and col for higher dim arrays.
   {
	d <- dim(x)
	id <- 1:length(d)
	id <- c(i, id[i!=id])
	y <- array(1:(d[i]), dim(x)[id])
	aperm(y, order(id))
   }
   indicate <-  (A==1.0)                # constants
   if (!is.null(Ac)) indicate <- indicate | Ac  # Ac is T for fixed entries
   plist$const <- c(plist$const,A[indicate])
   plist$const.location <- c(plist$const.location,rep(a,sum(indicate)))
   if (is.vector(A))# trend is a vector not an array (a=="t")
     {plist$const.l <- c(plist$const.l, rep(0,sum(indicate))) 
      plist$const.i <- c(plist$const.i, (1:length(A))[indicate]) 
      plist$const.j <- c(plist$const.j, rep(0,sum(indicate)))
     }
   else if (is.matrix(A))  # trend is a matrix
     {plist$const.l <- c(plist$const.l, rep(0,sum(indicate))) 
      plist$const.i <- c(plist$const.i, row(A)[indicate]) 
      plist$const.j <- c(plist$const.j, col(A)[indicate])
     }
   else if(3 == length(dim(A)))
     {plist$const.l <- c(plist$const.l, ind(A,1)[indicate]) 
      plist$const.i <- c(plist$const.i, ind(A,2)[indicate]) 
      plist$const.j <- c(plist$const.j, ind(A,3)[indicate])
     }
   else stop("The dimension of something in the ARMA model structure is bad.")
   indicate <- (A!=0.0) & (A!=1.0)
   if (!is.null(Ac)) indicate <- indicate & (!Ac)
   coef(plist) <- c(coef(plist), A[indicate])          # parameters
   plist$location <- c(plist$location,rep(a,sum(indicate)))
   if (is.vector(A))# trend is a vector not an array (a=="t")
     {plist$l <- c(plist$l, rep(0,sum(indicate))) 
      plist$i <- c(plist$i, (1:length(A))[indicate]) 
      plist$j <- c(plist$j, rep(0,sum(indicate)))
     }
   else if (is.matrix(A))  # trend is a matrix
     {plist$l <- c(plist$l, rep(0,sum(indicate))) 
      plist$i <- c(plist$i, row(A)[indicate]) 
      plist$j <- c(plist$j, col(A)[indicate])
     }
   else if(3 == length(dim(A)))
     {plist$l <- c(plist$l, ind(A,1)[indicate] )
      plist$i <- c(plist$i, ind(A,2)[indicate] )
      plist$j <- c(plist$j, ind(A,3)[indicate])
     }
   else stop("The dimension of something in the ARMA model structure is bad.")
   plist  
 }# end locateARMA

       m <-dim(model$C)[3]
       p <-dim(model$A)[2]
       a <-dim(model$A)[1]
       b <-dim(model$B)[1]
       cc <-dim(model$C)[1]
       if (p!= dim(model$A)[3]) stop("Model A array dim 2 and 3 should be equal.")
       if (p!= dim(model$B)[2]) stop("Model B array dim inconsistent with array A.")
       if (p!= dim(model$B)[3]) stop("Model B array dim inconsistent with array A.")
       if (!is.null(model$C))
          if (p!= dim(model$C)[2]) stop("Model C array dim inconsistent with array A.")

       plist <- list(coefficients=NULL,location=NULL,i=NULL,j=NULL,
                         const=NULL,const.location=NULL,
                         const.i=NULL,const.j=NULL,l=NULL,const.l=NULL)
       plist <- locateARMA(model$A, constants$A, "A",p,p,a,plist)
       plist <- locateARMA(model$B, constants$B, "B",p,p,b,plist)
       if(!is.null(cc)) plist <- 
                locateARMA(model$C, constants$C, "C",p,m,cc,plist)
       if(!is.null(model$TREND)) plist <- 
	        locateARMA(model$TREND, constants$TREND, "t",p,m,cc,plist)

#       list.add(model, names(plist) ) <- plist
       model[ names(plist) ] <- plist
       model
}


setArrays <- function(model, coefficients=NULL, constants=NULL)  
  # complete representaion info. based on parameter info.  
  UseMethod("setArrays")
   
setArrays.TSestModel <- function(model, coefficients=NULL, constants=NULL)  
 setArrays(TSmodel(model), coefficients=coefficients, constants=constants) 
    
setArrays.SS <- function(model, coefficients=NULL, constants=NULL){
	# N.B. Dimension and class (innov/ nonInnov) info. is assumed accurate
    if (is.null(coefficients)) coefficients   <- coef(model)
                        else   model$coefficients <- coefficients
    a.pos  <- model$location
    i.pos  <- model$i
    j.pos  <- model$j
    const  <-if (is.null(constants)) model$const else constants #untested
    ca.pos <- model$const.location
    ci.pos <- model$const.i
    cj.pos <- model$const.j  
    m <-dim(model$G)[2]
    n <-dim(model$F)[1]
    p <-dim(model$H)[1]
    f       <-  matrix(0,n,n)    # F     
    if(!is.null(m)) G        <-  matrix(0,n,m)    # control 
    H       <-  matrix(0,p,n)    # H       
    K        <-  matrix(0,n,p)    # Kalman gain
    if (is.null(model$Q)) Q <-  matrix(0,n,n)        # eventually discarded
    else                  Q <- array(0,dim(model$Q)) # system noise might not be nxn
    R        <-  diag(0,p,p)     # measurement noise
    z       <-  rep(0,n)        # initial state
    P       <-  diag(0,n)       # initial tracking error
    rP      <-  diag(0,n)       # root initial tracking error
    if (length(coefficients)>0) 
       {i <- a.pos == "f"
        f[cbind(i.pos[i],j.pos[i])] <- coefficients[i]
        if(!is.null(m)) 
          {i <- a.pos == "G"
           G[cbind(i.pos[i],j.pos[i])] <- coefficients[i]
          }
        i <- a.pos == "H"
        H[cbind(i.pos[i],j.pos[i])] <- coefficients[i]
        i <- a.pos == "K"
        K[cbind(i.pos[i],j.pos[i])] <- coefficients[i]
        i <- a.pos == "Q"
        Q[cbind(i.pos[i],j.pos[i])] <- coefficients[i]
        i <- a.pos == "R"
        R[cbind(i.pos[i],j.pos[i])] <- coefficients[i]
        i <- a.pos == "z"
        z[i.pos[i]] <- coefficients[i]
        i <- a.pos == "P"
        P[cbind(i.pos[i],j.pos[i])] <- coefficients[i]
        i <- a.pos == "rP"
        rP[cbind(i.pos[i],j.pos[i])] <- coefficients[i]
       }
    if (length(const)>0) 
       {i <- ca.pos == "f"
        f[cbind(ci.pos[i],cj.pos[i])] <- const[i]
        if(!is.null(m)) 
          {i <- ca.pos == "G"
           G[cbind(ci.pos[i],cj.pos[i])] <- const[i]
          }
        i <- ca.pos == "H"
        H[cbind(ci.pos[i],cj.pos[i])] <- const[i]
        i <- ca.pos == "K"
        K[cbind(ci.pos[i],cj.pos[i])] <- const[i]
        i <- ca.pos == "Q"
        Q[cbind(ci.pos[i],cj.pos[i])] <- const[i]
        i <- ca.pos == "R"
        R[cbind(ci.pos[i],cj.pos[i])] <- const[i]
        i <- ca.pos == "z"
        z[ci.pos[i]] <- const[i]
        i <- ca.pos == "P"
        P[cbind(ci.pos[i],cj.pos[i])] <- const[i]
        i <- ca.pos == "rP"
        rP[cbind(ci.pos[i],cj.pos[i])] <- const[i]
       }
    model$F <- f 
    if(!is.null(m)) model$G <- G
    model$H <- H
    if (is.innov.SS(model))
      model$K <- K  
    else
      {model$Q <-Q
       model$R <-R
      }
    if(!is.null(model$z0)) model$z0 <- z
    if(!is.null(model$P0)) model$P0 <- P
    if(!is.null(model$rootP0)) model$rootP0 <- rP
    model
} #end setArrays.SS

setArrays.ARMA <- function(model, coefficients=NULL, constants=NULL) { 
	# N.B. Dimension and class info. is assumed accurate
       if (is.null(coefficients)) coefficients    <- coef(model)
                        else   model$coefficients <- coefficients
       a.pos  <- model$location
       i.pos  <- model$i
       j.pos  <- model$j
       const  <-if (is.null(constants)) model$const else constants #untested
       ca.pos <- model$const.location
       ci.pos <- model$const.i
       cj.pos <- model$const.j  
       l.pos  <- model$l
       cl.pos <- model$const.l
       A  <-  array(0,dim(model$A))
       B  <-  array(0,dim(model$B))
       m <- dim(model$C)[3]
       p <- dim(model$A)[2]
       TREND <- if(is.null(model$TREND)) NULL else rep(0,length(model$TREND))
       if(!is.null(m)) C  <-  array(0,dim(model$C))
       if (length(coefficients)>0) 
          {i <- a.pos == "A"
           A[cbind(l.pos[i],i.pos[i],j.pos[i])] <- coefficients[i]
           i <- a.pos == "B"
           B[cbind(l.pos[i],i.pos[i],j.pos[i])] <- coefficients[i]
           if(!is.null(m))
             {i <- a.pos == "C"
              C[cbind(l.pos[i],i.pos[i],j.pos[i])] <- coefficients[i]
             }
           i <- a.pos == "t"
           if(!is.null(TREND)) TREND[i.pos[i]] <- coefficients[i]
          }
    if (length(const)>0) 
          {i <- ca.pos == "A"
           A[cbind(cl.pos[i],ci.pos[i],cj.pos[i])] <- const[i]
           i <- ca.pos == "B"
           B[cbind(cl.pos[i],ci.pos[i],cj.pos[i])] <- const[i]
           if(!is.null(m))
             {i <- ca.pos == "C"
              C[cbind(cl.pos[i],ci.pos[i],cj.pos[i])] <- const[i]
             }
           i <- ca.pos == "t"
           if(!is.null(TREND)) TREND[ci.pos[i]] <- const[i]
          }
      model$A <- A
      model$B <- B
      if(!is.null(m)) model$C <- C 
      model$TREND <- TREND
      model 
} #end setArrays.ARMA


############################################################

# following is a workaround for ar in ts library

DSE.ar <- function(data, ...) {
  #fix for ar in R ts library (so that univariate case also gives array result)
  # before R 1.9.0 required ts not stats
  res <- stats::ar(ts(data), ...)
  if (! is.array(res$ar)) res$ar <- array(res$ar, c(length(res$ar),1,1))
  res
  }

printTestValue <- function (x, digits = 16){
    cat("c( ")
    if (all(is.na(x)))       cat("NAs")
    else if (is.null(x))     cat("NULL")
    else if (is.logical(x))  cat(x)
    else if (!is.R())        print(x, digits = digits)
    else if (is.matrix(x)) {
       for (i in 1:nrow(x)) {
	 cat("\n      ")
         for (j in 1:ncol(x)) 
	   cat(formatC(x[i,j], digits=digits, format="g"), ", ")
	 }
       cat("), ", ncol(x), ", ", nrow(x), ")\n")
       }
    else for (i in 1:length(x)) cat(", ", formatC(x[i], digits=digits, format="g"))
    cat(")\n")
    invisible()
    }


############################################################

#  Model simulation functions (to generate data)   <<<<<<<<<<

############################################################



simulate <- function(model, ...)UseMethod("simulate")

simulate.TSestModel <- function(model, input=inputData(model),
			sd=NULL, Cov=NULL, ...)
  {if (is.null(sd) & is.null(Cov)) Cov <- model$estimates$cov
   simulate(TSmodel(model), input=input, sd=sd, Cov=Cov, ...)
  }

simulate.SS <- function(model, input=NULL,
                 start=NULL, freq=NULL, sampleT=100, 
                 noise=NULL, sd=1, Cov=NULL, rng=NULL, 
                 compiled=.DSEflags()$COMPILED, ...)
{#  (... further arguments, currently disregarded)
if (is.TSestModel(model)) model <- TSmodel(model)
if (!is.SS(model)) stop("simulate.SS expecting an SS TSmodel.")
if (!checkConsistentDimensions(model)) stop("The SS model is not correct.")

 FF<-    model$F
 G <-    model$G
 H <-    model$H
 m <- dim(G)[2]
 if(is.null(m)) m <-0
 n <- dim(FF)[2]
 p <- dim(H)[1]
 if (m!=0)
   {if( is.null(input)) stop("input series must be supplied for this model.")
    if (sampleT != Tobs(input) ) input <- tfTruncate(input, end=sampleT)
   }
 if (is.innov.SS(model))
   {K <-    model$K}
 else
   {Q <-    model$Q
    if (ncol(Q)<n) Q <- cbind(Q,matrix(0,n,n-ncol(Q))) # Q assumed square
    R <-    model$R
   }
 
 if(is.null(rng)) rng <- setRNG() # returns setting so don't skip if NULL
 else        {old.rng <- setRNG(rng);  on.exit(setRNG(old.rng))  }
 
if (!is.null(noise) && is.matrix(noise)) noise <- list(w=noise)

set.ts <- TRUE             
if (!is.null(start))
  {if (!is.null(freq))   tf <- list(start=start, frequency=freq)
   else
      {warning("start is specified but not freq. Using freq=1.") 
       tf <- list(start=start, frequency=1)
      }
  }  
else if( (!is.null(input))   && is.tframed(input))   tf <- tframe(input)
else if ((!is.null(noise$w)) && is.tframed(noise$w)) tf <- tframe(noise$w)
else set.ts <-  FALSE



# It would be better to use makeTSnoise here (lag=1 and w0<-c(noise$w0) ) and
#  possibly two calls to get e for non-innov models. However, this will affect
#  historical comparisons so it needs to be done carefully!
 e <-NULL
 if (is.null(noise))
   {if (!is.innov.SS(model)) 
      {w0 <- rnorm(n)
       e <- matrix(rnorm(sampleT*n),sampleT,n)
       w <- matrix(rnorm(sampleT*p),sampleT,p)
      }
    else 
      {w <- matrix(NA,sampleT,p)
       w0 <- rep(NA,p)
       if (!is.null(Cov))
         {if(length(Cov) == 1) Cov <- diag(Cov, p)
	  W <- t(chol(Cov))
	  w.help <- t(W %*% matrix(rnorm((sampleT+1)*p),nrow=p,ncol=sampleT+1))
	  w0 <- w.help[1,]
	  w <- w.help[-1,]
         }
       else
         {if (length(sd)==1) sd <-rep(sd,p)
          for (i in 1:p)
            {w0[i] <- rnorm(1,sd=sd[i])
             w[,i] <- rnorm(sampleT,sd=sd[i])
            }
         }
      }
   }
 else
   {#  noise is not null
   # already done if (is.matrix(noise)) noise <- list(w=noise)
   if (is.null(noise$w))
       stop("supplied noise structure is not correct. w must be specified.")
   if (is.null(noise$w0)) noise$w0 <- rep(0,p)
   if ((!is.innov.SS(model)) &&  is.null(noise$e))
       stop("supplied noise structure is not correct. e must be specified for non-innovations form models.")

    w0 <- noise$w0
    if (is.matrix(w0)) w0 <- rep(0,p) # see note in man re VAR compatiblity
    w <- noise$w
    e <- noise$e
    sampleT <- Tobs(w)
   }

 y <- matrix(0,sampleT,p)
 state <- matrix(0,sampleT,n)
 if(is.null(model$z0)) z <- rep(0,n) # initial state
 else  z <- model$z0        

 if (is.innov.SS(model)) 
   {z <-  c(FF%*% z) + c(K %*% w0)
    if (m !=0) z <-  z + c(G%*%input[1,]) 
    y[1,]  <-  c(H %*% z) + c(w[1,])
   }   
 else      
   {z <-  c(FF%*% z) + c(Q %*% w0)
    if (m !=0) z <-  z + c(G %*% input[1,])  
    y[1,]  <-  c(H %*% z) + c(R %*% w[1,])
   }
 state[1,]<-z
 # Note: the first period is done above for both compiled and S versions so
 #    that initial conditions do not have to be passed to compiled code.
 if (compiled)
   {if (is.innov.SS(model))
      {Q <-matrix(0,n,n)
       R <-matrix(0,p,p)
       e <- matrix(0,sampleT,n)
      }
    else K <- matrix(0,n,p)
    if (m==0) 
      {input <- matrix(0,sampleT,1)
       G <- matrix(0,n,1)
      }
 #   storage.mode(y)     <- "double"
 #   storage.mode(state) <- "double"
 #   storage.mode(input) <- "double"
 #   storage.mode(w) <- "double"
 #   storage.mode(e) <- "double"
 #   storage.mode(FF) <- "double"
 #   storage.mode(G) <- "double"
 #   storage.mode(H) <- "double"
 #   storage.mode(K) <- "double"
 #   storage.mode(Q) <- "double"
 #   storage.mode(R) <- "double"
   #r <- .Fortran(simss_sym, 
   r <- .Fortran("simss", 
                         y=if(is.double(y)) y else as.double(y), 
                         state=if(is.double(state)) state else as.double(state), 
                         as.integer(m),
                         as.integer(n),
                         as.integer(p), 
                         as.integer(sampleT),  
                         if(is.double(input)) input else as.double(input),
                         if(is.double(w)) w else as.double(w),
                         if(is.double(e)) e else as.double(e),
                         if(is.double(FF)) FF else as.double(FF),
                         if(is.double(G)) G else as.double(G),   
                         if(is.double(H)) H else as.double(H),
                         if(is.double(K)) K else as.double(K), 
                         if(is.double(Q)) Q else as.double(Q),	 
                         if(is.double(R)) R else as.double(R),    
                         as.integer(is.innov.SS(model)), 
			 PACKAGE="dse"
			 )[c("y","state")]
    y <- r$y
    state <- r$state
    if (m==0) input <- NULL
   }
 else
   {for (Time in 2:sampleT)  
      {z <-  c(FF%*% z) 
       if (m !=0) z <-  z + c(G%*%input[Time,]) 
       if (is.innov.SS(model))  z <-  z + c(K%*%w[Time-1,])
       else       z <-  z + c(Q%*%e[Time-1,])
       state[Time,] <- z
       if (is.innov.SS(model)) y[Time,]  <-  c(H %*% z) + c(w[Time,])   
       else      y[Time,]  <-  c(H %*% z) + c(R %*% w[Time,])
   }  }
 if (set.ts)
   { y     <- tframed(y,     tf=tf, names=seriesNamesOutput(model)) 
     state <- tframed(state, tf=tf )
   }
 else  seriesNames(y) <- seriesNamesOutput(model)
 TSdata(list(input=input,output=y, state=state, version=version, 
   model=model, description="data generated by simulate.ss",
   noise=list(w0=w0,w=w, e=e, rng=rng, Cov=Cov, sd=sd)))
}

simulate.ARMA <- function(model, y0=NULL, input=NULL, input0=NULL,
                start=NULL, freq=NULL, sampleT=100,
                noise=NULL, sd=1, Cov=NULL,
                rng=NULL, noise.model=NULL, 
                compiled=.DSEflags()$COMPILED, ...)
{# see details in help("ARMA") and help("simulate.ARMA")

if (is.TSestModel(model)) model <- TSmodel(model)
if (!is.ARMA(model)) stop("simulate.ARMA expecting an ARMA TSmodel.")
if (!checkConsistentDimensions(model)) stop("The ARMA model is not correct.")

A<-    model$A
B <-    model$B
C <-    model$C
m <- dim(C)[3]
if (is.null(m)) m <-0
p <- dim(A)[2]
a <-dim(A)[1]
b <-dim(B)[1]
cc <- if (is.null(C)) 0 else  dim(C)[1]

TREND <- model$TREND
if (p == length(TREND)) TREND <- t(matrix(TREND, p, sampleT))

if (m!=0)
   {if( is.null(input)) stop("input series must be supplied for this model.")
    if (sampleT != Tobs(input) ) input <- tfTruncate(input, end=sampleT)
   }
 
if (p==1) invA0 <- matrix(1/A[1,,],1,1)
else      invA0 <- solve(A[1,,])
for (l in 1:a) A[l,,] <- invA0 %*% A[l,,]      # set A(0) = I      
for (l in 1:b) B[l,,] <- invA0 %*% B[l,,] 
if (m!=0) for (l in 1:dim(C)[1]) C[l,,] <- invA0 %*% C[l,,]  
if(!is.null(TREND)) TREND <- t(invA0 %*% t(TREND))

set.ts <- TRUE             

if (!is.null(noise)) {
   if (is.matrix(noise)) noise <- list(w=noise)
   if (is.null(noise$w0)) noise$w0 <-matrix(0,b,p)
   }

if (!is.null(start))
  {if (!is.null(freq))   tf <- list(start=start, frequency=freq)
   else
      {warning("start is specified but not freq. Using freq=1.") 
       tf <- list(start=start, frequency=1)
      }
  }  
else if( (!is.null(input))   && is.tframed(input))   tf <- tframe(input)
else if ((!is.null(noise$w)) && is.tframed(noise$w)) tf <- tframe(noise$w)
else set.ts <-  FALSE


noise <- makeTSnoise(sampleT,p,b, noise=noise, rng=rng,
                        Cov=Cov, sd=sd, 
			noise.model=noise.model,
                        start=start, frequency=freq)
if (is.null(sampleT)) sampleT<-noise$sampleT
 
 if(is.null(y0)) y0<-matrix(0,a,p)
 if((m!=0) & is.null(input0)) input0 <- matrix(0,dim(C)[1],m)
 y <- matrix(0,sampleT,p)  
 if (compiled)
   {if (m==0) 
      {input <- matrix(0,sampleT,1)
       input0 <- matrix(0,1,1)
       C <- matrix(0,1,1)
      }
    if (is.null(TREND)) TREND<- matrix(0,sampleT, p)
#    yo<- list(y=y, y0,m,p, a, b, cc, sampleT,input,input0,w,w0,A,B, C,TREND)
#    storage.mode(y)     <- "double"
#    storage.mode(y0)     <- "double"
#    storage.mode(input0)     <- "double"
#    storage.mode(noise$w)     <- "double"
#    storage.mode(noise$w0)     <- "double"
#    storage.mode(A)     <- "double"
#    storage.mode(B)     <- "double"
#    storage.mode(C)     <- "double"
#    storage.mode(TREND)     <- "double"

    #y <- .Fortran(simrma_sym, 
     y <- .Fortran("simrma", 
                        y=if(is.double(y)) y else as.double(y), 
                         if(is.double(y0)) y0 else as.double(y0),
                         as.integer(m),
                         as.integer(p), 
                         as.integer(a), 
                         as.integer(b), 
                         as.integer(cc), 
                         as.integer(sampleT),  
                         as.double(input[1:sampleT,]),
                         if(is.double(input0)) input0 else as.double(input0),
                         if(is.double(noise$w)) noise$w else as.double(noise$w),
                         if(is.double(noise$w0)) noise$w0 else as.double(noise$w0),
                         if(is.double(A)) A else as.double(A),
                         if(is.double(B)) B else as.double(B),   
                         if(is.double(C)) C else as.double(C),
                         if(is.double(TREND)) TREND else as.double(TREND), 
			 PACKAGE="dse"
			 ) [["y"]]

    if (m==0) 
      {input  <- NULL
       input0 <- NULL
      }
   }
 else
  {w0 <- noise$w0
   w<-noise$w
   if(!is.null(TREND)) y <- TREND  
   for (Time in 1:sampleT)  
   {for (l in 2:a) 
       if(Time+1-l<=0)
          if (p==1) y[Time,] <- y[Time,]-c(A[l,,]  *  y0[l-Time,]) 
          else      y[Time,] <- y[Time,]-c(A[l,,] %*% y0[l-Time,])
       else                
          if (p==1) y[Time,] <- y[Time,]-c(A[l,,]  *  y[Time+1-l,])
          else      y[Time,] <- y[Time,]-c(A[l,,] %*% y[Time+1-l,])
    for (l in 1:b) 
       if (Time+1-l<=0)
          if (p==1) y[Time,]<- y[Time,] +c(B[l,,]  *  w0[l-Time,]) 
          else      y[Time,]<- y[Time,] +c(B[l,,] %*% w0[l-Time,])
       else
          if (p==1) y[Time,]<- y[Time,] +c(B[l,,]  *  w[Time+1-l,])
          else      y[Time,]<- y[Time,] +c(B[l,,] %*% w[Time+1-l,])
    if (m!=0) for (l in 1:cc) 
       if (Time+1-l<=0)
          if (m==1) y[Time,]<- y[Time,] + c(C[l,,]  *  input0[l-Time,])
          else      y[Time,]<- y[Time,] + c(C[l,,] %*% input0[l-Time,])
       else
          if (m==1) y[Time,]<-y[Time,] + c(C[l,,]  *  input[Time+1-l,])
          else      y[Time,]<-y[Time,] + c(C[l,,] %*% input[Time+1-l,])
  }}
 if (set.ts) y <- tframed(y, tf=tf, names=seriesNamesOutput(model)) 
 else seriesNames(y) <- seriesNamesOutput(model)
 TSdata(list(input=input,output=y, 
          model=model, input0=input0, 
          description="data generated by simulate.ARMA", 
#          prior.args=prior.args, post.args=post.args,
          noise=noise))
}

#######################################################################

#functions which work on data (and models)  <<<<<<<<<<

############################################################

#     likelihood and residual calculation functions  <<<<<<<<<<

############################################################

residualStats <- function(pred, data, sampleT=nrow(pred), warn=TRUE){  
   e <- if (is.null(pred))     -data[1:sampleT,,drop=FALSE]
        else if (is.null(data)) pred[1:sampleT,,drop=FALSE]
        else               pred[1:sampleT,,drop=FALSE] - data[1:sampleT,,drop=FALSE]
   p <- ncol(e)
   #Om <-t(e) %*% e /sampleT
   Om <-crossprod(e)/sampleT
   if (any(is.na(Om))) {like1 <- like2 <- 1e100}  
   else if (any(Om >1e100)) {like1 <- like2 <- 1e100}  
   else
     {v <- La.svd(Om) #eigenvalues are not robust to degenerate density.
#     i <- v$d!=0
     i <- v$d > (v$d[1]*sqrt(.Machine$double.eps))
      if (!all(i))
        {if(warn) warning("The cov. matrix is singular. Working on subspace.")
         v$d <- v$d[i]
      #   v$u <- v$u[i,i, drop=FALSE]
      #   v$vt <- v$vt[i,i, drop=FALSE]
      #   e <- e[,i, drop=FALSE]
        }
      like1 <- 0.5 * sampleT * log(prod(v$d)) # det
      # $vt is transposed in La.svd, but not v in svd
#      if (1 == length(v$d)) OmInv <-  v$u %*% (1/v$d) %*% v$vt 
#      else OmInv <-  v$u %*% diag(1/v$d) %*% v$vt # more robust than solve(Om)
#  following is equivalent
#      OmInv <-  v$u %*% sweep(v$vt,1,1/v$d, "*") 
#      like2 <- sum(e * (e %*% OmInv)) /2
#  but this works out to (sampleT*p/2) and fixing for degenerate distributions:
       like2 <- (sampleT*length(v$d))/2
     }
   const <- (sampleT * p * log(2 * pi))/2
   invisible(list(like=c(const+like1+like2, const, like1,like2),
                  cov =Om, pred=pred, sampleT=sampleT))
   }

sumSqerror <- function(coefficients, model=NULL, data=NULL, error.weights=NULL) {
 if ( is.null(model)) stop("model missing") 
 if ( is.null(data))  stop("data missing") 
 if ( is.null(error.weights)) stop("error.weights missing")  
 sum(l(setArrays(model,coefficients=coefficients), data,
       result="weighted.sqerror",error.weights=error.weights))
 }

l <- function(obj1, obj2, ...)UseMethod("l")
l.TSdata <- function(obj1, obj2, ...) {l(obj2, obj1, ...) }
l.TSestModel <- function(obj1, obj2, ...) {l(TSmodel(obj1), obj2, ...)}


l.ARMA <- function(obj1, obj2, sampleT=NULL, predictT=NULL,result=NULL,
                error.weights=0,  compiled=.DSEflags()$COMPILED, 
		warn=TRUE, return.debug.info=FALSE, ...)
{#  (... further arguments, currently disregarded)
  #  see help("l.ARMA")
model <- if (is.TSestModel(obj1)) TSmodel(obj1) else obj1
if (!is.ARMA(model)) stop("l.ARMA expecting an ARMA TSmodel.")

#dat <- freeze(obj2)
dat <- obj2
tf <- tfTruncate(tframe(outputData(dat)), end=predictT)

if(!checkConsistentDimensions(model,dat)) stop("dimension error")
if (is.null(sampleT))  sampleT  <- TobsOutput(dat)
if (is.null(predictT)) predictT <- sampleT
if (sampleT > predictT) stop("sampleT cannot exceed predictT.")
if (sampleT > TobsOutput(dat)) stop("sampleT cannot exceed length of data.")

if (0 != nseriesInput(dat))
  {if (TobsInput(dat) < predictT)
      stop("input data must be at least as long as requested prediction.")
   if (any(is.na(inputData(dat)))) stop("input data cannot contain NAs.")
   if(!all(startInput(dat) == startOutput(dat)))
	   stop("input and output data must start at the same time.")
  }
if (any(is.na(outputData(dat)))) stop("output data cannot contain NAs.")

u <- inputData(dat)
y <- outputData(dat)
A<-    model$A
B <-   model$B
C <-   model$C
m <- dim(C)[3]
if (is.null(m)) m <-0
p <- dim(A)[2]
a <-dim(A)[1]
b <-dim(B)[1]

TREND <- model$TREND
if (p == length(TREND)) TREND <- t(matrix(TREND, p, predictT))

if (m == 0) 
   {if(!is.null(u)) 
      stop("model parameters indicate an no input but input data exists!")
   }

if (compiled)
  {if (m==0)
     {C <- array(0,c(1,p,1))    # can't pass 0 length array to compiled
      u <- matrix(0,predictT,1)
     }
   else # since input and output are declared same dim in fortran
      if (Tobs(y) < Tobs(u)) y <- 
	   rbind(y, matrix(0, Tobs(u) - Tobs(y), nseries(y)))
   if (is.null(TREND)) TREND <- matrix(0,predictT, p)
   is  <- max(m,p)

#   storage.mode(error.weights)     <- "double"
#   storage.mode(u)     <- "double"
#   storage.mode(y)     <- "double"
#   storage.mode(A)     <- "double"
#   storage.mode(B)     <- "double"
#   storage.mode(C)     <- "double"
#   storage.mode(TREND)     <- "double"

   #r  <- .Fortran(arma_sym,
   r  <- .Fortran("arma",
                         pred=matrix(1e20,predictT,p), # pred,     
                         as.integer(length(error.weights)),
                         weighted.sqerror=matrix(0,sampleT,p),
                         error.weights= if(is.double(error.weights)) 
			           error.weights else as.double(error.weights),
                         as.integer( m), 
                         as.integer( p) ,      
                         as.integer( dim(A)[1]),  # 1+order of A  
                         as.integer( dim(B)[1]),  # 1+order of B  
                         as.integer( dim(C)[1]),  # 1+order of C  
                         sampleT=as.integer(sampleT),
                         predictT=as.integer(predictT),
                         as.integer(Tobs(y)), 
                         if(is.double(u)) u else as.double(u), # as.double() works ok with compiled but
                          #messes up the dim(u) returned in the list
                         if(is.double(y)) y else as.double(y),	     
                         if(is.double(A)) A else as.double(A),  
                         if(is.double(B)) B else as.double(B),   
                         if(is.double(C)) C else as.double(C),
                         if(is.double(TREND)) TREND else as.double(TREND),
                         as.integer(is),  # scratch array dim
                         matrix(double(1),is,is),  # scratch array
                         matrix(double(1),is,is),  # scratch array
                         double(is),         # scratch array
                         integer(is*is),         # scratch array IPIV
                         PACKAGE="dse"
			 ) [c("pred", "weighted.sqerror")]
   if (all(0==error.weights)) r$weighted.sqerror <- NULL
  }
else   # start S version
  {prederror <- matrix(0,sampleT,p)  
   wt.err <- NULL
   invB0 <- solve(B[1,,])
   for (l in 1:a) A[l,,] <- invB0 %*% A[l,,]      # set B(0) = I      
   for (l in 1:b) B[l,,] <- invB0 %*% B[l,,]  
   if (m!=0) for (l in 1:dim(C)[1]) C[l,,] <- invB0 %*% C[l,,]  
   if(!is.null(TREND)) TREND <- t(invB0 %*% t(TREND))
   if (1 < length(error.weights)) wt.err <- matrix(0,predictT,p)    
   for (Time in 1:sampleT)  
      {if(!is.null(TREND)) vt <- -TREND[Time,]
       else vt    <-  rep(0,p) 
       for (l in 1:a)
          if (l<=Time)  #this is cumbersome but drop=FALSE leaves A,B,C as 3 dim arrays
              if (p==1) vt <- vt + c(A[l,,]  *  y[Time+1-l,])
              else      vt <- vt + c(A[l,,] %*% y[Time+1-l,])  
       if (b >= 2) for (l in 2:b) 
          if (l<=Time) 
             if (p==1) vt <- vt - c(B[l,,]  *  prederror[Time+1-l,]) 
             else      vt <- vt - c(B[l,,] %*% prederror[Time+1-l,])
       if (m!=0) for (l in 1:dim(C)[1]) 
          if (l<=Time) 
             if (m==1) vt <- vt - c(C[l,,]  *  u[Time+1-l,])
             else      vt <- vt - c(C[l,,] %*% u[Time+1-l,])  
       prederror[Time,] <- vt #this is not really the pred error unless B0=I   

       if (any(0!=error.weights))          
        {wt.err[Time,] <- error.weights[1]*vt^2  # weighted sq prediction error
         if (length(error.weights)>1)
           {for (h in 2:length(error.weights))
            if ( (Time+h-1) <= sampleT)
              {if(!is.null(TREND)) vt <- -TREND[Time,]
               else vt    <-  rep(0,p) 
               for (l in 1:a)
                  if (l < Time+h) 
                     if (p==1) vt <- vt + c(A[l,,]  *  y[Time+h-l,])
                     else      vt <- vt + c(A[l,,] %*% y[Time+h-l,])  
               if (b >= 2) for (l in 2:b) 
                  if (l < Time+h) 
                     if (p==1) vt <- vt - c(B[l,,]  *  prederror[Time+h-l,]) 
                     else      vt <- vt - c(B[l,,] %*% prederror[Time+h-l,])
               if (m!=0) for (l in 1:dim(C)[1]) 
                   if (l < Time+h) 
                      if (m==1) vt <- vt - c(C[l,,]  *  u[Time+h-l,])
                      else      vt <- vt - c(C[l,,] %*% u[Time+h-l,]) 
               wt.err[Time,] <- wt.err[Time,] + error.weights[h]*(solve(invB0 )%*%vt)^2
     }  }  }
  }

   pred <- matrix(0,predictT,p)
   pred[1:sampleT,] <- y[1:sampleT,,drop=FALSE] - prederror[1:sampleT,,drop=FALSE] %*% t(solve(invB0)) 
   # now multi-step predictions to predictT
   if (predictT > sampleT)
    {for (Time in (sampleT+1):predictT)  
      {invA0 <- solve(A[1,,])
       for (l in 1:a) A[l,,] <- invA0 %*% A[l,,]      # set A(0) = I      
       for (l in 1:b) B[l,,] <- invA0 %*% B[l,,]  
       if (m!=0) for (l in 1:dim(C)[1]) C[l,,] <- invA0 %*% C[l,,]  
       if(!is.null(TREND)) TREND <- t(invA0 %*% t(TREND))
       if(!is.null(TREND)) pred[Time,] <- pred[Time,]+TREND[Time,]
       if (a >= 2) for (l in 2:a) 
          if(Time+1-l<=sampleT)
             if (p==1) pred[Time,] <- pred[Time,]-c(A[l,,]  *  y[Time+1-l,]) 
             else      pred[Time,] <- pred[Time,]-c(A[l,,] %*% y[Time+1-l,])
          else                
             if (p==1) pred[Time,] <- pred[Time,]-c(A[l,,]  *  pred[Time+1-l,])
             else      pred[Time,] <- pred[Time,]-c(A[l,,] %*% pred[Time+1-l,])
       if (b >= 2) for (l in 2:b) 
          if (Time+1-l<=sampleT)
             if (p==1)  pred[Time,] <- pred[Time,] +c(B[l,,]  *  prederror[Time+1-l,]) 
             else       pred[Time,] <- pred[Time,] +c(B[l,,] %*% prederror[Time+1-l,])
       if (m!=0) for (l in 1:dim(C)[1]) 
          if (m==1) pred[Time,] <- pred[Time,] + c(C[l,,]  *  u[Time+1-l,]) 
          else      pred[Time,] <- pred[Time,] + c(C[l,,] %*% u[Time+1-l,])
      }
    }
   r<-list(pred=pred,   weighted.sqerror=wt.err)
  } # end of S version

tframe(r$pred) <- tf
if (! is.null(r$weighted.sqerror))  tframe(r$weighted.sqerror) <- tf
if((!is.null(result)) && (result == "pred")) return(r$pred)
r <- append(residualStats(r$pred, y, sampleT, warn=warn), 
        list(error.weights=error.weights, weighted.sqerror=r$weighted.sqerror))

if (return.debug.info) 
    r$debug.info <-list(m=m, p=p, a=a,b=b,c=c, A=A,B=B,C=C,TREND=TREND, 
          u=u, y=y,
          prederror=prederror, pred=pred, invB0=invB0, wt.err=wt.err, 
          error.weights=error.weights, sampleT=sampleT)

if ( is.null(result)) return(classed( # TSestModel constructor (l.ARMA)
              list(estimates=r, data=dat, model=model), "TSestModel"))
else 
   {if (result =="like") return(r$like[1]) # neg.log.like. from residualStats
    else { return(r[[result]]) }
   }
stop("should never get to here in l.ARMA.")
}                       


l.SS <- function(obj1, obj2, sampleT=NULL, predictT=NULL, error.weights=0,
                 return.state=FALSE, return.track=FALSE, result=NULL, 
		 compiled=.DSEflags()$COMPILED,
                 warn=TRUE, return.debug.info=FALSE, ...)
{#  (... further arguments, currently disregarded)
 # see  help("l.SS") and help("SS")
# could check that Q is symmetric  or positive definite but ...

model <- if (is.TSestModel(obj1)) TSmodel(obj1) else obj1
if (!is.SS(model)) stop("l.SS expecting an SS TSmodel.")

#data <- freeze(obj2)
data <- obj2
if(!checkConsistentDimensions(model, data)) stop("dimension error.\n")
if (is.null(sampleT))  sampleT  <- TobsOutput(data)
if (is.null(predictT)) predictT <- sampleT
if (sampleT > predictT) stop("sampleT cannot exceed predictT.\n")
if (sampleT > TobsOutput(data))
    stop("sampleT cannot exceed length of data.\n")
if (0 != nseriesInput(data))
  {if (TobsInput(data) < predictT)
      stop("input data must be at least as long as requested prediction.\n")
   if (any(is.na(inputData(data)))) stop("input data cannot contain NAs.\n")
  }
if (any(is.na(outputData(data)))) stop("output data cannot contain NAs.\n")

Innov <- is.innov.SS(model)
if (Innov & return.track) 
   warning("Tracking error is zero for an innovations model. track will not be calculated.")

FF<-    model$F
H <-    model$H
n <- dim(FF)[2]
p <- dim(H)[1]
y <- outputData(data)
m <- if(is.null(model$G)) 0 else dim(model$G)[2]

if (m != 0)
  {G <-model$G
   u <- inputData(data)
  } 
if (Innov)            # K or Q,R can be NUll in model, which messes up compiled
   {K <-    model$K
    Q <-    matrix(0,1,1)      #not used
    R <-    matrix(0,1,1)      #not used
    track <-array(0,c(1,1,1))  #not used
   }
else
   {Q <-    model$Q
    if (ncol(Q)<n) Q <- cbind(Q,matrix(0,n,n-ncol(Q))) # Q assumed square in compiled
    R <-    model$R
    K <-    matrix(0,n,p)      # this is used
    if(return.track) track <-array(0,c(predictT,n,n))
    else             track <-array(0,c(1,1,1))  #not used
   }
if (return.state | return.debug.info) state <- matrix(double(1),predictT,n)
else              state <- matrix(double(1),1,1)    #not used
if(is.null(model$z0)) z <-rep(0,n)   # initial state
else  z <-model$z0

# initial state tracking error  --  not used in innov. models
 P <-  if(!is.null(model$rootP0)) t(model$rootP0) %*% model$rootP0
       else if(!is.null(model$P0))  model$P0
       else   diag(1,n)

if (compiled)
  {if (m == 0) {
      G<-matrix(0,n,1)       # can't call compiled with 0 length arrays
      u <- matrix(0,predictT,1)
      } 
   else {# since input and output are declared same dim in fortran
      if (Tobs(y) < Tobs(u)) y <- 
	   rbind(y, matrix(0, Tobs(u) - Tobs(y), nseries(y)))
      }
#   storage.mode(error.weights)     <- "double"
#   storage.mode(state) <- "double"
#   storage.mode(track) <- "double"
#   storage.mode(u)     <- "double"
#   storage.mode(outputData(data))     <- "double"
#   storage.mode(FF)     <- "double"
#   storage.mode(G)     <- "double"
#   storage.mode(H)     <- "double"
#   storage.mode(K)     <- "double"
#   storage.mode(Q)     <- "double"
#   storage.mode(R)     <- "double"
#   storage.mode(z)     <- "double"
#   storage.mode(P)     <- "double"
IS <- max(n,p)
   #r <- .Fortran(kf_sym,
   r <- .Fortran("kf",
                  pred=matrix(double(1),predictT,p), #note, as.double messes dim of pred    
                  as.integer(length(error.weights)), 
                  weighted.sqerror=matrix(0,sampleT,p),
                  error.weights=if(is.double(error.weights))
		                 error.weights else as.double(error.weights),   
                  as.integer(return.state),
                  state=state,         
                  as.integer(return.track & !Innov),
                  track=if(is.double(track)) track else as.double(track),                  
                  as.integer(m), 
                  as.integer(n), 
                  as.integer(p), 
                  sampleT=as.integer(sampleT), 
                  predictT=as.integer(predictT), 
                  as.integer(Tobs(y)),  
                  if(is.double(u)) u else as.double(u), 
                  if(is.double(y)) y else as.double(y),  
                  if(is.double(FF)) FF else as.double(FF),   
                  if(is.double(G)) G else as.double(G),	
                  if(is.double(H)) H else as.double(H),  
                  if(is.double(K)) K else as.double(K), 
                  if(is.double(Q)) Q else as.double(Q),	   
                  if(is.double(R)) R else as.double(R),	 
                  as.integer(Innov),
                  if(is.double(z)) z else as.double(z),
                  if(is.double(P)) P else as.double(P),
	          as.integer(IS),           # scratch arrays for KF, IS
	          matrix(double(1),IS,IS),  #A
	          matrix(double(1),IS,IS),  # AA
	          matrix(double(1),IS,IS),  # PP
	          matrix(double(1),n,n),  # QQ
	          matrix(double(1),p,p),  # RR 
	          rep(double(1),IS),  # Z
	          rep(double(1),IS), # ZZ
	          rep(double(1),IS), # WW		   
                  integer(IS*IS),         # scratch array IPIV
		  PACKAGE="dse"
		  ) [c("pred","state","track","weighted.sqerror")]
   if (all(0==error.weights)) r$weighted.sqerror <- NULL
  }
else                  #  S version
  {vt    <-  rep(0,p)     # initial prediction error
   pred  <- matrix(0,predictT,p) 
   wt.err <- NULL
   if (1 < length(error.weights)) wt.err <- matrix(0,predictT,p)

## This is the beginning of the  most important part of the algorithm ####

   if ( ! Innov ) 
     {RR <- R %*% t(R)  
      QQ <- Q %*% t(Q)  
     }                            
                                   
   for (Time in 1:sampleT)  {
       if ( ! Innov) 
         {PH  <-  P %*% t(H)
          ft    <- ( H %*% PH )  + RR         
          ft    <-  (ft + t(ft))/2   # force ft to be symmetric 
          K   <-  t(solve(ft,t(FF %*% PH)))  
          P   <-  (FF %*% P %*% t(FF) ) - ( K %*% H %*% P %*% t(FF) ) + QQ  # P(t|t-1)
          P   <-  (P + t(P))/2  # force symmetry (eliminate rounding error problems)
          if (return.track) track[Time,,] <- P   # P(t|t-1)
          # note P(t|t) = P-P%*%t(H)%*%solve(H%*%t(H)+RR)%*%P
         } 
         
       z <- c(FF%*%z) + c(K%*%vt)  # E[z(t)| t-1 ]
       if (m !=0) z <- z + c(G%*%u[Time,])
       if (return.state | return.debug.info) state[Time,] <- z
       pred[Time,] <- Ey  <-  c(H %*% z)       # predicted output     
       vt<-  y[Time,] - Ey                     # prediction error 

## This is the end of the most important part of the algorithm ####

       if (any(0!=error.weights))          
        {wt.err[Time,] <- error.weights[1]*vt^2  # weighted sq prediction error
         if (length(error.weights)>1)
          {zh <- z
           for (h in 2:length(error.weights))
            if ( (Time+h-1) <= sampleT)
              {zh <-  c(FF%*%zh)
               if (h==2)  zh <- zh + c(K%*%vt) # vt is 0 for h>2
               if (m !=0) zh <- zh + c(G%*%u[Time+h-1,])
               wt.err[Time,] <- wt.err[Time,] + 
                          error.weights[h]*(y[Time+h-1,] -  c(H %*% zh))^2
   }   }  }   }


#   prederror <- y[1:sampleT,,drop=FALSE]-pred[1:sampleT,,drop=FALSE]

   # now multi-step prediction to predictT 
   # This requires u but not y (y is ignored if it is supplied)
   if (predictT > sampleT)
    {for (Time in (sampleT+1):predictT)  
       {z <-  c(FF%*% z) 
        if (m !=0) z <-  z + c(G%*%u[Time,]) 
        if (Time==sampleT+1) z <- z + c(K%*%vt)
        if (return.state) state[Time,] <- z
        pred[Time,]  <-  c(H %*% z)                  # predicted output 
       }
     }
   r<- list(pred=pred, state=state, track=track, weighted.sqerror=wt.err) 
  }   # end of S version

tf <- tfTruncate(tframe(outputData(data)), end=predictT)
tframe(r$pred) <- tf
if (! is.null(r$weighted.sqerror))  tframe(r$weighted.sqerror) <- tf

filter <-NULL
if (return.track) filter$track <- if (Innov) NULL else tframed(r$track, tf)
if (return.state) filter$state <- tframed(r$state,tf)
 
if((!is.null(result)) && (result == "pred")) return(r$pred)
r <- append(residualStats(r$pred, outputData(data), sampleT, warn=warn), 
        list(error.weights=error.weights, weighted.sqerror=r$weighted.sqerror))

if (return.debug.info) 
   r$debug.info <- list(m=m, p=p, c=c, G=G,FF=FF,K=K,P=P, H=H, u=u,y=y,
        pred=pred, wt.err=wt.err, error.weights=error.weights, sampleT=sampleT)

if ( is.null(result)) return( 
    classed(list(estimates=r, data=data, model=model, filter=filter),
            "TSestModel")) # constructor (l.SS)
else if (result =="like") return(r$like[1]) # neg.log.like. from residualStats
else return(r[[result]]) 
  
"should never get to here in l.SS"
} # end of l.SS



smoother <- function(model, data, compiled=.DSEflags()$COMPILED) UseMethod("smoother")

smoother.TSestModel <- function(model, data=TSdata(model),
                                compiled=.DSEflags()$COMPILED){
    r <- smoother(TSmodel(model), data=data, compiled=compiled)
    r$estimates <- model$estimates # otherwise convergence info is lost
    r
    }
    
smoother.TSmodel <- function(model, data, compiled=.DSEflags()$COMPILED){
    if  (!is.nonInnov.SS(model))
       stop("smoothing needs an nonInnov state space model (class nonInnov.SS.TSmodel).")
    # should neot get to here. (smoother.nonInnov should have been called
    # instead, but ...)
    smoother.nonInnov(model, data, compiled=compiled)
    }

smoother.nonInnov <- function(model, data, compiled=.DSEflags()$COMPILED){
 # See help("smoother") and help("SS") for details of the model: 
 data <- TSdata(data)
 z <-l(model, data, return.state=TRUE,return.track=TRUE, compiled=compiled)
 filter    <- z$filter
 estimates <- z$estimates
 if  (is.null(filter$state) |  is.null(filter$track)) 
     filter <- l(model,data, return.state=TRUE,return.track=TRUE)$filter
 model     <- TSmodel(model)
 if  (!is.nonInnov.SS(model)) stop("smoother expecting an nonInnov.SS TSmodel.")
 if (is.null(model$G))
  {m<-0
   G<-matrix(0,dim(model$F)[2],1)   # can't call compiled with 0 length arrays
   u <- matrix(double(1),nrow(filter$state),1)
  }
 else
  {m <- dim(model$G)[2]
   G <-model$G
   u <- inputData(data)
   if (! is.double(u)) u <- array(as.double(u), dim(u))                  
  } 
 sampleT  <-min(nrow(u), nrow(filter$state), dim(filter$track)[1])
 QQ <- model$Q %*% t(model$Q)           
 RR <- model$R %*% t(model$R) 
 n <- dim(model$F)[2]          
 p <- dim(model$H)[1]
 if (compiled)
   {#storage.mode(filter$state)     <- "double"
    #storage.mode(filter$track)     <- "double"
    #storage.mode(u)     <- "double"
    #storage.mode(outputData(data))     <- "double"
    #storage.mode(model$F)     <- "double"
    #storage.mode(G)     <- "double"
    #storage.mode(model$H)     <- "double"
    #storage.mode(RR)     <- "double"
    
    #note, as.double messes dim needed in result, but need attributes removed
    IS <- max(n,p)

    #r <- .Fortran(smooth_sym, 
    r <- .Fortran("smooth", 
                         state=array(as.double(filter$state), dim(filter$state)),      
                         track=array(as.double(filter$track), dim(filter$track)),     
                         u,	       # input
                         if(is.double(outputData(data))) outputData(data) else 
			        as.double(outputData(data)),
			 as.integer(n),
                         as.integer(m),
                         as.integer(p),  
                         sampleT =as.integer(sampleT), 
                         as.double(model$F),   
                         as.double(G),	
                         as.double(model$H),  
                         as.double(RR),	 
                         as.integer(IS),           # scratch array dim IS
			 matrix(double(1),IS,IS),   # scratch array A
                         matrix(double(1),IS,IS),   # scratch array D
                         matrix(double(1),IS,IS),   # scratch array L
                         matrix(double(1),IS,IS),   # scratch array PT1
                         double(IS),         # scratch array ZT
                         integer(IS*IS),         # scratch array IPIV
                         #chk1=array(double(1), c(sampleT, IS,IS)),     
                         #chk2=array(double(1), c(sampleT, IS,IS)),     
                         PACKAGE="dse"
			 )[c("state","track")]
   }
 else   # S version
   {FF<-  model$F
    H <-  model$H
    #   filter$state is the one step ahead state estimate E{z(t)| y(t-1), u(t)}
    #   zt below is the filter state estimate E{z(t)| y(t), u(t+1)}
    #   filter$track is tracking error P(t|t-1)
    # A5-A7 just use the filter information , which was calculated 1:sampleT.
    # A8-A10 calculate smooth state and track, which is done in reverse sampleT:1 
    sm <- array(NA,dim(filter$state))  # smoother state estimate
    sm[sampleT,] <- filter$state[sampleT,]
    tr <- array(NA,dim(filter$track))  # smoother tracking error
    tr[sampleT,,] <- filter$track[sampleT,,]
    for (Time in (sampleT-1):1) 
      {#K <- filter$track[Time,,] %*% t(H) %*% 
       #       solve(H %*% filter$track[Time,,] %*% t(H) + RR) #(A5)
       K <- filter$track[Time,,] %*%  
              t(solve(t(H %*% filter$track[Time,,] %*% t(H) + RR), H)) #(A5)
       zt <- filter$state[Time,] + K %*% 
       		(outputData(data)[Time,] - H %*% filter$state[Time,]) #(A6)
       #if (m!=0) zt <- zt + G %*% u[Time,] 
       P <- filter$track[Time,,] - K %*% H %*% filter$track[Time,,] #P(t|t) (A7)
       P <- (P+t(P))/2 #force symmetry to avoid rounding problems
       #J <- P %*% t(FF) %*% solve(filter$track[Time+1,,])             #(A8) 
       J <- P %*%  t(solve(t(filter$track[Time+1,,]), FF))             #(A8) 
       if (m==0)                                          
          sm[Time,] <- zt + J %*% (sm[Time+1,] - FF %*% zt)              #(A9) 
       else 
          sm[Time,] <- zt + J %*% (sm[Time+1,] - FF %*% zt - G %*% u[Time+1,]) #check   
       P <- P + J %*% (tr[Time+1,,] - filter$track[Time+1,,]) %*% t(J)    #(A10)
       tr[Time,,]<- (P+t(P))/2   #force symmetry to avoid rounding problems
      }
     r <- list(state=sm, track=tr) 
    }  # end S version
 
  tf <- tframe(outputData(data))
  names <- seriesNames(filter$state) 
  classed(list(estimates=estimates, data=data, model=model, 
               filter=filter, 
	       smooth=list(state=tframed(r$state, tf, names=names),
	                   track=tframed(r$track, tf))),
       "TSestModel") # constructor (smoother)
} # end of smoother
  

state <- function(obj, smoother=FALSE, filter=!smoother) {
     if (!inherits(obj,"TSestModel"))
        stop("The argument needs to be a TSestModel from an SS model.")
     if(filter & smoother) stop("only one of filter and smoother can be specified.")
     if(filter) if(!is.null(obj$filter$state)) return(obj$filter$state)
     else return(l(TSmodel(obj), TSdata(obj), return.state=TRUE)$filter$state)
     if (!is.null(obj$smooth)) return(obj$smooth$state)
     if (!inherits(TSmodel(obj),"nonInnov"))
          stop("A smoother state cannot be calculated . The argument needs to be a TSestModel from an nonInnov SS model.")
     state(smoother(obj), smoother=TRUE)
     }
 
#trackingError <- function(obj, ... should do this as for state 



############################################################

#     parameter estimation functions   <<<<<<<<<<

############################################################


estVARXls <- function(data, subtract.means=FALSE, re.add.means=TRUE,
   standardize=FALSE, unstandardize=TRUE, max.lag=NULL, 
   trend=FALSE, lag.weight=1.0, warn=TRUE) 
{# Estimate VAR model with exogenous variable using lsfit(). 
 # Returns a TSestModel.
 # This is very similar to estimating Markov parameters (see estSSMittnik).
 # Residuals,etc, are calculated by evaluating the estimated model with ARMA.
 # Data should be of class TSdata.
 # lag.weight is an exponential weight applied to lags. It should be in (0,1].
   if (is.null(max.lag)) max.lag <- 6
   #data <- freeze(data)
   names <- seriesNames(data)
   m <-  nseriesInput(data)
   p <- nseriesOutput(data)
   if (0 == p) stop("estVARXls requires output data to estimate a model.")
   missing.data <- any(is.na(outputData(data)))
   if (0 != m) { # ensure same time window
	inputData(data) <- tfwindow(inputData(data),
              start=startOutput(data), end=endOutput(data), warn=FALSE)
	outputData(data) <- tfwindow(outputData(data),
              start=startInput(data), end=endInput(data), warn=FALSE)
	missing.data <- any(missing.data, is.na(inputData(data)))
   }
   N <- TobsOutput(data)
   if (standardize)
     {svd.cov <- La.svd(var(outputData(data)))
      scalefac <- svd.cov$u %*% diag(1/svd.cov$d^.5, ncol=p)
      data <- scale(data, scale=list(output=scalefac))
     }
   if (subtract.means)
    {if(m!=0)
       {input.means<-colMeans(inputData(data))
        inputData(data)<-inputData(data)-t(matrix( input.means, m,N))
       }
     output.means <- colMeans(outputData(data))
     outputData(data)  <- outputData(data) - t(matrix(output.means, p,N))
    }
      # shift input to give feedthrough in one period
   if (m != 0) {z <- cbind(inputData(data)[2:N,],outputData(data)[1:(N-1),])}
   else z <- outputData(data)

 # The matrix Past is blocks of data:
 #  [in | out-1 | in-1 | out-2 | ... | in-max.lag | out-max.lag-1 ]
 # so the coef. matrix M has a corresponding structure.
 
   Past <- matrix(NA,N-max.lag,(p+m)*(max.lag))
   for (i in 0:(max.lag-1)) 
      Past[,(1+(m+p)*i):((m+p)*(1+i))] <-z[(max.lag-i):(N-1-i),] /(lag.weight^i)
   fit <- lsfit(Past,outputData(data)[(max.lag+1):N,,drop=FALSE],intercept=trend)
   if(missing.data)
     fit.res <-TSdata(output=fit$residual) #used only in the case of missing data for scaling
   M <- t(fit$coef)
   # correct exogenous blocks for normalization:
   if (standardize && (m!=0)) 
     {Tinv <- diag(svd.cov$d^.5, ncol=p)%*%svd.cov$u
      for (i in 0:(max.lag-1)) 
         M[,(1+(m+p)*i):(m+(m+p)*i)] <- Tinv %*% M[,(1+(m+p)*i):(m+(m+p)*i)]
     }
   TREND <- NULL
   if (trend)
     {TREND <- M[,1,drop=FALSE]
      M<-M[,2:(dim(M)[2]),drop=FALSE]
     }
   if (subtract.means & re.add.means)
    {if(m!=0) 
       {inputData(data)<-inputData(data) + t(matrix( input.means, m,N))
        Past<-Past +
           t(matrix(c(input.means,output.means),(p+m)*(max.lag),N-max.lag))
       }
     else
        Past <-Past+t(matrix(output.means, (p+m)*(max.lag),N-max.lag))
     outputData(data)  <- outputData(data) + t(matrix(output.means, p,N))
     # and correct estimation for non-zero mean:
     M<-M+estVARXmean.correction(Past, 
                       outputData(data)[(max.lag+1):N,,drop=FALSE],M, warn=warn)
    }
   A <- array(NA, c(1 + max.lag, p, p))
   A[1,,] <-diag(1, p)
   for (i in 0:(max.lag-1)) 
          A[2+i,,] <-  -M[,(m+1+(m+p)*i):(m+p+(m+p)*i),drop=FALSE] %*% 
                                 diag(lag.weight^i,p)
   if (m==0) C <- NULL   # no exog. variable
   else               # NB. there is an implicit shift in inputData(data)
      {C <-array(NA,c(max.lag,p,m))
       for (i in 0:(max.lag-1)) 
         C[1+i,,] <- M[,(1+(m+p)*i):(m+(m+p)*i),drop=FALSE]  %*% 
                                 diag(lag.weight^i,m)
      }
   B <- array(diag(1,p),c(1,p,p))
   model <-ARMA( description="model estimated by estVARXls",
              A=A,B=B,C=C,TREND=TREND)
   seriesNames(model) <- seriesNames(data)

   if (standardize & unstandardize)
     {scalefac <- solve(scalefac)
      data <- scale(data,  scale=list(output=scalefac))
      model<- scale(model, scale=list(output=scalefac))
      if(missing.data) fit.res <- scale(fit.res, scale=list(output=scalefac))
     }

   if(missing.data)
       {if (warn) warning(
          "missing data kludge. Predictions are reconstructed from lsfit residuals")
        return(fake.TSestModel.missing.data(model,data, fit.res$output,max.lag))
       }
   else return(l(model, data, warn=warn))
}

fake.TSestModel.missing.data <- function(model,data, residual, max.lag, warn=TRUE)
  {pred <- rbind(matrix(NA,max.lag,nseriesOutput(data)),residual) + 
                                                           outputData(data)
   r <- list(estimates = residualStats(pred, outputData(data), warn=warn),
               data = data, model = model)
   classed(r, "TSestModel") # fake missing data constructor
  }

estVARXmean.correction <- function(X, y, bbar,
                     fuzz=sqrt(.Machine$double.eps), warn=TRUE)
{# correction for model estimated with means subtracted
 Xbar <- t(array(colMeans(X), rev(dim(X))))
 ybar <- t(array(colMeans(y), rev(dim(y))))
 v <- La.svd(t(X)%*%X) # this is more robust than solve()
 if (warn && any(abs(v$d[1]*fuzz) > abs(v$d) ) ) 
   warning("The covariance matrix is nearly singular. Check for linearly related data.")
# if(1 == length(v$d))OmInv <- v$u %*% (1/v$d) %*% v$vt
# else OmInv <- v$u %*% diag(1/v$d) %*% v$vt	
#  following is equivalent
# OmInv <-  t(v$vt) %*% sweep(t(v$u),1,1/v$d, "*") CHECK
 OmInv <-  crossprod(v$vt, sweep(t(v$u),1,1/v$d, "*")) 

# t(-OmInv %*% ( 
#        (2*t(Xbar)%*%X - t(Xbar)%*%Xbar) %*% t(bbar) 
#         - t(Xbar)%*%y - t(X)%*%ybar + t(Xbar)%*%ybar ))
 t(-OmInv %*% ( 
        (2*crossprod(Xbar, X) - crossprod(Xbar)) %*% t(bbar) 
         - crossprod(Xbar, y) - crossprod(X,ybar) + crossprod(Xbar, ybar)))
}


estVARXar <- function(data, subtract.means=FALSE,  re.add.means=TRUE,
      standardize=FALSE, unstandardize=TRUE, aic=TRUE, max.lag=NULL, 
      method="yule-walker", warn=TRUE) 
{
   #data <- freeze(data)
   m <-  nseriesInput(data)
   p <- nseriesOutput(data)
   if (0 == p) stop("estVARXar requires output data to estimate a model.")
   N <- TobsOutput(data)
   if (standardize)
     {svd.cov <- La.svd(var(outputData(data)))
      scalefac <- svd.cov$u %*% diag(1/svd.cov$d^.5, ncol=p)
      data <- scale(data, scale=list(output=scalefac))
     }
   if(m!=0) input.means <- colMeans(inputData(data))
   output.means <- colMeans(outputData(data))
   if (subtract.means)
    {if(m!=0) inputData(data) <- inputData(data)-t(matrix( input.means, m,N))
     outputData(data)  <- outputData(data) - t(matrix(output.means, p,N))
    }
   if (m==0)  zdata <- outputData(data)  # no exog. variable
   else       zdata <- cbind(inputData(data),outputData(data))
         # NB. there is an implicit shift in inputData(data) in the line above
         #   because ar estimates lag parameters,
         # so C[1,,] is for one lag in u.
   if (is.null(max.lag))  AC<- DSE.ar(zdata, aic=aic, method=method)
   else                   AC<- DSE.ar(zdata, aic=aic, method=method, order.max=max.lag)
   if (re.add.means)
     {if (subtract.means)
        {if(m!=0) inputData(data)<-inputData(data) + t(matrix( input.means, m,N))
         outputData(data)  <- outputData(data) + t(matrix(output.means, p,N))
        }
      # and correct estimation for non-zero mean:
      max.lag <- AC$order
      if (max.lag==0) stop("all lags eliminated by AIC order selection.")
      if (m==0)  zdata <- outputData(data)  # no exog. variable
      else       zdata <- cbind(inputData(data),outputData(data))
      Past <- matrix(NA,N-max.lag,(p+m)*(max.lag))
      for (i in 0:(max.lag-1)) 
         Past[,(1+(m+p)*i):((m+p)*(1+i))] <-zdata[(max.lag-i):(N-1-i),]
      M <-estVARXmean.correction(Past, zdata[(max.lag+1):N,,drop=FALSE],
                  matrix(aperm(AC$ar, c(2,3,1)),m+p,(m+p)*max.lag), warn=warn)
      AC$ar<-AC$ar + aperm(array(M,c(m+p,m+p,max.lag)), c(3,1,2))
     }
   A <- array(0, c(1 + AC$order, p, p))
   A[1,,] <-diag(1, p)
   if (0==AC$order)
      warning("lagged output variables eliminated by AIC order selection.")
   else 
      A[2:(1+AC$order),,] <-  -AC$ar[, (m+1):(m+p), (m+1):(m+p), drop=FALSE]
   if (m==0) C <- NULL
   else 
     {if (0==AC$order)
             {warning("input variables eliminated by AIC order selection.")
              C <-array(0, c(1,p,m))
             }
      else C <-array(AC$ar[,(m+1):(m+p),1:m],c(AC$order,p,m))
     }
   B <- array(diag(1,p),c(1,p,p))
   model <-ARMA(
            description="model estimated by estVARXar",
            A=A,B=B,C=C, 
            names = seriesNames(data) )
 
   if (standardize & unstandardize)
     {scalefac <- solve(scalefac)
      data <- scale(data,  scale=list(output=scalefac))
      model<- scale(model, scale=list(output=scalefac))
     }
   l(model, data, warn=warn)
}

estSSfromVARX <- function(data, warn=TRUE, ...){
	# estimate a VARX model and convert to state space
	#  estimate a nested-balanced state space model by svd a la Mittnik from
	#   least squares estimate of VAR coefficents.
	model <-estVARXls(data, warn=warn, ...)
	l(toSS(model$model), model$data, warn=warn)
}



############################################################################

#  functions for model estimation (see also VARX ) and reduction   <<<<<<<<<<<<<

############################################################################


estWtVariables <- function(data, variable.weights,
                        estimation="estVARXls", estimation.args=NULL)
{ if (is.matrix(variable.weights))
    {if (any(La.svd(variable.weights)$d == 0))  
       stop("variable.weights transformation must be non singular.")
    }
 else
   {if (any(variable.weights == 0))  stop("variable.weights must be non zero.")
    variable.weights <- diag(variable.weights)
   }
 inv.wts <- solve(variable.weights)
 dimnames(inv.wts)          <-list(NULL, seriesNamesOutput(data))
 dimnames(variable.weights) <-list(NULL, seriesNamesOutput(data))
 scaled.model <- do.call(estimation, append(list(
                  scale(data, scale=list(output=inv.wts))),  estimation.args))
#           freeze(scale(data, scale=list(output=inv.wts)))), estimation.args))
 model <-scale(TSmodel(scaled.model), scale=list(output=variable.weights))
 model$description <- 
    paste("model estimated by estWtVariables with", estimation)
 l(model, data)
}

estMaxLik <- function(obj1, obj2=NULL, ...) UseMethod("estMaxLik")

estMaxLik.TSestModel <- function(obj1, obj2=TSdata(obj1), ...) {
	 # if obj1 is result from a previous estMaxLik then the gradient
	 # hessian and other information should be extracted, but
	 estMaxLik(TSmodel(obj1), obj2, ...) }

estMaxLik.TSdata <- function(obj1, obj2, ...) 
	 estMaxLik(TSmodel(obj2), obj1, ...) 

estMaxLik.TSmodel <- function(obj1, obj2, 
	algorithm="optim",
	algorithm.args=list(method="BFGS", upper=Inf, lower=-Inf, hessian=TRUE),
	...)
{# maximum likelihood estimation
 # "nml" algorithm.args=list(hessian=T, iterlim=20, 
 #     dfunc=gradNumerical, line.search="nlm",ftol=1e-5, gtol=1e-3,)
 Shape <- obj1
 #data <- freeze(obj2)
 data <- obj2
 func.like <- function(coefficients, Shape,data)
      {l(setArrays(Shape, coefficients=coefficients), data, result="like",
         warn=FALSE) }

 if (algorithm=="optim")
    {results <- optim(coef(Shape), func.like, method=algorithm.args$method,
	gr=algorithm.args$gr, 
	lower=algorithm.args$lower, upper=algorithm.args$upper,
	control=algorithm.args$control, hessian=algorithm.args$hessian,
	Shape, data) 
     parms <- results$par
     converged <- results$convergence == 0
     convergenceCode <- results$convergence
     description <- paste("Estimated with max.like/optim")
    }
 else if (algorithm=="nlm")
   {warning("This has not been tested recently (and there have been changes which may affect it.")
    results <-nlm(func.like, coef(Shape), hessian=algorithm.args$hessian, 
    	iterlim=algorithm.args$iterlim)
    parms <- results$estimate
    converged <- results$code <= 2
    convergenceCode <- results$code
    description <- paste("Estimated with max.like/nlm")
   }
  else stop(paste("Minimization method ", algorithm, " not supported."))
# else if (algorithm=="nlmin")
#   {warning("This has not been tested recently (and there have been changes which may affect it.")
#     results <-nlmin(func.like, coef(Shape), max.iter=algorithm.args$max.iter, 
#     	max.fcal=5*algorithm.args$max.iter, ckfc=0.01)
#     parms <- results$parms
#     converged <- results$converged
#     convergenceCode <- results$converged
#     # emodel$est$converged should be improved with conv.type info
#     description <- paste("Estimated with max.like/nlmin")
#   }
# else if (algorithm=="dfpMin")
#    {stop("This optimization method is no longer supported.")
#     results <- dfpMin(func.like, coef(Shape), 
#	dfunc=algorithm.args$dfunc, 
#	max.iter=algorithm.args$max.iter, 
#	ftol=algorithm.args$ftol, gtol=algorithm.args$gtol, 
#	line.search=algorithm.args$line.search) 
#     parms <- results$parms
#     converged <- results$converged
#     convergenceCode <- results$converged
#     description <- paste("Estimated with max.like/dfpMin")
#    }
 emodel <- l(setTSmodelParameters(setArrays(Shape, coefficients=parms)),data)
 emodel$estimates$algorithm <- algorithm
 emodel$estimates$results   <- results
 emodel$estimates$converged <- converged
 emodel$estimates$convergenceCode <- convergenceCode
 emodel$model$description <- paste(description, 
       if(!converged) " (not " else " (", "converged",
       ") from initial model: ", Shape$description)
 emodel
}    


estBlackBox <- function(data,...)
{# call current best black box technique.
  estBlackBox4(data, ...)
}


estBlackBox1 <- function(data,estimation="estVARXls", 
     reduction="MittnikReduction", criterion="taic", trend=FALSE, 
     subtract.means=FALSE, verbose=TRUE, max.lag=6)
{if ((estimation!="estVARXls") && (trend) )
     {cat("Trend estimation only support with estVARXls.\n")
      cat("Proceeding using estVARXls.\n")
      estimation<-"estVARXls"
     }

 if(estimation=="estVARXls")
     model <- estVARXls(data,trend=trend, subtract.means=subtract.means, 
                          max.lag=max.lag)
 else if(estimation=="estVARXar")
     model <- estVARXar(data, subtract.means=subtract.means, max.lag=max.lag)
 else if(estimation=="estVARXls")
     model <- estVARXls(data,trend=trend, subtract.means=subtract.means, 
                        max.lag=max.lag)
 else if(estimation=="estSSMittnik")
     model <- estSSMittnik(data,max.lag=max.lag, 
                             subtract.means=subtract.means, normalize=FALSE)
 else
   stop("estimation technique not supported.")
 if (verbose) 
   cat("First VAR model,              lags= ", dim(model$model$A)[1]-1,
       ", -log likelihood = ", model$estimates$like[1], "\n")
 model <- l(toSS(model),data)
 n <- dim(model$model$F)[1]
 if (verbose) cat("Equivalent    state space model, n= ", n,
                  ", -log likelihood = ", model$estimates$like[1], "\n")
 if (1 < n)
   {model <- do.call(reduction,
                     list(model, criterion=criterion, verbose=verbose))
   #model <- eval(call(reduction,model,criterion=criterion, verbose=verbose))
    if (verbose) 
       cat("Final reduced state space model, n= ", nstates(model),
           ", -log likelihood = ", model$estimates$like[1], "\n")
   }
  if (verbose &&  dev.cur() != 1 ) checkResiduals(model)
 model
}


estSSMittnik <- function(data, max.lag=6, n=NULL, 
    subtract.means=FALSE, normalize=FALSE)
{ #data <- freeze(data)
  m <- ncol(inputData(data))
  if(is.null(m))  m <- 0
  p <- ncol(outputData(data))
  if (0 == p) stop("estSSMittnik requires output data to estimate a model.")
  N <- nrow(outputData(data))
  if (subtract.means)
    {if(m!=0)inputData(data)<-inputData(data)-t(matrix(colMeans(inputData(data)), m,N))
     outputData(data)<- outputData(data) - t(matrix(colMeans(outputData(data)), p,N))
    }
  if (normalize)
    {svd.cov <- La.svd(var(outputData(data)))
     outputData(data) <- outputData(data) %*% svd.cov$u %*% diag(1/svd.cov$d^.5)
    }
      # shift input to give feedthrough in one period
  if (m != 0) {z <- cbind(inputData(data)[2:N,],outputData(data)[1:(N-1),])}
  else z <- outputData(data)
  Past <- matrix(NA,N-1-max.lag,(p+m)*(1+max.lag))
  for (i in 0:max.lag) 
    Past[,(1+(m+p)*i):((m+p)*(1+i))] <-z[(1+max.lag-i):(N-1-i),]
  M <- t(lsfit(Past,outputData(data)[(max.lag+2):N,],intercept=FALSE)$coef)
  if (normalize && (m!=0))  # correct exogenous blocks for normalization
    {Tinv <- diag(svd.cov$d^.5)%*%svd.cov$u
     for (i in 0:max.lag) 
       M[,(1+(m+p)*i):(m+(m+p)*i)] <- Tinv %*% M[,(1+(m+p)*i):(m+(m+p)*i),drop=FALSE]
    }
  if (p==1) M <-matrix(M,1,length(M))
#  browser()
  model <-SVDbalanceMittnik( M, m=m, n=n )$model
  z <-"nested-balanced model from least sq. estimates of markov parameters a la Mittnik"
  if(subtract.means) z <-paste(z," - means subtracted")
  if(normalize)      z <-paste(z," - outputs normalized")
  model$description  <-z
  seriesNames(model)  <- seriesNames(data)
  l(model, data)
}


MittnikReduction <- function(model, data=NULL, criterion=NULL, 
      verbose=TRUE,warn=TRUE)
{if (is.TSestModel(model)) 
    {if (is.null(data)) data <-TSdata(model)
     model <- TSmodel(model)
    }
  if(!is.TSmodel(model)) stop("MittnikReduction expecting a TSmodel.")
  if(is.null(data))
    stop("Reduction requires data to calculate criteria (balancing does not).")
  nMax <- if(is.SS(model)) dim(model$F)[2] else NULL
  MittnikReduction.from.Hankel( markovParms(model), nMax=nMax, data=data, 
        criterion=criterion, verbose=verbose, warn=warn)
}

MittnikReduction.from.Hankel <- function(M, data=NULL, nMax=NULL, 
   criterion=NULL, verbose=TRUE, warn=TRUE) 
#   Spawn=if (exists(".SPAWN")) .SPAWN else FALSE)
{ # See the documentation for MittnikReduction.

   read.int <- function(prmt)
     {err <-TRUE
      while (err)
       {cat(prmt)
 	n <- as.integer(readline()) # crude. this truncates reals
 	if (is.na(n)) cat("value must be an integer\n")
 	else err <- FALSE
       }
      n
     }

   #data <- freeze(data)
   m <-ncol(inputData(data))      # dim of input series
   if(is.null(m))m<-0
   z <-SVDbalanceMittnik(M, m, nMax)
   largeModel <- z$model
   svd.crit    <-z$crit
   n <- dim(largeModel$F)[1]
#   if (!Spawn)
#     {
      values <- NULL 
      for (i in 1:n) 
        {if(m!=0) z <-largeModel$G[1:i,,drop=FALSE]
         else     z <-NULL
         z <-SS(F.=largeModel$F[1:i,1:i,drop=FALSE],G=z,
                  H=largeModel$H[,1:i,drop=FALSE],K= largeModel$K[1:i,,drop=FALSE])
         z <-informationTestsCalculations(l(z,data, warn=warn))
         values <-rbind(values, z)
         if (verbose) cat(".")
        }
 #     }
#   else 
#     {#  loop used below is to avoid S memory problems
#      if (verbose) cat("Spawning processes to calculate criteria test for state dimension 1 to ",n)
#      forloop <- function(largeModel, data, warn=TRUE)
#           {if(!is.null(largeModel$G)) z <-largeModel$G[1:forloop.i,,drop=FALSE]
#            else                       z <-NULL
#            z <-SS(F.=largeModel$F[1:forloop.i,1:forloop.i,drop=FALSE],G=z,
#                     H=largeModel$H[  , 1:forloop.i, drop=FALSE],
#                     K=largeModel$K[1:forloop.i,,drop=FALSE])
#            informationTestsCalculations(l(z,data, warn=warn))
#           }
#	assign("balance.forloop", forloop, where=1)
#	assign("forloop.n", n,   where=1 )
#	assign("forloop.values", matrix(NA,n,12),   where=1 )
#	assign("forloop.largeModel", largeModel, where=1)
#	assign("forloop.data", data,   where=1 )
#	assign("forloop.warn", warn,   where=1 )
#	on.exit(remove(c("balance.forloop", "forloop.i", "forloop.n", 
#	   "forloop.values", "forloop.largeModel", 
#	   "forloop.data", "forloop.warn"),where=1))
#	For (forloop.i=1:forloop.n,
#	  forloop.values[forloop.i,]<-balance.forloop(forloop.largeModel, 
#         forloop.data, forloop.warn), sync=TRUE)
#       values <-forloop.values
#      }
    dimnames(values) <- list(NULL,c("port","like","aic","bic", 
          "gvc","rice","fpe","taic","tbic","tgvc","trice","tfpe")) 
    if (verbose) cat("\n")
    opt <-apply(values,2,order)[1,]  # minimum
    if (verbose | is.null(criterion))
      {zz<-criteria.table.nheading()
       options(width=120)
       print(values,digits=4)
       cat("opt     ")
       for (i in 1:length(opt)) cat(opt[i],"    ")
       cat("\n")
       zz<-criteria.table.legend()
      }
    if (is.null(criterion))
      {n <- read.int("Enter the state dimension (enter 0 to stop): ")
       if( n<1) stop("TERMINATED! STATE DIMENSION MUST BE GREATER THAN 0!")
      }
    else { n <- opt[criterion == dimnames(values)[[2]]]  }
    if(m==0) z <-NULL 
      else   z <-largeModel$G[1:n,,drop=FALSE]
    model <- SS(description="nested model a la Mittnik",
          F.=largeModel$F[1:n,1:n,drop=FALSE],G=z,
          H=largeModel$H[,1:n,drop=FALSE],K= largeModel$K[1:n,,drop=FALSE], 
          names=seriesNames(data))
    l(model,data, warn=warn)          
}



############################################################

#     model balancing functions   <<<<<<<<<<

############################################################

balanceMittnik <- function(model, n=NULL){
  if (is.TSestModel(model)) model <- TSmodel(model)
  if  (!is.TSmodel(model)) stop("balanceMittnik expecting a TSmodel.")
  m <- nseriesInput(model)
  newmodel <- SVDbalanceMittnik(markovParms(model), m, n)$model
  newmodel$description <- paste(model$description,"converted to", newmodel$description)
  seriesNames(newmodel)  <-  seriesNames(model)
  newmodel
}


SVDbalanceMittnik <- function(M, m, n=NULL)
{#                   # Form k block Hankel Matrix from M.

   read.int <- function(prmt)
     {err <-TRUE
      while (err)
       {cat(prmt)
 	n <- as.integer(readline()) # crude. this truncates reals
 	if (is.na(n)) cat("value must be an integer\n")
 	else err <- FALSE
       }
      n
     }

   p <- nrow(M)     # dim of endo. series
   r <- m + p       # each sub-matrix is p x r   (= p x (m+p) )
   k<- dim(M)[2] / r
   Hkl <- matrix(0, p*k, r*k) 
   for(i in 1:k)             # Hankel with 0s in SE half  
     Hkl[(1+p*(i-1)):(p*i), 1:(r*(1+k-i))]<- M[,(1+r*(i-1)):(r*k),drop=FALSE]
# note: if last block of M is 0, which is often the case, then filling
#   with zeros, as above, is correct.
#   k <- k %/% 2
#   Hkl <- matrix(0, p*k, r*k) 
#   for (i in 1:k)              # (smaller) completely filled Hankel
#     {for (j in 1:k)       
#        {Hkl[(1+p*(i-1)):(p*i), (1+r*(j-1)):(r*j)] <-M[ ,(1+r*(i+j-2)):(r*(i+j-1))] }}

   svd.of.Hkl <- La.svd(Hkl)
   shifted.Hkl <- Hkl[,(r+1):dim(Hkl)[2]]
   shifted.Hkl <- cbind(shifted.Hkl,matrix(0,p*k,r))
   rtQ <- diag(sqrt(svd.of.Hkl$d))
   rtQ.inv <- diag(1/sqrt(svd.of.Hkl$d))   
   o.mtr <- svd.of.Hkl$u %*% rtQ 
   o.inv <- t(svd.of.Hkl$u %*% rtQ.inv )
   c.mtr <- rtQ %*% svd.of.Hkl$vt
   c.inv <- t(rtQ.inv %*% svd.of.Hkl$vt)
   svd.crit <- svd.of.Hkl$d
   c.mtr[svd.crit==0,] <-0     # this gets rid of NAs but the resulting
   c.inv[,svd.crit==0] <-0     #  model may not be minimal
   o.inv[svd.crit==0,] <-0     # eg. F may have zero rows and columns
   o.mtr[,svd.crit==0] <-0     #   (specifically if the model included a trend)
   if (is.null(n))
     {svd.crit <- svd.criteria(svd.of.Hkl$d)
      n <- read.int("Enter the number of singular values to use for balanced model: ")
     }
   H  <- o.mtr[1:p,1:n,drop=FALSE]
   FF <- o.inv[1:n,,drop=FALSE] %*% shifted.Hkl%*%c.inv[,1:n,drop=FALSE]
   if (m == 0) G <- NULL             # no exog.
   else        G <- c.mtr[1:n,1:m,drop=FALSE]  # exog.
   K  <- c.mtr[1:n,(m+1):(m+p),drop=FALSE]
   FF <- FF + K%*%H  #converts from system with lagged outputs as 
                    #   inputs see Mittnik p1190.
  #browser()
  #  good checks here are (with n=length(svd.of.Hkl$d)):
  #   0 == max(abs(shifted.Hkl - o.mtr[,1:n] %*% (FF-K%*%H) %*% c.mtr[1:n,]))
  #   0 != max(abs(shifted.Hkl))
  model<-SS(description="nested-balanced model a la Mittnik",
              F.=FF,G=G,H=H,K= K)
  list(crit=svd.crit,model=model)
}


MittnikReducedModels <- function(largeModel)
{# Return a list of models with all smaller state dimesions.
  largeModel <- TSmodel(toSS(largeModel))
  largeModel <- balanceMittnik(largeModel, n=dim(largeModel$F)[1])
  r <- vector("list", dim(largeModel$F)[1])
  for (j in 1:length(r))
    r[[j]] <- SS(F.=largeModel$F[1:j,1:j,drop=FALSE],
                G=if(is.null(largeModel$G)) NULL else largeModel$G[1:j,,drop=FALSE],
                H=largeModel$H[  , 1:j, drop=FALSE],   K= largeModel$K[1:j,,drop=FALSE])
  r
}


############################################################################
#
#       experimental estimation techniques    <<<<<<<<<<
#
############################################################################

estBlackBox2 <- function(data, estimation="estVARXls", 
          lag.weight=.9, 
          reduction="MittnikReduction", 
          criterion="taic", 
          trend=FALSE, 
          subtract.means=FALSE,  re.add.means=TRUE, 
          standardize=FALSE, verbose=TRUE, max.lag=12)
{if ((estimation!="estVARXls") && (trend) )
     {cat("Trend estimation only support with estVARXls.\n")
      cat("Proceeding using estVARXls.\n")
      estimation<-"estVARXls"
     }

 if(estimation=="estVARXls")
     model <- estVARXls(data,trend=trend, subtract.means=subtract.means, 
                         re.add.means=re.add.means, max.lag=max.lag, 
                         standardize=standardize, lag.weight=lag.weight)
 # else if(estimation=="estVARXar")
 #     model <-estVARXar(data, subtract.means=subtract.means, max.lag=max.lag)
 else
    stop("Only estVARXls estimation technique is supported to date.\n")
 if (verbose) cat("First VAR model,              lags= ", 
          dim(model$model$A)[1]-1, ", -log likelihood = ", model$estimates$like[1], "\n")
 model <- l(toSS(model),model$data) # data is standardized if standardize=T in estimation
 n <- dim(model$model$F)[1]
 if (verbose) cat("Equivalent    state space model, n= ", 
                  n, ", -log likelihood = ", model$estimates$like[1], "\n")
 if (1 < n)
   {model <- do.call(reduction,
                     list(model,criterion=criterion, verbose=verbose))
   #model <- eval(call(reduction,model,criterion=criterion, verbose=verbose))
    if (verbose) cat("Final reduced state space model, n= ",
              dim(model$model$F)[1], ", -log likelihood = ", model$estimates$like[1], "\n")
   }
  if (verbose &&  dev.cur() != 1 ) checkResiduals(model)
 model
}


bestTSestModel <- function(models, sample.start=10, sample.end=NULL,
    criterion="aic", verbose=TRUE){
  values <- NULL
  for (lst in models ) 
    {z <- informationTestsCalculations(lst, sample.start=sample.start, 
                  sample.end=sample.end)
     values <-rbind(values,z)
    }
  if (verbose)
     {cat("Criterion value for all models based on data starting in period: ",  
           sample.start, "\n")
      cat(values[,criterion], "\n")
     }
#Rbug rbind above looses dimnames
  dimnames(values) <- dimnames(z)
  opt <-order(values[,criterion])[1]  # minimum
  invisible(models[[opt]])
}


estBlackBox3 <- function(data, estimation="estVARXls", 
       lag.weight=1.0, 
       reduction="MittnikReduction", 
       criterion="aic", 
       trend=FALSE, subtract.means=FALSE,  re.add.means=TRUE, 
       standardize=FALSE, verbose=TRUE, max.lag=12, sample.start=10) {
    if ((estimation!="estVARXls") && (trend) )
     {cat("Trend estimation only support with estVARXls.\n")
      cat("Proceeding using estVARXls.\n")
      estimation<-"estVARXls"
     }
 models <- vector("list", max.lag)
 for (i in 1:max.lag)
   {if(estimation=="estVARXls")
      models[[i]] <- estVARXls(data,trend=trend, 
                          subtract.means=subtract.means,
                          re.add.means=re.add.means, max.lag=i, 
                          standardize=standardize, lag.weight=lag.weight)
    else
      stop("Only estVARXls estimation technique is supported to date.\n")
   }
 model <- bestTSestModel(models, criterion=criterion, sample.start=sample.start, verbose=verbose)
 if (verbose) cat("Selected VAR model,              lags= ", 
          dim(model$model$A)[1]-1, ", -log likelihood = ", model$estimates$like[1], "\n")
 model <- l(toSS(model),model$data) # data is standardized if standardize=T in estimation
 n <- dim(model$model$F)[1]
 if (verbose) cat("Equivalent    state space model,    n= ", 
                  n, ", -log likelihood = ", model$estimates$like[1], "\n")
 if (1 < n)
   {model <- do.call(reduction,
                      list(model,criterion=criterion, verbose=verbose))
   #model <- eval(call(reduction,model,criterion=criterion, verbose=verbose))# , sample.start=sample.start))
    if (verbose) cat("Final reduced state space model, n= ",
              dim(model$model$F)[1], ", -log likelihood = ", model$estimates$like[1], "\n")
   }
  if (verbose &&  dev.cur() != 1 ) checkResiduals(model)
 model
}


bft <- function(data, ...) estBlackBox4(data, ...)

estBlackBox4 <- function(data, estimation="estVARXls", 
                lag.weight=1.0,  variable.weights=1, 
                reduction="MittnikReduction", 
                criterion="taic", 
                trend=FALSE, subtract.means=FALSE,  re.add.means=TRUE, 
                standardize=FALSE, verbose=TRUE, max.lag=12, 
		sample.start=10, warn=TRUE)
{if ((estimation!="estVARXls") && (trend) )
     {cat("Trend estimation only support with estVARXls.\n")
      cat("Proceeding using estVARXls.\n")
      estimation<-"estVARXls"
     }
 models <- vector("list", max.lag)
 for (i in 1:max.lag)
   {if(estimation=="estVARXls")
      {model<- estVARXls(data,trend=trend, 
                          subtract.means=subtract.means,
                          re.add.means=re.add.means, max.lag=i, 
                          standardize=standardize, lag.weight=lag.weight,
                          warn=warn)
      }
    else if(estimation=="estWtVariables")
      {model<- estWtVariables(data, variable.weights,
                        estimation.args=list(trend=trend, 
                          subtract.means=subtract.means,
                          re.add.means=re.add.means, max.lag=i, 
                          standardize=standardize, lag.weight=lag.weight,
                          warn=warn) )
      }
    else
      stop("Estimation technique not yet is supported.\n")
    if (verbose) cat("Estimated  VAR   model       -log likelihood = ", 
        model$estimates$like[1],", lags= ",  dim(model$model$A)[1]-1,"\n")
    model <- l(toSS(model),model$data,warn=warn) # data is standardized if standardize=T in estimation
    n <- dim(model$model$F)[1]
    if (verbose) cat("Equivalent state space model -log likelihood = ", 
               model$estimates$like[1], ",   n = ", n, "\n")
    if (1 < n)
      {model <- do.call(reduction,
                    list(model,criterion=criterion, verbose=verbose, warn=warn))
      #model <- eval(call(reduction,model,criterion=criterion, 
      #              verbose=verbose, warn=warn))# , sample.start=sample.start))
       if (verbose) cat("Final reduced state space model, n= ",
           dim(model$model$F)[1], ", -log likelihood = ", 
                         model$estimates$like[1], "\n")
      }
     models[[i]] <- model
   }
 model <- bestTSestModel(models, criterion=criterion, sample.start=sample.start, verbose=verbose)
 if (verbose &&  dev.cur() != 1 ) checkResiduals(model)
 model
}

#  z<-estBlackBox4(eg.data,  max.lag=3 ) 



############################################################

#     statistical test functions  <<<<<<<<<<

############################################################

Portmanteau <- function(res){
  # Portmanteau statistic for residual
  # before R 1.9.0 required ts not stats
  ac <- stats::acf(as.ts(res),type="covariance", plot=FALSE)$acf
  p <- dim(ac)[1]
#  a0 <- solve(ac[1,,])  the following is more robust than solve for
#              degenerate densities
  v <- La.svd(ac[1,,])
#  if(1 == length(v$d)) a0 <- t(v$vt) %*%     (1/v$d) %*% t(v$u)
#	           else a0 <- t(v$vt) %*% diag(1/v$d) %*% t(v$u)	
#  following is equivalent
  a0 <-  t(v$vt) %*% sweep(t(v$u),1,1/v$d, "*") 
  P <-0
  for (i in 2:p) 
    {P <- P + sum(diag(t(ac[i,,]) %*% a0 %*% ac[i,,] %*% a0))
    }
  dim(res)[1]*P
}


checkResiduals <- function(obj, ...)  UseMethod("checkResiduals")

checkResiduals.TSestModel <- function(obj, ...){
   invisible(checkResiduals(obj$estimates$pred - outputData(obj), ...))}

checkResiduals.TSdata <- function(obj, ...){
	invisible(checkResiduals(outputData(obj), ...)) }

checkResiduals.default <- function(obj, ac=TRUE, pac=TRUE,
     select=seq(nseries(obj)), drop=NULL, plot.=TRUE, 
     graphs.per.page=5, verbose=FALSE, ...)
{#  (... further arguments, currently disregarded)
 #  select can indicate column indices of  residuals to use.
 # If drop is not null it should be a vector of the row indices of residuals
 #  which are not to be used. Typically this can be used to get rid
 #  of bad initial conditions (eg. drop=seq(10) ) or outliers.
  resid <- selectSeries(obj, series=select)
  if (!is.null(drop))  resid <- resid[-drop,, drop=FALSE]
  mn<-colMeans(resid)
  acr <-NULL
  pacr<-NULL
  cov <- var(resid)
  if (verbose) 
    {cat("residual means: \n") ; print(mn); cat("\n")
     cat("residual cov  : \n") ; print(cov) ; cat("\n")
    }
  p <- nseries(resid)
  resid0 <- sweep(resid, 2, mn, FUN="-")
#  resid0 <- resid - t(array(colMeans(resid),rev(dim(resid)))) # mean 0
  cusum <- apply(resid0,2,cumsum)/ t(array(diag(var(resid0)),rev(dim(resid0))))
  # before R 1.9.0 required ts not stats
  if(plot. &&  dev.cur() != 1 ) 
    {graphs.per.page <- min(p, graphs.per.page)
     names <- seriesNames(resid)
     old.par <-par(mfcol = c(3, graphs.per.page),
           mar = c(5.1, 4.1,3.1, 0.1), no.readonly=TRUE ) #c(5,4.1,5,0.1) c(2.1, 4.1,3.1, 0.1)
     on.exit(par(old.par))
     for (i in 1:p) 
          {tfOnePlot(resid[,i], xlab=names[i])
           if (i %% graphs.per.page == ceiling(graphs.per.page/2)) 
                 title(main ="Residuals ")
           if (i %% graphs.per.page == 0) 
                 title(main = paste("       page ", floor(i / graphs.per.page)))
           tfOnePlot(cusum[,i], xlab=names[i])
           if (i %% graphs.per.page == ceiling(graphs.per.page/2)) 
                 title(main = "Cusum")
           if (exists("density"))
                       rd <- density(resid[,i],       bw=var(resid[,i])^0.5)
           else if (exists("ksmooth") & ! is.R()) 
                       rd<-ksmooth(resid[,i],bandwidth=var(resid[,i])^0.5)
           else
        stop("Neither ksmooth nor density are available to calculate the plot.")
           plot(rd,type="l",xlab=names[i],ylab="")
           if (i %% graphs.per.page == ceiling(graphs.per.page/2))
              title(main = "kernel estimate of residual distribution")
          }
     if (ac)
       {#par(mfrow = c(1, 1), mar = c(2.1, 4.1,3.1, 0.1), oma=c(0,0,5,0) )
        #acr  <-acf(as.ts(resid), plot=TRUE)$acf
        #mtext("Autocorrelations", side=3, outer=TRUE, cex=1.5)
        acr  <-acf(as.ts(resid), plot=TRUE, mar=c(3,3,2,0.8), oma=c(1,1.2,4,1))$acf
        mtext("Autocorrelations", side=3, outer=TRUE, cex=1.5, line=3)
       }
     if (pac)
       {pacr <- acf(as.ts(resid), plot=TRUE, type= "partial",
	           mar=c(3,3,2,0.8), oma=c(1,1.2,4,1))$acf
        mtext("Partial Autocorrelations", side=3, outer=TRUE, cex=1.5, line=3)
       }
    }
  else 
    {if (ac)  acr  <- acf(as.ts(resid), plot=FALSE)$acf
     if (pac)	  pacr <- acf(as.ts(resid), plot=FALSE, type= "partial")$acf
    }
  if (ac & verbose)
    {cat("residual auto-correlations:\n")
     cat("      lag:   ")
     for (i in 1: dim(acr)[1])  cat(i,"        ")  
     cat("\n") 
     for (i in 1:p) { cat(i,": "); cat(acr[,i,i]); cat("\n")  }
    }
  if (pac & verbose)
    {cat("partial auto-correlations:\n")
     cat("      lag:   ")
     for (i in 1: dim(pacr)[1])  cat(i,"        ")  
     cat("\n") 
     for (i in 1:p) { cat(i,": "); cat(pacr[,i,i]); cat("\n")  }
    }
#  cat("residual normality tests:\n")
#  cat("hetros. tests:\n")

  skewness <- colMeans(sweep(resid,2,mn)^3) / diag(cov)^(3/2)
  kurtosis <- colMeans(sweep(resid,2,mn)^4) / diag(cov)^2
  invisible(list(residuals=resid, mean=mn, cov=cov, acf=acr, pacf=pacr, 
                 cusum=cusum, skewness=skewness, kurtosis=kurtosis))
}

informationTests <- function(..., sample.start=1,sample.end=NULL,
		 Print=TRUE, warn=TRUE){
  if (Print) criteria.table.heading()
  values <- NULL
  options(width=100)
  for (lst in list(...) ) 
    {z <-informationTestsCalculations(lst,
    		 sample.start=sample.start,sample.end=sample.end, warn=warn)
     values <-rbind(values, z)
     if (Print)
       {print(c(z),digits=4)
        cat("\n")
       }
    }
  if (Print & (1 <dim(values)[1]) )
    {cat("opt     ")
     opt <-apply(values,2,order)[1,]  # minimum
     for (i in 1:length(opt)) cat(opt[i],"      ")
     cat("\n")
    }
  if (Print) criteria.table.legend()
  invisible(values)
  }


informationTestsCalculations <- function(lst,
      sample.start=1,sample.end=NULL, warn=TRUE){
   resid <- residuals(lst)
    # the following line is just to work around a bug with old style time series
   if (ncol(outputData(lst$data))==1) dim(resid) <- dim(outputData(lst$data))
   if (is.null(sample.end)) sample.end <- nrow(resid)
   resid <- resid[sample.start:sample.end,,drop=FALSE]
   ml <-   residualStats(resid, NULL, warn=warn)$like[1] # neg.log(likelihood)
# previously   ml <-   L(resid)[1]     # neg. log( likelihood ).
   n  <-  length(coef(lst$model))      # No. of parameters.
   # nt is theorical dimension of parameter space    n(m+2p)
   if (is.ARMA(lst$model)) nt <- NA 
   if (is.SS(lst$model))  
      if (is.null(lst$model$G)) nt <- nrow(lst$model$F)*2*nrow(lst$model$H)
      else   nt <- nrow(lst$model$F)*(ncol(lst$model$G)+2*nrow(lst$model$H))
   r  <- nrow(lst$estimates$pred)*ncol(lst$estimates$pred) #No. of residuals.
   port <-Portmanteau(resid)
   aic <-  2*ml + 2*n             # AIC
   ops  <- options(warn=-1)
   on.exit(options(ops))
   bic <-  2*ml + n*log(r)        # BIC
   gvc <-  2*ml - 2*r*log(1-n/r)  # GCV
   rice<-  2*ml - r*log(1-2*n/r)  # RICE
   fpe <-  2*ml + r*(log(1+(n/r))-log(1-(n/r)))  # FPE
   taic <-  2*ml + 2*nt             # AIC
   tbic <-  2*ml + nt*log(r)        # BIC
   tgvc <-  2*ml - 2*r*log(1-nt/r)  # GCV
   trice<-  2*ml - r*log(1-2*nt/r)  # RICE
   tfpe <-  2*ml + r*(log(1+(nt/r))-log(1-(nt/r)))  # FPE
   z<- matrix(c(port,ml,aic,bic,gvc,rice,fpe,taic,tbic,tgvc,trice,tfpe),1,12)
   dimnames(z)<-list(NULL,c("port","like","aic","bic","gvc","rice","fpe","taic",
                   "tbic","tgvc","trice","tfpe"))
   z
}

criteria.table.heading <- function(){
        cat("                    based on no.of parameters     based on theoretical parameter dim.","\n")
        cat("       PORT  -ln(L)  AIC   BIC   GVC  RICE   FPE   AIC   BIC   GVC   RICE   FPE\n")
}

criteria.table.nheading <- function(){
        cat("                    based on no.of parameters     based on theoretical parameter dim.","\n")
        cat("dim.   PORT  -ln(L)  AIC   BIC   GVC  RICE   FPE   AIC   BIC   GVC   RICE   FPE\n")
}

criteria.table.legend <- function(){
        cat("  PORT  - Portmanteau test                     ") 
        cat("  -ln(L)- neg. log likelihood                  \n")
        cat("  AIC   - neg. Akaike Information Criterion    ")  
        cat("  BIC   - neg. Bayes  Information Criterion    \n")
        cat("  GVC   - Generalized Cross Validation         ") 
        cat("  RICE  - Rice Criterion                       \n")
        cat("  FPE   - Final Prediction Error               \n") 
        cat(" WARNING - These calculations do not account for trend parameters in ARMA models.\n") 
}


svd.criteria <- function(sv){
   cat("\n SINGULAR VALUES OF THE HANKEL MATRIX:\n")
   print(sv,digits=3)
   svsqr.vct <- sv^2/sum(sv^2)    
   cat("\n GELFAND & YAGLOM INFORMATION CRITERIA:\n")
   GY <-cumsum(log(1-svsqr.vct))/sum(log(1-svsqr.vct))
   print(GY,digits=3)
   cat("\n RATIO OF SINGULAR VALUES TO MAX SINGULAR VALUE:\n")
   print(sv/sv[1],digits=3)
   cbind(sv,GY,sv/sv[1])
}   

#######################################################################

tfwindow.TSdata <- function(x, tf=NULL, start=tfstart(tf), end=tfend(tf), warn=TRUE)
{# window a TSdata object
  if (0 != nseriesInput(x))  inputData(x) <- 
      tfwindow(inputData(x), start=start, end=end, tf=tf, warn=warn)
  if (0 != nseriesOutput(x)) outputData(x) <- 
      tfwindow(outputData(x), start=start, end=end, tf=tf, warn=warn)
  x
}

window.TSdata <- function(x, start=NULL, end=NULL, tf=NULL, warn=TRUE, ...)
   tfwindow.TSdata(x, start=start, end=end, tf=tf, warn=warn)
#  (... further arguments, currently disregarded)


combine <- function(e1,e2)UseMethod("combine")
combine.default <- function(e1,e2){list(e1,e2)}


combine.TSdata <- function(e1,e2)
{# make a new TSdata object with the two objects
 outputData(e1) <- tbind(outputData(e1),outputData(e2)) 
 if ((0 != (nseriesInput(e1))) & (0 != (nseriesInput(e2))))
       inputData(e1) <- tbind(inputData(e1),inputData(e2)) 
 else  {if (0 != (nseriesInput(e2)))  inputData(e1) <-inputData(e2)  }
 e1
}


trimNA.TSdata <- function(x, startNAs=TRUE, endNAs=TRUE){
 p <- nseriesOutput(x)
 m <- nseriesInput(x)
 if (m==0)
   mat <- trimNA(outputData(x))
 else
   mat <- trimNA(tbind(inputData(x),outputData(x)),startNAs=startNAs,endNAs=endNAs)
 tf <- tframe(mat)
 if (m!=0)
   {sn <- seriesNamesInput(x)
    inputData(x)  <- tframed(mat[,1:m, drop=FALSE], tf) 
    seriesNamesInput(x) <- sn
   }
 sn <- seriesNamesOutput(x)
 outputData(x) <- tframed(mat[,(m+1):(m+p), drop=FALSE], tf)
 seriesNamesOutput(x) <- sn
 x
}



percentChange.TSestModel <- function(obj, base=NULL, lag=1,
   cumulate=FALSE, e=FALSE, ...)
  {#  (... further arguments, currently disregarded)
      TSdata(output=percentChange(obj$estimates$pred))
  }

percentChange.TSdata <- function(obj, base=NULL, lag=1,
   cumulate=FALSE, e=FALSE, ...)
  {#  (... further arguments, currently disregarded)
   if (0 != (nseriesInput(obj)))  inputData(obj)  <- percentChange(inputData(obj))
   if (0 != (nseriesOutput(obj))) outputData(obj) <- percentChange(outputData(obj))
   obj
  }



############################################################

#     model and data scaling functions   <<<<<<<<<<

############################################################

scale.TSdata <- function(x, center=FALSE, scale=NULL) 
{# scale should be a list with two matrices or vectors, named input and output,
 # giving the multiplication factor for inputs and outputs.
 # vectors are treated as diagonal matrices.
 # If input or output are NULL then no transformation is applied.
 # The resulting data has inputs and outputs which are different from
 #  the original in proportion to scale. ie. if S and T are output and input
 #  scaling matrices then 
 #          y'(t) = S y(t) where y' is the new output
 #          u'(t) = S u(t) where u' is the new input
 sc <- scale$input # inputData(scale) causes problems if scale is a vector
 if (!is.null(sc))
   {if (! (is.matrix(sc) | is.vector(sc))) stop("input scale must be a vector or matrix")
    d <- inputData(x)
    tf <- tframe(d)
    names <- seriesNames(d)
    if(is.matrix(sc))      d <- d %*% t(sc)
    else if(1==length(sc)) d <- d * sc 
    else                   d <- d %*% diag(sc)                             
    tframe(d) <- tf
    seriesNames(d) <- names
    inputData(x) <- d 
   }
 sc <- scale$output # outputData(scale) causes problems if scale is a vector
 if (!is.null(sc))
   {if (! (is.matrix(sc) | is.vector(sc))) stop("output scale must be a vector or matrix")
    d <- outputData(x)
    tf <- tframe(d)
    names <- seriesNames(d)
    if(is.matrix(sc))      d <- d %*% t(sc)
    else if(1==length(sc)) d <- d * sc 
    else                   d <- d %*% diag(sc)                             
    tframe(d) <- tf
    seriesNames(d) <- names
    outputData(x) <- d 
   }
 x
}


scale.TSestModel <- function(x, center=FALSE, scale=NULL)
   {scale(TSmodel(x), center=center, scale=scale)}

checkScale <- function(x, scale) UseMethod("checkScale")

checkScale.TSestModel <- function(x, scale){checkScale(TSmodel(x), scale)}

checkScale.TSmodel <- function(x, scale) 
{# This function only checks for some error conditions.
 # inputData(scale) and outputData(scale) cause problems if scale is a vector
 if (!is.null(scale$input))
   {if (is.matrix(scale$input))
       {if (any(La.svd(scale$input)$d == 0))  
        stop("scale input transformations must be non singular.")
       }
    else if(any(scale$input == 0)) stop("input scale elements must be non zero.")
   }
 if (!is.null(scale$output))
   {if (is.matrix(scale$output))
      {if (any(La.svd(scale$output)$d == 0))  
       stop("output scale transformations must be non singular.")
      }
    else if(any(scale$output == 0)) stop("output scale elements must be non zero.")
   }
 invisible(TRUE)
}


scale.innov <- function(x, center=FALSE, scale=NULL)
{if (!checkScale(x, scale)) stop("scaling error.")
 if (!is.null(scale$input)) 
   {sc <- scale$input
    if(is.vector(sc)) sc <- diag(sc, nseriesInput(x))
    x$G <- x$G %*% solve(sc)
   }
 sc <- scale$output
 if(is.vector(sc)) sc <- diag(sc, nseriesOutput(x))          
 x$H <- sc %*% x$H
 x$K <- x$K %*% solve(sc)
 setTSmodelParameters(x)
}

scale.nonInnov <- function(x, center=FALSE, scale=NULL)
{if (!checkScale(x, scale)) stop("scaling error.")
 if (!is.null(scale$input)) 
   {sc <- scale$input
    if(is.vector(sc)) sc <- diag(sc, nseriesInput(x))
    x$G <- x$G %*% solve(sc)
   }
 sc <- scale$output
 if(is.vector(sc)) sc <- diag(sc, nseriesOutput(x))          
 x$H <- sc %*% x$H
 x$R <- sc %*% x$R %*% solve(sc)
 setTSmodelParameters(x)
}

scale.ARMA <- function(x, center=FALSE, scale=NULL)
{if (!checkScale(x, scale)) stop("scaling error.")
 sc <- scale$output
 if(is.vector(sc)) sc <- diag(sc, nseriesOutput(x))   
 x$A <- polyprod(sc, polyprod(x$A, solve(sc)))
 x$B <- polyprod(sc, polyprod(x$B, solve(sc)))
 if (!is.null(x$C)) 
   {sci <- scale$input
    if (!is.null(sci)) 
       {if(is.vector(sci)) sci <- diag(sci, nseriesInput(x))          
        x$C <- polyprod(x$C, solve(sci))
       }
    x$C <- polyprod(sc, x$C)
   }
 if (!is.null(x$TREND))  x$TREND <- t(sc %*% t(x$TREND))
 setTSmodelParameters(x)
}


############################################################

#     Methods which are generic for models and TSdata    <<<<<<<<<<

# generic methods for TS structures (ie with input and output) <<<<<<<<<<
#                   see also tframe.s                          <<<<<<<<<<

############################################################################
# Tobs and Tobs.default are defined in tfame.s

TobsInput <- function(x)UseMethod( "TobsInput")
TobsOutput <- function(x)UseMethod("TobsOutput")

startInput <- function(x)UseMethod("startInput")
startOutput <- function(x)UseMethod("startOutput")

endInput <- function(x)UseMethod("endInput")
endOutput <- function(x)UseMethod("endOutput")

frequencyInput <- function(x)UseMethod("frequencyInput")
frequencyOutput <- function(x)UseMethod("frequencyOutput")




inputData <- function(x,  series=seqN(nseriesInput(x)))UseMethod("inputData")
outputData <- function(x, series=seqN(nseriesOutput(x)))UseMethod("outputData")

inputData.default <- function(x, series=seqN(nseriesInput(x)))
   selectSeries(x$input, series=series)
outputData.default <- function(x, series=seqN(nseriesOutput(x)))
   selectSeries(x$output, series=series)

"inputData<-" <- function(x,  value)  UseMethod("inputData<-")
"outputData<-" <- function(x,  value) UseMethod("outputData<-")

"inputData<-.default" <- function(x,  value) {x$input <- value; x}
"outputData<-.default" <- function(x,  value){x$output <- value; x}

# The logic (revised as of May26, 1998) is that seriesNames can be an attribute
# of any object (like a matrix) and has been moved to tframe.



 seriesNamesInput <- function(x)UseMethod( "seriesNamesInput")
seriesNamesOutput <- function(x)UseMethod("seriesNamesOutput")

 "seriesNamesInput<-" <- function(x, value)UseMethod( "seriesNamesInput<-")
"seriesNamesOutput<-" <- function(x, value)UseMethod("seriesNamesOutput<-")


############################################################################

#    methods for TSmodels  <<<<<<<<<<
# and check elsewhere too

############################################################################

seriesNames.TSmodel <- function(x)
  {list(input=seriesNamesInput(x), output=seriesNamesOutput(x))}

seriesNamesOutput.TSmodel <- function(x)
 {# return output names if available in the object,
  # otherwise return "out" pasted with integers.
  # Note, there is no TSdata in this case, so tframe is not involved
  if (!is.null(attr(x, "seriesNamesOutput")))
                          return(attr(x, "seriesNamesOutput"))
  else if (0 != nseriesOutput(x)) 
                   return(paste("out", seq(nseriesOutput(x)), sep=""))  
  else return(NULL)
 }

seriesNamesInput.TSmodel <- function(x)
 {# return input names if available in the object,
  # otherwise return "in" pasted with integers.
  if (!is.null(attr(x, "seriesNamesInput")))
                          return(attr(x, "seriesNamesInput"))
  else if (0 != nseriesInput(x)) 
                   return(paste("in", seq(nseriesInput(x)), sep=""))  
  else return(NULL)
 }

"seriesNames<-.TSmodel" <- function(x, value)
   { seriesNamesInput(x) <-  value$input
    seriesNamesOutput(x) <-  value$output
    x
   }

"seriesNamesInput<-.TSmodel" <- function(x,  value)
   {if (!is.null(value) && length(value) != nseriesInput(x))
       stop("model input dimension and number of names do not match.")
    attr(x, "seriesNamesInput")  <- value
    x
   }

"seriesNamesOutput<-.TSmodel" <- function(x,  value)
   {if (!is.null(value) && length(value) != nseriesOutput(x))
       stop("model output dimension and number of names do not match.")
    attr(x, "seriesNamesOutput")  <- value
    x
   }


############################################################################

#    methods for TSestModels  <<<<<<<<<<
# and check elsewhere too   (esp. for start, end and frequency)

############################################################################
start.TSestModel <- function(x, ...)tfstart(x$data)
startInput.TSestModel <- function(x)startInput(x$data)
startOutput.TSestModel <- function(x)startOutput(x$data)

end.TSestModel <- function(x, ...)tfend(x$data)
endInput.TSestModel <- function(x)endInput(x$data)
endOutput.TSestModel <- function(x)endOutput(x$data)

frequency.TSestModel <- function(x, ...)tffrequency(x$data)
frequencyInput.TSestModel <- function(x)frequencyInput(x$data)
frequencyOutput.TSestModel <- function(x)frequencyOutput(x$data)

Tobs.TSestModel <- function(x)Tobs(x$data)
TobsInput.TSestModel <- function(x)TobsInput(x$data)
TobsOutput.TSestModel <- function(x)TobsOutput(x$data)

inputData.TSestModel <- function(x, series=seqN(nseriesInput(x)))
    inputData(x$data, series=series)
outputData.TSestModel <- function(x, series=seqN(nseriesOutput(x)))
    outputData(x$data, series=series)

nseriesInput.TSestModel <- function(x)  nseriesInput(x$data)
nseriesOutput.TSestModel <- function(x) nseriesOutput(x$data)

seriesNamesInput.TSestModel <- function(x)
 {m <- seriesNamesInput(x$model)
  d <- seriesNamesInput(x$data)
  if (!all(m==d)) 
    warning("data and model names do not correspond. Model names returned.")
  m
 }

seriesNamesOutput.TSestModel <- function(x)
 {m <- seriesNamesOutput(x$model)
  d <- seriesNamesOutput(x$data)
  if (!all(m==d)) 
    warning("data and model names do not correspond. Model names returned.")
  m
 }

seriesNames.TSestModel <- function(x)
  {list(input=seriesNamesInput(x), output=seriesNamesOutput(x))}

"seriesNames<-.TSestModel" <- function(x, value)
   { seriesNamesInput(x) <- value$input
    seriesNamesOutput(x) <- value$output
    x
   }


"seriesNamesInput<-.TSestModel" <- function(x, value)
   {seriesNamesInput(x$model) <- value;
    seriesNamesInput(x$data ) <- value;
    x
   }
"seriesNamesOutput<-.TSestModel" <- function(x, value) 
   {seriesNamesOutput(x$model) <- value;
    seriesNamesOutput(x$data ) <- value;
    x
   }


############################################################################

#    methods for TSdata  <<<<<<<<<<
#   (there are methods in other packages too)

############################################################################

is.TSdata <- function(obj) { inherits(obj, "TSdata")}

print.TSdata <- function(x, ...)
 {print.tframed <- function(x, ...)
     {nm <- seriesNames(x)
      tm <- time(x) 
      if (1 != tffrequency(x))
         tm <- paste( floor(tm), "[", 1+ (tm %% 1) * tffrequency(x), "]", sep="")
      dimnames(x) <- list(tm, nm) 
      attr(x,"tframe") <- NULL
      attr(x,"seriesNames") <- NULL
      class(x) <- if (all(class(x) == "tframed")) NULL else class(x)[class(x) != "tframed"]
      print(x,...)
      invisible(x)
     }
  if(0 != (nseriesInput(x)))
     {cat("input data:\n")
      print.tframed(inputData(x), ...)
      if(!is.null(x$input.transformations))
          {cat("input.transformations:\n")
           print(x$input.transformations, ...)
          }
      cat("\n")
     }
  if(0 != (nseriesOutput(x)))
     {cat("output data:\n")
      print.tframed(outputData(x),...)
      if(!is.null(x$output.transformations))
         {cat("output.transformations:\n")
          print(x$output.transformations, ...)
         }
     }
   cat("\n")
   if(!is.null(x$retrieval.date))
      cat("retrieval date: ", x$retrieval.date, "   ")
   if(!is.null(x$source)) 
     {cat("source:\n")
      print(x$source)
     }
  invisible(x)
}

summary.TSdata <- function(object, ...)
  {#  (... further arguments, currently disregarded)
   d <- outputData(object)
   if (!is.tframed(d)) d <- as.ts(d)
   st <- tfstart(d)
   en <- tfend(d)
   fr <- tffrequency(d)
   d <- cbind(d,inputData(object))

   classed(list(  # summary.TSdata constructor
      description=object$description,
      start=st,
      end=en,
      freq=fr,
      sampleT= nrow(outputData(object)),
      p=nseriesOutput(object),
      m=nseriesInput(object),
      ave=colMeans(d),
      max=apply(d,2,max),
      min=apply(d,2,min),
      retrieval.date=object$retrieval.date, 
      source=object$source), "summary.TSdata")
  }

print.summary.TSdata <- function(x, digits=options()$digits, ...)
{#  (... further arguments, currently disregarded)
   if (!is.null(x$description)) cat(x$description,"\n")
   cat("start: ", x$start, " end: ",  x$end," Frequency: ", x$freq,"\n")
   cat("sample length = ",x$sampleT,"\n")
   cat("number of outputs=",x$p, "   number of inputs=",x$m, "\n")
   cat("average :\n")
   print(x$ave)
   cat("max    :\n")
   print(x$max)
   cat("min    :\n")
   print(x$min)
   cat("\n")
   if(!is.null(x$retrieval.date)) 
      cat("retrieval date: ", x$retrieval.date, "   ")
   if(!is.null(x$source)) {cat("source:\n"); print(x$source)  }
   invisible(x)
}


tfplot.TSdata <- function(x, ..., 
        tf=NULL, start=tfstart(tf), end=tfend(tf),
        select.inputs  = seq(length=nseriesInput(x)),
        select.outputs = seq(length=nseriesOutput(x)),
        Title=NULL, xlab=NULL, ylab=NULL, 
	graphs.per.page=5, mar=par()$mar, reset.screen=TRUE)
{# plot input and output data.
 # ... is a list of objects of class TSdata (with similar input
 #       and output dimensions.
 #previously mar = if(is.R()) c(3.1,6.1,1.1,2.1) else c(5.1,6.1,4.1,2.1)
 # Note that using ... like this means it cannot be used to pass additional
 #   arguments to plot, so unfortunately all necessary plot arguments must be 
 #   explicit in the arguments to tfplot.TSdata
 # start is the starting point (date)  and end the ending point for
 # plotting. If not specified the whole sample is plotted.
# output graphics can be paused between pages by setting par(ask=T).
  #x <- freeze(x)
  Ngraphs <- length(select.outputs) + length(select.inputs)
  if(reset.screen)
    {Ngraphs <- min(Ngraphs, graphs.per.page)
     old.par <- par(mfcol = c(Ngraphs, 1), mar=mar, no.readonly=TRUE) # previously c(5.1,6.1,4.1,2.1)) 
     on.exit(par(old.par))
    }
  #if (is.null(Title)) Title <- ""

  if (!is.null(ylab))
       names <- rep(ylab, nseriesInput(x) + nseriesOutput(x))
  else names <-  c(seriesNamesInput(x),  seriesNamesOutput(x))

  if (1 == length(xlab)) xlab <- rep(xlab, length(names))

#  if (0 != length(select.inputs)) 
    {for (i in select.inputs) 
      {j <- 0
       z <- matrix(NA, TobsInput(x), length(append(list(x),list(...))))
       if(mode(i)=="character") i <- match(i, seriesNamesInput(x))
       for (d in append(list(x),list(...)) ) 
         {if (!is.TSdata(d))
            stop("Expecting TSdata objects. Do not truncate argument names as that can cause a problem here.")
          #d <- freeze(d)
          j <- j + 1
          z[, j] <- inputData(d, series = i)
         }
       tframe(z) <-tframe(inputData(x))
       tfOnePlot(z, xlab=xlab[i], ylab=names[i], start=start, end=end) 
       if(!is.null(Title) && (i==1) && (is.null(options()$PlotTitles)
                || options()$PlotTitles)) title(main = Title)
    } }
  for (i in select.outputs) 
    {j <-0
     if(mode(i)=="character") i <- match(i, seriesNamesOutput(x))
     z <- matrix(NA,TobsOutput(x),length(append(list(x),list(...))))
     #z <-NULL 
     for (d in append(list(x),list(...)) ) 
       {if (!is.TSdata(d))
            stop("Expecting TSdata objects. Do not truncate argument names as that can cause a problem here.")
        #d <- freeze(d)
        j <- j+1
        z[,j]<-outputData(d, series=i) 
       }
     tframe(z) <-tframe(outputData(x))
     tfOnePlot(z, xlab=xlab[i], ylab=names[nseriesInput(x) + i],start=start, end=end)
     if(!is.null(Title) && (0 == nseriesInput(x)) && (i==1) && 
      (is.null(options()$PlotTitles) || options()$PlotTitles)) title(main = Title)
    }
  invisible()
}
 
inputData.TSdata <- function(x, series=seqN(nseriesInput(x)))
  {if (is.null(x$input)) NULL  else selectSeries(x$input, series=series) }

outputData.TSdata <- function(x, series=seqN(nseriesOutput(x)))
  {if (is.null(x$output)) NULL  else selectSeries(x$output, series=series) }


"inputData<-.TSdata" <- function(x, value) 
   {cls <- class(x); x$input <- value; class(x) <- cls
    x
   }


"outputData<-.TSdata" <- function(x, value)
   {cls <- class(x); x$output <-value; class(x) <- cls
    x
   }

# Note: series names changed to an attribute of input and output for 
#     data but not for models!!!!

seriesNames.TSdata <- function(x)
 {list(input=seriesNamesInput(x), output=seriesNamesOutput(x))}

 seriesNamesInput.TSdata <- function(x) {seriesNames( inputData(x))}
seriesNamesOutput.TSdata <- function(x) {seriesNames(outputData(x))}

"seriesNames<-.TSdata" <- function(x, value){ 
    seriesNamesInput(x)  <-  value$input
    seriesNamesOutput(x) <-  value$output
    x
   }

"seriesNamesInput<-.TSdata"  <- function(x,  value)
   {seriesNames(inputData(x))  <- value; x }

"seriesNamesOutput<-.TSdata" <- function(x,  value) 
   {seriesNames(outputData(x)) <- value; x }

nseriesInput.TSdata <- function(x)
   {if (is.null(x$input)) 0 else nseries(x$input)}

nseriesOutput.TSdata <- function(x)
   {if (is.null(x$output)) 0 else nseries(x$output)}

start.TSdata <- function(x, ...)
{#  (... further arguments, currently disregarded)
 i  <-  startInput(x)
 o  <- startOutput(x)
 if (((!is.null(o)) & (!is.null(i))) && all(i==o)) return(o)
 else return(c(i,o))
}

startInput.TSdata <- function(x)
 {if (is.null(x$input))  return(NULL)
  else return(tfstart(x$input))
 }

startOutput.TSdata <- function(x)
 {if (is.null(x$output))  return(NULL)
  else return(tfstart(x$output))
 }

end.TSdata <- function(x, ...)
{#  (... further arguments, currently disregarded)
 i  <-  endInput(x)
 o  <- endOutput(x)
 if (((!is.null(o)) & (!is.null(i))) && all(i==o)) return(o)
 return(c(i,o))
}

endInput.TSdata <- function(x)
 {if (is.null(x$input))  return(NULL)
  else return(tfend(x$input))
 }

endOutput.TSdata <- function(x)
 {if (is.null(x$output))  return(NULL)
  else return(tfend(x$output))
 }

frequency.TSdata <- function(x, ...)
{#  (... further arguments, currently disregarded)
 i  <-  frequencyInput(x)
 o  <- frequencyOutput(x)
 if (((!is.null(o)) & (!is.null(i))) && all(i==o)) return(o)
 return(c(i,o))
}

frequencyInput.TSdata <- function(x)
 {if (is.null(x$input))  return(NULL)
  else return(tffrequency(x$input))
 }

frequencyOutput.TSdata <- function(x)
 {if (is.null(x$output))  return(NULL)
  else return(tffrequency(x$output))
 }


Tobs.TSdata <- function(x, ...) Tobs(outputData(x))
#  (... further arguments, currently disregarded)
TobsOutput.TSdata <- function(x)  Tobs(outputData(x))[1]
TobsInput.TSdata  <- function(x)  Tobs( inputData(x))[1]

tbind.TSdata <- function(x, d2, ..., pad.start=TRUE, pad.end=TRUE, warn=TRUE)
 {if( ! (is.TSdata(x) & is.TSdata(d2)))
     stop("tbind requires arguments to be of a similar type (ie. TSdata).")
  if ((0 != nseriesInput(x)) || (0 != nseriesInput(d2)) )
    inputData(x) <- tbind(inputData(x),inputData(d2),
                 ..., pad.start=pad.start, pad.end=pad.end, warn=warn)
  if ((0 != nseriesOutput(x)) || (0 != nseriesOutput(d2)) )
    outputData(x) <- tbind(outputData(x),outputData(d2),
                 ..., pad.start=pad.start, pad.end=pad.end, warn=warn)
  x
 }



testEqual.TSdata <- function(obj1, obj2, fuzz=1e-16)
  {r <- TRUE
   if (r & (!is.null(obj1$input)))
     {if(is.null(obj2$input)) r <- FALSE
      else  r <-testEqual(obj1$input, obj2$input, fuzz=fuzz)
     }
   if (r & (!is.null(obj1$output)))
     {if(is.null(obj2$output)) r <- FALSE
      else r <-testEqual(obj1$output, obj2$output, fuzz=fuzz)
     }
   r
  }




TSdata <- function (data=NULL, ...) UseMethod("TSdata") 
 

#TSdata.default <- function(data=NULL, input=NULL, output=NULL, ...)  
#{if (is.null(data) && (!is.null(input) | !is.null(output) ))
#    {if(!is.null(input) && is.vector(input)) input <- 
#	   tframed(matrix(input, length(input),1), tf=tframe(input),
#	           names=seriesNames(input))
#     if(!is.null(output) && is.vector(output)) output <- 
#	   tframed(matrix(output, length(output),1), tf=tframe(output),
#	           names=seriesNames(output))
#     data <- classed(list(input=input, output=output), "TSdata") # constructor
#  }else 
#     data <- classed(data, "TSdata")   # constructor keeps other list elements
#  
# if(!is.list(data)) stop("TSdata.default could not construct a TSdata format.")	    
# if   ( 0 == nseriesInput(data)   &    0 == nseriesOutput(data) 
#    |  (0 != nseriesOutput(data) && !is.matrix(outputData(data)))
#    |  (0 !=  nseriesInput(data) && !is.matrix( inputData(data))) )
#      stop("TSdata.default could not construct a TSdata format.")
# data
#}

TSdata.default <- function(data=NULL, input=NULL, output=NULL, ...)  
{if (is.null(data) & is.null(input)  & is.null(output) )
      stop("TSdata.default needs data, input, or output specified.")
 if (is.null(data))
    {if(!is.null(input)) input <- 
	   tframed(as.matrix(input), tf=tframe(input), names=seriesNames(input))
     if(!is.null(output)) output <- 
	   tframed(as.matrix(output), tf=tframe(output),names=seriesNames(output))
     data <- classed(list(input=input, output=output), "TSdata") # constructor
  }else 
     data <- classed(data, "TSdata")   # constructor keeps other list elements
  
 if(!is.list(data)) stop("TSdata.default could not construct a TSdata format.")	    
 if   ( 0 == nseriesInput(data)   &    0 == nseriesOutput(data) 
    |  (0 != nseriesOutput(data) && !is.matrix(outputData(data)))
    |  (0 !=  nseriesInput(data) && !is.matrix( inputData(data))) )
      stop("TSdata.default could not construct a TSdata format.")
 data
}

TSdata.TSdata <- function(data, ...)  {data} # already TSdata 

TSdata.TSestModel <- function(data, ...)  {data$data} # extract TSdata 

as.TSdata <- function(d) 
 {# Use whatever is actually $input and $output,
  #strip any other class, other parts of the list, and DO NOT use
  # inputData(d) and outputData(d) which may eg reconstitute something
  TSdata(input=d$input, output=d$output)
 }

"tframe<-.TSdata" <- function(x, value) {
 if (0 != nseriesInput(x))  tframe(inputData(x))  <- value
 if (0 != nseriesOutput(x)) tframe(outputData(x)) <- value
 x
}

tframed.TSdata <- function(x, tf=NULL, names=NULL, ...)  
{# switch to tframe representation
 if(0 != (nseriesOutput(x)))
       outputData(x) <-tframed(outputData(x), tf=tf, names=names$output, ...)
 if (0 != (nseriesInput(x)))
        inputData(x) <-tframed(inputData(x),  tf=tf, names=names$input, ...)
 x
}  


############################################################

#  Utility function for time series noise   <<<<<<<<<<

############################################################

makeTSnoise <- function(sampleT,p,lags,noise=NULL, rng=NULL,
                        Cov=NULL, sd=1, noise.model=NULL,
                        noise.baseline=0,
                        tf=NULL, start=NULL,frequency=NULL)
 {# CAUTION: changes here can affect historical comparisons.
  # noise.baseline is added to noise. It should be either a scalar, a matrix of
  #   the same dimension as noise (or noise generated by noise.model), or a
  #   vector of length equal to the dimension of the noise process (which will
  #   be replicated for all Tobs.)
 if(is.null(rng)) rng <- setRNG() # returns setting so don't skip if NULL
 else        {old.rng <- setRNG(rng);  on.exit(setRNG(old.rng))  }

  if ( (!is.null(noise)) & (!is.null(noise.model)) )
    stop("noise and noise.model cannot both be specified.")

  if(!is.null(noise.baseline) && is.matrix(noise.baseline) &&
    (sampleT < dim(noise.baseline)[1]))
      {warning("sampleT (and start date) for noise adjusted to match noise.baseline")
       sampleT <- dim(noise.baseline)[1]
      }

 # Note: noise is added to initial conditions.
 if (!is.null(noise.model))
   {noise <- outputData(simulate(noise.model, sampleT=sampleT+lags))
    noise <- list(w0=noise[1:lags,,drop=FALSE], w=noise[lags+seq(sampleT),,drop=FALSE])
   }
  if (is.null(noise)) {
    w0 <-matrix(NA,lags,p)
    w <- matrix(NA,sampleT,p)
    if (!is.null(Cov)) {
        if(length(Cov) == 1) Cov <- diag(Cov, p)
	W <- t(chol(Cov))
        w <- t(W %*% t(matrix(rnorm((lags+sampleT)*p),(lags+sampleT),p)))
#  above has unfortunate effect that w0 was not from the first p values below
#   would be better, but changes a lot? of tests
#        w <- t(W %*% matrix(rnorm((lags+sampleT)*p),p, (lags+sampleT)))
	w0 <- w[1:lags,]
	w  <- w[-c(1:lags),]
        }
    else {
       if (length(sd)==1) sd <-rep(sd,p)
       for (i in 1:p)
         {w0[,i] <- rnorm(lags,sd=sd[i])
          w[,i]  <- rnorm(sampleT,sd=sd[i])
         }
       }
    noise <- list(w=w, w0=w0)
    }
  else {
     if (is.null(noise$w0) || is.null(noise$w) )
       stop("supplied noise structure is not correct.")
     }
       
  if(!is.null(noise.baseline))
     {if (is.vector(noise.baseline))
        {if(length(noise.baseline)==1) noise$w <- noise$w + noise.baseline
         else if(length(noise.baseline)==1)
           noise$w <-noise$w + t(array(noise.baseline, rev(dim(noise$w))))
         else stop("noise.baseline vector is not correct length.")
        }
      else  noise$w <- noise$w + noise.baseline
     }

  if(!is.null(tf)) tframe(noise$w) <- tf
  else if(!is.null(start))
      {if (is.null(frequency))
         {frequency <- 1
          warning("start set but frequency not specified. Using frequency=1.")
         }
       else noise$w <-tframed(noise$w, list(start=start, frequency=frequency))
       if (is.tframed(noise.baseline) &&
           testEqual(tframe(noise.baseline),tframe(noise$w)))
           {warning("tframe of noise set to tframe of noise.baseline.")
            tframe(noise$w )<-tframe(noise.baseline)
            if(!all(dimnames(noise$w)[[2]] == dimnames(noise.baseline)[[2]]))
              warning("noise names and noise.baseline names do not correspond.")
           }
      }
  append(noise, list(sampleT=sampleT, rng=rng,
     Cov=Cov, sd=sd, noise.model=noise.model,version=as.vector(version)))
 }
