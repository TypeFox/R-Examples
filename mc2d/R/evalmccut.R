#<<BEGIN>>
evalmccut <- function(model, nsv = ndvar(), nsu = ndunc(), seed = NULL,  ind = "index")
#TITLE Evaluates a Two-Dimensional Monte Carlo Model in a Loop.
#NAME mccut
#DESCRIPTION
#\samp{evalmccut} evaluates a Two-Dimensional Monte Carlo model using a loop on the uncertainty dimension.
#Within each loop, it calculates statistics in the variability dimension and stores them 
#for further analysis.
#It allows to evaluate very high dimension models using (unlimited?) time instead of (limited) memory.</>
#\samp{mcmodelcut} builds a \samp{mcmodelcut} object that can be sent to \samp{evalmccut}.
#KEYWORDS methods
#INPUTS
#{model}<<a \samp{mcmodelcut} object obtained using \samp{mcmodelcut} function or (directly) a valid call including three blocks.
#See Details and Examples for the structure of the call.>>
#{x}<<a call or an expression (if \samp{is.expr=TRUE}) including three blocks. See Details and Examples for the structure of the call.>>
#[INPUTS]
#{nsv}<<The number of simulations for variability used in the evaluation.>>
#{nsu}<<The number of simulations for uncertainty used in the evaluation.>>
#{seed}<<The random seed used for the evaluation. If \samp{NULL} the \samp{seed} is unchanged.>>
#{ind}<<The variable name used in \samp{model} to refers to the uncertainty. see Details and Example.>>
#{is.expr}<< \samp{FALSE} to send a call,  \samp{TRUE} to send an expression (see \code{\link{mcmodel}} examples)>>
#{lim}<<A vector of values used for the quantile function (uncertainty dimension).>>
#{digits}<<Number of digits in the print.>>
#{\dots}<<Additional arguments to be passed in the final print function.>>
#VALUE
#An object of class \samp{mccut}. This is a list including statistics evaluated within the third block. Each list consists of all the \samp{nsu}
#values obtained. The \samp{print.mccut} method print the median, the mean, the \samp{lim} quantiles estimated on each statistics on the
#uncertainty dimension.
#NOTE
#The methods and functions available on the \samp{mccut} object is function of the statistics evaluated within the third block:
#{*}<<a \code{\link{print.mccut}} is available as soon as one statistic is evaluated within the third block;>>
#{*}<<a \code{\link{summary.mccut}} and a \code{\link{tornadounc.mccut}} are available if a \code{\link{summary.mc}}
#is evaluated within the third block;>>
#{*}<<\code{\link{converg}} may be used if a \code{\link{summary.mc}} is evaluated within the third block;>>
#{*}<<a \code{\link{plot.mccut}} is available if a \code{\link{plot.mc}} is evaluated within the third block.
#(Do not forget to use the argument \samp{draw = FALSE} in the third block);>>
#{*}<<a \code{\link{tornado}} is available if a \code{\link{tornado}} is evaluated within the third block.>>
#DETAILS
#This function should be used for high dimension Two-Dimensional Monte-Carlo simulations, when the memory limits of \R are attained.
#The use of a loop will take (lots of) time, but less memory.</>
#\samp{x} (or \samp{model} if a call is used directly in \samp{evalmccut}) should be built as three blocks, separated by \samp{\{}.
#{#}<<The first block is evaluated once (and only once) before the first loop (step 1).>>
#{#}<<The second block, which should lead to an \samp{mc} object, is evaluated using \samp{nsu = 1} (step 2).>>
#{#}<<The third block is evaluated on the \samp{mc} object. All resulting statistics are stored (step 3).>>
#{#}<<The steps 2 and 3 are repeated \samp{nsu} times. At each iteration, the values of the loop index (from 1 to \samp{nsu})
#is given to the variable specified in \samp{ind}.>>
#{#}<<Finally, the \samp{nsu} statistics are returned in an invisible object of class \samp{mccut}.>>
#Understanding this, the call should be built like this:
#\samp{{{block 1}{block 2}{block 3}}}
#{#}<<The first block (maybe empty) is an expression that will be evaluated only once.
#This block should evaluate all \samp{"V" mcnode} and \samp{"0" mcnode}s. It may evaluate and \samp{"U" mcnode}
#that will be sent in the second and third block by column, and, optionnaly, some other codes
#(even \samp{"VU" mcnode}, sent by columns) that can not be evaluated if \samp{ndunc=1}
#(e.g. sampling without replacement in the uncertainty dimension).>>
#{#}<<The second block is an expression that leads to the \samp{mc} object.
#It must end with an expression as \samp{mymc <- mc(...)}. The variable specified as \samp{ind}
#may be helpful to refer to the uncertainty dimension in this step >>
#{#}<<The last block should build a list of statistics refering to the \samp{mc}
#object. The function \samp{summary} should be used if a summary, a tornado on uncertainty (\code{\link{tornadounc.mccut}}) or a convergence diagnostic
#\code{\link{converg}} is needed,
#the function \code{\link{plot.mc}} should be used if a plot is needed, the function \code{\link{tornado}} should be used if a tornado is needed.
#Moreover, any other function that leads to a
#vector, a matrix, or a list of vector/matrix of statistics evaluated from the \samp{mc} object may be used. list are time consuming.>>
#IMPORTANT WARNING: do not forget to affect the results, since the print method provide only
#a summary of the results while all data may be stored in an \samp{mccut} object.
#NOTE
#The seed is set at the beginning of the evaluation. Thus, the complete similarity of two evaluations
#is not certain, depending of the structure of your model. Moreover, with a similar seed, the simulation will not be equal to
#the one obtained with \code{\link{evalmcmod}} since the random samples will not be obtained in the same order.</>
#In order to avoid conflicts between the \samp{model} evaluation and the function, the function uses upper case variables.
#Do not use upper case variables in your model.</>
#The function should be re-adapted if a new function to be applied on \samp{mc} objects is written.
#SEE ALSO
#\code{\link{evalmcmod}}
#EXAMPLE
#modEC3 <- mcmodelcut({
#
### First block:
###  Evaluates all the 0, V and U nodes.
#    { cook <- mcstoc(rempiricalD, type = "V", values = c(0, 1/5,
#        1/50), prob = c(0.027, 0.373, 0.6))
#      serving <- mcstoc(rgamma, type = "V", shape = 3.93, rate = 0.0806)
#      conc <- mcstoc(rnorm, type = "U", mean = 10, sd = 2)
#      r <- mcstoc(runif, type = "U", min = 5e-04, max = 0.0015)
#    }
### Second block:
###  Evaluates all the VU nodes
###  Leads to the mc object. 
#    {
#      expo <- conc * cook * serving
#      dose <- mcstoc(rpois, type = "VU", lambda = expo)
#      risk <- 1 - (1 - r)^dose
#      res <- mc(conc, cook, serving, expo, dose, r, risk)
#    }
### Third block:
###  Leads to a list of statistics: summary, plot, tornado
###  or any function leading to a vector (et), a list (minmax),
###  a matrix or a data.frame (summary)
#    {
#      list(
#        sum = summary(res),
#        plot = plot(res, draw=FALSE),
#        minmax = lapply(res,range)
#     )
#    }
#})
#
#x <- evalmccut(modEC3, nsv = 101, nsu = 101, seed = 666)
#summary(x)

#CREATED 08-04-16
#REVISED 08-04-16
#--------------------------------------------
#
{

  # DEFINE SOME FUNCTIONS
  # DEFINEOUT build the empty list that the results will fill 

  DEFINEOUT <- function(ARG){
    LAFUNC <- function(ARGJBD){
      if(inherits(ARGJBD,"tornado")) ARGJBD <- ARGJBD$value
      if(is.list(ARGJBD)) return(DEFINEOUT(ARGJBD))
      if(is.data.frame(ARGJBD)) ARGJBD <- as.matrix(ARGJBD)
      if(is.array(ARGJBD)){
         DIMS <- dim(ARGJBD)
         if(inherits(ARGJBD,"mcnode")){
            DIMS <- c(DIMS[1],nsu,DIMS[3])
            NOM <- vector(mode="list",length=3)}
         else {if(length(DIMS) > 2)                                              #Other arrays = matrix
                stop("Array > 2 dimensions are not supported in this function")
          DIMS <- c(DIMS[1],nsu,DIMS[2])
          NOM <- list(rownames(ARGJBD),NULL,colnames(ARGJBD))}
      }
      else {
        DIMS <- c(length(ARGJBD),nsu,1)
        NOM <- list(names(ARGJBD),NULL,NULL)
      }
      return(array(NA,dim=DIMS,dimnames=NOM))
      }
    OUT <- lapply(ARG,LAFUNC)
    names(OUT) <- names(ARG)
    class(OUT) <- class(ARG)
    return(OUT)
  }

  # FONCSPEC is a tricky "autocall" function to fill the results
  FONCSPEC <- function(STAT, PREC, LOOP){
    if(inherits(STAT,"tornado")) {STAT <- STAT$value}
    if(is.list(PREC)){
      PREC <- mapply(FONCSPEC,STAT,PREC,MoreArgs=list(LOOP=LOOP),SIMPLIFY=FALSE)
    }
    else PREC[,LOOP,] <- STAT
    return(PREC)
    }

OLDV <- ndvar()
OLDU <- ndunc()
RESULTAT <- try({                                                                # Debut du try
  ndvar(nsv)
  ndunc(nsu)
  if(!is.null(seed)) set.seed(seed)

  #Analyse expression
  NBEXPR <- length(model[[1]])
  if(NBEXPR != 4) stop("The expression should include three blocks")

  #Evaluate the first block and stock the U and VU mcnode
  LSBEG <- ls()
  eval(model[[1]][[2]])
  LSEND <- ls()
  NEW <- LSEND[!(LSEND %in% LSBEG)]                                             # New objects built
  NEW <- NEW[sapply(NEW,function(x) inherits(get(x),"mcnode"))]                 # New mcnode built
  TYPE <- sapply(NEW, function(x) attr(get(x),which="type"))
  OUTM <- sapply(NEW, function(x) attr(get(x),which="outm"))
    cat("'0' mcnode(s) built in the first block:",NEW[TYPE == "0"],"\n")
    cat("'V' mcnode(s) built in the first block:",NEW[TYPE == "V"],"\n")
    cat("'U' mcnode(s) built in the first block:",NEW[TYPE == "U"],"\n")
    cat("'VU' mcnode(s) built in the first block:",NEW[TYPE == "VU"],"\n")
    cat("The 'U' and 'VU' nodes will be sent column by column in the loop\n")

  QUEL <- (TYPE == "U" | TYPE == "VU")
  OUTM <- OUTM[QUEL]
  NN    <- sum(QUEL)
  TORIG <- TYPE[QUEL]                                                           # Type of Origin
  TORIG1D <- ifelse(TORIG %in% c("0","U"),"0","V")                              # Corresponding Type in 1D
  NORIG <- NEW[QUEL]                                                            # Name of the node to save
  NSAVED <- paste("SAVE",NORIG,sep="")                                          # New names

  for(i in 1:NN)  assign(NSAVED[i], get(NORIG[i]))                              # Save the nodes

  # Prepare the loop
  ndunc(1)                                                                      # nsu = 1
  for(i in 1:5) cat("---------|")
  cat("\n")
  quelcroix <- as.integer(seq(1,nsu,length=50))                                 # define the step for the progress bar
  cat(rep("*",sum(1 == quelcroix)),sep="")                                    # progress bar
  flush.console()

# First simulation
  assign(ind,1)

  # Get the first row of the stocked nodes
  NODE <- vector(mode="list",length=NN)
  for(i in 1:NN) {
    NODE[[i]] <- get(NSAVED[i])[,1,,drop=FALSE]                                 # Can not be put in an mapply function
    attr(NODE[[i]],which="type") <- TORIG[1]                                    # First simulation: keep the original "type"
    attr(NODE[[i]],which="outm") <- OUTM[i]
    class(NODE[[i]])  <- "mcnode"
    assign(NORIG[i],NODE[[i]]) }                                                 # And assign them to their original name


  # Evaluate the model
    eval(model[[1]][[3]])

    NOM <- as.character(model[[1]][[3]][[length(model[[1]][[3]])]][[2]])        # name of the mc second block
    MCOBJ <- get(NOM)
    if(!is.mc(MCOBJ)) stop("The second block result is not an mc object")

    TYPEORI <- sapply(MCOBJ,attr,which="type")                                  #Original type
    TYPENEW <- ifelse(TYPEORI %in% c("0","U"),"0","V")                          #New type

 #Evaluate the model with original type and keep the first results
    PREMS <- eval(model[[1]][[4]])

 #Evaluate the model with new types
    NODE <- mapply("attr<-",NODE,"type",TORIG1D,SIMPLIFY=FALSE)                 # Assign the new type
    for(i in 1:NN) assign(NORIG[i],NODE[[i]])                                   # And assign them to their original name

    MCOBJ[] <- mapply("attr<-",MCOBJ,"type",TYPENEW,SIMPLIFY=FALSE)
    assign(NOM,MCOBJ)

    PREMSPRIM <- eval(model[[1]][[4]])

 #Define the final structure and stock the first simu
    SORTIE <- DEFINEOUT(PREMSPRIM)
    SORTIE <- mapply(FONCSPEC,PREMSPRIM,SORTIE,MoreArgs=list(LOOP=1),SIMPLIFY=FALSE)

# Other Simulations
  for(LOOP in 2:nsu){

    assign(ind,LOOP)

    for(i in 1:NN) {
      NODE[[i]][] <- get(NSAVED[i])[,LOOP,,drop=FALSE]
      assign(NORIG[i],NODE[[i]])
    }

    eval(model[[1]][[3]])
    MCOBJ <- get(NOM)
    MCOBJ[] <- mapply("attr<-",MCOBJ,which="type",value=TYPENEW,SIMPLIFY=FALSE)
    assign(NOM,MCOBJ)

    STAT <- eval(model[[1]][[4]])
    SORTIE <- mapply(FONCSPEC,STAT,SORTIE,MoreArgs=list(LOOP=LOOP),SIMPLIFY=FALSE)

    cat(rep("*",sum(LOOP == quelcroix)),sep="")                                    # progress bar
    flush.console()
    }

  # Post Production : Affecte les classes + petites modifications pour qq fonctions connues

  for(JBD in 1:length(SORTIE)){
    if(inherits(PREMS[[JBD]],"summary.mc")) {
      LESTYPES <- sapply(PREMS[[JBD]],"attr","type")
      SORTIE[[JBD]] <- mapply("attr<-",SORTIE[[JBD]],"type",LESTYPES,SIMPLIFY=FALSE)
      class(SORTIE[[JBD]]) <- c("summary.mccut")
      }

    else if(inherits(PREMS[[JBD]],"mcnode")) {
      attr(SORTIE[[JBD]],"type") <- attr(PREMS[[JBD]],"type")
      attr(SORTIE[[JBD]],"outm") <- attr(PREMS[[JBD]],"outm")
      class(SORTIE[[JBD]]) <- c("mcnode")
      }

    else if(inherits(PREMS[[JBD]],"tornado")) {
      PREMS[[JBD]]$value <- SORTIE[[JBD]]
      SORTIE[[JBD]] <- PREMS[[JBD]]
      class(SORTIE[[JBD]]) <- c("tornado.mccut")
      }
      
    else if(inherits(PREMS[[JBD]],"plotmc")){
      TYPEPLOT <- sapply(PREMS[[JBD]],"attr",which="type")
      SORTIE[[JBD]] <- mapply("attr<-",SORTIE[[JBD]],"type",TYPEPLOT,SIMPLIFY=FALSE)
      class(SORTIE[[JBD]]) <- c("plot.mccut")}

  }

  class(SORTIE) <- "mccut"

SORTIE}, silent=TRUE)     # Fin du try
ndvar(OLDV)
ndunc(OLDU)
cat("\n") 
if(inherits(RESULTAT,"try-error")) stop(RESULTAT,call. = FALSE)
return(invisible(RESULTAT))
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

