#
#  multinomRob
#
#  Walter R. Mebane, Jr.
#  University of Michigan
#  http://www-personal.umich.edu/~wmebane
#  <wmebane@umich.edu>
#
#  Jasjeet Singh Sekhon 
#  UC Berkeley
#  http://sekhon.polisci.berkeley.edu
#  sekhon@berkeley.edu
#
#
###################################
#New Front End for Genoud, with tuned defaults
###################################

#sets genoud.parms defaults
genoudParms  <- function(genoud.parms)
  {
    #set user controlled defaults
    if (!is.null(genoud.parms$cluster)) {
      if (length(genoud.parms$cluster) > 1) {
        warning("cluster option cannot be used with 'multinomRob'")
      } else if (is.list(genoud.parms$cluster)) {
        warning("cluster option cannot be used with 'multinomRob'")
      } else if (genoud.parms$cluster!=FALSE) {
        warning("cluster option cannot be used with 'multinomRob'")
      }
    }
    genoud.parms$cluster  <- FALSE;

    if (is.null(genoud.parms$balance))
      genoud.parms$balance  <- FALSE  ;
    
    if (is.null(genoud.parms$pop.size))
      genoud.parms$pop.size  <- 1000;

    if (is.null(genoud.parms$max.generations))
      genoud.parms$max.generations  <- 100;
    
    if (is.null(genoud.parms$wait.generations))
      genoud.parms$wait.generations  <- 10;
    
    if (is.null(genoud.parms$hard.generation.limit))
      genoud.parms$hard.generation.limit  <- FALSE;

    #this is redundant, but maintains clarity  
    if (is.null(genoud.parms$MemoryMatrix))
      genoud.parms$MemoryMatrix  <- TRUE;
  
    if (is.null(genoud.parms$Debug))
      genoud.parms$Debug  <- FALSE ;

    #this is redundant, but maintains clarity   
    if (is.null(genoud.parms$Domains))
      genoud.parms$Domains  <- NULL;

    if (is.null(genoud.parms$scale.domains))
      genoud.parms$scale.domains  <- 10;
    
    if (is.null(genoud.parms$boundary.enforcement))
      genoud.parms$boundary.enforcement  <- 0;
    
    if (is.null(genoud.parms$solution.tolerance))
      genoud.parms$solution.tolerance  <- 0.0000001;
    
    if (is.null(genoud.parms$BFGS))
      genoud.parms$BFGS  <- TRUE;
  
    if (is.null(genoud.parms$unif.seed))
      genoud.parms$unif.seed  <- 812821;
    
    if (is.null(genoud.parms$int.seed))
      genoud.parms$int.seed  <- 53058;
    
    if (is.null(genoud.parms$print.level))
      genoud.parms$print.level  <- 0;
    
    if (is.null(genoud.parms$share.type))
      genoud.parms$share.type  <- 0;
    
    if (is.null(genoud.parms$instance.number))
      genoud.parms$instance.number  <- 0;
    
    if (is.null(genoud.parms$output.path))
      genoud.parms$output.path  <- "stdout";
    
    if (is.null(genoud.parms$output.append))
      genoud.parms$output.append  <- FALSE;
    
    if (is.null(genoud.parms$project.path))
      genoud.parms$project.path  <- "/dev/null";
    
    if (is.null(genoud.parms$P1))
      genoud.parms$P1  <- 50;
    
    if (is.null(genoud.parms$P2))
      genoud.parms$P2  <- 50;
    
    if (is.null(genoud.parms$P3))
      genoud.parms$P3  <- 50;
    
    if (is.null(genoud.parms$P4))
      genoud.parms$P4  <- 50;
    
    if (is.null(genoud.parms$P5))
      genoud.parms$P5  <- 50;
    
    if (is.null(genoud.parms$P6))
      genoud.parms$P6  <- 50;
    
    if (is.null(genoud.parms$P7))
      genoud.parms$P7  <- 50;
    
    if (is.null(genoud.parms$P8))
      genoud.parms$P8  <- 50;
    
    if (is.null(genoud.parms$P9))
      genoud.parms$P9  <- 0  ;

    return(genoud.parms);
  } #end genoudParms


genoudRob <- function(fn,nvars,starting.values,genoud.parms)
{
  #new options for genoud > 2.0 which are needed
  lexical=FALSE
  cluster  <- genoud.parms$cluster
  balance  <- genoud.parms$balance
  
  #set static defaults
  max  <- FALSE
  gradient.check  <- FALSE
  data.type.int  <- FALSE
  hessian  <- FALSE  

  #load up genoud.parms
  pop.size  <- genoud.parms$pop.size;
  max.generations  <- genoud.parms$max.generations;
  wait.generations  <- genoud.parms$wait.generations;
  hard.generation.limit  <- genoud.parms$hard.generation.limit;
  MemoryMatrix  <- genoud.parms$MemoryMatrix;
  Debug  <- genoud.parms$Debug;
  Domains  <- genoud.parms$Domains;
  scale.domains  <- genoud.parms$scale.domains;
  boundary.enforcement  <- genoud.parms$boundary.enforcement;
  solution.tolerance  <- genoud.parms$solution.tolerance;
  BFGS  <- genoud.parms$BFGS;
  unif.seed  <- genoud.parms$unif.seed;
  int.seed  <- genoud.parms$int.seed;
  print.level  <- genoud.parms$print.level;
  share.type  <- genoud.parms$share.type;
  instance.number  <- genoud.parms$instance.number;
  output.path  <- genoud.parms$output.path;
  output.append  <- genoud.parms$output.append;
  project.path  <- genoud.parms$project.path;
  P1  <- genoud.parms$P1;
  P2  <- genoud.parms$P2;
  P3  <- genoud.parms$P3;
  P4  <- genoud.parms$P4;
  P5  <- genoud.parms$P5;
  P6  <- genoud.parms$P6;
  P7  <- genoud.parms$P7;
  P8  <- genoud.parms$P8;
  P9  <- genoud.parms$P9;

  if (!(is.matrix(Domains)))
    {
      Domains <- matrix(nrow=nvars, ncol=2);
      for (i in 1:nvars)
        {
          Domains[i,1] <- starting.values[i] - abs(starting.values[i])*scale.domains;
          Domains[i,2] <- starting.values[i] + abs(starting.values[i])*scale.domains;
        } # end of for loop
    } # end of Domains if

  ret = genoud(fn, nvars=nvars, max=max, pop.size=pop.size, max.generations=max.generations, wait.generations=wait.generations,
               hard.generation.limit=hard.generation.limit, starting.values=starting.values, MemoryMatrix=MemoryMatrix, 
               Domains=Domains, solution.tolerance=solution.tolerance,
               gr=NULL, boundary.enforcement=boundary.enforcement, lexical=lexical, gradient.check=gradient.check, BFGS=BFGS, 
               data.type.int=data.type.int, hessian=hessian, unif.seed=unif.seed, int.seed=int.seed,
               print.level=print.level, share.type=share.type, instance.number=instance.number,
               output.path=output.path, output.append=output.append, project.path=project.path,
               P1=P1, P2=P2, P3=P3, P4=P4, P5=P5, P6=P6, P7=P7, P8=P8, P9=P9,
               cluster=cluster, balance=balance, debug=Debug)

  return(ret)
} #end of genoudRob()
