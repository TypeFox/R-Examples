#  Allan Strand 9/27/01
#
#

###
### coerce landscape tries to make sure that the landscape object is in a form that is palatable to 
### C++ code in rmetasim 
###

landscape.coerce <- function(rland,noind=F)
  {
    rland$intparam <- lapply(rland$intparam,as.integer)
    rland$switchparam <- lapply(rland$switchparam,as.integer)
    rland$floatparam <- lapply(rland$floatparam,as.numeric)
    
    rland$loci <- lapply(rland$loci,function(x)
                         {
                           x$type <- as.integer(x$type)
                           x$ploidy <- as.integer(x$ploidy)
                           x$trans <- as.integer(x$trans)
                           x$rate <- as.numeric(x$rate)
                           x$alleles <- lapply(x$alleles,function(y,typ)
                                               {
                                                 y$aindex <- as.integer(y$aindex)
                                                 y$birth <- as.integer(y$birth)
                                                 y$prop <- as.numeric(y$prop)
                                                 if (typ!=253)
                                                   {
                                                     y$state <- as.integer(y$state)
                                                   }
                                                 y
                                               },
                                               typ=x$type)
                           x
                         })
    if (!rland$switchparam$densdepdemo)
      {
        rland$demography$localdemK <- rland$demography$localdem
      }
    if (!noind)  {rland$individuals <- matrix(as.integer(rland$individuals),nrow=dim(rland$individuals)[1])}
    rland
  }
####
#default seed is same as calling environment.  Zero or any positive integer
#is assigned to R's random number seed.  The type of RNG is inherited from the
#calling environment
#
landscape.simulate <- function(Rland, numit, seed=-1, compress=FALSE, adj.lambda=0)
  {
    if (is.landscape(Rland))
      {
        if (!(seed<0))
          {
            set.seed(seed)
          }
        Rland <- landscape.coerce(Rland)
        .Call("iterate_landscape",as.integer(numit),Rland,as.integer(compress),as.integer(adj.lambda),PACKAGE = "rmetasim")
      }
    else
      {
        print("Rland not a landscape object...exiting")
      }
  }


landscape.survive <- function(Rland, seed=-1)
  {
    if (is.landscape(Rland))
      {
        if (!(seed<0))
          {
            set.seed(seed)
          }
        Rland <- landscape.coerce(Rland)
        .Call("survive_landscape",Rland,PACKAGE = "rmetasim")
      }
    else
      {
        print("Rland not a landscape object...exiting")
      }
  }

landscape.reproduce <- function(Rland, seed=-1)
  {
    if (is.landscape(Rland))
      {
        if (!(seed<0))
          {
            set.seed(seed)
          }
        Rland <- landscape.coerce(Rland)
        .Call("reproduce_landscape",Rland,PACKAGE = "rmetasim")
      }
    else
      {
        print("Rland not a landscape object...exiting")
      }
  }


landscape.carry <- function(Rland, seed=-1)
  {
    if (is.landscape(Rland))
      {
        if (!(seed<0))
          {
            set.seed(seed)
          }
        Rland <- landscape.coerce(Rland)
        .Call("carry_landscape",Rland,PACKAGE = "rmetasim")
      }
    else
      {
        print("Rland not a landscape object...exiting")
      }
  }

landscape.extinct <- function(Rland, seed=-1)
  {
    if (is.landscape(Rland))
      {
        if (!(seed<0))
          {
            set.seed(seed)
          }
        Rland <- landscape.coerce(Rland)
        .Call("extinct_landscape",Rland,PACKAGE = "rmetasim")
      }
    else
      {
        print("Rland not a landscape object...exiting")
      }
  }

landscape.advance <- function(Rland, seed=-1)
  {
    if (is.landscape(Rland))
      {
        if (!(seed<0))
          {
            set.seed(seed)
          }
        Rland <- landscape.coerce(Rland)
        .Call("advance_landscape",Rland,PACKAGE = "rmetasim")
      }
    else
      {
        print("Rland not a landscape object...exiting")
      }
  }

landscape.compress <- function(Rland)
    {
    if (is.landscape(Rland))
      {
        Rland <- landscape.coerce(Rland)
        .Call("compress_landscape",Rland,PACKAGE = "rmetasim")
      }
    else
      {
        print("Rland not a landscape object...exiting")
        Rland
      }
  }
