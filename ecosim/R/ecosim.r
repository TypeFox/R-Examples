# ==============================================================================
#
# ECOSIM - An R library for didactical ecological model simulations
# -----------------------------------------------------------------
#
# Class definitions                First version: Peter Reichert, Nov.  12, 2005
# -----------------                Last revision: Peter Reichert, Feb.  01, 2014
#
# ==============================================================================


# Load libraries:
# ===============

library(deSolve)


# ==============================================================================
# Class for biogeochemical transformation process
# ==============================================================================

# First version: Peter Reichert, Nov.  12, 2005
# Last revision: Peter Reichert, Jan.  07, 2007

# An object of this class defines a transformation process by a common
# transformation rate and substance-specific stoichiometric coefficients.

# A process is defined by 
#   name:     String defining the name of the process
#   rate:     Expression characterizing the transformation rate as a function
#             of substance or organism concentrations, model parameters,
#             environmental conditions, and time. All these variables are
#             identified by their name; the name "t" is reserved for time.
#             For the definition of parameters and environmental conditions
#             see member function "calcres" of the class "system";
#             concentration and time are calculated and provided internally.
#   stoich:   List of numerics or expressions for stoichiometric coefficients  
#             defining the relative transformation rates of different substances
#             or organisms affected by the process. The contribution of the 
#             process to the transformation rate of a substance is given as 
#             the product of the rate (see above) times the substance-specific 
#             stoichiometric coefficient. The stoichiometric coefficient can  
#             depend on the same variables as the rate (see above).
#   pervol:   Logical variable equal to TRUE if rate is per volume, FALSE if it
#             is per surface area.
# Note that surface to volume conversions are taken into account automatically
# in the member function "calc.rates.statevar.reactor" of the class "reactor"
# and must not be specified in the process stoichiometry. However, it must be
# specified if the provided rate is a rate per volume or per unit of 
# surface area.
#
# Available member functions for internal use:
#   calc.trans.rates(process,param,cond,conc,t):
#             Calculation of transformation rate contributions of a process
#             for given model parameters, environmental conditions, sub-
#             stance concentrations, and time.


# Class definition:
# =================

setClass(Class           = "process",
         representation 
          = list(name    = "character",   # name of process
                 rate    = "expression",  # process rate
                 stoich  = "list",        # list of stoichiometric
                                          # coefficients (identified by name)
                 pervol  = "logical"),    # type of rate (TRUE: per volume,
                                          # FALSE: per area)
         prototype
          = list(pervol  = TRUE))         # default rate type: per volume


# Method for calculating transformation rates:
# ============================================

# Calculation of transformation rate contributions of a process given 
#   param:    Model parameters.
#   cond:     Environmental conditions.
#   conc:     Substance concentrations.
#   t:        time.
# The process rate and stoichiometry are considered as specified in the 
# process definition.

calc.trans.rates <- function(proc,param,cond,conc,t) {return(NULL)}
setGeneric("calc.trans.rates")

setMethod(f          = "calc.trans.rates",
          signature  = "process",
          definition = function(proc,param,cond,conc,t)
                       {
                          # define environment:
                          # -------------------
                          
                          env = c(param,cond,as.list(conc),list(t=t))

                          # calculate rates as product of common process rate
                          # with substance specific stoichiometric coefficient:
                          # ---------------------------------------------------
                          
                          trans.rates <- rep(0,length(proc@stoich))
                          if ( length(proc@stoich) > 0 )
                          {
                             names(trans.rates) <- names(proc@stoich)
                             for ( i in 1:length(proc@stoich) )
                             {
                                trans.rates[i] <- eval(expr  = proc@stoich[[i]],
                                                       envir = env)
                             }
                             trans.rates <- trans.rates *
                                            eval(expr  = proc@rate,
                                                 envir = env)
                          }
                            
                          # return vector of calculated transformation rates:
                          # -------------------------------------------------
                          
                          return(trans.rates)
                       })



# ==============================================================================
# Class for mixed reactor
# ==============================================================================

# First version: Peter Reichert, Nov.  12, 2005
# Last revision: Peter Reichert, Jan.  07, 2006

# An object of this class defines a mixed reactor with substances or organisms
# dissolved or suspended in the water phase or attached to a surface in the
# reactor. The definition includes substance input, in- and outflow and the
# list of transformation processes active in the reactor.

# A mixed reactor of variable volume is defined by
#   name:             String defining the name of the reactor.
#   volume.ini:       Initial volume of the reactor. Starting from this 
#                     value, the volume is calculated dynamically according
#                     to the difference in inflow and outflow. Specify the
#                     same expression for outflow as for inflow if you want to
#                     keep the volume constant.
#   area:             Surface area of the reactor. This value is kept constant.
#                     It is required to convert surface densities of attached
#                     substances or organisms to total masses in the reactor
#                     and to convert transformation rates from rates per volume
#                     to rates per unit of surface area.
#   conc.pervol.ini:  Vector of initial concentrations of substances or 
#                     organisms dissolved or suspended in the reactor. The
#                     substances ae identified by their name. This vector
#                     determines the mass balance equations of which substances 
#                     are solved in the reactor. This means that initial
#                     concentrations of all substances to be calculated must
#                     be specified.
#   conc.perarea.ini: Vector of initial concentrations of substances or
#                     organisms attached to a surface in the reactor. The
#                     substances are identified by their name. This vector
#                     determines the mass balance equations of which substances 
#                     are solved in the reactor. This means that initial 
#                     concentrations of all substances to be calculated must 
#                     be specified. Attached substances must have different 
#                     names than dissolved or suspended substances in order
#                     to distinguish them in process rates.
#   input             Vector of net input fluxes of dissolved or suspended 
#                     substances into the reactor as mass per time. This  
#                     variable is used to describe exchange of substances with 
#                     the environment that are not associated with the  
#                     volumetric exchange (e.g. gas exchange across a surface).
#   inflow            Expression specifying the volumetric inflow into the 
#                     reactor.
#   inflow.conc       Vector of expressions specifying the concentrations of  
#                     dissolved or suspended substances in the inflow to the
#                     reactor. The substances are identified by their name.
#   outflow           Expression specifying the volumetric outflow of the 
#                     reactor. Specify the same expression as for inflow if you
#                     want to keep the volume constant.
#   cond              List of numerics or expressions specifying the local 
#                     environmental conditions to which the reactor is exposed
#   processes         List of processes (of class "process" defined above)
#                     that are active in the reactor.
#   
# Available member functions for internal use:
#   calc.rates.statevar.reactor(reactor,param,cond,conc,t):
#                     Calculation of rates of state variables in the reactor
#                     (reactor volume, masses of dissolved or suspended 
#                     substances or organisms, and masses of attached 
#                     substances or organisms) for given model parameters, 
#                     environmental conditions, substance concentrations, and
#                     time.


# Class definition:
# =================

setClass(Class                    = "reactor",
         representation
          = list(name             = "character",    # name of reactor
                 volume.ini       = "expression",   # initial volume
                 area             = "expression",   # surface area for attached
                                                    # substances or organisms
                 conc.pervol.ini  = "list",         # list of initial con-
                                                    # centrationstions of   
                                                    # dissolved or suspened  
                                                    # substances or organisms
                 conc.perarea.ini = "list",         # list of initial con-
                                                    # centrations of attached  
                                                    # substances or organisms
                 input            = "list",         # list of net input fluxs
                                                    # of dissolved or suspended
                                                    # substances or organisms 
                                                    # into the reactor
                 inflow           = "expression",   # inflow to reactor
                 inflow.conc      = "list",         # vector of concentrations
                                                    # of dissolved or suspended 
                                                    # substancees or organsims
                                                    # in the inflow
                 outflow          = "expression",   # outflow of reactor
                 cond             = "list",         # vector of environmental
                                                    # conditions
                 processes        = "list",         # list of active processes
                                                    # in the reactor (objects
                                                    # of the class "process"
                                                    # defined above)
                 a                = "numeric"),     # evaluated area
                                                    # (for internal use only)
         prototype
          = list(volume.ini       = expression(1),  # default initial volume: 1
                 area             = expression(1))) # default surface area:   1
                               

# Method for calculating rates of change of masses of substances in the reactor:
# ==============================================================================

# Calculation of rates of state variables in the reactor given 
#   param:    Model parameters.
#   cond:     Environmental conditions.
#   volume:   Current volume of the reactor.
#   conc:     Substance concentrations.
#   t:        Current time.
# The rates are calculated under consideration of input, inflow, outflow, and
# active processes specified in the reactor definition.

calc.rates.statevar.reactor <- function(reactor,param,cond.global,volume,conc,t) 
                               {return(NULL)}
setGeneric("calc.rates.statevar.reactor")

setMethod(f          = "calc.rates.statevar.reactor",
          signature  = "reactor",
          definition = function(reactor,param,cond.global,volume,conc,t)
                       {
                          # evaluate local environmental conditions:
                          # ----------------------------------------

                          env = c(param,cond.global,as.list(conc),list(t=t))
                          cond.local <- list()
                          if ( length(reactor@cond) > 0 )
                          {
                             for ( l in 1:length(reactor@cond) )
                             {
                                cond.l <- 
                                   eval(expr  = reactor@cond[[l]],
                                        envir = env)
                                cond.local <- c(cond.local,cond.l)
                             }
                             names(cond.local) <- names(reactor@cond)
                             n <- names(env)
                             env <- c(env,list(cond.l))
                             names(env) <- c(n,names(reactor@cond[[l]]))
                          }

                          # define environment:
                          # -------------------

                          env = c(param,cond.local,cond.global,
                                  as.list(conc),list(t=t))

                          # calculate inflow and outflow and change of volume:
                          # --------------------------------------------------

                          Q.in <- 0
                          if ( length(reactor@inflow) > 0 )
                          {                      
                             Q.in  <- eval(expr  = reactor@inflow,
                                           envir = env)
                          }
                          Q.out <- 0
                          if ( length(reactor@outflow) > 0 )
                          {                      
                             Q.out <- eval(expr  = reactor@outflow,
                                           envir = env)
                          }
                          rate.volume <- Q.in - Q.out

                          # initialize rate vector of mass changes to zero:
                          # -----------------------------------------------
                          
                          rates.mass <- rep(0,length(reactor@conc.pervol.ini)+
                                              length(reactor@conc.perarea.ini))
                          if ( length(rates.mass) > 0 )
                          {
                             names(rates.mass) <- 
                                c(names(reactor@conc.pervol.ini),
                                  names(reactor@conc.perarea.ini))
                                  
                             # add transformation rates:
                             # -------------------------

                             if ( length(reactor@processes) > 0 )
                             {
                                for ( i in 1:length(reactor@processes) )
                                {
                                   trans.rates <- 
                                      calc.trans.rates(reactor@processes[[i]],param,
                                                       c(cond.global,cond.local),
                                                       conc,t)
                                   if ( reactor@processes[[i]]@pervol )
                                   {
                                      trans.rates <- trans.rates*volume
                                   }
                                   else
                                   {
                                      trans.rates <- trans.rates*reactor@a
                                   }
                                   for ( j in 1:length(trans.rates) )
                                   {
                                      if ( !is.na(match(names(trans.rates)[j],
                                                        names(rates.mass))) )
                                      {
                                         rates.mass[names(trans.rates)[j]] <- 
                                            rates.mass[names(trans.rates)[j]] +
                                            trans.rates[j]
                                      }
                                   }
                                }
                             }

                             # add rates of change due to net input:
                             # -------------------------------------
                          
                             if ( length(reactor@input) > 0 )
                             {
                                for ( j in 1:length(reactor@input) )
                                {
                                   if ( !is.na(match(names(reactor@input)[[j]],
                                              names(reactor@conc.pervol.ini))) )
                                   {
                                      inp <- eval(expr  = reactor@input[[j]],
                                                  envir = env)
                                      rates.mass[names(reactor@input)[j]] <- 
                                         rates.mass[names(reactor@input)[j]]+inp
                                   }
                                }
                             }
                             
                             # add rates of change due to inflow:
                             # ----------------------------------
                          
                             if ( length(reactor@inflow.conc) > 0 )
                             {
                                for ( j in 1:length(reactor@inflow.conc) )
                                {
                                   if ( !is.na(match(
                                              names(reactor@inflow.conc)[j],
                                              names(reactor@conc.pervol.ini))) )
                                   {
                                     inp <- eval(expr  = reactor@inflow.conc[[j]],
                                                 envir = env) * Q.in
                                     rates.mass[names(reactor@inflow.conc)[j]]<- 
                                      rates.mass[names(reactor@inflow.conc)[j]]+
                                      inp
                                   }
                                }
                             }
                             
                             # subtract rates of change due to outflow:
                             # ----------------------------------------
                             
                             if ( length(reactor@conc.pervol.ini) > 0 )
                             {
                                for ( j in 1:length(reactor@conc.pervol.ini) )
                                {
                                   if ( is.na(match(
                                              names(reactor@conc.pervol.ini)[j],
                                              names(conc))) )
                                   {
                                      stop(
                                        paste("mass of substance",
                                          names(reactor@conc.pervol.ini)[j],
                                          "not found as reactor state variable",
                                          sep=" "))
                                   }
                                   else
                                   {
                                      rates.mass[j] = rates.mass[j] - 
                                      Q.out*
                                         conc[names(reactor@conc.pervol.ini)[j]]
                                   }
                                }
                             }
                          }
                          
                          # return vector of calculated rates of change: 
                          # --------------------------------------------

                          return(c(rate.volume,rates.mass))
                       })



# ==============================================================================
# Class for link between mixed reactors
# ==============================================================================

# First version: Peter Reichert, April 09, 2006
# Last revision: Peter Reichert, Feb.  13, 2013


setClass(Class             = "link",
         representation 
         = list(name       = "character",   # name of advective link
                from       = "character",   # name of reactor from which the 
                                            # link starts
                to         = "character",   # name of reactor at which the 
                                            # link ends
                flow       = "expression",  # water flow between reactors
                                            # carrying dissolved or suspended
                                            # substances; the following
                                            # specifications result in fluxes
                                            # not associated with water flow
                qadv.gen   = "expression",  # general transfer coefficient for
                                            # all dissolved or suspended 
                                            # substances
                                            # (F=qadv*C_from if qadv>0, 
                                            #  F=qadv*C_to   if qadv<0)
                qadv.spec  = "list",        # list of substance-specific
                                            # transfer coeff.
                qdiff.gen  = "expression",  # general exchange coefficient for
                                            # all dissolved or suspended 
                                            # substances
                                            # (F=qdiff*(C_to-C_from)
                qdiff.spec = "list"))       # list of substance-specific 
                                            # exchange coefficients


calc.rates.statevar.link <- function(link,param,cond.global,conc.from,conc.to,t)
                            {return(NULL)}
setGeneric("calc.rates.statevar.link")

setMethod(f          = "calc.rates.statevar.link",
          signature  = "link",
          definition = function(link,param,cond.global,conc.from,conc.to,t)
                       {
                          # define environment:
                          # -------------------

                          env = c(param,cond.global,as.list(conc.from),
                                  list(t=t))

                          # calculate flow:
                          # ---------------
                          
                          Q <- 0
                          if ( length(link@flow) > 0 )
                          {
                             Q <- eval(expr  = link@flow,
                                       envir = env)
                          }
                          rate.volume <- Q
                          names(rate.volume) <- "Q"

                          # initiate mass flux vector:
                          # --------------------------
                                                    
                          names <- unique(c(names(conc.from),names(conc.to)))
                          rates.mass <- rep(0,length(names))
                          names(rates.mass) <- names

                          # add fluxes associated with water flow:
                          # --------------------------------------
                          
                          if ( Q > 0 )
                          {
                             if ( length(conc.from) > 0 )
                             {
                                for ( i in 1:length(names(conc.from)) )
                                {
                                   rates.mass[names(conc.from)[i]] <-
                                      rates.mass[names(conc.from)[i]] +
                                                                Q * conc.from[i]
                                }
                             }
                          }
                          else
                          {
                             if ( Q < 0 )
                             {
                                if ( length(conc.to) > 0 )
                                {
                                   for ( i in 1:length(names(conc.to)) )
                                   {
                                      rates.mass[names(conc.to)[i]] <-
                                         rates.mass[names(conc.to)[i]] +
                                                                  Q * conc.to[i]
                                   }
                                }
                             }
                          }
                          
                          # add general advective fluxes:
                          # -----------------------------
                          
                          if ( length(link@qadv.gen) > 0 )
                          {
                             qadv <- eval(expr  = link@qadv.gen,
                                          envir = env)
                             if ( qadv > 0 )
                             {
                                if ( length(conc.from) > 0 )
                                {
                                   for ( i in 1:length(names(conc.from)) )
                                   {
                                      rates.mass[names(conc.from)[i]] <-
                                         rates.mass[names(conc.from)[i]] +
                                                             qadv * conc.from[i]
                                   }
                                }
                             }
                             else
                             {
                                if ( qadv < 0 )
                                {
                                   if ( length(conc.to) > 0 )
                                   {
                                      for ( i in 1:length(names(conc.to)) )
                                      {
                                         rates.mass[names(conc.to)[i]] <-
                                            rates.mass[names(conc.to)[i]] +
                                                               qadv * conc.to[i]
                                      }
                                   }
                                }
                             }
                          }
                          
                          # add substance-specific advective fluxes:
                          # ----------------------------------------

                          if ( length(link@qadv.spec) > 0 )
                          {
                             for ( i in 1:length(link@qadv.spec) )
                             {
                                name <- names(link@qadv.spec)[i]
                                qadv <- eval(expr  = link@qadv.spec[[i]],
                                             envir = env)
                                if ( qadv > 0 )
                                {
                                   if ( length(conc.from) > 0 )
                                   {
                                      if (!is.na(match(name,names(conc.from))) )
                                      {
                                         rates.mass[name] <-
                                            rates.mass[name] +
                                                          qadv * conc.from[name]
                                      }
                                   }
                                }
                                else
                                {
                                   if ( qadv < 0 )
                                   {
                                      if ( length(conc.to) > 0 )
                                      {
                                         if (!is.na(match(name,names(conc.to))))
                                         {
                                            rates.mass[name] <-
                                               rates.mass[name] + 
                                                            qadv * conc.to[name]
                                         }
                                      }
                                   }
                                }
                             }
                          }

                          # add general diffusive fluxes:
                          # -----------------------------
                          
                          if ( length(link@qdiff.gen) > 0 )
                          {
                             if ( length(conc.from) > 0 & length(conc.to) > 0 )
                             {
                                qdiff <- eval(expr  = link@qdiff.gen,
                                              envir = env)
                                for ( i in 1:length(names(conc.from)) )
                                {
                                   name <- names(conc.from)[i]
                                   if ( !is.na(match(name,names(conc.to))) )
                                   {
                                      rates.mass[name] <-
                                         rates.mass[name] +
                                            qdiff *
                                            (conc.from[name] - conc.to[name])
                                   }
                                }
                             }
                          }

                          # add substance-specific diffusive fluxes:
                          # ----------------------------------------

                          if ( length(link@qdiff.spec) > 0 )
                          {
                             if ( length(conc.from) > 0 & length(conc.to) > 0 )
                             {
                                for ( i in 1:length(link@qdiff.spec) )
                                {
                                   name <- names(link@qdiff.spec)[i]
                                   qdiff <- eval(expr  = link@qdiff.spec[[i]],
                                                 envir = env)
                                   if ( !is.na(match(name,names(conc.from))) &
                                        !is.na(match(name,names(conc.to))) )
                                   {
                                      rates.mass[name] <-
                                         rates.mass[name] +
                                            qdiff*
                                            ( conc.from[name] - conc.to[name] )
                                   }
                                }
                             }
                          }
                          
                          return(c(rate.volume,rates.mass))
                       })



# ==============================================================================
# Class for system of mixed reactors
# ==============================================================================

# First version: Peter Reichert, Nov.  12, 2005
# Last revision: Peter Reichert, Jan.  26, 2014

# An object of this class defines a system consisting of mixed reactors.
# Through the definition of the reactors, it includes the substances or
# organisms contained in the reactors, their external in- and outflows,
# and their transformation processes. In addition, links can be specified 
# to describe interactions between the reactors.
# This makes an object of this class define a complete model.
# Consequently, the class contains a member function to perform dynamic 
# simulations of the model.


# A system of mixed reactors is defined by
#   name:             String defining the name of the system.
#   cond:             Vector of expressions defining global environmental
#                     conditions. 
#   reactors:         List of rectors building the system (objects of the class
#                     "reactor" defined above). 
#   links:            List of links connecting the reactors (objects 
#                     of the class "link" defined above). 
#   
# Available member functions for external use:
#   calcres(system,method="lsoda",...): 
#                     Calculation of results (volume and concentration time
#                     series for the system) in the form of a matrix.


# Class definition:
# =================

setClass(Class             = "system",
         representation
          = list(name      = "character",  # name of system
                 reactors  = "list",       # list of reactors
                 links     = "list",       # list of links connecting the
                                           # reactors
                 cond      = "list",       # list of global environmental
                                           # conditions (numeric or expression)
                 param     = "list",       # list of model parameters (numeric)
                 t.out     = "numeric"),   # output times
         prototype
          = list(t.out     = seq(0,365,by=1)))  # default output times


# Method for calculating concentration time series:
# =================================================
       
# Calculation of model results  
#   system:   Object of type system defining the model.
#   method:   Integration algorithm passed to function ode from deSolve.
#   ...       Further arguments passed to ode.
# The rates are calculated under consideration of input, inflow, outflow, and
# active processes specified in the reactor definition.

calcres <- function(system,method="lsoda",...) {return(NULL)}
setGeneric("calcres")


# Define auxiliary function for getting index:
# --------------------------------------------

get.reactor.index <- function(name,reactors)
{
   if ( length(reactors) < 1 ) return(NA)
   {
      for ( i in 1:length(reactors) )
      {
         if ( identical(name,reactors[[i]]@name) ) return(i)
      }
   }
   return(NA)
}

# Define auxiliary function to interpolate time-dependent parameters:
# -------------------------------------------------------------------

get.param.val <- function(param,t)
{
   param.val <- param
   for ( i in 1:length(param) )
   {
     if ( length(param[[i]]) != 1 ) param.val[[i]] <- approx(x=param[[i]],xout=t,rule=2)$y
   }
   return(param.val)
}

      
# Define auxiliary function for right hand side of differential equation:
# -----------------------------------------------------------------------

system.rhs <- function(t,x,param,system)
{
   # state variables are volume, masses of dissolved or suspended substances
   # or organisms, and masses of attached substances or organisms.
   
   # interpolate time-dependent parameters:
   # --------------------------------------
  
   param.val <- get.param.val(system@param,t)
  
   # combine rates of all reactors:
   # ------------------------------
  
   num.react  <- length(system@reactors)
   offsets    <- numeric(num.react)
   num.vol    <- numeric(num.react)
   offsets[1] <- 0
   rhs <- numeric(0)
   conc <- list()
   for ( k in 1:num.react )
   {
      # calculate concentrations:
      # -------------------------
      
      reactor <- system@reactors[[k]]
      volume <- x[offsets[k]+1]
      num.vol[k]  <- length(reactor@conc.pervol.ini)
      num.area <- length(reactor@conc.perarea.ini)
      conc.pervol <- numeric(0)
      if ( num.vol[k] > 0 )
      {
         conc.pervol  <- x[offsets[k]+1+(1:num.vol[k])]/x[offsets[k]+1]
      }
      conc.perarea <- numeric(0)
      if ( num.area > 0 )
      {
         ind <- seq(offsets[k]+1+num.vol[k]+1,
                    offsets[k]+1+num.vol[k]+num.area,
                    length.out = num.area)
         conc.perarea <- x[ind]/reactor@a
      }
      conc[[k]] <- c(conc.pervol,conc.perarea)

      # evaluate environmental conditions
      # ---------------------------------

      cond.global <- list()
      if ( length(system@cond) > 0 )
      {
         env <- c(param.val,as.list(conc[[k]]),list(t=t))
         for ( l in 1:length(system@cond) )
         {
            cond.l <- eval(expr  = system@cond[[l]],
                           envir = env)
            cond.global <- c(cond.global,cond.l)
            n <- names(env)
            env <- c(env,list(cond.l))
            names(env) <- c(n,names(system@cond[[l]]))
         }
         names(cond.global) <- names(system@cond)
      }
      
      # calculate rates within individual reactors:
      # -------------------------------------------
      
      rhs.current <- calc.rates.statevar.reactor(reactor,param.val,cond.global,
                                                 volume,conc[[k]],t)
      
      # combine rates of current reactor with rates of previous reactors:
      # -----------------------------------------------------------------
      
      rhs <- c(rhs,rhs.current)
      
      # calculate new index offset:
      # ---------------------------
      
      if ( k < num.react ) { offsets[k+1] <- offsets[k]+1+num.vol[k]+num.area }
   }
   
   # add rates from links:
   # ---------------------

   if ( length(system@links) > 0 )
   {
      for ( j in 1:length(system@links) )
      {
         # calculate rate contributions:
         # -----------------------------

         link <- system@links[[j]]
         ind.from <- get.reactor.index(link@from,system@reactors)
         ind.to   <- get.reactor.index(link@to,system@reactors)
         offset.from <- offsets[ind.from]
         offset.to   <- offsets[ind.to]
         rhs.contrib <- 
            calc.rates.statevar.link(link,param.val,cond.global,
                                     conc[[ind.from]][1:num.vol[ind.from]],
                                     conc[[ind.to]][1:num.vol[ind.to]],t)

         if ( length(rhs.contrib) > 0 )
         {
            # add water flow:
            # ---------------

            rhs[offset.from+1] <- rhs[offset.from+1] - rhs.contrib[1]
            rhs[offset.to+1]   <- rhs[offset.to+1]   + rhs.contrib[1]
         
            # add mass fluxes:
            # ----------------
         
            if ( length(rhs.contrib) > 1 )
            {
               for ( i in 2:length(rhs.contrib) )
               {
                  ind <- offset.from + 1 +
                         match(names(rhs.contrib)[i],
                               names(rhs)[(offset.from+2):
                                          (offset.from+1+num.vol[ind.from])])
                  if ( !is.na(ind) )
                  {
                     rhs[ind] = rhs[ind] - rhs.contrib[i]
                  }
                  else
                  {
                     stop(paste("Substance",
                                names(rhs.contrib)[i],
                                "not found in reactor",
                                link@from,
                                "from which the link starts"))
                  }
                  ind <- offset.to + 1 +
                         match(names(rhs.contrib)[i],
                               names(rhs)[(offset.to+2):
                                          (offset.to+1+num.vol[ind.to])])

                  if ( !is.na(ind) )
                  {
                     rhs[ind] = rhs[ind] + rhs.contrib[i]
                  }
                  else
                  {
                     stop(paste("Substance",
                                names(rhs.contrib)[i],
                                "not found in reactor",
                                link@to,
                                "at which the link ends"))
                  }
               }
            }
         }
      }
   }
      
   # return aggregate rates:
   # -----------------------

   return(list(rhs))
}


setMethod(f          = "calcres",
          signature  = "system",
          definition = function(system,method="lsoda",...)
                       {
                          # interpolate time-dependent parameters:
                          # --------------------------------------
            
                          param.val <- get.param.val(system@param,system@t.out[1])
            
                          # evaluate global environmental conditions:
                          # -----------------------------------------

                          cond.global <- list()
                          if ( length(system@cond) > 0 )
                          {
                             env <- c(param.val,list(t=system@t.out[1]))
                             for ( l in 1:length(system@cond) )
                             {
                                cond.l <- eval(expr  = system@cond[[l]],
                                               envir = env)
                                cond.global <- c(cond.global,cond.l)
                                n <- names(env)
                                env <- c(env,list(cond.l))
                                names(env) <- c(n,names(system@cond[[l]]))
                             }
                             names(cond.global) <- names(system@cond)
                          }
                          env.global = c(param.val,cond.global,
                                         list(t=system@t.out[1]))

                          # initialize masses:
                          # ------------------
                          
                          num.react  <- length(system@reactors)
                          offsets    <- numeric(num.react)
                          offsets[1] <- 0
                          x.ini <- numeric(0)
                          for ( k in 1:num.react )
                          {
                             # evaluate local environmental conditions:
                             # ----------------------------------------

                             reactor  <- system@reactors[[k]]
                             cond.local <- list()
                             if ( length(reactor@cond) > 0 )
                             {
                                env <- env.global
                                for ( l in 1:length(reactor@cond) )
                                {
                                   cond.l <- 
                                      eval(expr  = reactor@cond[[l]],
                                           envir = env)
                                   cond.local <- c(cond.local,cond.l)
                                   n <- names(env)
                                   env <- c(env,list(cond.l))
                                   names(env) <- c(n,names(reactor@cond[[l]]))
                                }
                                names(cond.local) <- names(reactor@cond)
                             }
                             env = c(env.global,cond.local)

                             # calculate masses:
                             # -----------------

                             volume   <- eval(expr  = reactor@volume.ini,
                                              envir = env)
                             area     <- eval(expr  = reactor@area,
                                              envir = env)
                             # store evaluated area:
                             system@reactors[[k]]@a <- area
                             num.vol  <- length(reactor@conc.pervol.ini)
                             num.area <- length(reactor@conc.perarea.ini)
                             conc.pervol <- numeric(num.vol)
                             if ( num.vol > 0 )
                             {
                                names(conc.pervol) <- 
                                   names(reactor@conc.pervol.ini)
                                for ( i in 1:num.vol )
                                {
                                   conc.pervol[i] <- 
                                      eval(expr  = reactor@conc.pervol.ini[[i]],
                                           envir = env)
                                }
                             }
                             conc.perarea <- numeric(num.area)
                             if ( num.area > 0 )
                             {
                                names(conc.perarea) <- 
                                   names(reactor@conc.perarea.ini)
                                for ( i in 1:num.area )
                                {
                                   conc.perarea[i] <- 
                                      eval(expr  = reactor@conc.perarea.ini[[i]],
                                           envir = env)
                                }
                             }
                             x.pervol  <- conc.pervol*volume
                             x.perarea <- conc.perarea*area
                             x.current <- c(volume=volume,x.pervol,x.perarea)
      
                             # combine masses with masses from previous react.:
                             # ------------------------------------------------
                             
                             x.ini <- c(x.ini,x.current)
      
                             # calculate new index offset:
                             # ---------------------------
      
                             if ( k < num.react ) 
                             { 
                                offsets[k+1] <- offsets[k] + 1+num.vol+num.area 
                             }
                          }

                          # calculate mass time series by numerical integration:
                          # ----------------------------------------------------

                          x <- ode(y      = x.ini,
                                   times  = system@t.out,
                                   func   = system.rhs,
                                   parms  = NA,
                                   method = method,
                                   system = system,
                                   ...)
                          rownames(x) <- x[,1]
                          x <- x[,-1]
                          
                          # convert masses back to concentrations:
                          # --------------------------------------
                          
                          for ( k in 1:num.react )
                          {
                             reactor <- system@reactors[[k]]
                             num.vol  <- length(reactor@conc.pervol.ini)
                             num.area <- length(reactor@conc.perarea.ini)
                             if ( num.vol > 0 )
                             {
                                for ( i in 1:num.vol )
                                {
                                   x[,offsets[k]+1+i] <- x[,offsets[k]+1+i]/
                                                         x[,offsets[k]+1]
                                }
                             }
                             if ( num.area > 0 )
                             {
                                for ( i in 1:num.area )
                                {
                                   x[,offsets[k]+1+num.vol+i] <- 
                                      x[,offsets[k]+1+num.vol+i]/reactor@a
                                }
                             }
                             if ( num.react > 1 )
                             {
                                names <- colnames(x)
                                names[offsets[k]+1:(1+num.vol+num.area)] <-
                                 paste(names[offsets[k]+1:(1+num.vol+num.area)],
                                       ".",reactor@name,sep="")
                                colnames(x) <- names
                             }
                          }
                                                    
                          # return concentration time series:
                          # ---------------------------------
                          
                          return(x)
                       })


# ==============================================================================
# Draw from Normal distribution
# ==============================================================================

# First version: Peter Reichert, Jan. 25, 2014
# Last revision: Peter Reichert, Jan. 25, 2014

# This function draws from a Normal or Lognormal random variable with parameters
# valid in original units.

# Arguments: 
#   mean:     Mean of the random variable.
#   sd:       Standard deviation of the random variable.
#   log:      Indicator whether the log of the variable should be
#             normally distributed (log=TRUE) rather than the
#             variable itself.
#             (Note: mean and sd are interpreted in original units
#             also for log=TRUE.)
#   n:        Sample size.


randnorm <- function(mean=0,sd=1,log=FALSE,n=1)
{
  # consistency checks:
  
  if ( n < 1 )
  {
    warning("n must a positive integer")
    return(NA)
  }
  if ( sd < 0 )
  {
    warning("sd must be non-negative")
    return(rep(NA,n))
  }
  
  # draw and return sample:
  
  if ( !log )
  {
    return(rnorm(n=n,mean=mean,sd=sd))
  }
  else
  {
    if ( mean <= 0 )
    {
      warning("if log=TRUE, mean must be positive")
      return(rep(NA,n))
    }
    meanlog <- log(mean/(sqrt(1+sd*sd/(mean*mean))))
    sdlog   <- sqrt(log(1+sd*sd/(mean*mean)))
    return(exp(rnorm(n=n,mean=meanlog,sd=sdlog)))
  }
}


# ==============================================================================
# Draw from Ornstein-Uhlenbeck process
# ==============================================================================

# First version: Peter Reichert, Jan. 25, 2014
# Last revision: Peter Reichert, Jan. 25, 2014

# This function draws a realization of an Ornstein-Uhlenbeck process.

# Arguments: 
#   mean:     Asymptotic mean of the process.
#   sd:       Asymptotic standard deviation of the process.
#   tau:      Correlation time of the process.
#   y0:       Starting value of the process.
#   t:        Time points at which the process should be sampled.
#             (Note: the value at t[1] will be the starting value y0.)
#   log:      Indicator whether the log of the variable should be an
#             Ornstein-Uhlenbeck process (log=TRUE) rather than the
#             variable itself.
#             (Note: mean and sd are interpreted in original units
#             also for log=TRUE.)

randou <- function(mean=0,sd=1,tau=0.1,y0=NA,t=0:1000/1000,log=FALSE)
{
  # consistency checks:
  
  if ( min(diff(t)) <= 0 )
  {
     warning("t must be strictly increasing")
     return(NA)
  }
  if ( sd < 0 )
  {
    warning("sd must be non-negative")
    return(NA)
  }
  if ( tau <= 0 )
  {
    warning("tau must be positive")
    return(NA)
  }
    
  # calculate exp of log if log=TRUE (note: mean and sd are in original units):
  
  if ( log ) 
  {
    if ( mean <= 0 )
    {
      warning("if log=TRUE, mean must be positive")
      return(NA)
    }
    meanlog <- log(mean/(sqrt(1+sd*sd/(mean*mean))))
    sdlog   <- sqrt(log(1+sd*sd/(mean*mean)))
    res <- randou(mean = meanlog,
                  sd   = sdlog,
                  tau  = tau,
                  y0   = log(y0),
                  t    = t,
                  log  = FALSE)
    res$y <- exp(res$y)
    return(res)
  }
  else
  {  
    # draw and return sample:
  
    y <- rep(NA,length(t))
    y[1] <- y0
    if ( is.na(y0) ) y[1] <- rnorm(n=1,mean=mean,sd=sd)
    for ( i in 2:length(t) )
    {
      y[i] <- rnorm(n    = 1,
                    mean = mean+(y[i-1]-mean)*exp(-(t[i]-t[i-1])/tau),
                    sd   = sd*sqrt(1-exp(-2*(t[i]-t[i-1])/tau)))
    }
    return(list(x=t,y=y))
  }
}


# ==============================================================================
# Plot simulation results
# ==============================================================================

# First version: Peter Reichert, April 09, 2006
# Last revision: Peter Reichert, Feb.  01, 2014


# Plot all columns of a matrix using the row labels as the common x-axis
# ======================================================================

plotres <- function(res,colnames=list(),file=NA,...)
{
   r <- res
   if ( !is.list(r) ) r <- list(r)
   cols <- colnames
   if ( !is.list(cols) )   cols <- list(cols)
   if ( length(cols) < 1 ) cols <- as.list(colnames(r[[1]]))
   num.col <- as.integer(sqrt(length(cols))+0.9999)
   num.row <- as.integer(length(cols)/num.col+0.9999)
   if ( !is.na(file) ) { pdf(file=file,...) }
   par.def <- par(no.readonly=TRUE)
   par(mfrow=c(num.row,num.col),
       xaxs="i",yaxs="i",
       mar=c(5,4.5,2,2))  # bottom,left,top,right
   t <- as.numeric(row.names(r[[1]]))
   for ( i in 1:length(cols) )
   {
      ymax <- 0
      for ( j in 1:length(cols[[i]]) )
      {
         ymax <- max(ymax,r[[1]][is.finite(r[[1]][,cols[[i]][j]]),cols[[i]][j]])
      }
      plot(numeric(0),numeric(0),
           xlim=c(min(t,na.rm=TRUE),max(t,na.rm=TRUE)),
           ylim=c(0,1.4*ymax),
           main=paste(cols[[i]],collapse=", "),
           xlab="t",ylab=paste(cols[[i]],collapse=", "))
      for ( j in 1:length(cols[[i]]) )
      {
         for ( k in 1:length(r) )
         {
            lines(t,r[[k]][,cols[[i]][j]],lty=j,col=j)
         }
      }
      legend("topright",legend=cols[[i]],
             lty=1:length(cols[[i]]),col=1:length(cols[[i]]))
   }
   par(par.def)
   if ( !is.na(file) ) { dev.off() }
}


