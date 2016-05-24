########################################################################################################
###################################### ORGANISM CLASS ##################################################
########################################################################################################

#' Structure of the S4 class "Organism"
#' 
#' Structure of the S4 class \code{Organism} representing the organisms present in the environment.
#' @export Organism
#' @exportClass Organism
#' @import sybil
#' @importFrom stats na.omit
#'
#' @slot lbnd A numeric vector containing the lower bounds of the model structure.
#' @slot ubnd A numeric vector containing the upper bounds of the model structure.
#' @slot type A character vector containing the description of the organism.
#' @slot medium A character vector containing all exchange reactions of the organism.
#' @slot lpobj A sybil optimization object containing the linear programing problem.
#' @slot fbasol A list with the solutions of the flux balance analysis.
#' @slot lyse A boolean variable indicating if the organism should lyse after death.
#' @slot feat A list containing conditional features for the object (contains at the momement only biomass components for lysis).
#' @slot deathrate A numeric value giving the factor by which the growth should be reduced in every iteration (unit: fg)
#' @slot growthlimit A numeric value giving the growth limit at which the organism dies.
#' @slot growtype A character vector giving the functional type for growth (linear or exponential).
#' @slot kinetics A List containing Km and v_max values for each reactions.
#' @slot speed A integer vector representing the speed by which bacterium is moving (given by cell per iteration).
#' @slot cellarea A numeric value indicating the surface that one organism occupies (unit: mu cm^2)
#' @slot cellweight A numeric value giving the maximal dry weight of single organism (unit: fg)
setClass("Organism",
         representation(
           lbnd="numeric",
           ubnd="numeric",
           type="character",
           medium="character",
           lpobj="sysBiolAlg",
           fbasol="list",
           lyse="logical",
           feat="list",
           deathrate="numeric",
           growthlimit="numeric",
           growtype="character",
           kinetics="list",
           cellarea="numeric",
           cellweight="numeric",
           speed="numeric"
         ),
         prototype(
           deathrate = 0.21,
           growthlimit = 0.083,
           growtype = "exponential",
           kinetics = list(),
           cellarea = 4.42,
           cellweight = 1.172,
           speed = 2       
         )
)

########################################################################################################
###################################### CONSTRUCTOR #####################################################
########################################################################################################

#' Constructor of the S4 class \code{Organism}
#' 
#' The constructor to get a new object of class \code{Organism}
#' @export
#' @name Organism-constructor
#' 
#' @param model model
#' @param algo A single character string giving the name of the algorithm to use. See \link[sybil]{SYBIL_SETTINGS}
#' @param ex Identifier for exchange reactions
#' @param ex_comp ex_comp
#' @param csuffix csuffix
#' @param esuffix esuffix
#' @param feat A list containing conditional features for the object (contains at the momement only biomass components for lysis).
#' @param typename A string defining the name (set to model name in default case)
#' @param lyse A boolean variable indicating if the organism should lyse after death.
#' @param ... Arguments of \code{\link{Organism-class}}
#' @return Object of class Organism
Organism <- function(model, algo="fba", ex="EX_", ex_comp=NA, csuffix="\\[c\\]", esuffix="\\[e\\]", lyse=F, feat=list(), 
                     typename=NA, ...){
  #the constructor requires the model, after that it is not stored anymore  
  if(is.na(typename)) typename <- sybil::mod_desc(model)
  rxname = sybil::react_id(model)
  lpobject <- sybil::sysBiolAlg(model, algorithm=algo)
  fbasol <- sybil::optimizeProb(lpobject)
  names(fbasol$fluxes) = rxname
  lobnd = sybil::lowbnd(model)
  names(lobnd) = rxname
  upbnd = sybil::uppbnd(model)
  names(upbnd) = rxname
  if(is.na(ex)){
    medc <- sybil::react_id(sybil::findExchReact(model))
  }else{
    medc <- sybil::react_id(sybil::findExchReact(model))
    medc <- medc[grep(ex, medc)]
  }
  if(!is.na(ex_comp)){
    medc <- medc[grep(ex_comp, medc)]
  }
  if(lyse){
    stochmat <- as.matrix(sybil::S(model))
    colnames(stochmat) <- rxname
    rownames(stochmat) <- sybil::met_id(model)
    stoch <- stochmat[,which(model@obj_coef==1)] #find stochiometry of biomass components
    biomets <- stoch[-which(stoch==0)]
    exs <- sybil::findExchReact(model)
    extrans <- sybil::react_id(exs)
    names(extrans) <- sybil::met_id(exs)
    extrans <- extrans[which(extrans %in% medc)]
    names(biomets) <- gsub(csuffix,esuffix,names(biomets))
    biomets <- biomets[which(names(biomets) %in% names(extrans))]
    names(biomets) <- extrans[names(biomets)]
    feat[["biomass"]] <- biomets
  }
  new("Organism", lbnd=lobnd, ubnd=upbnd, type=typename, medium=medc, lpobj=lpobject,
      fbasol=fbasol, feat=feat, lyse=lyse, ...)
}

########################################################################################################
###################################### GET METHODS FOR ATTRIBUTES ######################################
########################################################################################################

setGeneric("lbnd", function(object){standardGeneric("lbnd")})
setMethod("lbnd", "Organism", function(object){return(object@lbnd)})
setGeneric("ubnd", function(object){standardGeneric("ubnd")})
setMethod("ubnd", "Organism", function(object){return(object@ubnd)})
setGeneric("type", function(object){standardGeneric("type")})
setMethod("type", "Organism", function(object){return(object@type)})
setGeneric("medium", function(object){standardGeneric("medium")})
setMethod("medium", "Organism", function(object){return(object@medium)})
setGeneric("lpobj", function(object){standardGeneric("lpobj")})
setMethod("lpobj", "Organism", function(object){return(object@lpobj)})
setGeneric("fbasol", function(object){standardGeneric("fbasol")})
setMethod("fbasol", "Organism", function(object){return(object@fbasol)})
setGeneric("lyse", function(object){standardGeneric("lyse")})
setMethod("lyse", "Organism", function(object){return(object@lyse)})
setGeneric("feat", function(object){standardGeneric("feat")})
setMethod("feat", "Organism", function(object){return(object@feat)})
setGeneric("deathrate", function(object){standardGeneric("deathrate")})
setMethod("deathrate", "Organism", function(object){return(object@deathrate)})
setGeneric("growthlimit", function(object){standardGeneric("growthlimit")})
setMethod("growthlimit", "Organism", function(object){return(object@growthlimit)})
setGeneric("growtype", function(object){standardGeneric("growtype")})
setMethod("growtype", "Organism", function(object){return(object@feat)})
setGeneric("kinetics", function(object){standardGeneric("kinetics")})
setMethod("kinetics", "Organism", function(object){return(object@kinetics)})
setGeneric("speed", function(object){standardGeneric("speed")})
setMethod("speed", "Organism", function(object){return(object@speed)})

########################################################################################################
###################################### METHODS #########################################################
########################################################################################################

#' @title Function for constraining the models based on metabolite concentrations
#'
#' @description The generic function \code{constrain} changes the constraints of the model representation of an organism.
#' @export
#' @rdname constrain
#'
#' @param object An object of class Organisms.
#' @param reacts A character vector giving the names of reactions which should be constrained.
#' @param lb A numeric vector giving the constraint values of lower bounds (e.g. avaible metabolite concentrations
#' @param dryweight A number giving the current dryweight of the organism.
#' @param time A number giving the time intervals for each simulation step.
#' @param scale A numeric defining the scaling (units for linear programming has to be in certain range)
#' @return Returns the lower bounds, which carry the constraints and names of relevant reactions.
#' @details The constraints are calculated according to the flux definition as mmol/(gDW*hr) with the parameters \code{dryweight} and \code{time}.
#' @seealso \code{\link{Organism-class}}
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' org <- Organism(Ec_core,deathrate=0.05,
#'            growthlimit=0.05,growtype="exponential") #initialize an organism
#' lobnds <- constrain(org,org@medium,org@lbnd[org@medium],1,1)
setGeneric("constrain", function(object, reacts, lb, dryweight, time, scale){standardGeneric("constrain")})
#' @export
#' @rdname constrain
setMethod("constrain", "Organism", function(object, reacts, lb, dryweight, time, scale){
  lobnd <- object@lbnd*dryweight*time #costrain according to flux definition: mmol/(gDW*hr)
  #lobnd[reacts] <- ifelse(lb<=lobnd[reacts], ifelse(lobnd[reacts]==0, lb, lobnd[reacts]), lb) #check if lower bounds in biological relevant range
  lobnd[reacts] <- ifelse(lb<=lobnd[reacts], lobnd[reacts], lb) #check if lower bounds in biological relevant range
  if(length(object@kinetics) != 0){
    lobnd[names(object@kinetics)] <- unlist(lapply(names(object@kinetics), function(name){
      Km  <- (object@kinetics[[name]][["Km"]]*0.01*scale)*10^12 #scale mM to fmol/gridcell
      vmax <- object@kinetics[[name]][["vmax"]]*10^12*(dryweight*10^(-12)) #scale to fmol/h
      s   <- -lb[name] # change sign to get concentrations
      lnew = -(vmax*s/(Km + s))*time
      if(abs(lnew)>s){if(-s>=lobnd[name]){lnew=-s}else{lnew=lobnd[name]}}
      return(lnew)
    }))
  }
  return(lobnd)
})


setGeneric("setKinetics", function(object, exchangeR, Km, vmax){standardGeneric("setKinetics")})
setMethod("setKinetics", "Organism", function(object, exchangeR, Km, vmax){
  if( !(exchangeR %in% names(object@lbnd)) ){
    stop("Incorrect exchange reaction setup for kinetics")
  }
  nKinetics <- object@kinetics
  if(exchangeR %in% names(nKinetics)){
    warning("Kinetics for this exchange reaction already given, set to new value")
    nKinetics[[exchangeR]] <- list()
  }
  nKinetics[[exchangeR]] <- c(nKinetics[[exchangeR]], "Km"=Km, "vmax"=vmax)
  eval.parent(substitute(object@kinetics <- nKinetics))
})


#' @title Function for computing the linear programming according to the model structure 
#'
#' @description The generic function \code{optimizeLP} implements a linear programming based on the problem structure and refined constraints.
#' @export
#' @rdname optimizeLP
#'
#' @param object An object of class Organisms.
#' @param lpob A linear programing object encoding the problem to solve.
#' @param lb A numeric vector giving the constraint values of lower bounds.
#' @param ub A numeric vector giving the constraint values of upper bounds.
#' @details The problem object \code{lpob} is modified according to the constraints and then solved with \code{optimizeProb}.
#' @seealso \code{\link{Organism-class}}, \code{\link{optimizeProb}} and \code{\link{sysBiolAlg}}
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' org <- Organism(Ec_core,deathrate=0.05,
#'            growthlimit=0.05,growtype="exponential") #initialize a organism
#' optimizeLP(org)
setGeneric("optimizeLP", function(object, lpob=object@lpobj, lb=object@lbnd, ub=object@ubnd){standardGeneric("optimizeLP")})
#' @export
#' @rdname optimizeLP
setMethod("optimizeLP", "Organism", function(object, lpob=object@lpobj, lb=object@lbnd, ub=object@ubnd){ 
  fbasl <- sybil::optimizeProb(lpob, react=1:length(lb), ub=ub, lb=lb)
  names(fbasl$fluxes) <- names(object@lbnd)
  eval.parent(substitute(object@fbasol <- fbasl))
})

#' @title Function to account for the consumption and production of substances
#'
#' @description The generic function \code{consume} implements the consumption and production of substances based on the flux of exchange reactions of organisms
#' @export
#' @rdname consume
#'
#' @param object An object of class Organisms.
#' @param sublb A vector containing the substance concentrations in the current position of the individual of interest.
#' @param cutoff A number giving the cutoff value by which value of objective function is considered greater than 0.
#' @param bacnum integer indicating the number of bacteria individuals per gridcell
#' @return Returns the updated vector containing the substance concentrations in the current position of the individual of interest.
#' @details The consumption is implemented by adding the flux of the exchange reactions to the current substance concentrations.
#' @seealso \code{\link{Organism-class}}
#' @examples
#' NULL
setGeneric("consume", function(object, sublb, cutoff=1e-6, bacnum){standardGeneric("consume")})
#' @export
#' @rdname consume
setMethod("consume", "Organism", function(object, sublb, cutoff=1e-6, bacnum){
  if(object@fbasol$obj>=cutoff && !is.na(object@fbasol$obj)){
    flux = object@fbasol$fluxes[object@medium]*bacnum #scale flux to whole population size
    flux = na.omit(ifelse(abs(flux)<=cutoff,NA,flux))
    sublb[names(flux)] = round(sublb[names(flux)]+flux, 6)
  }
  return(sublb)
})



#' @title Function to extract the phenotype of an organism object
#'
#' @description The generic function \code{getPhenotype} implements an identification of organism phenotypes.
#' @export
#' @rdname getPhenotype
#'
#' @param object An object of class Organisms.
#' @param cutoff A number giving the cutoff value by which value of objective function is considered greater than 0.
#' @return Returns the phenotype of the organisms where the uptake of substances is indicated by a negative and production of substances by a positive number
#' @details The phenotypes are defined by flux through exchange reactions, which indicate potential differential substrate usages. Uptake of substances is indicated by a negative and production of substances by a positive number.
#' @seealso \code{\link{Organism-class}}, \code{\link{checkPhen}} and \code{\link{minePheno}}
#' @examples
#' \dontrun{
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' org <- Organism(Ec_core,deathrate=0.05,
#'            growthlimit=0.05,growtype="exponential") #initialize a organism
#' getPhenotype(org)
#' }
setGeneric("getPhenotype", function(object, cutoff=1e-6){standardGeneric("getPhenotype")})
#' @export
#' @rdname getPhenotype
setMethod("getPhenotype", "Organism", function(object, cutoff=1e-6){
  exflux=object@fbasol$fluxes[object@medium]
  exflux=ifelse(abs(exflux)<cutoff,0,1)*exflux
  exflux=ifelse(exflux>0,1,exflux)
  exflux=ifelse(exflux<0,2,exflux)
  return(exflux[which(exflux!=0)])
})

#' @title Function for letting organisms grow linearly
#'
#' @description The generic function \code{growLin} implements a growth model of organisms in their environment.
#' @export
#' @rdname growLin
#'
#' @param object An object of class Organisms.
#' @param growth A number indicating the current biomass, which has to be updated. 
#' @return Returns the updated biomass of the organisms of interest.
#' @details Linear growth of organisms is implemented by adding the calculated growthrate by \code{optimizeLP} to the already present growth value.
#' @seealso \code{\link{Organism-class}} and \code{\link{optimizeLP}}
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' org <- Organism(Ec_core,deathrate=0.05,
#'            growthlimit=0.05,growtype="exponential") #initialize a organism
#' growLin(org,1)
setGeneric("growLin", function(object, growth){standardGeneric("growLin")})
#' @export
#' @rdname growLin
setMethod("growLin", "Organism", function(object, growth){
  if(object@fbasol$obj > 0) grow_accum <- object@fbasol$obj + growth
  else grow_accum <- growth - object@deathrate
  return(grow_accum)
})

#' @title Function for letting organisms grow exponentially
#'
#' @description The generic function \code{growExp} implements a growth model of organisms in their environment.
#' @export
#' @rdname growExp
#'
#' @param object An object of class Organisms.
#' @param growth A number indicating the current biomass, which has to be updated. 
#' @return Returns the updated biomass of the organisms of interest.
#' @details Exponential growth of organisms is implemented by adding the calculated growthrate multiplied with the current growth calculated by \code{optimizeLP} plus to the already present growth value
#' @seealso \code{\link{Organism-class}} and \code{\link{optimizeLP}}
#' @examples
#' \dontrun{
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' org <- Organism(Ec_core,deathrate=0.05,
#'            growthlimit=0.05,growtype="exponential") #initialize a organism
#' growExp(org,1)
#' }
setGeneric("growExp", function(object, growth){standardGeneric("growExp")})
#' @export
#' @rdname growExp
setMethod("growExp", "Organism", function(object, growth){
  if(object@fbasol$obj > 0) grow_accum <- (object@fbasol$obj * growth + growth)
  else grow_accum <- growth - object@deathrate
  return(grow_accum)
  #return(object@fbasol$obj * growth + growth - object@deathrate)
})

#' @title Lysis function of organismal cells by adding biomass_compounds to the medium
#'
#' @description The generic function \code{lysis} implements cell lysis by the stochiometric concentration of the biomass compounds of organisms to the concentration of substances in the environment
#' @export
#' @rdname lysis
#'
#' @param object An object of class Organisms.
#' @param sublb A vector containing the substance concentrations in the current position of the individual of interest.
#' @param factor A number given the factor with which the biomass compound concentrations are multiplied to achieve the final concentration which is added to the environment
#' @return Returns the updated vector containing the substance concentrations in the current position of the dead individual of interest.
#' @details Lysis is implemented by taking the intersect between biomass compounds and the substances in the environment and adding the normalized stochiometric concentrations of the biomass compounds to the medium.
#' @seealso \code{\link{Organism-class}} and \code{\link{optimizeLP}}
#' @examples
#' NULL
setGeneric("lysis", function(object, sublb, factor=object@growthlimit){standardGeneric("lysis")})
#' @export
#' @rdname lysis
setMethod("lysis", "Organism", function(object, sublb, factor=object@growthlimit){
  stoch = object@feat[["biomass"]]
  lysate = round(abs(stoch)*factor, 6)
  sublb[names(lysate)] = sublb[names(lysate)] + lysate
  return(sublb)
})

#' @title Function to check if the there is a free place in the Moore neighbourhood
#'
#' @description The generic function \code{emptyHood} gives a free space which is present in the Moore neighbourhood of an individual of interest.
#' @export
#' @rdname emptyHood
#'
#' @param object An object of class Organisms.
#' @param x A number giving the x position of the individual of interest in its environment.
#' @param y A number giving the y position of the individual of interest in its environment.
#' @param n A number giving the horizontal size of the environment.
#' @param m A number giving the vertical size of the environment.
#' @param pos A dataframe with all occupied x and y positions 
#' @return Returns the free position in the Moore neighbourhood, which is not occupied by other individuals. If there is no free space \code{NULL} is returned.
#' @seealso \code{\link{Organism-class}}
#' @examples
#' NULL
setGeneric("emptyHood", function(object, pos, n, m, x, y){standardGeneric("emptyHood")})
#' @export
#' @rdname emptyHood
setMethod("emptyHood", "Organism", function(object, pos, n, m, x, y){
  xp = c(x-1,x,x+1)
  yp = c(y-1,y,y+1)
  xp=ifelse(xp<=0,NA,xp)
  xp=na.omit(ifelse(xp>n,NA,xp))
  yp=ifelse(yp<=0,NA,yp)
  yp=na.omit(ifelse(yp>m,NA,yp))
  #xp = xp[xp>0 & xp<=n]
  #xp = xp[yp>0 & yp<=m]
  nb=sapply(xp,function(x,y){return(paste(x,y,sep='_'))},y=yp)
  pos = pos[which(pos$x %in% xp),]
  pos = pos[which(pos$y %in% yp),]
  freenb=setdiff(nb,paste(pos$x,pos$y,sep='_'))
  if(length(freenb)==0){return(NULL)}else{return(freenb)}
})

#' @title Function to check if the there is a free place in the Moore neighbourhood
#'
#' @description The generic function \code{NemptyHood} gives a free space which is present in the Moore neighbourhood of an individual of interest.
#' @export
#' @rdname NemptyHood
#'
#' @param object An object of class Organisms.
#' @param x A number giving the x position of the individual of interest in its environment.
#' @param y A number giving the y position of the individual of interest in its environment.
#' @param n A number giving the horizontal size of the environment.
#' @param m A number giving the vertical size of the environment.
#' @param pos A dataframe with all occupied x and y positions 
#' @return Returns the free position in the Moore neighbourhood, which is not occupied by other individuals. If there is no free space \code{NULL} is returned.
#' @seealso \code{\link{Organism-class}}
#' @examples
#' NULL
setGeneric("NemptyHood", function(object, pos, n, m, x, y){standardGeneric("NemptyHood")})
#' @export
#' @rdname NemptyHood
setMethod("NemptyHood", "Organism", function(object, pos, n, m, x, y){
  xp = c(x-1,x,x+1)
  yp = c(y-1,y,y+1)
  for(i in 2:object@speed){
    xp = c(xp,x-i,x+i)
    yp = c(yp,y-i,y+i)
  }
  xp=ifelse(xp<=0,NA,xp)
  xp=na.omit(ifelse(xp>n,NA,xp))
  yp=ifelse(yp<=0,NA,yp)
  yp=na.omit(ifelse(yp>m,NA,yp))
  nb=sapply(xp,function(x,y){return(paste(x,y,sep='_'))},y=yp)
  pos = pos[which(pos$x %in% xp),]
  pos = pos[which(pos$y %in% yp),]
  freenb=setdiff(nb,paste(pos$x,pos$y,sep='_'))
  if(length(freenb)==0){return(NULL)}else{return(freenb)}
})

#' @title Function for random movement of organisms
#'
#' @description The generic function \code{move} implements a random movement in the Moore neighbourhood of an individual.
#' @export
#' @rdname move
#'
#' @param object An object of class Organism.
#' @param j The number of the iteration of interest.
#' @param n A number giving the horizontal size of the environment.
#' @param m A number giving the vertical size of the environment.
#' @param pos A dataframe with all occupied x and y positions 
#' @details Organisms move in a random position the Moore neighbourhood, which is not occupied by other individuals. If there is no free space the individuals stays in the same position.
#' @seealso \code{\link{Organism-class}}, \code{\link{emptyHood}}
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            growthlimit=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' addOrg(arena,bac,amount=10) #add 10 organisms
#' addSubs(arena,40) #add all possible substances
#' move(bac,n=20,m=20,j=1,pos=arena@orgdat[,c('x','y')])
setGeneric("move", function(object, pos, n, m, j){standardGeneric("move")})
#' @export
#' @rdname move
setMethod("move", "Organism", function(object, pos, n, m, j){
  if(object@speed == 1){
    freenb <- emptyHood(object, pos, n, m, pos[j,1], pos[j,2])
  }else{
    freenb <- NemptyHood(object, pos, n, m, pos[j,1], pos[j,2])
  }
  if(length(freenb) != 0){
    npos = freenb[sample(length(freenb),1)]
    npos = as.numeric(unlist(strsplit(npos,'_')))
    pos[j,] = npos
  }
  return(pos)
})

#show function for class Organism

setMethod(show, signature(object="Organism"), function(object){
  print(paste('Organism ',object@type,' of class Organism.',sep=''))
})

# Bac is a subclass of Organism containing bacteria specific features

########################################################################################################
###################################### BAC CLASS #######################################################
########################################################################################################

#' Structure of the S4 class "Bac"
#' 
#' Structure of the S4 class \code{Bac} inheriting from class \code{\link{Organism-class}} representing bacterial cells.
#' @export Bac
#' @exportClass Bac
#' @rdname Bac
#'
#' @slot chem A character vector indicating name of substance which is the chemotaxis attractant. Empty character vector if no chemotaxis.
setClass("Bac",
         contains="Organism",
         representation(
           chem="character" # name of substance which is the chemotaxis attractant
         ),
         prototype(
           deathrate = 0.21,
           growthlimit = 0.083,
           cellarea = 4.42,
           cellweight = 1.172
         )
)

########################################################################################################
###################################### CONSTRUCTOR #####################################################
########################################################################################################

#' Constructor of the S4 class \code{\link{Bac-class}}
#'
#' @name Bac-Constructor
#' @export
#' 
#' @return Object of class \code{\link{Bac-class}}
#' @param model model
#' @param ... Arguments of \code{\link{Organism-class}}
#' @param chem A character vector indicating name of substance which is the chemotaxis attractant. Empty character vector if no chemotaxis.
Bac <- function(model, chem='', ...){
  new("Bac", Organism(model=model, ...), chem=chem)
}

########################################################################################################
###################################### GET METHODS FOR ATTRIBUTES ######################################
########################################################################################################

setGeneric("chem", function(object){standardGeneric("chem")})
setMethod("chem", "Bac", function(object){return(object@chem)})

########################################################################################################
###################################### METHODS #########################################################
########################################################################################################

#' @title Function implementing a growth model of a bacterium
#'
#' @description The generic function \code{growth} implements different growth models for an object of class Bac.
#' @export
#' @rdname growth
#'
#' @param object An object of class Bac.
#' @param population An object of class Arena.
#' @param j The number of the iteration of interest.
#' @return Boolean variable of the \code{j}th individual indicating if individual died.
#' @details Linear growth of organisms is implemented by adding the calculated growthrate by \code{optimizeLP} to the already present growth value. Exponential growth of organisms is implemented by adding the calculated growthrate multiplied with the current growth calculated by \code{optimizeLP} plus to the already present growth value
#' @seealso \code{\link{Bac-class}}, \code{\link{growLin}} and \code{\link{growExp}}
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            growthlimit=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' addOrg(arena,bac,amount=10) #add 10 organisms
#' addSubs(arena,40) #add all possible substances
#' growth(bac,arena,1)
setGeneric("growth", function(object, population, j){standardGeneric("growth")})
#' @export
#' @rdname growth
setMethod("growth", "Bac", function(object, population, j){
  neworgdat <- population@orgdat
  popvec <- neworgdat[j,]
  switch(object@growtype,
         "linear"= {popvec$growth <- growLin(object, popvec$growth)},
         "exponential"= {popvec$growth <- growExp(object, popvec$growth)},
         stop("Growth type must be either linear or exponential"))
  dead <- F
  neworgdat[j,'growth'] <- popvec$growth
  if(popvec$growth > object@cellweight){
    freenb <- emptyHood(object, population@orgdat[,c('x','y')],
              population@n, population@m, popvec$x, popvec$y)
    if(length(freenb) != 0){
      npos = freenb[sample(length(freenb),1)]
      npos = as.numeric(unlist(strsplit(npos,'_')))
      daughter <- popvec
      daughter$growth <- popvec$growth/2
      daughter$x <- npos[1]
      daughter$y <- npos[2]
      neworgdat[nrow(neworgdat)+1,] <- daughter
      neworgdat[j,'growth'] <- popvec$growth/2
    }
  }
  else if(popvec$growth < object@growthlimit){
    neworgdat[j,'growth'] <- NA
    dead <- T
  }
  eval.parent(substitute(population@orgdat <- neworgdat))
  return(dead)
})

# chemotaxis still needs to be edited

#' @title Function for chemotaxis of bacteria to their prefered substrate
#'
#' @description The generic function \code{chemotaxis} implements a bacterial movement in the Moore neighbourhood to the highest substrate concentration.
#' @export
#' @rdname chemotaxis
#'
#' @param object An object of class Bac.
#' @param population An object of class Arena.
#' @param j The number of the iteration of interest.
#' @details Bacteria move to a position in the Moore neighbourhood which has the highest concentration of the prefered substrate, which is not occupied by other individuals. The prefered substance is given by slot \code{chem} in the \code{Bac} object. If there is no free space the individuals stays in the same position. If the concentration in the Moore neighbourhood has the same concentration in every position, then random movement is implemented.
#' @seealso \code{\link{Bac-class}} and \code{\link{emptyHood}}
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05, chem = "EX_o2(e)",
#'            growthlimit=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' addOrg(arena,bac,amount=10) #add 10 organisms
#' addSubs(arena,40) #add all possible substances
#' chemotaxis(bac,arena,1)
setGeneric("chemotaxis", function(object, population, j){standardGeneric("chemotaxis")})
#' @export
#' @rdname chemotaxis
setMethod("chemotaxis", "Bac", function(object, population, j){
  popvec <- population@orgdat[j,]
  attract <- population@media[[object@chem]]@diffmat
  freenb <- emptyHood(object, population@orgdat[,c('x','y')],
                      population@n, population@m, popvec$x, popvec$y)
  if(length(freenb) != 0){
    conc <- sapply(freenb, function(x, attract){
      npos = as.numeric(unlist(strsplit(x,'_')))
      return(attract[npos[1],npos[2]])
    }, attract=attract)
    abs <- freenb[which(conc==max(conc))]
    if(length(abs)!=1){
      abs <- abs[sample(1:length(abs),1)]
    }
    npos = as.numeric(unlist(strsplit(abs,'_')))
    eval.parent(substitute(population@orgdat[j,]$x <- npos[1]))
    eval.parent(substitute(population@orgdat[j,]$y <- npos[2]))
  }
})

#function for one iteration for Bac class
#' @title Function for one simulation iteration for objects of Bac class
#'
#' @description The generic function \code{simBac} implements all neccessary functions for the individuals to update the complete environment. 
#' @export
#' @rdname simBac
#'
#' @param object An object of class Bac.
#' @param arena An object of class Arena defining the environment.
#' @param j The number of the iteration of interest.
#' @param bacnum integer indicating the number of bacteria individuals per gridcell
#' @param sublb A vector containing the substance concentrations in the current position of the individual of interest.
#' @return Returns the updated enivironment of the \code{population} parameter with all new positions of individuals on the grid and all new substrate concentrations.
#' @details Bacterial individuals undergo step by step the following procedures: First the individuals are constrained with \code{constrain} to the substrate environment, then flux balance analysis is computed with \code{optimizeLP}, after this the substrate concentrations are updated with \code{consume}, then the bacterial growth is implemented with \code{growth}, the potential new phenotypes are added with \code{checkPhen}, finally the additional and conditional functions \code{lysis}, \code{move} or \code{chemotaxis} are performed. Can be used as a wrapper for all important bacterial functions in a function similar to \code{simEnv}.
#' @seealso \code{\link{Bac-class}}, \code{\link{Arena-class}}, \code{\link{simEnv}}, \code{constrain}, \code{optimizeLP}, \code{consume}, \code{growth}, \code{checkPhen}, \code{lysis}, \code{move} and \code{chemotaxis}
#' @examples
#' NULL
setGeneric("simBac", function(object, arena, j, sublb, bacnum){standardGeneric("simBac")})
#' @export
#' @rdname simBac
setMethod("simBac", "Bac", function(object, arena, j, sublb, bacnum){
  lobnd <- constrain(object, object@medium, lb=-sublb[j,object@medium]/bacnum, #scale to population size
                     dryweight=arena@orgdat[j,"growth"], time=arena@tstep, scale=arena@scale)
  optimizeLP(object, lb=lobnd)
  
  eval.parent(substitute(sublb[j,] <- consume(object, sublb[j,], bacnum=bacnum))) #scale consumption to the number of cells?

  dead <- growth(object, arena, j)
  arena@orgdat[j,'phenotype'] <- as.integer(checkPhen(arena, object))
  
    type <- object@type
  arena@mflux[[type]] <- arena@mflux[[type]] + object@fbasol$fluxes # remember active fluxes

  if(dead && object@lyse){
    eval.parent(substitute(sublb[j,] <- lysis(object, sublb[j,])))
  }
  pos <- arena@orgdat[,c('x','y')]
  if(!dead && !arena@stir && object@speed != 0){
    if(object@chem == ''){
      pos <- move(object, pos, arena@n, arena@m, j)
    }else{
      chemotaxis(object, arena, j)
    }
  }
  arena@orgdat[,c('x','y')] <- pos
  return(arena)
})

#show function for class Bac

setMethod(show, signature(object="Bac"), function(object){
  print(paste('Bacterium ',object@type,' of class Bac.',sep=''))
})

# Human is a subclass of Organism containing human specific features

########################################################################################################
###################################### HUMAN CLASS #####################################################
########################################################################################################

#' Structure of the S4 class "Human"
#' 
#' Structure of the S4 class \code{Human} inheriting from class \code{\link{Organism-class}} representing human cells.
#' @export Human
#' @exportClass Human
#' @rdname Human
#'
#' @slot objective A character vector representing the current reaction which should be used as an objective function for the flux balance analysis.
setClass("Human",
         contains="Organism",
         representation(
           objective="character"
         )
)

########################################################################################################
###################################### CONSTRUCTOR #####################################################
########################################################################################################

#' Constructor of the S4 class \code{\link{Human-class}}
#' 
#' @name Human-constructor
#' @export
#' 
#' @param model model
#' @param objective A character vector representing the current reaction which should be used as an objective function for the flux balance analysis.
#' @param speed A integer vector representing the speed by which bacterium is moving (given by cell per iteration).
#' @param ... Arguments of \code{\link{Organism-class}}
#' @return Object of class \code{\link{Human-class}}
Human <- function(model, objective=model@react_id[which(model@obj_coef==1)], speed=0, ...){
  model <- sybil::changeObjFunc(model, objective)
  new("Human", Organism(model=model, speed=speed, ...), objective=objective)
}

########################################################################################################
###################################### GET METHODS FOR ATTRIBUTES ######################################
########################################################################################################

setGeneric("objective", function(object){standardGeneric("objective")})
setMethod("objective", "Human", function(object){return(object@objective)})

########################################################################################################
###################################### METHODS #########################################################
########################################################################################################

#' @title Function for changing the objective function of the model
#'
#' @description The generic function \code{changeFobj} changes the objective function, which is used for the linear programming in \code{optimizeLP}.
#' @export
#' @rdname changeFobj
#'
#' @param object An object of class Human.
#' @param new_fobj A character vector giving the reaction name of the new objective function.
#' @param model The original model structure which is converted into a problem object used for the next optimization.
#' @param alg A character vector giving the algorithm which should be used for the optimization (default is flux balance analysis).
#' @details To avoid the bias to just one particular objective function, the objective can be changed dynamically in this function. 
#' @seealso \code{\link{Human-class}} and \code{\link{optimizeLP}}
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' human <- Human(Ec_core,deathrate=0.05,
#'            growthlimit=0.05,growtype="exponential") #initialize a bacterium
#' changeFobj(human,'EX_glc(e)',Ec_core)
setGeneric("changeFobj", function(object, new_fobj, model, alg="fba"){standardGeneric("changeFobj")})
#' @export
#' @rdname changeFobj
setMethod("changeFobj", "Human", function(object, new_fobj, model, alg="fba"){
  eval.parent(substitute(object@objective <- new_fobj)) #(pseudo) call by reference implementation
  model <- sybil::changeObjFunc(model, new_fobj)
  eval.parent(substitute(object@lpobj <- sysBiolAlg(model, algorithm=alg))) #the lp object has to be updated according to the new objective
})

#' @title Function implementing a growth model of a human cell
#'
#' @description The generic function \code{cellgrowth} implements different growth models for an object of class Human.
#' @export
#' @rdname cellgrowth
#'
#' @param object An object of class Human.
#' @param population An object of class Arena.
#' @param j The number of the iteration of interest.
#' @return Boolean variable of the \code{j}th individual indicating if individual died.
#' @details Linear growth of organisms is implemented by adding the calculated growthrate by \code{optimizeLP} to the already present growth value. Exponential growth of organisms is implemented by adding the calculated growthrate multiplied with the current growth calculated by \code{optimizeLP} plus to the already present growth value.
#' @seealso \code{\link{Human-class}}, \code{\link{growLin}} and \code{\link{growExp}}
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' human <- Human(Ec_core,deathrate=0.05,
#'            growthlimit=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' addOrg(arena,human,amount=10) #add 10 organisms
#' addSubs(arena,40) #add all possible substances
#' cellgrowth(human,arena,1)
setGeneric("cellgrowth", function(object, population, j){standardGeneric("cellgrowth")})
#' @export
#' @rdname cellgrowth
setMethod("cellgrowth", "Human", function(object, population, j){
  neworgdat <- population@orgdat
  popvec <- neworgdat[j,]
  switch(object@growtype,
         "linear"= {popvec$growth <- growLin(object, popvec$growth)},
         "exponential"= {popvec$growth <- growExp(object, popvec$growth)},
         stop("Growth type must be either linear or exponential"))
  dead <- F
  neworgdat[j,'growth'] <- popvec$growth
  if(popvec$growth > object@cellweight){
    freenb <- emptyHood(object, population@orgdat[,c('x','y')],
                        population@n, population@m, popvec$x, popvec$y)
    if(length(freenb) != 0){
      npos = freenb[sample(length(freenb),1)]
      npos = as.numeric(unlist(strsplit(npos,'_')))
      daughter <- popvec
      daughter$growth <- popvec$growth/2
      daughter$x <- npos[1]
      daughter$y <- npos[2]
      neworgdat[nrow(neworgdat)+1,] <- daughter
      neworgdat[j,'growth'] <- popvec$growth/2
    }
  }
  else if(popvec$growth < object@growthlimit){
    neworgdat[j,'growth'] <- NA
    dead <- T
  }
  eval.parent(substitute(population@orgdat <- neworgdat))
  return(dead)
})

#' @title Function for one simulation iteration for objects of Human class
#'
#' @description The generic function \code{simHum} implements all neccessary functions for the individuals to update the complete environment. 
#' @export
#' @rdname simHum
#'
#' @param object An object of class Human.
#' @param j The number of the iteration of interest.
#' @param arena An object of class Arena defining the environment.
#' @param bacnum integer indicating the number of bacteria individuals per gridcell
#' @param sublb A vector containing the substance concentrations in the current position of the individual of interest.
#' @return Returns the updated enivironment of the \code{arena} parameter with all new positions of individuals on the grid and all new substrate concentrations.
#' @details Human cell individuals undergo the step by step the following procedures: First the individuals are constrained with \code{constrain} to the substrate environment, then flux balance analysis is computed with \code{optimizeLP}, after this the substrate concentrations are updated with \code{consume}, then the cell growth is implemented with \code{cellgrowth}, the potential new phenotypes are added with \code{checkPhen}, finally the conditional function \code{lysis} is performed. Can be used as a wrapper for all important cell functions in a function similar to \code{simEnv}.
#' @seealso \code{\link{Human-class}}, \code{\link{Arena-class}}, \code{\link{simEnv}}, \code{constrain}, \code{optimizeLP}, \code{consume}, \code{cellgrowth}, \code{checkPhen} and \code{lysis}
#' @examples
#' NULL
setGeneric("simHum", function(object, arena, j, sublb, bacnum){standardGeneric("simHum")})
#' @export
#' @rdname simHum
setMethod("simHum", "Human", function(object, arena, j, sublb, bacnum){
  lobnd <- constrain(object, object@medium, lb=-sublb[j,object@medium]/bacnum, #scale to population size
                     dryweight=arena@orgdat[j,"growth"], time=arena@tstep, scale=arena@scale)
  optimizeLP(object, lb=lobnd)
  eval.parent(substitute(sublb[j,] <- consume(object, sublb[j,], bacnum=bacnum))) #rescale from population size
  dead <- cellgrowth(object, arena, j)
  arena@orgdat[j,'phenotype'] <- as.integer(checkPhen(arena, object))
  if(dead && object@lyse){
    eval.parent(substitute(sublb[j,] <- lysis(object, names(arena@media), sublb[j,])))
  }
  return(arena)
})

#show function for class Human

setMethod(show, signature(object="Human"), function(object){
  print(paste('Cell culture ',object@type,' of class Human.',sep=''))
})
