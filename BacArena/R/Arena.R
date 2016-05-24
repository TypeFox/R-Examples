globalVariables(c("diffuseNaiveCpp","diffuseSteveCpp"))

########################################################################################################
###################################### Arena CLASS ################################################
########################################################################################################

#' Structure of the S4 class "Arena"
#' 
#' Structure of the S4 class \code{Arena} to represent the environment in which Organisms and Substances interact.
#' @import Matrix Rcpp
#' @export Arena
#' @exportClass Arena
#'
#' @slot orgdat A data frame collecting information about the accumulated growth, type, phenotype, x and y position for each individual in the environment.
#' @slot specs A list of organism types and their associated parameters.
#' @slot media A list of objects of class \code{\link{Substance-class}} for each compound in the environment.
#' @slot phenotypes A list of unique phenotypes (metabolites consumed and produced), which occurred in the environment.
#' @slot mediac A character vector containing the names of all substances in the environment.
#' @slot tstep A number giving the time (in h) per iteration.
#' @slot stir A boolean variable indicating if environment should be stirred.
#' @slot mflux A vector containing highly used metabolic reactions within the arena
#' @slot n A number giving the horizontal size of the environment.
#' @slot m A number giving the vertical size of the environment.
#' @slot Lx A number giving the horizontal grid size in cm.
#' @slot Ly A number giving the vertical grid size in cm.
#' @slot gridgeometry A list containing grid geometry parameter 
#' @slot seed An integer refering to the random number seed used to be reproducible
#' @slot scale A numeric defining the scale factor used for intern unit conversion.
setClass("Arena",
         representation(
           orgdat="data.frame",
           specs="list",
           media="list",
           phenotypes="character",
           mediac="character",
           tstep="numeric",
           stir="logical",
           mflux="list",
           n="numeric",
           m="numeric",
           gridgeometry="list",
           Lx="numeric",
           Ly="numeric",
           seed="numeric",
           scale="numeric"
        ),
        prototype(
          orgdat = data.frame(growth=numeric(0),type=integer(0),phenotype=integer(0),x=integer(0),y=integer(0)),
          specs = list(),
          media = list(),
          phenotypes = character(),
          mediac = character(),
          tstep = 1,
          stir = F,
          mflux = list(),
          seed=sample(1:10000,1)
        )
)

########################################################################################################
###################################### CONSTRUCTOR #####################################################
########################################################################################################

#' Constructor of the S4 class \code{\link{Arena-class}}
#' 
#' @export
#' @name Arena-constructor
#' @param n A number giving the horizontal size of the environment.
#' @param m A number giving the vertical size of the environment.
#' @param Lx A number giving the horizontal grid size in cm.
#' @param Ly A number giving the vertical grid size in cm.
#' @param ... Arguments of \code{\link{Arena-class}}
Arena <- function(Lx=0.05, Ly=0.05, n=100, m=100, ...){
  gridgeometry = list(grid2D=ReacTran::setup.grid.2D(ReacTran::setup.grid.1D(x.up = 0, L = Lx, N = n), 
                                                     ReacTran::setup.grid.1D(x.up = 0, L = Ly, N = m)))
  scale=(Lx*Ly)/(n*m)
  new("Arena", Lx=Lx, Ly=Ly, n=n, m=m, scale=scale, gridgeometry=gridgeometry, ...)
}


########################################################################################################
###################################### GET METHODS FOR ATTRIBUTES ######################################
########################################################################################################

setGeneric("orgdat", function(object){standardGeneric("orgdat")})
setMethod("orgdat", "Arena", function(object){return(object@orgdat)})
setGeneric("specs", function(object){standardGeneric("specs")})
setMethod("specs", "Arena", function(object){return(object@specs)})
setGeneric("media", function(object){standardGeneric("media")})
setMethod("media", "Arena", function(object){return(object@media)})
setGeneric("phenotypes", function(object){standardGeneric("phenotypes")})
setMethod("phenotypes", "Arena", function(object){return(object@phenotypes)})
setGeneric("mediac", function(object){standardGeneric("mediac")})
setMethod("mediac", "Arena", function(object){return(object@mediac)})
setGeneric("tstep", function(object){standardGeneric("tstep")})
setMethod("tstep", "Arena", function(object){return(object@tstep)})
setGeneric("stir", function(object){standardGeneric("stir")})
setMethod("stir", "Arena", function(object){return(object@stir)})
setGeneric("mflux", function(object){standardGeneric("mflux")})
setMethod("mflux", "Arena", function(object){return(object@mflux)})
setGeneric("n", function(object){standardGeneric("n")})
setMethod("n", "Arena", function(object){return(object@n)})
setGeneric("m", function(object){standardGeneric("m")})
setMethod("m", "Arena", function(object){return(object@m)})
setGeneric("gridgeometry", function(object){standardGeneric("gridgeometry")})
setMethod("gridgeometry", "Arena", function(object){return(object@gridgeometry)})
setGeneric("scale", function(object){standardGeneric("scale")})
setMethod("scale", "Arena", function(object){return(object@scale)})
setGeneric("seed", function(object){standardGeneric("seed")})
setMethod("seed", "Arena", function(object){return(object@seed)})

########################################################################################################
###################################### METHODS #########################################################
########################################################################################################


#' @title Add individuals to the environment
#'
#' @description The generic function \code{addOrg} adds individuals to the environment.
#' @export
#' @rdname addOrg
#'
#' @param object An object of class Arena.
#' @param specI An object of class Organism.
#' @param amount A numeric number giving the number of individuals to add.
#' @param x A numeric vector giving the x positions of individuals on the grid.
#' @param y A numeric vector giving the y positions of individuals on the grid.
#' @param growth A numeric vector giving the starting biomass of the individuals.
#' @param mean A numeric giving the mean of starting biomass (used for normal distribution) if growth is not defined
#' @param sd A numeric giving the standard derivation of starting biomass (used for normal distribution) if growth is not defined
#' @details The arguments \code{x} and \code{y} should be in the same length as the number of organisms added (given by the argument \code{amount}).
#' @seealso \code{\link{Arena-class}} and \code{\link{Bac-class}} 
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            growthlimit=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' addOrg(arena,bac,amount=10) #add 10 organisms
setGeneric("addOrg", function(object, specI, amount, x=NULL, y=NULL, growth=NA, mean=0.489, sd=0.132){standardGeneric("addOrg")})
#' @export
#' @rdname addOrg
setMethod("addOrg", "Arena", function(object, specI, amount, x=NULL, y=NULL, growth=NA, mean=0.489, sd=0.132){
  if(amount+nrow(object@orgdat) > object@n*object@m){
    stop("More individuals than space on the grid")
  }
  n <- object@n
  m <- object@m
  spectype <- specI@type
  neworgdat <- object@orgdat
  newspecs <- object@specs
  newspecs[[spectype]] <- specI
  type <- which(names(newspecs)==spectype)
  newmflux <- object@mflux
  
  # mflux
  newmflux[[spectype]] <- numeric(length(specI@lbnd))
  names(newmflux[[spectype]]) <- names(specI@lbnd)
  
  type <- which(names(newspecs)==spectype) 
  lastind <- nrow(object@orgdat)
  if(length(x*y)==0){
    cmbs = expand.grid(1:n,1:m)
    rownames(cmbs) = paste(cmbs[,1],cmbs[,2],sep='_')
    taken <- paste(object@orgdat$x,object@orgdat$y,sep='_')
    if(length(taken)!=0){
      cmbs <- cmbs[-which(rownames(cmbs) %in% taken),]
    }
    sel <- sample(1:nrow(cmbs),amount)
    xp = cmbs[sel,1]
    yp = cmbs[sel,2]
    neworgdat[(lastind+1):(amount+lastind),'x']=xp
    neworgdat[(lastind+1):(amount+lastind),'y']=yp
    if(is.numeric(growth)) neworgdat[(lastind+1):(amount+lastind),'growth'] = rep(growth, amount)
    else neworgdat[(lastind+1):(amount+lastind),'growth'] = abs(rnorm(amount, mean=mean, sd=sd))
    neworgdat[(lastind+1):(amount+lastind),'type']=rep(type, amount)
    neworgdat[(lastind+1):(amount+lastind),'phenotype']=rep(NA, amount)
  }else{
    neworgdat[(lastind+1):(amount+lastind),'x']=x
    neworgdat[(lastind+1):(amount+lastind),'y']=y
    if(is.numeric(growth)) neworgdat[(lastind+1):(amount+lastind),'growth'] = rep(growth, amount)
    else neworgdat[(lastind+1):(amount+lastind),'growth'] = abs(rnorm(amount, mean=mean, sd=sd))
    neworgdat[(lastind+1):(amount+lastind),'type']=rep(type, amount)
    neworgdat[(lastind+1):(amount+lastind),'phenotype']=rep(NA, amount)
  }
  if(sum(duplicated(paste(neworgdat$x,neworgdat$y,sep="_")))!=0){
    stop("You have multiple individuals in the same position! Make sure that your x an y positions are unique")
  }
  eval.parent(substitute(object@orgdat <- neworgdat))
  eval.parent(substitute(object@specs <- newspecs))
  #eval.parent(substitute(object@phenotypes[[spectype]] <- newphens))
  eval.parent(substitute(object@mediac <- union(object@mediac, specI@medium)))
  eval.parent(substitute(object@mflux <- newmflux))
})

#' @title Add substances to the environment
#'
#' @description The generic function \code{addSubs} adds specific substances to the environment.
#' @export
#' @rdname addSubs
#'
#' @param object An object of class Arena.
#' @param mediac A character vector giving the names of substances, which should be added to the environment (the default takes all possible substances).
#' @param smax A numeric vector indicating the maximum substance concentration per grid cell.
#' @param unit A character used as chemical unit to set the amount of the substances to be added (valid values are: mmol/cell, mmol/cm2, mmol/arena, mM)
#' @param difunc A character vector ("pde","cpp" or "r") describing the function for diffusion.
#' @param difspeed A number indicating the diffusion speed (given by number of cells per iteration).
#' @param add A boolean variable defining whether the amount of substance should be summed or replaced
#' @details If nothing but \code{object} is given, then all possible substrates are initilized with a concentration of 0. Afterwards, \code{\link{changeSub} can be used to modify the concentrations of specific substances.} 
#' @seealso \code{\link{Arena-class}} and \code{\link{changeSub}} 
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            growthlimit=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' addOrg(arena,bac,amount=10) #add 10 organisms
#' addSubs(arena,20,c("EX_glc(e)","EX_o2(e)","EX_pi(e)")) #add substances glucose, oxygen and phosphate
setGeneric("addSubs", function(object, smax=0, mediac=object@mediac, difunc="pde", difspeed=6.7e-6, unit="mmol/cell", add=T){standardGeneric("addSubs")})
#' @rdname addSubs
#' @export
setMethod("addSubs", "Arena", function(object, smax=0, mediac=object@mediac, difunc="pde", difspeed=6.7e-6, unit="mmol/cell", add=T){
  if(length(smax) != length(mediac) && length(smax) != 1){
    stop("The parameter smax should be of the same size of mediac or equal to 1.")
  }
  if(sum(mediac %in% object@mediac) != length(mediac)){stop("Substance does not exist in medium.")}
  if(length(intersect(unit,c("mmol/cell","mM","mmol/arena","mmol/cm2")))==0){stop("Wrong unit for concentration.")}
  if(length(smax) == 1){
    smax = rep(smax,length(mediac))
  }
  if(unit=="mM"){smax <- (smax*0.01)*object@scale}  # conversion of mMol in mmol/grid_cell
  if(unit=="mmol/cm2"){smax <- smax*object@scale}  # conversion of mmol/arena in mmol/grid_cell
  if(unit=="mmol/arena"){smax <- smax/(object@n*object@m)}  # conversion of mmol/arena in mmol/grid_cell
  if(length(difspeed)!=length(mediac)){difspeed = rep(difspeed,length(mediac))}
  if(length(object@media) == 0){
    newmedia <- list()
    for(i in 1:length(mediac)){
      newmedia[[mediac[i]]] <- Substance(object@n, object@m, 0, name=mediac[i], difunc=difunc, difspeed=difspeed[i], gridgeometry=object@gridgeometry)
    }
  }else{newmedia <- object@media}
  if(length(object@mediac) > length(newmedia)){
    appendlist <- list()
    for(i in setdiff(object@mediac, names(newmedia))){
      appendlist[[i]] <- Substance(object@n, object@m, 0, name=i, difunc=difunc, difspeed=6.7e-6, gridgeometry=object@gridgeometry)
    }
    newmedia = c(newmedia,appendlist)
  }
  if(add){
    for(i in 1:length(mediac)){
      newdmat = newmedia[[mediac[i]]]@diffmat + Matrix::Matrix(smax[i], nrow=object@n, ncol=object@m, sparse=TRUE)
      newmedia[[mediac[i]]]@diffmat <- newdmat
      #newmedia[[mediac[i]]]@difspeed = difspeed[i]
    }
  }else{
    for(i in 1:length(mediac)){
      newmedia[[mediac[i]]]@diffmat <- Matrix::Matrix(smax[i], nrow=object@n, ncol=object@m, sparse=TRUE)
      newmedia[[mediac[i]]]@difspeed = difspeed[i]
    }
  }
  eval.parent(substitute(object@media <- newmedia))
})

#' @title Change substances in the environment
#'
#' @description The generic function \code{changeSub} changes specific substances in the environment.
#' @export
#' @rdname changeSub
#'
#' @param object An object of class Arena.
#' @param smax A number indicating the maximum substance concentration per grid cell.
#' @param mediac A character vector giving the names of substances, which should be added to the environment (the default takes all possible substances).
#' @param unit A character used as chemical unit to set the amount of the substances to be added (valid values are: mmol/cell, mmol/cm2, mmol/arena, mM)
#' @details If nothing but \code{object} is given, then all possible substrates are initilized with a concentration of 0. Afterwards, \code{\link{changeSub}} can be used to modify the concentrations of specific substances.
#' @seealso \code{\link{Arena-class}} and \code{\link{addSubs}} 
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            growthlimit=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' addOrg(arena,bac,amount=10) #add 10 organisms
#' addSubs(arena) #add all substances with no concentrations.
#' changeSub(arena,20,c("EX_glc(e)","EX_o2(e)","EX_pi(e)")) 
#' #add substances glucose, oxygen and phosphate
setGeneric("changeSub", function(object, smax, mediac, unit="mmol/cell"){standardGeneric("changeSub")})
#' @rdname changeSub
#' @export
setMethod("changeSub", "Arena", function(object, smax, mediac, unit="mmol/cell"){
  if(sum(mediac %in% names(object@media))==length(mediac)){
    if(unit=="mM"){smax <- (smax*0.01)*object@scale}  # conversion of mMol in mmol/grid_cell
    if(unit=="mmol/cm2"){smax <- smax*object@scale}  # conversion of mmol/arena in mmol/grid_cell
    if(unit=="mmol/arena"){smax <- smax/(object@n*object@m)}  # conversion of mmol/arena in mmol/grid_cell
    for(i in 1:length(mediac)){
      eval.parent(substitute(object@media[mediac[i]] <- Substance(object@n, object@m, smax=smax, name=mediac[i],
                                                                  difunc=object@media[[mediac[i]]]@difunc,
                                                                  difspeed=object@media[[mediac[i]]]@difspeed, gridgeometry=object@gridgeometry)))
    }
  }else stop("Substance does not exist in medium.")
})

#' @title Remove all substances in the environment
#'
#' @description The generic function \code{flushSubs} removes specific substances in the environment.
#' @export
#' @rdname flushSubs
#'
#' @param object An object of class Arena.
#' @seealso \code{\link{Arena-class}} and \code{\link{addSubs}} 
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            growthlimit=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' addOrg(arena,bac,amount=10) #add 10 organisms
#' addSubs(arena, smax=40) #add all substances with no concentrations.
#' changeSub(arena,20,c("EX_glc(e)","EX_o2(e)","EX_pi(e)")) 
#' #add substances glucose, oxygen and phosphate
#' flushSubs(arena) #remove all created substance concentrations
setGeneric("flushSubs", function(object){standardGeneric("flushSubs")})
#' @export
#' @rdname flushSubs
setMethod("flushSubs", "Arena", function(object){
  eval.parent(substitute(object@media <- list()))
})

#' @title Change substance concentration patterns in the environment
#'
#' @description The generic function \code{changeDiff} changes specific substance concentration patterns in the environment.
#' @export
#' @rdname changeDiff
#'
#' @param object An object of class Arena.
#' @param newdiffmat A matrix giving the new gradient matrix of the specific substances in the environment.
#' @param mediac A character vector giving the names of substances, which should be added to the environment (the default takes all possible substances).
#' @details This function can be used to add gradients of specific substances in the environment. The default conditions in \code{changeSubs} assumes an equal concentration in every grid cell of the environment. 
#' @seealso \code{\link{Arena-class}} and \code{\link{changeSub}} 
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            growthlimit=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' addOrg(arena,bac,amount=10) #add 10 organisms
#' addSubs(arena,30) #add all substances with no concentrations.
#' gradient <- matrix(1:200,20,20)
#' changeDiff(arena,gradient,c("EX_glc(e)","EX_o2(e)","EX_pi(e)"))
#' # add substances glucose, oxygen and phosphate
setGeneric("changeDiff", function(object, newdiffmat, mediac){standardGeneric("changeDiff")})
#' @export
#' @rdname changeDiff
setMethod("changeDiff", "Arena", function(object, newdiffmat, mediac){
  if(nrow(newdiffmat)==object@n && ncol(newdiffmat)==object@m){
    for(i in 1:length(mediac)){
      eval.parent(substitute(object@media[[mediac[i]]]@diffmat <- Matrix::Matrix(newdiffmat, sparse=TRUE)))
    }
  }else stop("Given matrix is not compatible in dimensions with the environment.")
})

#' @title Change substance concentration patterns in the environment according to a gradient
#'
#' @description The generic function \code{createGradient} changes specific substance concentration patterns in the environment.
#' @export
#' @rdname createGradient
#'
#' @param object An object of class Arena.
#' @param mediac A character vector giving the names of substances, which should be added to the environment (the default takes all possible substances).
#' @param position A character vector giving the position (top, bottom, right and left) of the gradient.
#' @param smax A number giving the maximum concentration of the substance.
#' @param add A boolean variable defining whether the amount of substance should be summed or replaced
#' @param steep A number between 0 and 1 giving the steepness of the gradient (concentration relative to the arena size).
#' @param unit A character used as chemical unit to set the amount of the substances to be added (valid values are: mmol/cell, mmol/cm2, mmol/arena, mM)
#' @details This function can be used to add gradients of specific substances in the environment. 
#' @seealso \code{\link{Arena-class}} and \code{\link{changeSub}} 
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            growthlimit=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' addOrg(arena,bac,amount=10) #add 10 organisms
#' addSubs(arena,30) #add all substances with no concentrations.
#' createGradient(arena,smax=50,mediac=c("EX_glc(e)","EX_o2(e)","EX_pi(e)"),
#'              position='top',steep=0.5, add=FALSE)
setGeneric("createGradient", function(object, mediac, position, smax, steep, add=FALSE, unit='mmol/cell'){standardGeneric("createGradient")})
#' @export
#' @rdname createGradient
setMethod("createGradient", "Arena", function(object, mediac, position, smax, steep, add=FALSE, unit='mmol/cell'){
  if(steep<=0 || steep>=1){stop("Steepness must be in between 0 and 1.")}
  if(length(intersect(unit,c("mmol/cell","mM","mmol/cm2","mmol/arena")))==0){stop("Wrong unit for concentration.")}
  mediac = intersect(mediac,object@mediac)
  if(unit=="mM"){smax <- (smax*0.01)*object@scale}  # conversion of mMol in mmol/grid_cell
  if(unit=="mmol/cm2"){smax <- smax*object@scale}  # conversion of mmol/arena in mmol/grid_cell
  if(unit=="mmol/arena"){smax <- smax*object@scale}  # conversion of mmol/arena in mmol/grid_cell
  newdiffmat <- matrix(0,nrow=object@n,ncol=object@m)
  gradn = floor(object@n*steep)
  gradm = floor(object@m*steep)
  switch(position,
         'top'={for(i in 1:object@m){newdiffmat[0:gradm+1,i]=rev(seq(0,smax,length.out=gradm+1))}},
         'bottom'={for(i in 1:object@m){newdiffmat[object@m:(object@m-gradm),i]=rev(seq(0,smax,length.out=gradm+1))}},
         'right'={for(i in 1:object@n){newdiffmat[i,gradn:object@n]=seq(0,smax,length.out=gradn+1)}},
         'left'={for(i in 1:object@n){newdiffmat[i,0:gradn+1]=seq(smax,0,length.out=gradn+1)}},
         stop("Positions must be top, bottom, right, or left."))
  for(i in 1:length(mediac)){
    if(add){
      eval.parent(substitute(object@media[[mediac[i]]]@diffmat <- Matrix::Matrix(as.matrix(object@media[[mediac[i]]]@diffmat)+newdiffmat, sparse=TRUE)))
    }else{
      eval.parent(substitute(object@media[[mediac[i]]]@diffmat <- Matrix::Matrix(newdiffmat, sparse=TRUE)))
    }
  }
})

#' @title Change organisms in the environment
#'
#' @description The generic function \code{changeOrg} changes organisms in the environment.
#' @export
#' @rdname changeOrg
#'
#' @param object An object of class Arena.
#' @param neworgdat A data frame with new information about the accumulated growth, type, phenotype, x and y position for each individual in the environment.
#' @details The argument \code{neworgdat} contains the same information as the \code{orgdat} slot of \code{\link{Arena-class}}. The \code{orgdat} slot of an \code{Arena} object can be used to create \code{neworgdat}.
#' @seealso \code{\link{Arena-class}} and \code{\link{addOrg}}
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            growthlimit=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' addOrg(arena,bac,amount=10) #add 10 organisms
#' neworgdat <- arena@orgdat #get the current orgdat
#' neworgdat <- neworgdat[-1,] #remove the first individual
#' changeOrg(arena,neworgdat)
setGeneric("changeOrg", function(object, neworgdat){standardGeneric("changeOrg")})
#' @export
#' @rdname changeOrg
setMethod("changeOrg", "Arena", function(object, neworgdat){
  eval.parent(substitute(object@orgdat <- neworgdat))
})

#' @title Function for checking phenotypes in the environment
#'
#' @description The generic function \code{checkPhen} checks and adds the phenotypes of organisms in the environment.
#' @export
#' @rdname checkPhen
#'
#' @param object An object of class Arena.
#' @param org An object of class Organism.
#' @param cutoff A number giving the cutoff for values of the objective function and fluxes of exchange reactions.
#' @return Returns a number indicating the number of the phenotype in the phenotype list.
#' @details The phenotypes are defined by flux through exchange reactions, which indicate potential differential substrate usages. Uptake of substances are indicated by a negative and production of substances by a positive number.
#' @seealso \code{\link{Arena-class}} and \code{\link{getPhenotype}}
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            growthlimit=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' checkPhen(arena,bac) #returns 1 as the index of the current phenotype in the list.
setGeneric("checkPhen", function(object, org, cutoff=1e-6){standardGeneric("checkPhen")})
#' @export
#' @rdname checkPhen
setMethod("checkPhen", "Arena", function(object, org, cutoff=1e-6){
  pind <- 0
  if(org@fbasol$obj>=cutoff){
    test = getPhenotype(org, cutoff=1e-6)
    tspec = org@type
    pvec = rep(0,length(object@mediac))
    names(pvec) = object@mediac
    pvec[names(test)] = test
    pvec <- paste(pvec,collapse='')
    phenc <- object@phenotypes
    phensel <- phenc[which(names(phenc)==tspec)]
    pind <- which(phensel==pvec)
    if(length(pind)==0){
      pind = length(phensel)+1
      names(pvec) = tspec
      eval.parent(substitute(object@phenotypes <- c(phenc,pvec)))
    }
  }
  return(pind)
})

#' @title Main function for simulating all processes in the environment
#'
#' @description The generic function \code{simEnv} for a simple simulation of the environment.
#' @export
#' @rdname simEnv
#'
#' @param object An object of class Arena or Eval.
#' @param time A number giving the number of iterations to perform for the simulation
#' @param lrw A numeric value needed by solver to estimate array size (by default lwr is estimated in the simEnv() by the function estimate_lrw())
#' @param continue A boolean indicating whether the simulation should be continued or restarted.
#' @return Returns an object of class \code{Eval} which can be used for subsequent analysis steps.
#' @details The returned object itself can be used for a subsequent simulation, due to the inheritance between \code{Eval} and \code{Arena}.
#' @seealso \code{\link{Arena-class}} and \code{\link{Eval-class}}
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            growthlimit=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' addOrg(arena,bac,amount=10) #add 10 organisms
#' addSubs(arena,40) #add all possible substances
#' eval <- simEnv(arena,10)
setGeneric("simEnv", function(object, time, lrw=NULL, continue=F){standardGeneric("simEnv")})
#' @export
#' @rdname simEnv
setMethod("simEnv", "Arena", function(object, time, lrw=NULL, continue=F){
  switch(class(object),
         "Arena"={arena <- object; evaluation <- Eval(arena)},
         "Eval"={arena <- getArena(object); evaluation <- object},
         stop("Please supply an object of class Arena."))
  if(is.null(lrw)){lrw=estimate_lrw(arena@n,arena@m)}
  for(i in names(arena@specs)){
    phensel <- arena@phenotypes[which(names(arena@phenotypes)==i)]
    if(length(phensel)==0){
      test = getPhenotype(arena@specs[[i]], cutoff=1e-6)
      pvec = rep(0,length(arena@mediac))
      names(pvec) = arena@mediac
      pvec[names(test)] = test
      pvec <- paste(pvec,collapse='')
      names(pvec) = i
      arena@phenotypes <- c(arena@phenotypes,pvec)
    }
  }
  if(class(object)!="Eval"){addEval(evaluation, arena)}
  sublb <- getSublb(arena)
  for(i in 1:time){
    cat("iter:", i, "Organisms:",nrow(arena@orgdat),"\n")
    arena@mflux <- lapply(arena@mflux, function(x){numeric(length(x))}) # empty mflux pool
    if(nrow(arena@orgdat) > 0){ # if there are organisms left
      sublb[,arena@mediac] = sublb[,arena@mediac]*(10^12) #convert to fmol per gridcell
      for(j in 1:nrow(arena@orgdat)){ # for each organism in arena
        org <- arena@specs[[arena@orgdat[j,'type']]]
        bacnum = round((arena@scale/(org@cellarea*10^(-8)))) #calculate the number of bacteria individuals per gridcell
        switch(class(org),
               "Bac"= {arena = simBac(org, arena, j, sublb, bacnum)}, #the sublb matrix will be modified within this function
               "Human"= {arena = simHum(org, arena, j, sublb, bacnum)}, #the sublb matrix will be modified within this function
               stop("Simulation function for Organism object not defined yet."))
      }
      sublb[,arena@mediac] = sublb[,arena@mediac]/(10^12) #convert again to mmol per gridcell
      test <- is.na(arena@orgdat$growth)
      if(sum(test)!=0) arena@orgdat <- arena@orgdat[-which(test),]
      rm("test")
    }
    if(!arena@stir){
      sublb_tmp <- matrix(0,nrow=nrow(arena@orgdat),ncol=(length(arena@mediac)))
      sublb <- as.data.frame(sublb) #convert to data.frame for faster processing in apply
      for(j in seq_along(arena@media)){ #get information from sublb matrix to media list
        submat <- as.matrix(arena@media[[j]]@diffmat)
        if(nrow(sublb) != sum(sublb[,j+2]==mean(submat))){
          apply(sublb[,c('x','y',arena@media[[j]]@name)],1,function(x){submat[x[1],x[2]] <<- x[3]})
        }
        #skip diffusion if already homogenous (attention in case of boundary/source influx in pde!)
        homogenous = arena@n*arena@m != sum(submat==mean(submat))
        diffspeed  = arena@media[[j]]@difspeed!=0
        diff2d     = arena@media[[j]]@pde=="Diff2d"
        if( diffspeed && ( diff2d&&homogenous || !diff2d ) ){  
          switch(arena@media[[j]]@difunc,
                 "pde"  = {submat <- diffusePDE(arena@media[[j]], submat, gridgeometry=arena@gridgeometry, lrw, tstep=object@tstep)},
                 "pde2" = {diffuseSteveCpp(submat, D=arena@media[[j]]@difspeed, h=1, tstep=arena@tstep)},
                 "naive"= {diffuseNaiveCpp(submat, donut=FALSE)},
                 "r"    = {for(k in 1:arena@media[[j]]@difspeed){diffuseR(arena@media[[j]])}},
                 stop("Diffusion function not defined yet.")) 
          arena@media[[j]]@diffmat <- Matrix::Matrix(submat, sparse=TRUE)
        }
        sublb_tmp[,j] <- apply(arena@orgdat, 1, function(x,sub){return(sub[x[4],x[5]])},sub=submat)
      }
      sublb <- cbind(as.matrix(arena@orgdat[,c(4,5)]),sublb_tmp)
      colnames(sublb) <- c('x','y',arena@mediac)
      rm("sublb_tmp")
      rm("submat")
    }else{
      sublb <- stirEnv(arena, sublb)
    }
    addEval(evaluation, arena)
    if(nrow(arena@orgdat)==0 & !continue){
      print("All organisms died!")
      break
    }
  }
  return(evaluation)
})

#' @title Function for calculated the substrate concentration for every organism
#'
#' @description The generic function \code{getSublb} calculates the substrate concentration for every individual in the environment based on their x and y position.
#' @export
#' @rdname getSublb
#'
#' @param object An object of class Arena.
#' @return Returns the substrate concentration for every individual in the environment with substrates as well as x and y positions as columns and rows for each organism.
#' @seealso \code{\link{Arena-class}}
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            growthlimit=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' addOrg(arena,bac,amount=10) #add 10 organisms
#' addSubs(arena,40) #add all possible substances
#' sublb <- getSublb(arena)
setGeneric("getSublb", function(object){standardGeneric("getSublb")})
#' @export
#' @rdname getSublb
setMethod("getSublb", "Arena", function(object){
  sublb <- matrix(0,nrow=nrow(object@orgdat),ncol=(length(object@mediac)))
  for(j in seq_along(object@media)){
    submat <- as.matrix(object@media[[j]]@diffmat)
    sublb[,j] <- apply(object@orgdat, 1, function(x,sub){return(sub[x[4],x[5]])},sub=submat)
  }
  sublb <- cbind(as.matrix(object@orgdat[,c(4,5)]),sublb)
  colnames(sublb) <- c('x','y',object@mediac)
  return(sublb)
})

#' @title Function for stirring/mixing the complete evironment
#'
#' @description The generic function \code{stirEnv} simulates the event of mixing all substrates and organisms in the environment.
#' @export
#' @rdname stirEnv
#'
#' @param object An object of class Arena.
#' @param sublb A matrix with the substrate concentration for every individual in the environment based on their x and y position.
#' @return Returns the substrate concentration for every individual in the environment with substrates as well as x and y positions as columns and rows for each organism.
#' @details The stirring is implemented as a random permutation of organism positions and the equalization of of all substrate concentrations.
#' @seealso \code{\link{Arena-class}} and \code{\link{getSublb}}
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            growthlimit=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' addOrg(arena,bac,amount=10) #add 10 organisms
#' addSubs(arena,40) #add all possible substances
#' sublb <- getSublb(arena)
#' stirEnv(arena,sublb)
setGeneric("stirEnv", function(object, sublb){standardGeneric("stirEnv")})
#' @export
#' @rdname stirEnv
setMethod("stirEnv", "Arena", function(object, sublb){
  #stir all the organism except the ones which are not moving
  neworgdat <- object@orgdat
  cmbs = expand.grid(1:object@n,1:object@m)
  sit <- which(unlist(lapply(object@specs,function(x)(return(x@speed))))==0)
  if(length(sit) != 0){
    siti <- which(neworgdat$type %in% sit)
    sitorgdat <- neworgdat[siti,]
    neworgdat <- neworgdat[-siti,]
    rownames(cmbs) <- paste(cmbs$Var1,cmbs$Var2,sep="_")
    cmbs <- cmbs[setdiff(rownames(cmbs),c(paste(sitorgdat$x,sitorgdat$y,sep="_"))),]
  }
  if(nrow(neworgdat) > nrow(cmbs)){ #not so nice -> there is a problem
    selength <- nrow(cmbs)
  }else{
    selength <- nrow(neworgdat)
  }
  sel <- sample(1:nrow(cmbs),selength)
  neworgdat[,'x'] <- cmbs[sel,1]
  neworgdat[,'y'] <- cmbs[sel,2]
  if(length(sit) != 0){neworgdat <- rbind(neworgdat,sitorgdat)}
  eval.parent(substitute(object@orgdat <- neworgdat))
  #stir all the substrates + modify substrates
  sublb_tmp <- matrix(0,nrow=nrow(object@orgdat),ncol=(length(object@mediac)))
  sublb <- as.data.frame(sublb) #convert to data.frame for faster processing in apply
  for(j in seq_along(object@media)){ #get information from sublb matrix to media list
    sval <- sum(sublb[,object@media[[j]]@name])/nrow(sublb)
    submat <- matrix(sval,object@n,object@m)
    eval.parent(substitute(object@media[[j]]@diffmat <- Matrix::Matrix(submat, sparse=TRUE)))
    sublb_tmp[,j] <- apply(object@orgdat, 1, function(x,sub){return(sub[x[4],x[5]])},sub=submat)
  }
  sublb <- cbind(as.matrix(object@orgdat[,c(4,5)]),sublb_tmp)
  colnames(sublb) <- c('x','y',object@mediac)
  return(sublb)
})

#' @title Function for transforming the organism data frame to a presence/absence matrix of organisms
#'
#' @description The generic function \code{dat2mat} simulates the event of mixing all substrates and organisms in the environment.
#' @export
#' @rdname dat2mat
#'
#' @param object An object of class Arena.
#' @return Returns the presence/absence matrix of organisms on the grid based on the \code{orgdat} slot of the \code{Arena} class.
#' @seealso \code{\link{Arena-class}} and \code{\link{getSublb}}
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            growthlimit=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' addOrg(arena,bac,amount=10) #add 10 organisms
#' occmat <- dat2mat(arena)
#' image(occmat)
setGeneric("dat2mat", function(object){standardGeneric("dat2mat")})
#' @export
#' @rdname dat2mat
setMethod("dat2mat", "Arena", function(object){
  newoccmat <- matrix(0,object@n,object@m)
  for(i in 1:nrow(object@orgdat)){
    newoccmat[object@orgdat[i,'y'],object@orgdat[i,'x']] = object@orgdat[i,'type']
  }
  return(newoccmat)
})

#show function for class Arena
setMethod(show, "Arena", function(object){
  ecoli_cellarea = 4.42
  ecoli_cellweight = 1172
  print(paste('Arena of size ',object@n,'x',object@m,' with ',nrow(object@orgdat),
              ' organisms of ',length(object@specs),' species.',sep=''))
  print(paste("flux unit:","mmol/(h*g_dw)"))
  print(paste("arena grid cells:",object@n,"x",object@m))
  print(paste("arena grid size [cm]:",object@Lx,"x",object@Ly))
  print(paste("area of one grid cell [cm^2]:", (object@Lx*object@Ly)/(object@n*object@m)))
  print(paste("maximal amount of E. coli cells in one grid cell:", round((object@Lx*object@Ly)/(object@n*object@m)/(ecoli_cellarea*10^(-8)),2) ))
  dwpgc <- (object@Lx*object@Ly)/(object@n*object@m)/(ecoli_cellarea*10^(-8))*ecoli_cellweight
  print(paste("maximal amount of E.coli dry weight in one grid cell [fg]:", round(dwpgc,1)))
  scaleF <- 15-3*floor(log10(dwpgc)/3)
  print(paste(3*floor(log10(dwpgc)/3)))
  scaleN <- c("Milli", "Mikro", "Nano", "Piko", "Femto", "Atto", "Zepto", "Yokto")
  print(paste("scale factor [10^-x]:", scaleF, scaleN[scaleF/3]))
  conc  <- 20 #mM
  apspgc <- conc * object@Lx * object@Ly / 1000 / (object@n * object@m)
  print(paste("concentration [mM]:", conc, "thus amount of substance per grid cell [mmol]:", apspgc ))
  print(paste("scaled: amount of substance per grid cell [", scaleN[scaleF/3+1], "mol]:", apspgc*10^(3*floor(log10(dwpgc)/3+3))))
  print(paste("scaled: maximal amount of E.coli dry weight in one grid cell [", scaleN[scaleF/3],"gram]:", dwpgc/10^(3*floor(log10(dwpgc)/3))))
})




# Eval is a subclass of Arena containing function to reduce the size of simulations and evalution of results

########################################################################################################
###################################### EVAL CLASS ######################################################
########################################################################################################

#' Structure of the S4 class "Eval"
#' 
#' Structure of the S4 class \code{Eval} inheriting from class \code{\link{Arena-class}} for the analysis of simulations.
#' @export Eval
#' @exportClass Eval
#' @importFrom  graphics barplot legend lines par points segments text
#' @importFrom stats cor dist hclust na.omit prcomp rnorm
#' @importFrom utils combn
#'
#' @slot medlist A list of compressed medium concentrations (only changes of concentrations are stored) per time step.
#' @slot simlist A list of the organism features per time step.
#' @slot mfluxlist A list of containing highly used metabolic reactions per time step. 
#' @slot subchange A vector of all substrates with numbers indicating the degree of change in the overall simulation.
setClass("Eval",
         contains="Arena",
         representation(
           medlist="list",
           simlist="list",
           mfluxlist="list",
           subchange="numeric"
         ),
         prototype(
           medlist = list(),
           simlist = list(),
           mfluxlist = list()
         )
)

########################################################################################################
###################################### CONSTRUCTOR #####################################################
########################################################################################################

#' Constructor of the S4 class \code{\link{Eval-class}}
#' 
#' @name Eval-constructor
#' @export
#' @param arena An object of class Arena.
Eval <- function(arena){
  subc = rep(0, length(arena@mediac))
  names(subc) <- arena@mediac
  new("Eval", n=arena@n, m=arena@m, tstep=arena@tstep, specs=arena@specs, mediac=arena@mediac, subchange=subc,
      phenotypes=arena@phenotypes, media=arena@media, orgdat=arena@orgdat, medlist=list(), simlist=list(), stir=arena@stir, mfluxlist=list())
}

########################################################################################################
###################################### GET METHODS FOR ATTRIBUTES ######################################
########################################################################################################

setGeneric("medlist", function(object){standardGeneric("medlist")})
setMethod("medlist", "Eval", function(object){return(object@medlist)})
setGeneric("simlist", function(object){standardGeneric("simlist")})
setMethod("simlist", "Eval", function(object){return(object@simlist)})
setGeneric("mfluxlist", function(object){standardGeneric("mfluxlist")})
setMethod("mfluxlist", "Eval", function(object){return(object@mfluxlist)})
setGeneric("subchange", function(object){standardGeneric("subchange")})
setMethod("subchange", "Eval", function(object){return(object@subchange)})

########################################################################################################
###################################### METHODS #########################################################
########################################################################################################

#' @title Function for adding a simulation step
#'
#' @description The generic function \code{addEval} adds results of a simulation step to an \code{Eval} object.
#' @export
#' @rdname addEval
#'
#' @param object An object of class Eval.
#' @param arena An object of class Arena.
#' @param replace A boolean variable indicating if the last simulation step should be replaced by the new simulation step \code{arena}.
#' @details The function \code{addEval} can be used in iterations to manipulate an \code{Arena} object and store the results in an \code{Eval} object.
#' @seealso \code{\link{Eval-class}} and \code{\link{Arena-class}}
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            growthlimit=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' addOrg(arena,bac,amount=10) #add 10 organisms
#' addSubs(arena,40) #add all possible substances
#' eval <- simEnv(arena,10)
#' addEval(eval,arena)
setGeneric("addEval", function(object, arena, replace=F){standardGeneric("addEval")})
#' @export
#' @rdname addEval
setMethod("addEval", "Eval", function(object, arena, replace=F){
  if(!replace){
    subch = rep(0, length(arena@mediac))
    if(length(object@medlist)!=0){
      names(subch) <- arena@mediac
      sapply(names(subch), function(x, oldmed, newmed){
        subch[x] <<- subch[x]+sum(abs(oldmed[[x]]-as.vector(newmed[[x]]@diffmat)))
      },oldmed=extractMed(object), newmed=arena@media)
      subch[names(object@subchange)] <- object@subchange + subch[names(object@subchange)]
      eval.parent(substitute(object@subchange <- subch))
    }
    if(sum(subch)!=0){
      eval.parent(substitute(object@medlist[[length(object@medlist)+1]] <- lapply(arena@media, function(x, subc){
        if(subc[x@name]!=0){
          return(as.vector(x@diffmat))
        }else{return(vector())}
      }, subc=subch)))
    }else{
      eval.parent(substitute(object@medlist[[length(object@medlist)+1]] <- lapply(arena@media, function(x){
        return(as.vector(x@diffmat))
      })))
    }
    eval.parent(substitute(object@simlist[[length(object@simlist)+1]] <- arena@orgdat))
    eval.parent(substitute(object@phenotypes <- arena@phenotypes))
    eval.parent(substitute(object@mfluxlist[[length(object@mfluxlist)+1]] <- arena@mflux))
  }else{
    eval.parent(substitute(object@medlist[[length(object@medlist)]] <- lapply(arena@media, function(x){
      return(as.vector(x@diffmat))
    })))
    eval.parent(substitute(object@simlist[[length(object@simlist)]] <- arena@orgdat))
    eval.parent(substitute(object@phenotypes <- arena@phenotypes))
    eval.parent(substitute(object@specs <- arena@specs)) 
    eval.parent(substitute(object@mediac <- arena@mediac))
    eval.parent(substitute(object@media <- arena@media))
  }
})

#' @title Function for re-constructing an Arena object from a simulation step
#'
#' @description The generic function \code{getArena} re-constructs an \code{Arena} object from a simulation step within an \code{Eval} object.
#' @export
#' @rdname getArena
#'
#' @param object An object of class Eval.
#' @param time A number giving the simulation step of interest.
#' @return Returns an object of class \code{Arena} containing the organisms and substance conditions in simulation step \code{time}.
#' @details The function \code{addEval} can be used to manipulate an \code{Arena} object from a simulation step to modify the subsequent simulation steps.
#' @seealso \code{\link{Eval-class}} and \code{\link{Arena-class}}
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            growthlimit=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' addOrg(arena,bac,amount=10) #add 10 organisms
#' addSubs(arena,40) #add all possible substances
#' eval <- simEnv(arena,10)
#' arena5 <- getArena(eval,5)
setGeneric("getArena", function(object, time=(length(object@medlist)-1)){standardGeneric("getArena")})
#' @export
#' @rdname getArena
setMethod("getArena", "Eval", function(object, time=(length(object@medlist)-1)){ #index in R start at 1, but the first state is 0
  time = time+1 #index in R start at 1, but the first state is 0
  newmedia <- lapply(object@media, function(x, meds, n, m){
    x@diffmat <- Matrix::Matrix(meds[[x@name]],nrow=n,ncol=m,sparse=TRUE)
    return(x)
  },meds=extractMed(object,time), n=object@n, m=object@m)
  occdat <- object@simlist[[time]]
  
  arena <- Arena(n=object@n, m=object@m, tstep=object@tstep, specs=object@specs, mediac=object@mediac,
                 phenotypes=object@phenotypes , media=newmedia, orgdat=occdat, stir=object@stir)
  return(arena)
})

#' @title Function for re-constructing a medium concentrations from simulations
#'
#' @description The generic function \code{extractMed} re-constructs a list of vectors of medium concentrations from a simulation step in an \code{Eval} object.
#' @export
#' @rdname extractMed
#'
#' @param object An object of class Eval.
#' @param time A number giving the simulation step of interest.
#' @return Returns a list containing concentration vectors of all medium substances.
#' @details Medium concentrations in slot \code{medlist} of an object of class \code{Eval} store only the changes of concentrations in the simulation process. The function \code{extractMed} reconstructs the original and uncompressed version of medium concentrations.
#' @seealso \code{\link{Eval-class}} and \code{\link{Arena-class}}
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            growthlimit=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' addOrg(arena,bac,amount=10) #add 10 organisms
#' addSubs(arena,40) #add all possible substances
#' eval <- simEnv(arena,10)
#' med5 <- extractMed(eval,5)
setGeneric("extractMed", function(object, time=length(object@medlist)){standardGeneric("extractMed")})
#' @export
#' @rdname extractMed
setMethod("extractMed", "Eval", function(object, time=length(object@medlist)){
  medl <- object@medlist
  medlind <- medl[[time]]
  for(i in 1:length(object@mediac)){
    if(length(medl[[time]][[i]])==0){
      j <- time
      while(length(medl[[j]][[i]])==0){j <- j-1}
      medlind[[i]] <- medl[[j]][[i]]
    }
  }
  return(medlind)
})

#' @title Function for plotting spatial and temporal change of populations and/or concentrations
#'
#' @description The generic function \code{evalArena} plots heatmaps from the simulation steps in an \code{Eval} object.
#' @export
#' @rdname evalArena
#'
#' @param object An object of class Eval.
#' @param plot_items A character vector giving the items, which should be plotted.
#' @param phencol A boolean variable indicating if the phenotypes of the organisms in the environment should be integrated as different colors in the population plot.
#' @param retdata A boolean variable indicating if the data used to generate the plots should be returned.
#' @param time A numeric vector giving the simulation steps which should be plotted.
#' @return Returns several plots of the chosen plot items. Optional the data to generate the original plots can be returned.
#' @details If \code{phencol} is \code{TRUE} then different phenotypes of the same organism are visualized by varying colors, otherwise different organism types are represented by varying colors. The parameter \code{retdata} can be used to access the data used for the returned plots to create own custom plots. 
#' @seealso \code{\link{Eval-class}} and \code{\link{Arena-class}}
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            growthlimit=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' addOrg(arena,bac,amount=10) #add 10 organisms
#' addSubs(arena,40) #add all possible substances
#' eval <- simEnv(arena,10)
#' evalArena(eval)
#'\dontrun{
#' ## if animation package is installed a movie of the simulation can be stored:
#' library(animation)
#' saveVideo({evalArena(eval)},video.name="Ecoli_sim.mp4")
#' }
setGeneric("evalArena", function(object, plot_items='Population', phencol=F, retdata=F, time=(seq_along(object@simlist)-1)){standardGeneric("evalArena")})
#' @export
#' @rdname evalArena
setMethod("evalArena", "Eval", function(object, plot_items='Population', phencol=F, retdata=F, time=(seq_along(object@simlist)-1)){ #index in R start at 1, but the first state is 0
  time = time+1
  #old.par <- par(no.readonly = TRUE)
  if(retdata){
    retlist = list()
    for(i in 1:length(plot_items)){
      retlist[[i]] = list()
    }
    names(retlist) = plot_items
  }
  for(i in time){
    subnam <- names(object@medlist[[i]])
    inds <- which(subnam %in% plot_items)
    meds <- extractMed(object, i)
    if(length(plot_items)==1){
      if(length(inds)!=0){
        for(j in 1:length(inds)){
          if(retdata){
            retlist[[subnam[inds[j]]]][[paste0("time",(i-1))]] = matrix(meds[[subnam[inds[j]]]],nrow=object@n,ncol=object@m)
          }
          image(matrix(meds[[subnam[inds[j]]]],nrow=object@n,ncol=object@m),axes=F,main=paste(subnam[inds[j]], ": #", i),
                zlim=c(0,max(unlist(lapply(object@medlist,function(x, snam){return(x[[snam]])},snam=subnam[inds[j]])))))
        }
      }
    }else if(length(plot_items)<=6){
      par(mfrow=c(2,ceiling(length(plot_items)/2)))
      for(j in 1:length(inds)){
        if(retdata){
          retlist[[subnam[inds[j]]]][[paste0("time",(i-1))]] = matrix(meds[[subnam[inds[j]]]],nrow=object@n,ncol=object@m)
        }
        image(matrix(meds[[subnam[inds[j]]]],nrow=object@n,ncol=object@m),axes=F,main=paste(subnam[inds[j]], ": #", i),
              zlim=c(0,max(unlist(lapply(object@medlist,function(x, snam){return(x[[snam]])},snam=subnam[inds[j]])))))
      }
    }else{
      par(mfrow=c(3,ceiling(length(plot_items)/3)))
      for(j in 1:length(inds)){
        if(retdata){
          retlist[[subnam[inds[j]]]][[paste0("time",(i-1))]] = matrix(meds[[subnam[inds[j]]]],nrow=object@n,ncol=object@m)
        }
        image(matrix(meds[[inds[j]]],nrow=object@n,ncol=object@m),axes=F,main=paste(subnam[inds[j]], ": #", i),
              zlim=c(0,max(unlist(lapply(object@medlist,function(x, snam){return(x[[snam]])},snam=subnam[inds[j]])))))
      }
    }
    if(plot_items[1]=='Population'){
      if(retdata){
        retlist[['Population']][[paste0("time",(i-1))]] = object@simlist[[i]]
      }
      if(phencol){
        plot(object@simlist[[i]][,c('x','y')],xlim=c(0,object@n),ylim=c(0,object@m),xlab='',ylab='',
             pch=object@simlist[[i]]$type-1,axes=FALSE,cex=1,main=paste('Population', ": #", i), col=object@simlist[[i]]$phenotype+1)
      }else{
        plot(object@simlist[[i]][,c('x','y')],xlim=c(0,object@n),ylim=c(0,object@m),xlab='',ylab='',
             pch=object@simlist[[i]]$type-1,axes=FALSE,cex=1,main=paste('Population', ": #", i), col=object@simlist[[i]]$type)

      }
    }
  }
  #title(main=paste("step",i), outer=TRUE)
  #par(old.par)
  if(retdata){
    return(retlist)
  }
})

#' @title Function for plotting the overall change as curves
#'
#' @description The generic function \code{plotCurves} plots the growth curves and concentration changes of substances from simulation steps in an \code{Eval} object.
#' @export
#' @rdname plotCurves
#'
#' @param object An object of class Eval.
#' @param medplot A character vector giving the name of substances which should be plotted.
#' @param retdata A boolean variable indicating if the data used to generate the plots should be returned.
#' @param remove A boolean variable indicating if substances, which don't change in their concentration should be removed from the plot.
#' @param legend Boolean variable indicating if legend should be plotted
#' @return Returns two graphs in one plot: the growth curves and the curves of concentration changes. Optional the data to generate the original plots can be returned.
#' @details The parameter \code{retdata} can be used to access the data used for the returned plots to create own custom plots. 
#' @seealso \code{\link{Eval-class}} and \code{\link{Arena-class}}
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            growthlimit=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' addOrg(arena,bac,amount=10) #add 10 organisms
#' addSubs(arena,40) #add all possible substances
#' eval <- simEnv(arena,10)
#' plotCurves(eval)
setGeneric("plotCurves", function(object, medplot=object@mediac, retdata=F, remove=F, legend=F){standardGeneric("plotCurves")})
#' @export
#' @rdname plotCurves
setMethod("plotCurves", "Eval", function(object, medplot=object@mediac, retdata=F, remove=F, legend=F){
  old.par <- par(no.readonly = TRUE)
  growths <- matrix(0, nrow=length(object@specs), ncol=length(object@simlist))
  rownames(growths) = names(object@specs)
  subs <- matrix(0, nrow=length(medplot), ncol=length(object@simlist))
  rownames(subs) = medplot
  for(i in seq_along(object@simlist)){
    simdat <- object@simlist[[i]]
    abun <- table(simdat$type)
    for(j in 1:length(abun)){
      growths[as.numeric(names(abun[j])),i] = abun[j]
    }
    subdat <- extractMed(object, i)
    for(j in 1:length(medplot)){
      subs[medplot[j],i] <- sum(subdat[[medplot[j]]])/(object@n*object@m)
    }
  }
  if(remove){
    test <- which(apply(subs,1,function(x){return(sum(ifelse(x==x[1],T,F))==length(x))}))
    if(length(test)!=0){
      subs <- subs[-test,]
    }
  }
  par(mfrow=c(2,1))
  times = c((seq_along(object@simlist)-1)*object@tstep)
  plot(times, times, xlim=c(0,max(times)), ylim=c(0,max(growths)),
       type='n', xlab='time in h', ylab='number of individuals on grid',
       main='Population')
  for(i in 1:nrow(growths)){
    lines(times, growths[i,], col=i, type='b', pch=i-1)
  }
  if(legend){legend('topleft',legend=rownames(growths),col=1:nrow(growths), cex=ifelse(length(object@specs)==1,1,0.4/log10(nrow(growths)+1)),pch=(0:nrow(growths)-1),lwd=1, bty="n")}
  plot(times, times, xlim=c(0,max(times)), ylim=c(0,max(subs)),
       type='n', xlab='time in h', ylab='concentration in mmol per gridcell',
       main='Substance concentrations')
  for(i in 1:nrow(subs)){
    lines(times, subs[i,], col=i)
  }
  par(old.par)
  if(legend){legend('right',legend=rownames(subs),col=1:nrow(subs),cex=0.4/log10(nrow(subs)+1),lwd=1)}
  if(retdata){
    return(list('Population'=growths,'Substances'=subs))
  }
})

#' @title Function for plotting the overall change as curves with maximally distinct colors
#'
#' @description The generic function \code{plotCurves2} plots the growth curves and concentration changes of the most changing substances from simulation steps in an \code{Eval} object using maximally distinct colors.
#' @export
#' @rdname plotCurves2
#'
#' @param object An object of class Eval.
#' @param legendpos A character variable declaring the position of the legend
#' @param ignore A list of character variables with substance names that sould be omitted in the plot
#' @param num An integer defining the number of substrates to be plot
#' @param phencol Boolean variable indicating whether phenotypes should be higlighted
#' @param dict List defining new substance names. List entries are intepreted as old names and the list names as the new ones.
#' @return Returns two graphs in one plot: the growth curves and the curves of concentration changes
#' @details The parameter \code{retdata} can be used to access the data used for the returned plots to create own custom plots. 
#' @seealso \code{\link{Eval-class}} and \code{\link{Arena-class}}
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            growthlimit=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' addOrg(arena,bac,amount=10) #add 10 organisms
#' addSubs(arena,40) #add all possible substances
#' eval <- simEnv(arena,10)
#' plotCurves2(eval)
setGeneric("plotCurves2", function(object, legendpos="topleft", ignore=c("EX_h(e)","EX_pi(e)", "EX_h2o(e)"),
                                   num=10, phencol=F, dict=NULL){standardGeneric("plotCurves2")})
#' @export
#' @rdname plotCurves2
setMethod("plotCurves2", "Eval", function(object, legendpos="topright", ignore=c("EX_h(e)","EX_pi(e)", "EX_h2o(e)"), 
                                          num=10, phencol=F, dict=NULL){
  if(num>length(object@mediac) || num<1) stop("Number of substances invalid")
  # first get the correct (ie. complete) medlist
  prelist <- lapply(seq_along(object@medlist), function(i){extractMed(object, i)})
  list <- lapply(prelist, function(x){lapply(x, sum)})
  mat <- matrix(unlist(list), nrow=length(object@media), ncol=length(object@medlist))
  #remove substances that should be ignored
  ignore_subs <- which(object@mediac %in% ignore | gsub("\\(e\\)","", gsub("EX_","",object@mediac)) %in% ignore)
  if(length(ignore_subs) != 0){
    mat <- mat[-ignore_subs,]
    mediac <- object@mediac[-ignore_subs]
  } else mediac <- object@mediac
  
  rownames(mat) <- gsub("\\(e\\)","", gsub("EX_","",mediac))
  mat_var  <- rowSums((mat - rowMeans(mat))^2)/(dim(mat)[2] - 1)
  mat_nice <- tail(mat[order(mat_var),], num)
  
  if(num>length(colpal3)) cols <- colpal1[1:num] else cols <- colpal3[1:num]
  matplot(t(mat_nice), type='l', col=cols, pch=1, lty=1, lwd=5,
          xlab=paste0('time in ', ifelse(object@tstep==1, "", object@tstep), 'h'), ylab='amount of substance in mmol',
          main='Strongly changing substances')
  if(length(dict) > 0){
    new_names = unlist(lapply(rownames(mat_nice), function(x){dict[[x]]}))
    legend(legendpos, new_names, col=cols, cex=0.7, fill=cols)
  } else legend(legendpos, rownames(mat_nice), col=cols, cex=0.7, fill=cols)
  
  # get bacs
  list <- lapply(object@simlist, function(x){
    occ <- table(x$type)
    unlist(lapply(seq_along(object@specs), function(i){ifelse(i %in% names(occ),occ[paste(i)], 0)})) # ugly ;P
  })
  mat_bac  <- do.call(cbind, list)
  rownames(mat_bac) <- names(object@specs)
  
  
  # biomass
  list <- lapply(object@simlist, function(x){
    sum(x$growth)
  })
  mat_biom  <- do.call(cbind, list)
  rownames(mat_biom) <- "biomass [fg]"
  
  # get pheno
  if(phencol){
    pheno_nr <- table(names(object@phenotypes))
    list <- lapply(object@simlist, function(x){ # time step
      unlist(lapply(seq_along(object@specs), function(j){ # bac type
        occ <- table(x[which(x$type==j),]$phenotype)
        #p <- unlist(lapply(seq_along(object@phenotypes[[j]]), function(i){ifelse(i %in% names(occ),occ[paste(i)], 0)})) # ugly ;P
        p <- unlist(lapply(seq(0,pheno_nr[[names(object@specs[j])]]), function(i){ifelse(i %in% names(occ),occ[paste(i)], 0)})) # ugly ;P
        names(p) <- paste0(names(object@specs)[j], "_pheno", seq(0,pheno_nr[[names(object@specs[j])]]))
        p
      }))})
    mat_phen  <- do.call(cbind, list)
    #mat_with_phen <- rbind(mat_bac, mat_phen, mat_biom)
    mat_with_phen <- rbind(mat_bac, mat_phen)
  } else{
    #mat_with_phen <- rbind(mat_bac, mat_biom)
    mat_with_phen <- mat_bac
  }

  
  
  
  len <- dim(mat_with_phen)[1]
  if(len>length(colpal3)) cols <- colpal1[1:len] else cols <- colpal3[1:len]
  matplot(t(mat_with_phen), type='b', col=cols, pch=1, lty=1, lwd=5,
          xlab=paste0('time in ', ifelse(object@tstep==1, "", object@tstep), 'h'), ylab='amount of organisms',
          main='Growth curve')
  legend(legendpos, rownames(mat_with_phen), col=cols, cex=0.7, fill=cols)
})


#' @title Function for plotting the overall change in reaction activity
#'
#' @description The generic function \code{plotTotFlux} plots the time course of reactions with high variation in activity for an \code{Eval} object.
#' @export
#' @rdname plotTotFlux
#'
#' @param object An object of class Eval.
#' @param legendpos A character variable declaring the position of the legend
#' @param num An integer defining the number of substrates to be plot
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            growthlimit=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' addOrg(arena,bac,amount=10) #add 10 organisms
#' addSubs(arena,40) #add all possible substances
#' eval <- simEnv(arena,10)
#' plotTotFlux(eval)
setGeneric("plotTotFlux", function(object, legendpos="topright", num=20){standardGeneric("plotTotFlux")})
#' @export
#' @rdname plotTotFlux
setMethod("plotTotFlux", "Eval", function(object, legendpos="topright", num=20){
  if(num<1) stop("Number of reactions invalid")
  list <- lapply(object@mfluxlist, function(x){
    unlist(x)
  })
  mat  <- do.call(cbind, list)
  mat_var  <- rowSums((mat - rowMeans(mat))^2)/(dim(mat)[2] - 1)
  mat_nice <- tail(mat[order(mat_var),], num)
  
  if(num>length(colpal3)) cols <- colpal1[1:num] else cols <- colpal3[1:num]
  matplot(t(mat_nice), type='l', col=cols, pch=1, lty=1, lwd=3,
          xlab='time in h', ylab='reaction activity in mmol/(h * g_DW)',
          main='Highly active reactions')
  legend(legendpos, rownames(mat_nice), col=cols, cex=0.6, fill=cols)
  return(rownames(mat_nice))
})


#' @title Function for getting a matrix of phenotypes from the dataset
#'
#' @description The generic function \code{getPhenoMat} reconstructs a matrix with the usage of exchange reactions of the different organisms in the environment.
#' @export
#' @rdname getPhenoMat
#'
#' @param object An object of class Eval.
#' @param time An integer indicating the time step to be used (default value is character "total")
#' @param sparse A boolean indicating whether zero entries should be removed from return matrix
#' @return Returns a matrix with different phenotypes of the organism as rows and all possible exchange reactions as columns. A value of 1 means secretion, 2 means uptake and 0 means no usage of the substance of interest.
#' @details The phenotypes are defined by flux through exchange reactions, which indicate potential differential substrate usages.
#' @seealso \code{\link{Eval-class}} and \code{\link{getPhenotype}}
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            growthlimit=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' addOrg(arena,bac,amount=10) #add 10 organisms
#' addSubs(arena,40) #add all possible substances
#' eval <- simEnv(arena,10)
#' phenmat <- getPhenoMat(eval)
setGeneric("getPhenoMat", function(object, time="total", sparse=F){standardGeneric("getPhenoMat")})
#' @export
#' @rdname getPhenoMat
setMethod("getPhenoMat", "Eval", function(object, time="total", sparse=F){
  phenspec = list()
  for(i in names(object@specs)){
    phenspec[[i]] = object@phenotypes[which(names(object@phenotypes)==i)]
    names(phenspec[[i]]) = seq_along(phenspec[[i]])
  }
  phens = unlist(phenspec)
  if(time != "total"){
    time = time+1 #index in R start at 1, but the first state is 0
    tdat = object@simlist[[time]][,c('type','phenotype')]
    #tdat = tdat[-which(is.na(tdat$phenotype)),]
    tphen = factor(paste(names(object@specs)[tdat$type],tdat$phenotype,sep="."))
    pinds = which(names(phens) %in% levels(tphen))
    if(length(pinds)!=0){phens = phens[pinds]}
  }
  phenmat <- matrix(0, nrow=length(phens), ncol=length(object@mediac))
  colnames(phenmat) <- object@mediac
  rownames(phenmat) <- names(phens)
  for(i in 1:nrow(phenmat)){
    phenmat[i,] = as.numeric(unlist(strsplit(phens[i],split={})))
  }
  if(sparse){
    phenmat <- phenmat[,which(colSums(abs(phenmat))!=0)]
  }
  return(phenmat)
})

#' @title Function for mining/analyzing phenotypes which occured on the arena
#'
#' @description The generic function \code{minePheno} mines the similarity and differences of phenotypes reconstructed by \code{getPhenoMat} for each simulation step in an \code{Eval} object.
#' @export
#' @rdname minePheno
#'
#' @param object An object of class Eval.
#' @param plot_type A character vector giving the plot which should be returned (either "pca" for a principle coordinate analysis or "hclust" for hierarchical clustering).
#' @param legend Boolean variable indicating if legend should be plotted
#' @return Returns a plot for each simulation step representing the similarity of phenotypes of organisms within the environment. 
#' @param time An integer indicating the time step to be used (default value is character "total")
#' @details The phenotypes are defined by flux through exchange reactions, which indicate potential differential substrate usages.
#' @seealso \code{\link{Eval-class}} and \code{\link{getPhenoMat}}
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            growthlimit=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' addOrg(arena,bac,amount=10) #add 10 organisms
#' addSubs(arena,40) #add all possible substances
#' eval <- simEnv(arena,10)
#' minePheno(eval)
setGeneric("minePheno", function(object, plot_type="pca", legend=F, time="total"){standardGeneric("minePheno")})
#' @export
#' @rdname minePheno
setMethod("minePheno", "Eval", function(object, plot_type="pca", legend=F, time="total"){
  phenmat <- getPhenoMat(object, time)
  splitr = strsplit(rownames(phenmat),'\\.')
  rownames(phenmat) = unlist(lapply(splitr, function(x){return(paste(x[1:(length(x)-1)],collapse='.'))}))
  if(nrow(phenmat)<=1){
    stop('not enough phenotypes to analyze.')
  }
  old.par <- par(no.readonly = TRUE)
  pcount <- as.vector(table(rownames(phenmat)))
  plabs <- vector()
  plabs2 <- vector()
  for(i in 1:length(pcount)){
    plabs <- c(plabs,paste(rep(levels(as.factor(rownames(phenmat)))[i],pcount[i]),1:pcount[i],sep='_'))
    plabs2 <- c(plabs2,paste(i,1:pcount[i],sep='_'))
  }
  par(mfrow=c(1,1))
  if(plot_type=="pca"){
    phenpca <- prcomp(phenmat)
    plot(phenpca$x[,1:2], xlab='PC1', ylab='PC2', pch=20, cex=0.8, col=as.numeric(as.factor(rownames(phenpca$x))))
    text(phenpca$x[,1:2], labels=plabs2, col=as.numeric(as.factor(rownames(phenpca$x))), cex=0.7)
    typ <- names(object@specs)
    for(i in 1:length(object@specs)){
      segments(phenpca$x[which(rownames(phenpca$x)==typ[i]),1],phenpca$x[which(rownames(phenpca$x)==typ[i]),2],
               mean(phenpca$x[which(rownames(phenpca$x)==typ[i]),1]),mean(phenpca$x[which(rownames(phenpca$x)==typ[i]),2]),col=i)
      points(mean(phenpca$x[which(rownames(phenpca$x)==typ[i]),1]),mean(phenpca$x[which(rownames(phenpca$x)==typ[i]),2]),col=i,pch=i-1,cex=1.5)
    }
    if(legend){legend('topright',legend=names(object@specs),col=1:length(object@specs),pch=(0:length(object@specs)-1),cex=0.9,lwd=4)}
  }
  if(plot_type=="hclust"){
    rownames(phenmat) <- plabs
    plot(hclust(dist(phenmat)))
    par(old.par)
  }
})

#' @title Function for selecting phenotypes which occured on the arena from specific iterations and species
#'
#' @description The generic function \code{selPheno} selects phenotypes from specific simulation step in an \code{Eval} object.
#' @export
#' @rdname selPheno
#'
#' @param object An object of class Eval.
#' @param time A numeric vector giving the simulation steps which should be plotted. 
#' @param type A names indicating the species of interest in the arena.
#' @param reduce A boolean variable indicating if the resulting matrix should be reduced.
#' @return Returns a matrix with the substrate usage and the number of individuals using the phenotype. 
#' @details The phenotypes are defined by flux through exchange reactions, which indicate potential differential substrate usages.
#' @seealso \code{\link{Eval-class}} and \code{\link{getPhenoMat}}
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            growthlimit=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' addOrg(arena,bac,amount=10) #add 10 organisms
#' addSubs(arena,40) #add all possible substances
#' eval <- simEnv(arena,10)
#' selPheno(eval,time=10,type='ecoli_core_model',reduce=TRUE)
setGeneric("selPheno", function(object, time, type, reduce=F){standardGeneric("selPheno")})
#' @export
#' @rdname selPheno
setMethod("selPheno", "Eval", function(object, time, type, reduce=F){
  arena = getArena(object, time)
  type_num = which(names(arena@specs)==type)
  arena@orgdat[which(is.na(arena@orgdat$phenotype)),] = 0
  pabund = as.matrix(table(arena@orgdat[which(arena@orgdat$type == type_num),'phenotype']))
  rownames(pabund) = paste(type,'phen',rownames(pabund),sep='_')
  rownames(pabund)[which(rownames(pabund)==paste(type,'phen_0',sep='_'))] = 'inactive'
  colnames(pabund) = 'individuals'
  pmat = getPhenoMat(object, time)
  splitr = strsplit(rownames(pmat),'\\.')
  rownames(pmat) = unlist(lapply(splitr, function(x){return(paste(x[1:(length(x)-1)],collapse='.'))}))
  pmatsp = pmat[which(rownames(pmat) == type),]
  if(is.vector(pmatsp)){
    pmatsp = t(as.matrix(pmatsp))
  }else{
    pmatsp = as.matrix(pmatsp)
  }
  if(reduce){
    if(nrow(pmatsp)==1){
      pmatsp = t(as.matrix(pmatsp[,-which(apply(pmatsp,2,sum)==0)]))
    }else{
      pmatsp = pmatsp[,-which(apply(pmatsp,2,sum)==0)]
    }
  }
  if(length(grep('inactive',rownames(pabund)))!=0){
    rownames(pmatsp) = rownames(pabund)[-which(rownames(pabund)=="inactive")]
    pmatsp = as.data.frame(pmatsp)
    pmatsp['inactive',]=rep(0,ncol(pmatsp))
  }else{
    rownames(pmatsp) = rownames(pabund)
    pmatsp = as.data.frame(pmatsp)
  }
  pmatsp[,'individuals']=rep(NA,nrow(pmatsp))
  pmatsp[rownames(pabund),'individuals'] = pabund[,'individuals']
  pmatsp
  return(as.matrix(pmatsp))
})

#show function for class Eval
setMethod(show, signature(object="Eval"), function(object){
  print(paste('Evaluation results of ',length(object@medlist)-1,' simulation steps.',sep=''))
})



#' @title Function for investigating a specific phenotype of an organism
#'
#' @description The generic function \code{statPheno} provides statistical and visual information about a certain phenotype.
#' @export
#' @rdname statPheno
#'
#' @param object An object of class Eval.
#' @param type_nr A number indicating the Organism type of the phenotype to be investigated (from orgdat)
#' @param phenotype_nr A number indicating the phenotype to be investigated (from orgdat)
#' @param dict A character vector of all substance IDs with names that should be used instead of possibly cryptic IDs
#' @details The phenotypes are defined by flux through exchange reactions, which indicate potential differential substrate usages.
#' @seealso \code{\link{Eval-class}}
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            growthlimit=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' addOrg(arena,bac,amount=10) #add 10 organisms
#' addSubs(arena,40) #add all possible substances
#' eval <- simEnv(arena,10)
#' statPheno(eval, type_nr=1, phenotype_nr=2)
setGeneric("statPheno", function(object, type_nr=1, phenotype_nr, dict=NULL){standardGeneric("statPheno")})
#' @export
#' @rdname statPheno
setMethod("statPheno", "Eval", function(object, type_nr=1, phenotype_nr, dict=NULL){
  spec <- names(object@specs[type_nr])
  all  <- object@phenotypes[ names(object@phenotypes) == spec ]
  phen <- all[phenotype_nr]
  mediac <- gsub("\\(e\\)","", gsub("EX_","",object@mediac))
  if(length(dict) > 0) mediac <- unlist(lapply(mediac, function(x){dict[[x]]}))
  
  v <- as.numeric(unlist(strsplit(phen, split={})))
  names(v) <- mediac
  
  cat(paste("\nproduced by", spec,"\n"))
  print(names(v[which(v==1)]))
  cat(paste("\nconsumed by", spec,"\n"))
  print(names(v[which(v==2)]))
  
  # substrates
  prelist <- lapply(seq_along(object@medlist), function(i){extractMed(object, i)})
  list <- lapply(prelist, function(x){lapply(x, sum)})
  mat_sub <- matrix(unlist(list), nrow=length(object@media), ncol=length(object@medlist))
  rownames(mat_sub) <- mediac
  
  occ <- unlist(lapply(seq_along(object@simlist), function(t){dim(object@simlist[[t]][which(object@simlist[[t]]$type==type_nr & object@simlist[[t]]$phenotype==phenotype_nr),])[1]}))
  if(sum(occ)==0){
    cat("\n occurence of phenotype over time\n")
    print(occ)
    stop("This phenotype has never lived?!")
  }
  t_lb <- min(which(occ>0))
  t_ub <- max(which(occ>0))
  cat(paste("\n", spec, "occured between time steps"))
  print(paste0(t_lb, "-", t_ub))
  
  # get substances whose dynamic correlate with investigated phenotype
  mat <- rbind(mat_sub, phenotype=occ)[,t_lb:t_ub]
  corr <- cor(t(mat))
  corr[is.na(corr)] <- 0
  sorted <- sort(corr[dim(mat)[1],])
  high_corr <- c(tail(sorted, n=5), head(sorted, n=5))
  cat("\nsubstance with highest correlation\n")
  print(high_corr)
  
  if( length(high_corr) > 0){
    par(mfrow=c(2,1))
    mat_nice <- mat_sub[names(high_corr)[names(high_corr)!="phenotype"],][,t_lb:t_ub]
    num <- dim(mat_nice)[1]
    if(num>length(colpal3)) cols <- colpal1[1:num] else cols <- colpal3[1:num]
    matplot(x=seq(t_lb, t_ub), y=t(mat_nice), type='l', col=cols, pch=1, lty=1, lwd=5,
            xlab=paste0('time in ', ifelse(object@tstep==1, "", object@tstep), 'h'), ylab='amount of substance in mmol',
            main='correlated substances')
    legend("topleft", rownames(mat_nice), col=cols, cex=0.4, fill=cols)
    
    mat_pheno <- mat["phenotype",]
    plot(y=mat_pheno, x=seq(t_lb, t_ub), type='l', col="black", pch=1, lty=1, lwd=5,
            xlab=paste0('time in ', ifelse(object@tstep==1, "", object@tstep), 'h'), ylab='amount of organisms',
            main=paste('Growth curve', spec, "phenotype", phenotype_nr))
    
    par(mfrow=c(1,1))
  }
  
})


#' @title Function for investigation of feeding between phenotypes
#'
#' @description The generic function \code{findFeeding} 
#' @export
#' @rdname findFeeding
#' @importFrom igraph graph.empty add.edges delete.edges delete.vertices V E degree vcount layout_with_fr
#' 
#' @param object An object of class Eval.
#' @param tcut Integer giving the minimal mutual occurence ot be considered (dismiss very seldom feedings)
#' @param scut List of substance names which should be ignored
#' @param legendpos A character variable declaring the position of the legend
#' @param dict List defining new substance names. List entries are intepreted as old names and the list names as the new ones.
#' @return Graph (igraph)
#' 
setGeneric("findFeeding", function(object, dict=NULL, tcut=5, scut=list(), legendpos="topleft"){standardGeneric("findFeeding")})
#' @export
#' @rdname findFeeding
setMethod("findFeeding", "Eval", function(object, dict=NULL, tcut=5, scut=list(), legendpos="topleft"){

  # possible problem inactive phenotype is not mentioned in object@phenotypes...

  # 1) Time: get occupation matrix for all phenotypes (occ_phen)
  pheno_nr <- table(names(object@phenotypes))
  pheno_index <- vector()
  list <- lapply(object@simlist, function(x){ # time step
    unlist(lapply(seq_along(object@specs), function(j){ # bac type
      occ <- table(x[which(x$type==j),]$phenotype)
      p <- unlist(lapply(seq(0,pheno_nr[[names(object@specs[j])]]), function(i){ifelse(i %in% names(occ),occ[paste(i)], 0)})) # ugly ;P
      names(p) <- paste0(names(object@specs)[j], "_pheno", seq(0,pheno_nr[[names(object@specs[j])]]))
      p
    }))})
  
  mat_phen  <- do.call(cbind, list)
  mat_tmp <- mat_phen
  mat_tmp[mat_tmp>0] <- 1
  cutoff <- which(rowSums(mat_tmp)>tcut)
  if(length(cutoff) < 2) {
    stop("tcut too high, no phenotypes found")
  } else mat_phen <- mat_phen[which(rowSums(mat_tmp)>tcut),] # reduce considered phenotypes (should exists >= tcut time steps)
  occ_phen  <- mat_phen != 0

  # 2) Substances: get matrix of substrates that could consumed and produced by phenotypes in principle
  mediac <- gsub("\\(e\\)","", gsub("EX_","",object@mediac))
  if(length(dict) > 0) mediac <- unlist(lapply(mediac, function(x){dict[[x]]}))
  phens <- object@phenotypes
  phenmat <- matrix(0, nrow=length(phens), ncol=length(object@mediac))
  colnames(phenmat) <- mediac
  counter = vector("numeric", length(object@specs))
  names(counter) <- names(object@specs)
  new_names = unlist(lapply(names(phens), function(x){
    counter[x] <<- counter[x] + 1
    paste0(x, "_pheno", counter[x])
  }))
  rownames(phenmat) <- new_names
  for(i in 1:nrow(phenmat)){
    phenmat[i,] = as.numeric(unlist(strsplit(phens[i],split={})))
  }
  phenmat_bin <- replace(phenmat, phenmat==2, -1)
  phenmat_abs <- abs(phenmat_bin)
  res <- phenmat_bin[,which(abs(colSums(phenmat_bin)) != colSums(phenmat_abs))]
  if(length(scut)>0) res <- res[,-which(colnames(res) %in% scut)] # reduce substrates

  # graph
  pindex <- rownames(mat_phen)# phenotype index
  cindex <- colnames(res) # substance color index
  g <-igraph::graph.empty(n=length(pindex), directed=TRUE)
  igraph::V(g)$name <- gsub("pheno","",pindex)
  

  # 3) Combinatorics: check for all pairs of phenotypes if they 
  #    i) occured in the same time steps and 
  #   ii) exchange at least one substance
  #combi <- combn(rownames(phenmat), 2)
  combi <- combn(rownames(mat_phen), 2)
  for(i in 1:ncol(combi)){
    if(length(grep("pheno0", combi[,i])) > 0) next # do not consider cuples with phenotype0 (i.d. metabolic inactive)
    co_occ <- which(occ_phen[combi[,i][1],] & occ_phen[combi[,i][2],]==T)
    if( length(co_occ) > 0 ){
      ex_both     <- res[c(combi[,i][1], combi[,i][2]),]
      feeding_index <- which(colSums(ex_both)==0 & colSums(abs(ex_both))!=0)
      if(length(feeding_index)>0){
        # if only one substance is exchanged some hack to get name of substance into returned data structure of feeding
        if(length(feeding_index)==1){
          feeding <- unlist(list(colnames(ex_both)[feeding_index], ex_both[,feeding_index]))
        } else {
          feeding <- ex_both[,feeding_index]
          lapply(seq(dim(feeding)[2]), function(x){
            if(feeding[1,x] == -1){
              new_edge <- c(which(pindex==combi[,i][1]), which(pindex==combi[,i][2]))
            }else new_edge <- c(which(pindex==combi[,i][2]), which(pindex==combi[,i][1]))
            col <- colpal3[which(cindex == colnames(feeding)[x])]
            g <<- igraph::add.edges(g, new_edge, color=col, weight=length(co_occ))
          })
        }
        #cat("\npossible cross feeding at time steps\n")
        #print(co_occ)
        #print(feeding)
      }
    }
  }
  
  g <- igraph::delete.edges(g, which(igraph::E(g)$weight<tcut)) # delete seldom feedings
  g <- igraph::delete.vertices(g, which(igraph::degree(g, mode="all") == 0)) # delete unconnected
  
  if(igraph::vcount(g) >= 2){
    plot(g, layout=igraph::layout_with_fr, vertex.size=5,
         edge.arrow.size=0.3, edge.width=E(g)$weight/10)
    legend(legendpos,legend=cindex, col=colpal3, pch=19, cex=0.7)
  }
  return(g)
})




setGeneric("statSpec", function(object, type_nr=1, dict=NULL,
                                legend_show=TRUE, legend_pos="center", legend_cex=0.75){standardGeneric("statSpec")})
setMethod("statSpec", "Eval", function(object, type_nr=1, dict=NULL, 
                                       legend_show=TRUE, legend_pos="center", legend_cex=0.75){
  if(type_nr <= 0 | type_nr > length(object@specs)){
    stop("Invalid type number, should be number indicating a species of arena@specs")
  }
  sname <- names(object@specs[type_nr])
  pheno_nr <- table(names(object@phenotypes))[[sname]]
  print(paste(sname, "with phenotypes:", pheno_nr))
  
  # zero-phenotyp is inactive!
  phens <- rep(0, length(object@phenotypes)+1) # total number of occuring cells for each phenotyp
  ptimes <- rep(0, length(object@phenotypes)+1) # number of times steps a phenotyp occured
  availTimes <- rep(FALSE, length(object@simlist)) # true/false if cell is living for each timestep
  occ <- lapply(seq(from=2, to=length(object@simlist)), function(t){
    res <- table(object@simlist[[t]][which(object@simlist[[t]]$type==type_nr),]$phenotype)
    #print(phens)
    phens[as.numeric(names(res))+1]   <<- phens[as.numeric(names(res))+1] + res # -1 because of zero-phenotyp
    ptimes[as.numeric(names(res))+1]  <<- ptimes[as.numeric(names(res))+1] + 1
    if(length(res) > 0) availTimes[t] <<- TRUE
    p <- unlist(lapply(seq(0,pheno_nr), function(i){ifelse(i %in% names(res),res[paste(i)], 0)}))
    names(p) <- paste0("pheno", seq(0,pheno_nr))
    p
  })
  
  # species is available in the following time steps
  availTimesteps <- which(availTimes==TRUE)
  print(paste("alive between time steps:", min(availTimesteps), "-", max(availTimesteps)))
  
  len <- length(object@phenotypes)+1
  if(len>length(colpal3)) cols <- colpal1[1:len] else cols <- colpal3[1:len]

  # plot growth curve with all phenotypes
  mat_phen  <- do.call(cbind, occ)
  matplot(t(mat_phen), type='b', col=cols, pch=1, lty=1, lwd=5,
          xlab=paste0('time in ', ifelse(object@tstep==1, "", object@tstep), 'h'), ylab='amount of organisms',
          main='Growth curve')
  legend(legend_pos, rownames(mat_phen), col=cols, cex=legend_cex, fill=cols)  
  
  # analyze phenotypes: plot total cell number vs. time steps of occurence
  plot(ptimes, phens, bg=cols,type="p", pch = 23, log="y",
       xlab="Number of occured time steps ", ylab="Total number of cells", main=paste(sname,"phenotypes: cells vs. occurrence"))
  if(legend_show){
    legend(legend_pos, pch = 23, pt.bg=cols, cex=legend_cex,
           legend=seq(0, length(object@phenotypes)))
  }
})



#' @title Function to compute and return correlation matrix
#'
#' @description The generic function \code{getCorrM} returns the correlation matrix of several objects.
#' @export
#' @rdname getCorrM
#'
#' @param object An object of class Eval.
#' @param reactions A boolean indicating whether reactions should be included in correlation matrix
#' @param bacs A boolean indicating whether bacteria should be included in correlation matrix
#' @param substrates A boolean indicating whether substrates should be included in correlation matrix
#' @return correlation matrix
#' @details Returns correlation matrix which can be used for statistical analysis
#' @seealso \code{\link{Eval-class}}
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            growthlimit=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' addOrg(arena,bac,amount=10) #add 10 organisms
#' addSubs(arena,40) #add all possible substances
#' eval <- simEnv(arena,10)
#' getCorrM(eval)
setGeneric("getCorrM", function(object, reactions=TRUE, bacs=TRUE, substrates=TRUE){standardGeneric("getCorrM")})
#' @export
#' @rdname getCorrM
setMethod("getCorrM", "Eval", function(object, reactions=TRUE, bacs=TRUE, substrates=TRUE){
  mat <- matrix(0,0,length(object@medlist))
  if(substrates){
    prelist <- lapply(seq_along(object@medlist), function(i){extractMed(object, i)})
    list <- lapply(prelist, function(x){lapply(x, sum)})
    mat_sub <- matrix(unlist(list), nrow=length(object@media), ncol=length(object@medlist))
    rownames(mat_sub) <- gsub("\\(e\\)","", gsub("EX_","",object@mediac))
    mat <- rbind(mat, mat_sub)
  }
  
  if(reactions){
    list <- lapply(object@mfluxlist, function(x){
      unlist(x)
    })
    mat_rea  <- do.call(cbind, list)
    mat <- rbind(mat, mat_rea)
  }
  
  if(bacs){
    list <- lapply(object@simlist, function(x){
      occ <- table(x$type)
      unlist(lapply(seq_along(object@specs), function(i){ifelse(i %in% names(occ),occ[paste(i)], 0)})) # ugly ;P
    })
    mat_bac  <- do.call(cbind, list)
    rownames(mat_bac) <- names(object@specs)
    mat <- rbind(mat, mat_bac)
  }
  
  corr <- cor(t(mat))
  corr[is.na(corr)] <- 0
  return(corr)
})
  
#' @title Function to show correlations of a simulated organism or substrate
#'
#' @description The generic function \code{checkCorr} returns the correlation matrix of several objects.
#' @export
#' @rdname checkCorr
#'
#' @param object An object of class Eval.
#' @param corr A correlation matrix (\code{\link{getCorrM}})
#' @param tocheck A list with substrate, reactions or organism names whose correlations should be shown
#' @details Returns correlation matrix which can be used for statistical analysis
#' @seealso \code{\link{Eval-class}} and \code{\link{getCorrM}}
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            growthlimit=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' addOrg(arena,bac,amount=10) #add 10 organisms
#' addSubs(arena,40) #add all possible substances
#' eval <- simEnv(arena,10)
#' checkCorr(eval, tocheck="o2")
setGeneric("checkCorr", function(object, corr=NULL, tocheck=list()){standardGeneric("checkCorr")})
#' @export
#' @rdname checkCorr
setMethod("checkCorr", "Eval", function(object, corr=NULL, tocheck=list()){
  if(is.null(corr)) corr <- getCorrM(object)

  lapply(tocheck, function(feature){
    dat <- corr[feature,]
    dat <- dat[-which(names(dat)==feature)]
    dat <- dat[order(dat)]
    dat_interest <- c(head(dat), tail(dat))
    barplot(names.arg=names(dat_interest), height=dat_interest, las=2, main=paste("Highest correlations of", feature))
  })
})


