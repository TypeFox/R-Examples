GPCE.quad <- function(InputDim, # =c(3), # random variable definition
                      PCSpace, # ="Uniform",
                      InputDistrib, # =rep('Uniform',InputDim),
                      ParamDistrib, # =NULL,
                      Model=NULL, # =NULL, # definition of numerical model, i.e. the blackbox
                      ModelParam=NULL,
                      Output, #=NULL,
                      DesignInput, #=NULL,
                      # PCE expansion definition
                      p, #=c(2),                             # polynomial expansion (cardinal order)
                      ExpPoly, #=rep("LEGENDRE",InputDim),   # polynomial to use
                      # quadrature definitions
                      QuadType, # =c("FULL"), # or "SPARSE"
                      QuadPoly, # =rep("LEGENDRE",InputDim),   # polynomial to use
                      QuadLevel# =c(3)#rep(3,InputDim)
                      ){
  ### initialization Output
  ResultObject =list()

  ### Check Distribution
  OkDistrib=c("Uniform","Gaussian","Gamma","Beta")
  if (PCSpace!="FromData"){if(prod(InputDistrib%in%OkDistrib)==0){stop("the given distributions cannot be treated.")}}
    
  ### Get Design
  Args=list(InputDim=InputDim,
            PCSpace=PCSpace,
            InputDistrib=InputDistrib,
            ParamDistrib=ParamDistrib,
            Model=Model,
            Output=Output,
            DesignInput=DesignInput,
            p=p,
            ExpPoly=ExpPoly, # rep("LEGENDRE",InputDim),   # polynomial to use
            QuadType=QuadType,
            QuadPoly=QuadPoly, # rep("LEGENDRE",InputDim),   # polynomial to use
            QuadLevel=QuadLevel)
  
  ResultObject$Design <- CreateQuadrature(InputDim,QuadLevel,QuadPoly,ExpPoly,QuadType,ParamDistrib) # quadrature
  
  # NEED FUTURE CORRECTION
  Physical <- list(mu=0,sd=1)
  ResultObject$Design$PhysicalNodes <- Physical$mu + ResultObject$Design$QuadNodes*Physical$sd

  x= c(ResultObject,list(Design2Eval=ResultObject$Design$PhysicalNodes,Args=Args))
  class(x) <- "GPCE.quad"

  ### Case Model=NULL return Input design
  if(is.null(Output)){ # there is not output results 
      if(is.null(Model)){ # no output results or model -> return DEO
          return(x)
      } else { # not output but there is model -> run model
          ### Get Output when Model is given, function sampling
          if (is.null(ModelParam)){
            # print('no Output, no ModelParam, yes Model')
            Output <- Model(ResultObject$Design$PhysicalNodes)}
          else{
            # print('no Output, yes ModelParam, yes Model')
            Output <- Model(ResultObject$Design$PhysicalNodes,ModelParam)}
      }       
  } else { # we are given output determined from some model a priori
    if (identical(DesignInput,ResultObject$Design$QuadNodes)){ # check if output evaluated at right pts
      warning("Outputs evalauted at the correct quadrature points. \n")
    } else {
      stop("The input design at which the outputs are evaluated are different from the quadrature requested.")
    }
  }
  ResultObject <- tell.GPCE.quad(x,Output)
  return(ResultObject)          
}