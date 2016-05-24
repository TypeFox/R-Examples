setMethod("initialize", "cogarch.info",
          function(.Object,
                   p=numeric(),
                   q=numeric(),
                   ar.par=character(),
                   ma.par=character(),
                   loc.par=character(),
                   Cogarch.var=character(),
                   V.var=character(),
                   Latent.var=character(),
                   XinExpr=logical(),
                   measure=list(),
                   measure.type=character()){
            .Object@p <- p
            .Object@q <- q
            .Object@ar.par <- ar.par
            .Object@ma.par <- ma.par
            .Object@loc.par <- loc.par
            .Object@Cogarch.var <- Cogarch.var
            .Object@V.var <- V.var
            .Object@Latent.var <- Latent.var
            .Object@XinExpr <- XinExpr
            .Object@measure<-measure
            .Object@measure.type<-measure.type  
            return(.Object)
          })

setMethod("initialize", "yuima.cogarch",
          function(.Object,
                   info = new("cogarch.info"),
                   drift = expression() ,
                   diffusion = list() ,
                   hurst = 0.5,
                   jump.coeff = list(),
                   measure=list(),
                   measure.type=character(),
                   parameter = new("model.parameter"),
                   state.variable = "x",
                   jump.variable = "z",
                   time.variable = "t",
                   noise.number = numeric(),
                   equation.number = numeric(),
                   dimension = numeric(),
                   solve.variable = character(),
                   xinit = expression(),
                   J.flag = logical()
                  # quadr.var=new("quadratic.variation")
                   ){
            .Object@info <- info
            .Object@drift <- drift
            .Object@diffusion <- diffusion
            .Object@hurst <- hurst  	   
            .Object@jump.coeff <- jump.coeff
            .Object@measure <- measure
            .Object@measure.type <- measure.type
            .Object@parameter <- parameter
            .Object@state.variable <- state.variable
            .Object@jump.variable <- jump.variable
            .Object@time.variable <- time.variable
            .Object@noise.number <- noise.number
            .Object@equation.number <- equation.number
            .Object@dimension <- dimension
            .Object@solve.variable <- solve.variable
            .Object@xinit <- xinit
            .Object@J.flag <- J.flag
            #.Object@quadr.var <- quadr.var
            return(.Object)
          })


setCogarch<-function(p,
                   q,
                   ar.par="b",
                   ma.par="a",
                   loc.par="a0",
                   Cogarch.var="g",
                   V.var="v",
                   Latent.var="y",
                   jump.variable="z",                  
                   time.variable="t",
                   measure=NULL,
                   measure.type=NULL,
                   XinExpr=FALSE,
                   startCogarch=0,
                   work=FALSE,...){
# We use the same parametrization as in Brockwell (2000) 
  call <- match.call()
  mydots <- as.list(call)[-(1:2)]
  
  
  if (work==TRUE){
    Model_Cogarch=NULL
    return(Model_Cogarch)
  }
  if(is.null(measure)){
    yuima.warn("The Levy measure must be specified. The Brownian 
                 is not allowed as underlying process for Cogarch model.")
    return(NULL)
  }
  if(is.null(measure.type)){
    yuima.warn("The measure.type must be specified. See yuima documentation.")
    return(NULL)
  }
  # We need to build the auxiliar carma model that is the variance process
  # If we consider a Cogarch(p,q) the variance process is a Carma(q,p-1) model
  if(is.null(mydots$xinit)){
    aux.Carma<-setCarma(p=q,
                        q=p-1,
                        ar.par=ar.par,
                        ma.par=ma.par,
                        loc.par=loc.par,
                        lin.par=ma.par,
                        Carma.var=V.var,
                        Latent.var=Latent.var,
                        XinExpr=XinExpr,
                        Cogarch=TRUE)
    #  In order to have a representation of a Cogarch(p,q) model coherent with the 
    # Chaadra Brockwell and Davis we need to modify the slot xinit and drift[1]
    
  }else{
    if(!is.null(mydots$xinit)){
      aux.Carma<-setCarma(p=q,
                          q=p-1,
                          loc.par=loc.par,
                          ar.par=ar.par,
                          ma.par=ma.par,
                          lin.par=ma.par,
                          Carma.var=V.var,
                          Latent.var=Latent.var,
                          xinit=mydots$xinit,
                          Cogarch=TRUE)
      #  In order to have a representation of a Cogarch(p,q) model coherent with the 
      # Chaadra Brockwell and Davis we need to modify the slot xinit and drift[1]
    }
  }   

  newdrift<-c(0,as.character(aux.Carma@drift))
  newdiffusion<-c(0,as.character(eval(aux.Carma@diffusion)))
  line1<-c(paste0("sqrt(",V.var,")"),as.character(matrix(0,nrow=(q+1),ncol=1)))
  
  dumm<-as.character(eval(aux.Carma@diffusion))
  len.dumm<-length(dumm)
  line2<-character(length = (len.dumm+1))
  line2[1]<-"0"
  for(i in 1:length(dumm)){
    dumm.tmp<-nchar(x=dumm[i])
    line2[i+1]<-substring(dumm[i],13,dumm.tmp-2)  
  }
  
  state<-c(Cogarch.var,aux.Carma@state.variable)
  # We need now to modify the setModel in order to introduce a new component that represents
  # the discrete parts of the quadratic variation.
  # A possible solution is to write a new class that extends the yuima.model class and gives the possibility 
  # to add the additional components.
  # The way to introduce the differentiation of the discrete part of the quadratic variation have to follow the same
  # structure used for the jump component.
  Lev.coeff<-matrix(c(line1,line2),ncol=2)
  resdummy1<-setModel(drift=newdrift,
                 jump.coeff=Lev.coeff,
                 #quadr.coeff=line2,
                 measure=measure,
                 measure.type=measure.type,
                 #quadr.measure=measure,
                 #quadr.measure.type=measure.type,
                 state.variable=state,
                 jump.variable=jump.variable,
                 #quadr.variable=quadr.variable,
                 time.variable=time.variable,
                 xinit=c(startCogarch,as.character(aux.Carma@xinit)))
  
  resdummy2 <- new("cogarch.info",
                   p=p,
                   q=q,
                   ar.par=ar.par,
                   ma.par=ma.par,
                   loc.par=loc.par,
                   Cogarch.var=Cogarch.var,
                   V.var=V.var,
                   Latent.var=Latent.var,
                   XinExpr=XinExpr,
                   measure=measure,
                   measure.type=measure.type
                   )
  
#   res.quadr.var<-new("quadratic.variation",
#                      quadr.coeff=resdummy1@quadr.var@quadr.coeff,
#                      measure=resdummy1@quadr.var@measure,
#                      measure.type=resdummy1@quadr.var@measure.type,
#                      parms.quadr=resdummy1@quadr.var@parms.quadr,
#                      parms.quadr.meas = resdummy1@quadr.var@parms.quadr.meas,
#                      quadr.variable = resdummy1@quadr.var@quadr.variable,
#                      Q.flag = resdummy1@quadr.var@Q.flag)
  
  res<-new("yuima.cogarch",
           info=resdummy2,
           drift= resdummy1@drift,
           diffusion= resdummy1@diffusion,
           hurst=resdummy1@hurst,
           jump.coeff=resdummy1@jump.coeff,
           measure= resdummy1@measure,
           measure.type= resdummy1@measure.type,
           parameter= resdummy1@parameter, 
           state.variable= resdummy1@state.variable,
           jump.variable= resdummy1@jump.variable,
           time.variable= resdummy1@time.variable,
           noise.number= resdummy1@noise.number,
           equation.number= resdummy1@equation.number,
           dimension= resdummy1@dimension,
           solve.variable= resdummy1@solve.variable,
           xinit= resdummy1@xinit,
           J.flag = resdummy1@J.flag
          # quadr.var=res.quadr.var
           )
           
  return(res)
}
