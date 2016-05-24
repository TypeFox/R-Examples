

setMethod("initialize", "yuima.poisson",
          function(.Object,
                   drift = expression() ,
                   diffusion = list() ,
                   hurst = 0.5,
                   jump.coeff = expression(),
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
                   J.flag = logical()){
            .Object@drift <- drift
            .Object@diffusion <- diffusion
            .Object@hurst <- 0.5
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
            return(.Object)
          })

setPoisson <- function(intensity=1, df=NULL, scale=1, dimension=1, ...){
    
    poisson.model <- setModel(drift="0", diffusion="0",
    jump.coeff=as.expression(scale),
    measure=list(intensity=intensity, df=df),
    measure.type="CP",...)
    
 
    CPoisson <- new("yuima.poisson",
    drift=poisson.model@drift,
    diffusion = poisson.model@diffusion,
    hurst=poisson.model@hurst,
    jump.coeff=poisson.model@jump.coeff,
    measure=poisson.model@measure,
    measure.type=poisson.model@measure.type,
    parameter=poisson.model@parameter,
    state.variable=poisson.model@state.variable,
    jump.variable=poisson.model@jump.variable,
    time.variable=poisson.model@time.variable,
    noise.number =  poisson.model@noise.number,
    equation.number = dimension,
    dimension = dimension,
    solve.variable=poisson.model@solve.variable,
    xinit=poisson.model@xinit,
    J.flag = poisson.model@J.flag
    )
    
    return(CPoisson)
    
}

dconst <- function(x,k=1) k*(x==k)
rconst <- function(n,k=1) rep(k,n)
