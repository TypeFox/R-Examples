setMethod("initialize", "carma.info",
          function(.Object,
                   p=numeric(),
                   q=numeric(),
                   loc.par=character(),
                   scale.par=character(),
                   ar.par=character(),
                   ma.par=character(),
                   lin.par=character(),
                   Carma.var=character(),
                   Latent.var=character(),
                   XinExpr=logical()){
            .Object@p <- p
            .Object@q <- q
            .Object@loc.par <- loc.par
            .Object@scale.par <- scale.par
            .Object@ar.par <- ar.par
            .Object@ma.par <- ma.par
            .Object@lin.par <- lin.par
            .Object@Carma.var <- Carma.var
            .Object@Latent.var <- Latent.var
            .Object@XinExpr <- XinExpr
            return(.Object)
          })

setMethod("initialize", "yuima.carma",
          function(.Object,
                   info = new("carma.info"),
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
            return(.Object)
          })


setCarma<-function(p,
                   q,
                   loc.par=NULL,
                   scale.par=NULL,
                   ar.par="a",
                   ma.par="b",
                   lin.par=NULL,
                   Carma.var="v",
                   Latent.var="x",
                   XinExpr=FALSE,
                   Cogarch=FALSE,
                   ...){
# We use the same parametrization as in Brockwell (2000) 
  
# mydots$Carma.var= V
# mydots$Latent.var= X ?????  
  # The Carma model is defined by two equations:  
  # \begin{eqnarray} 
  # V_{t}&=&\alpha_{0}+a'Y_{t-}\\
  # dY_{t}&=& BY_{t-}dt+ e dL\\
  # \end{eqnarray}  
  # p is the number of the moving average parameters \alpha
  # q  is the number of the autoregressive coefficient \beta
                     
#Default parameters
  if (is.null(scale.par)){
    ma.par1<-ma.par
  } else{
    ma.par1<-paste(scale.par,ma.par,sep="*")
  }
  
  
  call <- match.call()
  quadratic_variation<-FALSE
    
  mydots <- as.list(call)[-(1:2)]
  

#   hurst<-0.5 
#   jump.coeff <- character()
#   measure <- list() 
#   measure.type <- character() 
#   state.variable <- "Y"
#   jump.variable <- "z" 
#   time.variable <- "t"
#   mydots$xinit<- NULL
  
  if (is.null(mydots$hurst)){
    mydots$hurst<-0.5
  }
  
  if(is.null(mydots$time.variable)){
    mydots$time.variable<-"t"
  }

  if(is.null(mydots$jump.variable)){
    mydots$jump.variable<-"z"
  }
    
#   if(is.null(mydots$xinit)){ 
#     if(is.null(mydots$XinExpr)){
#       mydots$xinit<-as.character(0*c(1:p))
#     }else{  
#       if(mydots$XinExpr==TRUE){
#         Int.Var<-paste(Latent.var,"0",sep="")
#         mydots$xinit<-paste(Int.Var,c(0:(p-1)),sep="")
#       }
#     }
#   } else{
#     dummy<-as.character(mydots$xinit)
#     mydots$xinit<-dummy[-1]
#   }
  
  
  if(is.null(mydots$xinit)){ 
    if(XinExpr==FALSE){
      mydots$xinit<-as.character(0*c(1:p))
    }else{  
      if(XinExpr==TRUE){
        Int.Var<-paste(Latent.var,"0",sep="")
        if(Cogarch==FALSE){
          mydots$xinit<-paste(Int.Var,c(0:(p-1)),sep="")
        }else{
          mydots$xinit<-paste(Int.Var,c(1:p),sep="")
        }
      }
    }
  } else{
    dummy<-as.character(mydots$xinit)
    mydots$xinit<-dummy[-1]
  }
  
  if(p<q){
    yuima.stop("order of AR must be larger than MA order")
  }
    
  beta_coeff0<-paste("-",ar.par,sep="")
  beta_coeff<-paste(beta_coeff0,p:1,sep="")
  if(Cogarch==FALSE){
    coeff_alpha<-c(paste(ma.par1,0:q,sep=""),as.character(matrix(0,1,p-q-1)))
  }else{
    coeff_alpha<-c(paste(ma.par1,1:(q+1),sep=""),as.character(matrix(0,1,p-q-1)))
  }
  fin_alp<-length(coeff_alpha)
  if(Cogarch==FALSE){
    Y_coeff<-paste(Latent.var,0:(p-1),sep="")
  }else{
    Y_coeff<-paste(Latent.var,1:p,sep="")
  }
  fin_Y<-length(Y_coeff)
  V1<-paste(coeff_alpha,Y_coeff,sep="*")
  V2<-paste(V1,collapse="+")
#   alpha0<-paste(ma.par1,0,sep="")
  
  if(is.null(loc.par)){
    Vt<-V2
    V<-paste0("(",V2,")",collapse="")
  } else {
    Vt<-paste(loc.par,V2,sep="+") 
    V<-paste0("(",Vt,")",collapse="")
  }
  drift_last_cond<-paste(paste(beta_coeff,Y_coeff,sep="*"),collapse="")
  # Drift condition for the dV_{t}   
  
  drift_first_cond_1<-c(paste(coeff_alpha[-fin_alp],Y_coeff[-1],sep="*"))
  drift_first_cond_2<-paste(drift_first_cond_1,collapse="+")
  drift_first_cond_a<-paste("(",drift_last_cond,")",sep="")
  drift_first_cond_b<-paste(coeff_alpha[fin_alp],drift_first_cond_a,sep="*")
  drift_first_cond<-paste(drift_first_cond_2,drift_first_cond_b,sep="+")
  
  if(length(Y_coeff)>1)
  {drift_Carma<-c(drift_first_cond,Y_coeff[2:length(Y_coeff)],drift_last_cond)}else 
  {drift_Carma<-c(drift_first_cond,drift_last_cond)}
  # We need to consider three different situations   
  
  if(is.null(mydots$jump.coeff) & is.null(mydots$measure)  &  
        is.null(mydots$measure.type) & quadratic_variation==FALSE){
    # The Carma model is driven by a Brwonian motion
    if (is.null(lin.par)){ 
      diffusion_Carma<-matrix(c(coeff_alpha[fin_alp],as.character(matrix(0,(p-1),1)),"1"),(p+1),1)
      #     Latent.var<-Y_coeff
      Model_Carma<-setModel(drift=drift_Carma, 
                             diffusion=diffusion_Carma,
                             hurst=mydots$hurst, 
                             state.variable=c(Carma.var,Y_coeff),  
                             solve.variable=c(Carma.var,Y_coeff),
                             xinit=c(Vt,mydots$xinit))
      #25/11
#       
#       carma.info<-new("carma.info",
#                       p=p,
#                       q=q,
#                       loc.par="character",
#                       scale.par="character",
#                       ar.par=ar.par,
#                       ma.par=ma.par,
#                       Carma.var=Carma.var,
#                       Latent.var=Latent.var,
#                       XinExpr=XinExpr)
      if(length(Model_Carma)==0){
        stop("Yuima model was not built") 
      } else { 
      #     return(Model_Carma1)
      }
    } else{
      if(ma.par==lin.par){
        first_term<-paste(coeff_alpha[fin_alp],V,sep="*")
        diffusion_Carma<-matrix(c(first_term,as.character(matrix(0,(p-1),1)),V),(p+1),1)  
        if(Cogarch==FALSE){
          Model_Carma<-setModel(drift=drift_Carma, 
                               diffusion=diffusion_Carma,
                               hurst=mydots$hurst, 
                               state.variable=c(Carma.var,Y_coeff),  
                               solve.variable=c(Carma.var,Y_coeff),
                               xinit=c(Vt,mydots$xinit))
        }else{# We add this part to have as initial condition for the COGARCH model V_0=a_0+a'X_0 LM 13/02/2015
           V01<-paste(coeff_alpha,mydots$xinit,sep="*")
           V02<-paste(V01,collapse="+")
           V0t<-paste(loc.par,V02,sep="+") 
           
           
           Model_Carma<-setModel(drift=drift_Carma, 
                              diffusion=diffusion_Carma,
                              hurst=mydots$hurst, 
                              state.variable=c(Carma.var,Y_coeff),  
                              solve.variable=c(Carma.var,Y_coeff),
                              xinit=c(V0t,mydots$xinit))
        }
#         return(Model_Carma1)
      }else{ 
#         coeff_gamma<-c(paste(lin.par,1:p,sep=""),as.character(matrix(0,1,p-q)))
        coeff_gamma<-c(paste(lin.par,1:p,sep=""))
        Gamma1<-paste(coeff_gamma,Y_coeff,sep="*")
        Gamma2<-paste(Gamma1,collapse="+")
        gamma0<-paste(lin.par,0,sep="")
        Gammat<-paste(gamma0,Gamma2,sep="+") 
        Gamma<-paste("(",Gammat,")",collapse="")
        first_term<-paste(coeff_alpha[fin_alp],Gamma,sep="*")
        
        diffusion_Carma<-matrix(c(first_term,as.character(matrix(0,(p-1),1)),Gamma),(p+1),1)  
#         Model_Carma1<-setModel(drift=drift_Carma, 
#                                diffusion=diffusion_Carma,
#                                hurst=mydots$hurst, 
#                                state.variable=c(Carma.var,Y_coeff),  
#                                solve.variable=c(Carma.var,Y_coeff),
#                                xinit=c(V,mydots$xinit))

        Model_Carma<-setModel(drift=drift_Carma, 
                               diffusion=diffusion_Carma,
                               hurst=mydots$hurst, 
                               state.variable=c(Carma.var,Y_coeff),  
                               solve.variable=c(Carma.var,Y_coeff),
                               xinit=c(Vt,mydots$xinit))
        
#         return(Model_Carma1)
      }
      
    }                      
  } else {
    if( is.null(mydots$jump.coeff) & is.null(mydots$measure)  &  
          is.null(mydots$measure.type) & is.null(lin.par) & 
          quadratic_variation==TRUE){
      
      stop("The CoGarch model requires a Carma process driven by the discrete part of the quadratic covariation:
           You Must specify the Levy Measure")
      
    } else {
      
      if(quadratic_variation==FALSE & is.null(lin.par)){
        #         warning("At the moment, we need specify the underlying L?vy directly")
        #         diffusion_Carma<-matrix(c(coeff_alpha[fin_alp],as.character(matrix(0,(q-1),1)),"1"),(q+1),1)
        #         Model_Carma1<-setModel(drift=drift_Carma, diffusion=diffusion_Carma,
        #                                hurst=hurst, state.variable=c(V,Y_coeff),  solve.variable=c(V,Y_coeff))
        #         under_Lev1<-setModel(drift="0",diffusion="0",jump.coeff="1" ,
        #                              measure=measure ,measure.type=measure.type , 
        #                              jump.variable=jump.variable , time.variable=time.variable)
        #         if(length(Model_Carma1)==0){
        #           stop("Yuima model was not built") 
        #         } else {
        #           Model_Carma<-Carma_Model()
        #           Model_Carma@model <- Model_Carma1
        #           Model_Carma@Cogarch_Model_Log <- Cogarch_Model 
        #           Model_Carma@Under_Lev <-under_Lev1
        #           return(Model_Carma)
        #         }
        
        # LM 27/09 We use a modified 
        # setModel that allows us to build a sde where the slot model@jump.coeff is an vector
        
        # jump_Carma<-matrix(c(coeff_alpha[fin_alp],as.character(matrix(0,(q-1),1)),"1"),(q+1),1)
        jump_Carma<-c(coeff_alpha[fin_alp],as.character(matrix(0,(p-1),1)),"1")
#         Model_Carma<-setModel(drift=drift_Carma, 
#                               diffusion = NULL, 
#                               hurst=mydots$hurst, 
#                               jump.coeff=jump_Carma,
#                               measure=eval(mydots$measure),
#                               measure.type=mydots$measure.type, 
#                               jump.variable=mydots$jump.variable, 
#                               time.variable=mydots$time.variable,
#                               state.variable=c(Carma.var,Y_coeff),  
#                               solve.variable=c(Carma.var,Y_coeff),
#                               xinit=c(V,mydots$xinit))
#         
        Model_Carma<-setModel(drift=drift_Carma, 
                              diffusion = NULL, 
                              hurst=mydots$hurst, 
                              jump.coeff=jump_Carma,
                              measure=eval(mydots$measure),
                              measure.type=mydots$measure.type, 
                              jump.variable=mydots$jump.variable, 
                              time.variable=mydots$time.variable,
                              state.variable=c(Carma.var,Y_coeff),  
                              solve.variable=c(Carma.var,Y_coeff),
                              xinit=c(Vt,mydots$xinit))
  #      return(Model_Carma)
      } else {
        if (quadratic_variation==FALSE ){
        # Selecting Quadratic_Variation==FALSE and specifying the Heteroskedatic.param in the model, 
        # The coefficient of the error term is a vector in which the last element is an affine linear function 
        # of the vector space Y               
        
        # We have to consider two different cases:
        # a) The last component of the error term is  $V_{t-}=\alpha_{0}+a'Y_{t-}$. Usually 
        #    this condition is used for the definition of the COGARCH(p,q) introduced in Brockwell and Davis     and 
        # b) The last component of the error term is a linear function of the state variable $Y_{t}$
        #     different of the observable variable V.  
        if(ma.par==lin.par){
          jump_first_comp<-paste(coeff_alpha[fin_alp],V,sep="*")
          jump_Carma<-c(jump_first_comp,as.character(matrix(0,(p-1),1)),V)
        }else{
#           coeff_gamma<-c(paste(lin.par,1:p,sep=""),as.character(matrix(0,1,q-p)))
          coeff_gamma<-c(paste(lin.par,1:p,sep=""))
          Gamma1<-paste(coeff_gamma,Y_coeff,sep="*")
          Gamma2<-paste(Gamma1,collapse="+")
          gamma0<-paste(lin.par,0,sep="")
          Gammat<-paste(gamma0,Gamma2,sep="+") 
          Gamma<-paste("(",Gammat,")",collapse="")
          jump_first_comp<-paste(coeff_alpha[fin_alp],Gamma,sep="*")
          jump_Carma<-c(jump_first_comp,as.character(matrix(0,(p-1),1)),Gamma)
        }
        
#         Model_Carma<-setModel(drift=drift_Carma,
#                               diffusion =NULL,
#                               hurst=0.5,
#                               jump.coeff=jump_Carma,
#                               measure=eval(mydots$measure),
#                               measure.type=mydots$measure.type,
#                               jump.variable=mydots$jump.variable,
#                               time.variable=mydots$time.variable,
#                               state.variable=c(Carma.var,Y_coeff),  
#                               solve.variable=c(Carma.var,Y_coeff),
#                               c(V,mydots$xinit))
#         return(Model_Carma)
        }
        
        Model_Carma<-setModel(drift=drift_Carma,
                              diffusion =NULL,
                              hurst=mydots$hurst,
                              jump.coeff=jump_Carma,
                              measure=eval(mydots$measure),
                              measure.type=mydots$measure.type,
                              jump.variable=mydots$jump.variable,
                              time.variable=mydots$time.variable,
                              state.variable=c(Carma.var,Y_coeff),  
                              solve.variable=c(Carma.var,Y_coeff),
                              xinit=c(Vt,mydots$xinit))
   #     return(Model_Carma)
         if(quadratic_variation==TRUE){
#           
                stop("Work in Progress: Implementation of CARMA model for CoGarch. 
                     We need the increments of Quadratic Variation")
#           
#           diffusion_first_cond<-paste(coeff_alpha[fin_alp],V,sep="*")
#           diffusion_Carma<-matrix(c(diffusion_first_cond,as.character(matrix(0,(q-1),1)),Vt),(q+1),1)
#           #    At the present time, Yuima does not support Multi - Jumps 
#           Model_Carma1<-setModel(drift=drift_Carma, diffusion=diffusion_Carma,
#                                  hurst=hurst, state.variable=c(V,Y_coeff),  solve.variable=c(V,Y_coeff))
#           under_Lev1<-setModel(drift="0",diffusion="0",jump.coeff="1" ,
#                                measure=measure ,measure.type=measure.type , 
#                                jump.variable=jump.variable , time.variable=time.variable)
#           if(length(Model_Carma1)==0){
#             stop("Yuima model was not built") 
#           } else {
#             Model_Carma<-Carma_Model()
#             Model_Carma@model <- Model_Carma1
#             Model_Carma@Cogarch_Model_Log <- Cogarch_Model 
#             Model_Carma@Under_Lev <-under_Lev1
#             return(Model_Carma)
#           }
#           
         } 
      } 
    }
    
  }
  # 25/11
  if(is.null(loc.par)){loc.par<-character()}
  if(is.null(scale.par)){scale.par<-character()}
  if(is.null(lin.par)){lin.par<-character()}
  
  
  carmainfo<-new("carma.info",
                  p=p,
                  q=q,
                  loc.par=loc.par,
                  scale.par=scale.par,
                  ar.par=ar.par,
                  ma.par=ma.par,
                  lin.par=lin.par,
                  Carma.var=Carma.var,
                  Latent.var=Latent.var,
                  XinExpr=XinExpr)
  
  Model_Carma1<-new("yuima.carma",
                    info=carmainfo,
                    drift=Model_Carma@drift,            
                    diffusion =Model_Carma@diffusion,
                    hurst=Model_Carma@hurst,
                    jump.coeff=Model_Carma@jump.coeff,
                    measure=Model_Carma@measure,
                    measure.type=Model_Carma@measure.type,
                    parameter=Model_Carma@parameter,
                    state.variable=Model_Carma@state.variable,
                    jump.variable=Model_Carma@jump.variable,
                    time.variable=Model_Carma@time.variable,
                    noise.number =  Model_Carma@noise.number,
                    equation.number = Model_Carma@equation.number,
                    dimension = Model_Carma@dimension, 
                    solve.variable=Model_Carma@solve.variable,
                    xinit=Model_Carma@xinit,
                    J.flag = Model_Carma@J.flag
                    )
  
  return(Model_Carma1)
}
