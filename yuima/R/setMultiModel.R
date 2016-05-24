setMethod("initialize", "yuima.multimodel",
          function(.Object,
                   #model = new("yuima.model")
                   drift = expression() ,
                   diffusion = list() ,
                   hurst = 0.5,
                   jump.coeff = list(),
                   #jump.coeff = expression(),
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
                   ){
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



# We need a function that construct a Multivariate Model
setMultiModel <- function(drift=NULL,
                     diffusion=NULL,
                     hurst=0.5,
                     jump.coeff=NULL,
                     intensity = NULL,
                     df = NULL,
#                     jump.dimens = NULL,
#                     measure=list(),
                     measure.type=character(),
                     state.variable="x",
                     jump.variable="z",
                     time.variable="t",
                     solve.variable,
                     xinit=NULL){
  ## we need a temp env for simplifications

  yuimaENV <- new.env()
  ##::measure and jump term #####################################

  ##::initialize objects ########
  MEASURE <- list()

  ##::end initialize objects ########

  ##::make type function list ####
  CPlist <- c("dnorm", "dgamma", "dexp", "dconst")
  codelist <- c("rIG", "rNIG", "rgamma", "rbgamma", "rngamma", "rstable")
  ##::end make type function list ####
  jump.dimens <- dim(jump.coeff)[2]
  numbMeasure <- length(measure.type)
  if(numbMeasure>0){
    if(numbMeasure!=1){
      if(jump.dimens!=numbMeasure){
        yuima.stop("dimension of jump is not coherent")
      }
    }


#   if(!length(measure.type)){
#     if( length(jump.coeff) || length(measure) ){
#       yuima.warn("measure type does not match with jump term.")
#       return(NULL)
#     }
#     jump.variable <- character()

    measure.par <- character()
      if(any(measure.type=="CP")){
        tmp.measure <- list(df=list(func=vector(mode="list", length=1),
        expr=as.expression(rep("0",1))),
        intensity=as.expression(rep("0",length(intensity))))
      }else{

        tmp.measure <- list(df=list(func=vector(mode="list", length=length(measure.type)),
                                  expr=as.expression(rep("0",1))))

      }
#     measure.par <- character()
#   }else{
#     tmp.measure <- list(df=list(func=vector(mode="list",
#       length=length(measure.type)),
#       expr=as.expression(rep("0",length(measure.type)))),
#       intensity=as.expression(rep("0",sum(measure.type=="CP"))))
#     if( !length(jump.coeff) || !length(measure) ){
#       yuima.warn("measure type isn't matched with jump term.")
#       return(NULL)
#       # }else
#       #       if(length(jump.coeff)!=1){
#       #        yuima.warn("multi dimentional jump term is not supported yet.")
#       #
#       #         return(NULL)
#       #     }
#
#     }
    if("CP" %in% measure.type){
      condCP <- (measure.type%in%"CP")
      numbCP<-sum(condCP)
      h <- 0
      for(i in c(1:length(measure.type))){
        if(length(measure.type[condCP[i]])!=0){
          if(!is.na(intensity[i-h])){
            tmp.measure$intensity[(i-h)] <- parse(text = intensity[i-h])
          }
          measure.par <- c(measure.par,all.vars(tmp.measure$intensity[(i-h)]))
        }else{
          h<-h+1
        }
      }
    }
    tmp.measure$df$expr <- parse(text=df[[1]])
    measure.par <- c(measure.par, all.vars(tmp.measure$df$expr))
    tmp <- regexpr("\\(", df[[1]])[1]
    measurefunc <- substring(df[[1]], 1, tmp-1)
    if(existsFunction(measurefunc)){
      tmp.measure$df$func[[1]] <- eval(parse(text=measurefunc))
    }else{
      yuima.stop("density function for jump must be specified")
    }


#     if("CP"%in%measure.type){
#       condCP <- (measure.type%in%"CP")
#       numbCP<-sum(condCP)
#       h <- 0
#       for(i in c(1:length(measure.type))){
#         if(length(measure[condCP[i]])!=0){
#           tmp.measure$df$expr[i] <- parse(text=measure$df[[i]])
#           tmp.measure$intensity[(i-h)] <- parse(text = measure$intensity[i-h])
#
#           tmp <- regexpr("\\(", measure$df[[i]])[1]
#           measurefunc <- substring(measure$df[[i]], 1, tmp-1)
#           if(!is.na(match(measurefunc, codelist))){
#             yuima.warn(paste("distribution function", measurefunc, "should be defined as type code."))
#             return(NULL)
#           }
#           tmp.measure$df$func[[i]] <- eval(parse(text=measurefunc))
#         }
#         h<-h+1
#       }
#
#
       MEASURE$df$func<-tmp.measure$df$func
       MEASURE$df$expr<-tmp.measure$df$expr
       if("CP" %in% measure.type){
         MEASURE$intensity<-tmp.measure$intensity
       }
       measure.par<-unique(measure.par)

       n.eqn3 <- dim(jump.coeff)[1]
       n.jump <- dim(jump.coeff)[2]
  }
#
#
#       ##measure.par$intensity <- unique(all.vars(MEASURE$intensity))
#       ##::end check df name ####################
#       ##::end CP
#
#     }
#     if("code"%in%measure.type){ ##::code
#       if(!is.list(measure)){
#         measure <- list(df=measure)
#       }
# #      else{
# #         if(length(measure[[1]])!=1){
# #           yuima.warn("multi dimentional jump term is considered.")
# #         }
#         ##::naming measure list #############
# #         if(is.null(names(measure)) || names(measure)=="df"){
# #           names(measure) <- "df"
# #         }else{
# #           yuima.warn("name of measure is incorrect.")
# #           return(NULL)
# #         }
#         ##::end naming measure list #############
# #      }
#
#       condcode <- (measure.type%in%"code")
#       numbcode<-sum(condcode)
#       h <- 0
#       for(i in c(1:length(measure.type))){
#         if(condcode[i]){
#           tmp.measure$df$expr[i]<-parse(text=measure$df[[i]])
#
#           tmp <- regexpr("\\(", measure$df[[i]])[1]
#           measurefunc <- substring(measure$df[[i]], 1, tmp-1)
#           if(!is.na(match(measurefunc, CPlist))){
#             yuima.warn(paste("\ndistribution function", measurefunc, "should be defined as type CP."))
#             return(NULL)
#           }else if(is.na(match(measurefunc, codelist))){
#             warning(paste("\ndistribution function", measurefunc, "is not officialy supported as type code.\n"))
#           }
#           MEASURE$df$func[[i]] <- eval(parse(text=measurefunc))
#           MEASURE$df$expr[i] <- tmp.measure$df$expr[i]
#
#
# #           measure.par <- unique(all.vars(MEASURE$df$expr))
# #
# #           tmp.measure$func
#         }
#       }
#
#       if( "CP" %in% measure.type){
#         measure.par <- unique( c( all.vars(MEASURE$intensity), all.vars(MEASURE$df$expr) ) )
#       }else if("code" %in% measure.type){
#         measure.par <- unique( c( all.vars(MEASURE$intensity), all.vars(MEASURE$df$expr) ) )
#       }
#
#
# #       ##::check df name ####################
# #       tmp <- regexpr("\\(", measure$df)[1]
# #       measurefunc <- substring(measure$df, 1, tmp-1)
# #       if(!is.na(match(measurefunc, CPlist))){
# #         yuima.warn(paste("\ndistribution function", measurefunc, "should be defined as type CP."))
# #         return(NULL)
# #       }else if(is.na(match(measurefunc, codelist))){
# #         warning(paste("\ndistribution function", measurefunc, "is not officialy supported as type code.\n"))
# #       }
# #       ##MEASURE$df$func <- eval(parse(text=measurefunc))
# #       MEASURE$df$expr <- parse(text=measure$df)
# #
# #       measure.par <- unique(all.vars(MEASURE$df$expr))
# #       ##::end check df name ####################
# #       ##::end code
#     }else if(measure.type=="density"){ ##::density
#       if(length(measure)!=1){
#         yuima.warn(paste("length of measure must be one on type", measure.type, "."))
#         return(NULL)
#       }
#
#       if(!is.list(measure)){
#         measure <- list(df=measure)
#       }else{
#         if(length(measure[[1]])!=1){
#           yuima.warn("multi dimentional jump term is not supported yet.")
#           return(NULL)
#         }
#
#         ##::naming measure list #############
#         if(is.null(names(measure))){
#           names(measure) <- "df"
#         }else if(names(measure)!="density" && names(measure)!="df"){
#           yuima.warn("name of measure is incorrect.")
#           return(NULL)
#         }
#         ##::end naming measure list #############
#       }
#
#       ##::check df name ####################
#       tmp <- regexpr("\\(", measure[[names(measure)]])[1]
#       measurefunc <- substring(measure[[names(measure)]], 1, tmp-1)
#       if(!is.na(match(measurefunc, CPlist))){
#         yuima.warn(paste("distribution function", measurefunc, "should be defined as type CP."))
#         return(NULL)
#       }else if(!is.na(match(measurefunc, codelist))){
#         yuima.warn(paste("distribution function", measurefunc, "should be defined as type code."))
#         return(NULL)
#       }
#       MEASURE[[names(measure)]]$func <- eval(parse(text=measurefunc))
#       MEASURE[[names(measure)]]$expr <- parse(text=measure[[names(measure)]])
#
#       measure.par <- unique(all.vars(MEASURE[[names(measure)]]$expr))
#       ##::end check df name ####################
#       ##::end density
#     }else{ ##::else
#       yuima.warn(paste("measure type", measure.type, "isn't supported."))
#       return(NULL)
#     }
#     n.eqn3 <- dim(jump.coeff)[1]
#     n.jump <- length(measure.type)
#   }

  ##::end measure and jump term #####################################

  ##:: check for errors and reform values
  if(any(time.variable %in% state.variable)){
    yuima.warn("time and state(s) variable must be different.")
    return(NULL)
  }
  if(is.null(dim(drift))){ # this is a vector
    n.eqn1 <- length(drift)
    n.drf <- 1
  }else{ # it is a matrix
    n.eqn1 <- dim(drift)[1]
    n.drf <- dim(drift)[2]
  }

  if(is.null(dim(diffusion))){ # this is a vector
    n.eqn2 <- length(diffusion)
    n.noise <- 1
  }else{ # it is a matrix
    n.eqn2 <- dim(diffusion)[1]
    n.noise <- dim(diffusion)[2]
  }

  if(is.null(diffusion)){
    diffusion <- rep("0", n.eqn1)
    n.eqn2 <- n.eqn1
    n.noise <- 1
  }

  ## TBC
  n.eqn3 <- n.eqn1

#   if(!length(measure)){
#     n.eqn3 <- n.eqn1
#   }

  if(n.eqn1 != n.eqn2 || n.eqn1 != n.eqn3){
    yuima.warn("Malformed model, number of equations in the drift and diffusion do not match.")
    return(NULL)
  }
  n.eqn <- n.eqn1

  if(is.null(xinit)){
    # xinit <- numeric(n.eqn)
    xinit <- character(n.eqn)
  }else if(length(xinit) != n.eqn){
    if(length(xinit)==1){
      xinit <- rep(xinit, n.eqn)
    }else{
      yuima.warn("Dimension of xinit variables missmatch.")
      return(NULL)
    }
  }

  if(missing(solve.variable)){
    yuima.warn("Solution variable (lhs) not specified. Trying to use state variables.")
    solve.variable <- state.variable
  }
  if(n.eqn != length(solve.variable)){
    yuima.warn("Malformed model, number of solution variables (lhs) do no match number of equations (rhs).")
    return(NULL)
  }

  loc.drift <- matrix(drift, n.eqn, n.drf)
  loc.diffusion <- matrix(diffusion, n.eqn, n.noise)
  # Modification starting point 6/11
  loc.xinit<-matrix(xinit,n.eqn,n.drf)

  ##:: allocate vectors
  DRIFT <- vector(n.eqn, mode="expression")
  DIFFUSION <- vector(n.eqn, mode="list")
  # Modification starting point 6/11
  XINIT<-vector(n.eqn, mode = "expression")

  ##:: function to make expression from drift characters
  pre.proc <- function(x){
    for(i in 1:length(x)){
      if(length(parse(text=x[i]))==0){
        x[i] <- "0"
      }
    }
    parse(text=paste(sprintf("(%s)", x), collapse="+"))
  }
  ##22/11:: function to simplify expression in drift, diffusion, jump and xinit characters
  yuima.Simplifyobj<-function(x){
    dummy<-yuima.Simplify(x, yuima.env=yuimaENV)
    dummy1<-yuima.Simplify(dummy, yuima.env=yuimaENV)
    dummy2<-as.character(dummy1)
    res<-parse(text=paste0("(",dummy2,")",collapse=NULL))
    return(res)
  }


  ##:: make expressions of drifts and diffusions and jump
  for(i in 1:n.eqn){
    DRIFT[i] <- pre.proc(loc.drift[i,])
    # 22/11 Simplify expressions
    DRIFT[i] <- yuima.Simplifyobj(DRIFT[i])
    # Modification starting point 6/11
    XINIT[i]<-pre.proc(loc.xinit[i, ])
    XINIT[i]<- yuima.Simplifyobj(XINIT[i])
    for(j in 1:n.noise){
      expr <- parse(text=loc.diffusion[i,j])
      if(length(expr)==0){
        expr <- expression(0)  # expr must have something
      }
      #       DIFFUSION[[i]][j] <- expr
      #22/11
      DIFFUSION[[i]][j] <- yuima.Simplifyobj(expr)
    }
    #22/11

    #if (length(JUMP)>0){
    #    JUMP[i] <- parse(text=jump.coeff[i])
    #    JUMP[i] <- yuima.Simplifyobj(JUMP[i])
    #}

  }



  #print(length(jump.coeff))
  #if (length(jump.coeff)==0){
  #    JUMP <- list(parse(text=jump.coeff))
  #}else{
  #    #    JUMP <- vector(n.eqn, mode="expression")
  #    JUMP <- vector(n.eqn, mode="list")
  #}

  if(length(jump.coeff)==0){
    JUMP <- list()
  } else {
    if(length(jump.coeff)==1 & !is.matrix(jump.coeff)){ # is a scalar
      expr <- parse(text=jump.coeff)
      if(length(expr)==0){
        expr <- expression(0)  # expr must have something
      }
      JUMP <- list(yuima.Simplifyobj(expr))
    } else { # must be matrix, n.col = dimension of Levy noise
      jump.coeff <- as.matrix(jump.coeff)
      c.j <- ncol(jump.coeff)
      r.j <- nrow(jump.coeff)
      #print(c.j)
      #print(r.j)
      #print(jump.coeff)
      JUMP <- vector(r.j, mode="list")
      for(i in 1:r.j){
        for(j in 1:c.j){
          #cat(sprintf("\ni=%d,j=%d\n",i,j))
          expr <- parse(text=jump.coeff[i,j])
          if(length(expr)==0){
            expr <- expression(0)  # expr must have something
          }
          JUMP[[i]][j] <- yuima.Simplifyobj(expr)
        }
      }
    }
  }
  #print(str(JUMP))

  #

  ##:: get parameters in drift expression
  drift.par <- unique(all.vars(DRIFT))
  # Modification starting point 6/11
  xinit.par <- unique(all.vars(XINIT))

  drift.idx <- as.numeric(na.omit(match(c(state.variable, time.variable, jump.variable, solve.variable), drift.par)))
  if(length(drift.idx)>0){
    drift.par <- drift.par[-drift.idx]
  }

  ##:: get parameters in diffusion expression
  diff.par <- unique(unlist(lapply(DIFFUSION, all.vars)))
  diff.idx <- as.numeric(na.omit(match(c(state.variable, time.variable, jump.variable, solve.variable), diff.par)))
  if(length(diff.idx)>0){
    diff.par <- diff.par[-diff.idx]
  }

  ##:: get parameters in jump expression
  J.flag <- FALSE
  #  jump.par <- unique(all.vars(JUMP))
  jump.par <- unlist(lapply(JUMP,all.vars))
  if(is.null(jump.par))
    jump.par <- character()
  if(length(na.omit(match(jump.par, jump.variable)))){
    J.flag <- TRUE
  }
  jump.idx <- as.numeric(na.omit(match(c(state.variable, time.variable, jump.variable, solve.variable), jump.par)))
  if(length(jump.idx)>0){
    jump.par <- jump.par[-jump.idx]
  }

  ##:: get parameters in measure expression
  measure.idx <- as.numeric(na.omit(match(c(state.variable, time.variable, jump.variable, solve.variable), measure.par)))
  if(length(measure.idx)>0){
    measure.par <- measure.par[-measure.idx]
  }

  ##:: order parameters for 'yuima.pars'
  ##id1 <- which(diff.par %in% drift.par)
  ##id2 <- which(drift.par %in% diff.par)
  ##common <- unique(c(diff.par[id1], drift.par[id2]))
  common <- c(drift.par, diff.par)
  common <- common[duplicated(common)]

  common1<-common
  # modification 06/11 common1 contains only
  # parameters that appear in both drift and diffusion terms.

  # Modification 06/11 common contains only parameters that appear
  # in drift, diff, Jump and xinit
  if (length(xinit)) {
    common <- c(common, xinit.par)
    common <- common[duplicated(common)]
    common <- c(common, xinit.par)
    common <- common[duplicated(common)]
  }


  if(length(measure.type)){
    common <- c(common, jump.par)
    common <- common[duplicated(common)]
    common <- c(common, measure.par)
    common <- common[duplicated(common)]
  }
  #   all.par <- unique(c(drift.par, diff.par, jump.par, measure.par))
  all.par <- unique(c(drift.par, diff.par, jump.par, measure.par, xinit.par))

  ##:: instanciate class
  tmppar <- new("model.parameter",
                all= all.par,
                #                 common= common,
                common= common1,
                diffusion= diff.par,
                drift= drift.par,
                jump= jump.par,
                measure= measure.par,
                xinit=xinit.par)

  tmp <- new("yuima.multimodel",
             drift= DRIFT,
             diffusion= DIFFUSION,
             hurst=as.numeric(hurst),
             jump.coeff=JUMP,
             measure= MEASURE,
             measure.type= measure.type,
             parameter= tmppar,
             state.variable= state.variable,
             jump.variable= jump.variable,
             time.variable= time.variable,
             noise.number= n.noise,
             equation.number= n.eqn,
             dimension= c(
               length(tmppar@all),
               length(tmppar@common),
               length(tmppar@diffusion),
               length(tmppar@drift),
               length(tmppar@jump),
               length(tmppar@measure)
             ),
             solve.variable= solve.variable,
             xinit= XINIT,
             J.flag <- J.flag)
  return(tmp)
}
