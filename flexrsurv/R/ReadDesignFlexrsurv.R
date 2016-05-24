ReadDesignFlexrsurv<-function(modframe,
                              Terms,
                              data=parent.frame(),
                              rate=NULL,
                              Spline_t0,
                              na.action=options("na.action"),
                                contrasts=NULL){
  # extract X0, X and Z componenets of a flexrsurv formula
  # X0 : linear and nonlinear designe matrix
  # X  : NPH design matrix
  # Z  : NPHNLL design matrix for multiplicative NLLNPH terms 

  # dans cette version, on ne gère pas les variables de type facteur.
  # ainsi, un variable LIN peut être une indicatrice
  
  mt <- attr(modframe, "terms")
  intercept <- attr(mt, "intercept")

  AllVariables <- as.character(attr(Terms,"variables"))[-1]
  # get the name of the ordered diminsion of de design matrix
  # replace multiplicative by additive to get the correct number of df for each NPHNLL effect

  tempf <- as.character(Terms)
  tempf <- as.formula(paste(tempf[2], tempf[1], gsub("multiplicative", "additive", tempf[3], fixed=TRUE), collapse=""))
  XZ<-model.matrix(Terms, data)[, drop = FALSE]
  # names of the coef in the formula

  names_coef_XZ <- dimnames(XZ)[[2]][-1]

  
  # remove the first T basis of each multiplicative NPHNLL effects
  # 
  ismultiplicative <- grep("multiplicative", names_coef_XZ)
  isfirstTbasis <- grep(":TtTtTTtTtT 1", names_coef_XZ)
  rmvname <- ismultiplicative[ismultiplicative %in% isfirstTbasis]
  if(length(rmvname)){
    names_coef_XZ <- names_coef_XZ[-rmvname]
  }
  
    NPHNLLVars <- attr(Terms, "specials")$NPHNLL
    # get additif and multiplicatif NPHNLL terms
  NPHNLLVarsAdd <- NPHNLLVarsMult <- NULL
      if( !is.null( NPHNLLVars )) {
        indxvarNPHNLL <- NPHNLLVars
        nvarsNPHNLL <- length(indxvarNPHNLL)        
        var_list <- NULL
        for( i in 1:nvarsNPHNLL){
          fun <- mget("NPHNLL",
                      mode = "function",
                      envir = parent.frame(), inherits=TRUE,
                      ifnotfound=list(NULL))[[1]]
          thecall <- match.call(fun,  attr(Terms,"variables")[[indxvarNPHNLL[i]+1]])
          if( thecall[["model"]] == "additive" ){
            NPHNLLVarsAdd <- c( NPHNLLVarsAdd ,  NPHNLLVars[[i]]) 
          }
          else {
            NPHNLLVarsMult <- c( NPHNLLVarsMult ,  NPHNLLVars[[i]]) 
          }
        }
        NPHNLLTerms <-  variable2term(NPHNLLVars, Terms)
        NPHNLLTermsAdd <-  variable2term(NPHNLLVarsAdd, Terms)
        NPHNLLTermsMult <-  variable2term(NPHNLLVarsMult, Terms)
      } else {
        NPHNLLVarsMult <- NULL
        NPHNLLVarsAdd <- NULL
        NPHNLLTerms <- NULL
        NPHNLLTermsAdd <-  NULL
        NPHNLLTermsMult <-  NULL
      }
  # merge additif NPHNLL termes with NLL and NPH terms

  NLLVars <- attr(Terms, "specials")$NLL
  if( !is.null( NLLVars )) {
    # get the indexes of the terms
    NLLTerms <- variable2term(NLLVars, Terms)
  } else {
    NLLTerms <-  NULL
  }
  
  NPHVars <- attr(Terms, "specials")$NPH
  if( !is.null( NPHVars )) {
    # get the indexes of the terms
    NPHTerms <-  variable2term(NPHVars, Terms)
  } else {
    NPHTerms <-  NULL
  }

  if(length(c(NLLTerms, NPHTerms, NPHNLLTerms))>0){
    LinVars <- (1:(length(AllVariables)))[-c(1,NLLVars, NPHVars, NPHNLLVars,attr(Terms, "offset"))  ]
    LinTerms <- (1:(length(attr(Terms, "term.labels"))))[-c(NLLTerms, NPHTerms, NPHNLLTerms)]
  } else if(length(attr(Terms, "term.labels")) > 0){
    LinVars <- (1:(length(AllVariables)))[-c(1,attr(Terms, "offset"))]
    LinTerms <- (1:(length(attr(Terms, "term.labels"))))
  } else {
    # model is Surv() ~ Intercept
    LinVars <- NULL
    LinTerms <- NULL
  }

  # get the indexes of the terms

  
  is.LIN <- is.NLL <- is.NPH <- is.NPHNLL <- is.NPHNLLMult <- is.NPHNLLAdd <- rep(FALSE, length( attr(Terms, "term.labels")))
  names(is.LIN) <- names(is.NLL) <- names(is.NPH) <- names(is.NPHNLL) <- names(is.NPHNLLMult) <- names(is.NPHNLLAdd) <- attr(Terms, "term.labels")

  if(length(LinTerms)>0) {
    is.LIN[LinTerms] <- TRUE
  }
  if(length(NLLTerms)>0) {
    is.NLL[NLLTerms] <- TRUE
  }
  if(length(NPHTerms)>0) {
    is.NPH[NPHTerms] <- TRUE
  }
  if(length(NPHNLLTerms)>0) {
    is.NPHNLL[NPHNLLTerms] <- TRUE
  }
  if(length(NPHNLLTermsMult)>0) {
    is.NPHNLLMult[NPHNLLTermsMult] <- TRUE
  }
  if(length(NPHNLLTermsAdd)>0) {
    is.NPHNLLAdd[NPHNLLTermsAdd] <- TRUE
  }
    

  is.X0 <- is.LIN + is.NLL + is.NPHNLLAdd
  is.X <- is.NPH   + is.NPHNLLAdd
  is.Z <- is.NPHNLLMult


  
  # index of the term in the 3 parts of the variable matrices 
  num.TermsX0 <- cumsum(is.X0)
  num.TermsX <- cumsum(is.X)
  num.TermsZ <- cumsum(is.Z)
    
  
  # in TermsNLL : NLLVars and NPHNLLVarsAdd
  # keep.intercept=TRUE but Intercept will be removed later
  if(length(c( NLLTerms, NPHNLLTermsAdd))){
    # in keep.terms, the second arg is the vector of terms that are kept 
    TermsNLL <- keep.terms(Terms, sort(c(NLLTerms, NPHNLLTermsAdd)), keep.intercept = TRUE, keep.response=FALSE)
        # replace NPHNLL by NLL
    for(k in  attr(TermsNLL,"specials")$NPHNLL ){
      attr(TermsNLL,"term.labels")[k] <- paste("NLL", substring(attr(TermsNLL,"term.labels")[k], first=7)) 
    }
    # recast the termsobj
    TermsNLL <- keep.terms(TermsNLL, NULL, keep.intercept = TRUE, keep.response=FALSE)
  } else {
    TermsNLL <- NULL
  }
  
  if(length(LinTerms)){
    # in keep.terms, the second arg is the vector of terms that are kept 
    TermsLin <- keep.terms(Terms, LinTerms, keep.intercept = TRUE, keep.response=FALSE)
  } else {
    TermsLin <- NULL
  }

  if(length(c( LinTerms, NLLTerms, NPHNLLTermsAdd))){
    # in keep.terms, the second arg is the vector of terms that are kept 
    TermsX0i <- keep.terms(Terms, sort(c(LinTerms, NLLTerms, NPHNLLTermsAdd)), keep.intercept = TRUE)
        # replace NPHNLL by NLL
    for(k in  attr(TermsX0i,"specials")$NPHNLL ){
      oneterm <- terms(as.formula(paste("~",attr(TermsX0i, "term.labels")[k])))
      thecall <-  match.call(NPHNLL, attr(oneterm,"variables")[[2]], expand.dots = TRUE)
      # rename NLL instead og NPHNLL
      thecall[[1]] <- as.name("NLL")
      # force intercept = FALSE
      thecall[["Intercept"]] <- FALSE
      thecall <-  match.call(NLL, thecall, expand.dots = FALSE)
      if( length(thecall[["..."]])){
        thecall[["..."]] <- NULL
      }
      attr(TermsX0i,"term.labels")[k] <- deparse(thecall, width.cutoff = 500L)
    }
# with intercept
    TermsX0i <- terms(reformulate(attr(TermsX0i,"term.labels"), response = NULL, intercept = TRUE),
                      specials=c("NPH","NLL", "NPHNLL") )
# without intercept
    TermsX0 <- terms(reformulate(attr(TermsX0i,"term.labels"), response = NULL, intercept = FALSE),
                      specials=c("NPH","NLL", "NPHNLL") )
  } else {
    TermsX0 <- NULL
  }
  
  # in X, , NPHTerms and NPHNLLTermsAdd
  # keep.intercept=FALSE because all X terms are "continuous"

  if(length(c( NPHTerms, NPHNLLTermsAdd))){
    TermsX <- Terms
# replace NPHNLL by NPH of additives NPHNLL terms
    for(k in NPHNLLTermsAdd){
      oneterm <- terms(as.formula(paste("~",attr(TermsX,"term.labels")[k])))
      thecall <-  match.call(NPHNLL, attr(oneterm,"variables")[[2]], expand.dots = TRUE)
      # rename NPH instead og NPHNLL
      thecall[[1]] <- as.name("NPH")
      # force intercept.T = FALSE
      thecall[["Intercept.t"]] <- FALSE
      thecall <-  match.call(NPH, thecall, expand.dots = FALSE)
      if( length(thecall[["..."]])){
        thecall[["..."]] <- NULL
      }
      attr(TermsX,"term.labels")[k] <- deparse(thecall, width.cutoff = 500L)
    }
    TermsX <- keep.terms(TermsX, sort(c( NPHTerms, NPHNLLTermsAdd)), keep.intercept = FALSE)
  } else {
    TermsX <- NULL
  }
  TermsNPH <- TermsX 
  
  # in Z  NPHNLLTermsMult
  # keep.intercept=FALSE because all Z terms are "continuous"
  if(length(NPHNLLTermsMult)){
    TermsZ <- keep.terms(Terms, NPHNLLTermsMult, keep.intercept = FALSE)
    NamesZ <- dimnames(model.matrix(TermsZ, modframe)[, drop = FALSE])[[2]]
  } else {
    TermsZ <- NULL
    NamesZ <- NULL
  }
  TermsNPHNLLMult <- TermsZ 

   
# vector of variable names involed in the terms in the original formula
  # in LINVars, there are no interaction!

  NamesLinVars <- all_LIN_vars(Terms)
  
  NamesNLLVars <- all_specials_vars( TermsNLL,
                                    specials=c("NLL"),
                                    unique = TRUE,
                                    order="formula")
  NamesNPHVars<- all_specials_vars( Terms,
                                   specials=c("NPH", "NPHNLL"),
                                   unique = TRUE,
                                   order="formula")
  
  NamesNPHNLLVars<- all_specials_vars( Terms,
                                      specials="NPHNLL",
                                      unique = TRUE,
                                      order="formula")
  
  NPHNLLVarsAdd<- all_specials_vars( TermsNLL,
                                    specials="NPHNLL",
                                    unique = TRUE,
                                    order="formula")
  
  NPHNLLVarsMult<- all_specials_vars( TermsZ,
                                     specials="NPHNLL",
                                     unique = TRUE,
                                     order="formula")
  


# in X0, LIN vars + basis of NLIN vars 
  if (is.null(TermsX0)){
    X0 <- NULL
    X0Vars <- NULL
    df.TermsX0 <- NULL
    paramX0.fin <- NULL
    paramX0.deb <- NULL
  }
  else {
    # matrix of pure linear and pure nonlinear effects
    # splines are expanded
    X0i <- model.matrix(TermsX0i,data=data, contrasts.args=contrasts)
    col.intercept <- attr(X0i, "assign")==0 

    X0 <- X0i[, dimnames(X0i)[[2]]!="(Intercept)", drop=FALSE]
    X0Vars <- all.vars(TermsX0)
# index of the term in X0term corresponding to colons in X0 : names(coef)<-attr(TermsXO, "term.labels")[assignX0]
    assignX0 <-  attr(X0i, "assign")[-1]
# number of degree of freedom of each effect
    df.TermsX0 <- table(attr(X0i, "assign")[!col.intercept])
    paramX0.fin <- cumsum(df.TermsX0) 
    paramX0.deb <- c(1, paramX0.fin[-length(paramX0.fin)]+1)

  }


# in X NPH vars
  if (is.null(TermsX)){
    X <- NULL
    XVars <- NULL
    Spline_XT <- NULL
    df.TermsX <- NULL
    paramX.fin <- NULL
    paramX.deb <- NULL
  }
  else {
    # matrix of true NPH effets (X matrix in the fitter)
    # raw variables, splines bases are not expanded
    XVars <- all_specials_vars(TermsX,
                               specials=c("NPH", "NPHNLL"),
                               unique=TRUE, order="formula")
    Spline_XT <- get_TimeSplinebasis(TermsX,
                                     data=data,
                                     specials=c("NPH", "NPHNLL"),
                                     all.vars.func=all_specials_vars, 
                                     unique=TRUE, order="formula")
    inpX <- parse(text = paste("list(", paste(XVars, collapse = ","), ")"))
    X <- as.data.frame(eval(inpX, envir=data))
    Xi <- model.matrix(TermsX,data=data, contrasts.args=contrasts)
    
    df.TermsX <- table(attr(Xi, "assign"))
    paramX.fin <- cumsum(df.TermsX) 
    paramX.deb <- c(1, paramX.fin[-length(paramX.fin)]+1)

    if(length(XVars)==1){
      X <- as.matrix(X)
    }
    dimnames(X)[[2]] <- XVars 

    
#    X <- get_special_vars(TermsX,data=data,
#                          specials=c("NPH", "NPHNLL"),
#                          all.vars.func=all_specials_vars, 
#                          unique=TRUE, order="formula")[, , drop = FALSE]
  }

# in Z raw multiplicative NPHNLL vars 
  if (is.null(TermsZ)){
    Z <- NULL
    ZVars <- NULL
    Spline_Z <- NULL
    Spline_ZT <- NULL
    df.TermsZ <- NULL
    df.TermsZT <- NULL
    paramZ.fin <- NULL
    paramZ.deb <- NULL

    paramZT.fin <- NULL
    paramZT.deb <- NULL

  }
  else {
    ZVars <- all_specials_vars(TermsZ,
                               specials=c("NPH", "NPHNLL"),
                               unique=TRUE, order="formula")

    # matrix of multiplicative NPHNLL variables (Z matrix in the fitter)
    # raw variables, splines bases are not expanded
    inp <- parse(text = paste("list(", paste(ZVars, collapse = ","), ")"))
    Z <- as.data.frame(eval(inp, envir=data))
    Zi <- model.matrix(TermsZ,data=data, contrasts.args=contrasts)
    # with multiplicative NPHNLL, the model matrix of each variable ase dfZ + dfT + 1 degree of freedom (Intercept.T=TRUE)

#
    
    
    df.TermsZZT <- table(attr(Zi, "assign")) -1


    
    
    
#Z <- get_all_vars(TermsZ,data=environment(modframe)) 
#    Z <- get_specials_vars(TermsZ,
#                           data=data,
#                           specials="NPHNLL",
#                           all.vars.func=all_specials_vars, 
#                           unique=TRUE, order="formula")[, , drop = FALSE]

# get Spline parameters of each Z
    Spline_Z <- get_Splinebasis(TermsZ,
                                 data=data,
                                 specials="NPHNLL",
                                 all.vars.func=all_specials_vars, 
                                 unique=TRUE, order="formula")
    Spline_ZT <- get_TimeSplinebasis(TermsZ,
                                     data=data,
                                     specials="NPHNLL",
                                     all.vars.func=all_specials_vars, 
                                     unique=TRUE, order="formula")
    names(Z) <- ZVars

    
    # df of all the beta(t)
    df.TermsZT <- rep(0, length(ZVars))
    for( s in 1:length(Spline_ZT)){ 
      df.TermsZT[s] <- getNBases(Spline_ZT[[s]])-1  # because first coef of Tbasis is constraints to 1 in NPHNLL()
    }

    df.TermsZ <-     df.TermsZZT - df.TermsZT 
    
    
    paramZ.fin <- cumsum(df.TermsZ) 
    paramZ.deb <- c(1, paramZ.fin[-length(paramZ.fin)]+1)

    paramZT.fin <- cumsum(df.TermsZT) 
    paramZT.deb <- c(1, paramZT.fin[-length(paramZT.fin)]+1)
    
  }

  df.T0 <- getNBases(Spline_t0)
# nb of coefficients
  ncoef   <- df.T0 + sum(df.TermsX0) + sum(df.TermsX) + sum(df.TermsZ) + sum(df.TermsZT)
# nb of terms in the formula
  nterms <- length(attr(Terms,"term.labels"))

  ngamma0 <- df.T0
  nalpha0 <- sum(df.TermsX0)
  nbeta0 <- sum(df.TermsX)
  nalpha <- sum(df.TermsZ)
  nbeta <- sum(df.TermsZT)
#  nterms0 <- length(df.TermsX0) + length(df.TermsX) + length(df.TermsZ) 

  
# set the transition vectors between the estimated parameters by FlexRsurv.LL and the terms of the formula 
# assign is in the compact R form, which is a vector (0, 0, 0, 0, 1, 2, 3, 3, 3); that can be
#   read as "the 4 first column of the X matrix (baseline hazard) goes with none of
#   the terms', 'the 5th column goes with term 1', etc.
# assignList is the Splus style assign component, which is a list
#      $(Intercept)     1 2 3 4
#      $age             5
#      $sex             6
#      $factor(other) 7 8 9   
#
#  see atrrassign() function in package survival
  
# df by terms
  param2term <- rep(0, ncoef)
  assign <- rep(0, ncoef)
  df.Terms <- rep(0, nterms)
  coef.deb <- rep(0, nterms)
  coef.fin <- rep(0, nterms)
# baseline hazard
  param2term[1:df.T0] <- 1:df.T0
  assign[1:df.T0] <- 0
  assignList <- list(intercept=1:df.T0)
  names(assignList)[1] <- "(Intercept)"
# list to pass from vector of coefficients to the list used by flexrsurv.**.fit
# param$xx = coef[term2param$xx]
  term2param <- list(gamma0=1:df.T0, alpha0=NULL, beta0=NULL, alpha=NULL, beta=NULL)

  
  count <- df.T0


  if(nterms > 0){
    for(it in 1: nterms){
      count0 <- count
      
      if( is.X0[it] ){
      # term it is LIN or NLIN
        df.Terms[it] <- df.TermsX0[num.TermsX0[it]]
        param2term[count0 + (1:df.Terms[it])] <- ngamma0 + paramX0.deb[num.TermsX0[it]]:paramX0.fin[num.TermsX0[it]]
        term2param$alpha0=c(term2param$alpha0, count0 + (1:df.Terms[it]))
        count0=count+df.Terms[it] 
      }
      if( is.X[it] ){
      # term it is NPH
        df.Terms[it] <- df.Terms[it] + df.TermsX[num.TermsX[it]]
        param2term[count0 + (1:df.TermsX[num.TermsX[it]])] <- ngamma0 + nalpha0 + paramX.deb[num.TermsX[it]]:paramX.fin[num.TermsX[it]]
        term2param$beta0=c(term2param$beta0, count0 + (1:df.TermsX[num.TermsX[it]]))
      }
      if( is.Z[it] ){
      # term it is multiplicative NPHNLL
        df.Terms[it]   <- df.TermsZ[num.TermsZ[it]] + df.TermsZT[num.TermsZ[it]]
        ndlZ  <- df.TermsZ[num.TermsZ[it]]
        ndlZT <- df.TermsZT[num.TermsZ[it]]
# alpha(X) term
        param2term[count+(1:ndlZ)] <- ngamma0 + nalpha0 + nbeta0 + paramZ.deb[num.TermsZ[it]]:paramZ.fin[num.TermsZ[it]]
        term2param$alpha=c(term2param$alpha, count+(1:ndlZ))
# beta(T) term
        param2term[count+ndlZ+(1:ndlZT)] <- ngamma0 + nalpha0 + nbeta0 + nalpha + paramZT.deb[num.TermsZ[it]]:paramZT.fin[num.TermsZ[it]]
        term2param$beta=c(term2param$beta, count+ndlZ+(1:ndlZT))
      }
      coef.deb[it] <- count + 1
      coef.fin[it] <- count + df.Terms[it] 
      assign[coef.deb[it]:coef.deb[it]] <- it
      assignList <- c(assignList, list(coef.deb[it]:coef.deb[it]))
      count <- count+df.Terms[it]
    }
  } 

  
return(list(X0=X0,
            X=X,  
            Z=Z,
            rate = rate,
            TermsX0=TermsX0,
            TermsX=TermsX,  
            TermsZ=TermsZ,
            LinVars=LinVars,
            NLLVars=NLLVars,
            NPHVars=NPHVars,
            NPHNLLVars=NPHNLLVars,
            NPHNLLVarsAdd=NPHNLLVarsAdd,
            NPHNLLVarsMult=NPHNLLVarsMult,

            X0Vars=X0Vars,
            XVars=XVars,
            Spline_XT=Spline_XT,
            ZVars=ZVars,
            Spline_Z=Spline_Z,
            Spline_ZT=Spline_ZT,
                                             # these are the number of DF of each term
            df.T0=df.T0,
            df.TermsX0=df.TermsX0,
            df.TermsX=df.TermsX,
            df.TermsZ=df.TermsZ,
            df.TermsZT=df.TermsZT,

            paramX0.deb=paramX0.deb,
            paramX.deb=paramX.deb,
            paramZ.deb=paramZ.deb,
            paramZT.deb=paramZT.deb,
            
            paramX0.fin=paramX0.fin,
            paramX.fin=paramX.fin,
            paramZ.fin=paramZ.fin,
            paramZT.fin=paramZT.fin,
            
# vector of length nterms =1 if Terms[i] is in X0 reps X, restp Z
            is.X0=is.X0,
            is.X=is.X,
            is.Z=is.Z,
            
                                             # 
            num.TermsX0=num.TermsX0,        # index of the LL/NLL effects in the formula (vectors of lenghth lenght(TermsX0)
            num.TermsX=num.TermsX,          #              NPH    effects in the formula (vectors of lenghth lenght(TermsX)
            num.TermsZ=num.TermsZ,          #              NPHNLL effects in the formula (vectors of lenghth lenght(TermsZ)
                                             # the following vectors are of lenghth nb_of_effects/terms in the formula
                                             # for example if terms[[i]] is NPH then num.TermsX[i]=the order of that effect in TermsX, else it is 0 
#            num.alltermsX0=num.alltermsX0,  # index of the LL/NLL effects in the LL/NLL effects 
#            num.alltermsX=num.alltermsX,    #              NPH    effects in the NPH effects    (vectors of lenghth nb_of_effects/terms in the formula)
#            num.alltermsZ=num.alltermsZ,    #              NPHNLL effects in the NPHNLL effects (vectors of lenghth nb_of_effects/terms in the formula)

            ncoef=ncoef,
            nterms=nterms,
            df.Terms=df.Terms,
            coef.deb = coef.deb,
            coef.fin = coef.fin,
            assign = assign,
            assignList = assignList,
            param2coef = param2term,
            coef2param = term2param,
            names_coef_XZ = names_coef_XZ

            ))
}



