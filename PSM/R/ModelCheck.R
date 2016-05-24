ModelCheck <- function(Model , Data , Par, DataHasY=TRUE) {
  # Check the dimensions of the Model and Data
  #
  # Data is data for a single subject
  # Data$Time ; Data$Y ; Data$U

  errmsg <- ""
  ModelHasInput <- TRUE
  if ( is.null(Data$U) ) {
    ModelHasInput <- FALSE
  } else if(any(is.na(Data$U))) {
    errmsg <- (paste("Input U contains NA.")); return(list(errmsg = errmsg, ok=FALSE))
  } else if(!is.matrix(Data$U)) {
    errmsg <- (paste("Input U is not a matrix.")); return(list(errmsg = errmsg, ok=FALSE))
  } 

  if (ModelHasInput) {
    Uk <- Data$U[,1,drop=FALSE]
  } else {
    Uk <- NULL
  }

  if("Matrices" %in% names(Model)) {
    Linear = TRUE
  } else if ("Functions" %in% names(Model)) {
    Linear = FALSE
  } else {
    errmsg <- ("Model must contain an object named either Functions or Matrices") ; return(list(errmsg = errmsg, ok=FALSE)) }
      
  # Model contains the correct object
  if( any(!( c("X0","SIG","S","ModelPar")   %in% names(Model))) ) {
    errmsg <- ("Model does not contain $X0, $SIG, $S and $ModelPar objects. They must all be defined.") ; return(list(errmsg = errmsg, ok=FALSE)) }

  
  # Model objects are functions
  for (name in c("X0","SIG","S","ModelPar")) 
    if(!is.function(Model[[name]])) {
      errmsg <- (paste("Model object",name,"is not a function")); return(list(errmsg = errmsg, ok=FALSE)) }
  if(Linear) {
    if(!is.function(Model$Matrices)) {
      errmsg <- (paste("Model object Matrices is not a function")); return(list(errmsg = errmsg, ok=FALSE)) }
  }
       
#  if ( any(!(unlist(lapply(Model , function(x) {is.function(x)})))) ) {
#     errmsg <- ("Model objects are not all functions") ; return(list(errmsg = errmsg, ok=FALSE)) }

  # Check  Par
  if( !("Init" %in% names(Par)) ) {
    errmsg <- ("Init missing in Par") ; return(list(errmsg = errmsg, ok=FALSE))
  }
  if( "LB" %in% names(Par) ) {
    #The parameter has bounds
    if( any( Par$LB > Par$Init) ) {
      errmsg <- ("Parameter LB is greater than Init") ; return(list(errmsg = errmsg, ok=FALSE)) }
    if( any( Par$UB < Par$Init) ) {
      errmsg <- ("Parameter UB is lower than Init") ; return(list(errmsg = errmsg, ok=FALSE)) }
  }

  # Check data
  if( any(!( "Time"  %in% names(Data))) ) {
    errmsg <- ("Individual data does not contain Time.") ; return(list(errmsg = errmsg, ok=FALSE)) }
  if(DataHasY) {
    if( any(!( "Y"  %in% names(Data))) ) {
      errmsg <- ("Individual data does not contain Y.") ; return(list(errmsg = errmsg, ok=FALSE)) }
    if( any(is.nan(Data$Y)) ) {
      errmsg <- ("Data$Y contains NaN. Use NA instead.");return(list(errmsg = errmsg, ok=FALSE)) }
  }
  if( any(c(is.nan(Data$Time),is.na(Data$Time))) ){
    errmsg <- ("Data$Time contains NA or NaN.");return(list(errmsg = errmsg, ok=FALSE)) }
  if (ModelHasInput && any(c(is.nan(Data$U),is.na(Data$U))) ){
    errmsg <- ("Data$U contains NA or NaN.");return(list(errmsg = errmsg, ok=FALSE)) }

  
  # Calculate parameter phi
  Parlist <- Model$ModelPar(THETA=Par$Init)
  
  if( !( "h" %in% names(Model)) ) {
      errmsg <- ("Model is missing function h") ; return(list(errmsg = errmsg, ok=FALSE)) }

  if(!is.null(Parlist$OMEGA)) {
    if(!is.matrix(Parlist$OMEGA)) {
      errmsg <- ("Model$OMEGA is not a matrix") ; return(list(errmsg = errmsg, ok=FALSE)) }
    dimEta <- nrow(Parlist$OMEGA)
    phi <- Model$h(eta=rep(0,dimEta) , theta=Parlist$theta ,covar=Data$covar)
  } else {
    # OMEGA is NULL  
      phi <- Model$h(eta=NULL,theta=Parlist$theta,covar=Data$covar)      
  }
  if( prod(unlist(lapply(phi,length)))==0) {
    errmsg <- ("An element in phi returned from $h has length zero. The error could be caused by using
  eta in $h while using OMEGA=NULL.") ; return(list(errmsg = errmsg, ok=FALSE)) }
  

  # Matrices
  if(Linear) {
    tmp   <- Model$Matrices(phi=phi)
    matA  <- tmp$matA
    matB  <- tmp$matB
    matC  <- tmp$matC
    matD  <- tmp$matD
    
    X0 <- Model$X0( Time=Data$Time[1], phi=phi, U=Uk)
    SIG <-  Model$SIG(  phi=phi)
    S <- Model$S(  phi=phi )
  } else {
    X0 <- Model$X0( Time=Data$Time[1], phi=phi, U=Uk)
    f  <- Model$Functions$f(x=X0,u=Uk ,time=Data$Time[1],phi=phi)
    df <- Model$Functions$df(x=X0,u=Uk,time=Data$Time[1],phi=phi)
    g  <- Model$Functions$g(x=X0,u=Uk ,time=Data$Time[1],phi=phi)
    dg <- Model$Functions$dg(x=X0,u=Uk,time=Data$Time[1],phi=phi)
    SIG<- Model$SIG(u=Uk,time=Data$Time[1],phi=phi)
    S  <- Model$S(u=Uk,time=Data$Time[1],phi=phi)
  }

  #type check!
  if(!is.matrix(S)) {
    errmsg <- ("Model$S is not a matrix") ; return(list(errmsg = errmsg, ok=FALSE)) }
  if(Linear) {
    if(!is.matrix(matA)) {
      errmsg <- ("A is not a matrix") ; return(list(errmsg = errmsg, ok=FALSE)) }
    if(!is.matrix(matC)) {
      errmsg <- ("C is not a matrix") ; return(list(errmsg = errmsg, ok=FALSE)) }
    if(ModelHasInput && !is.matrix(matB)) {
      errmsg <- ("B is not a matrix") ; return(list(errmsg = errmsg, ok=FALSE)) }
    if(ModelHasInput && !is.matrix(matD)) {
      errmsg <- ("D is not a matrix") ; return(list(errmsg = errmsg, ok=FALSE)) }
  }


  all.equal2 <- function(target,current,...) {
    out <- all.equal(target,current,...)
    b <- isTRUE(out)
    attributes(b) <- list(msg=as.character(out))
    b
  }

  #f, df, g, dg
  if(!Linear) {
    dimX <- length(X0)
    if(!(is.matrix(f) & is.matrix(df) & is.matrix(g) & is.matrix(dg))) {
      errmsg <- ("Either f, df, g or df is not a matrix.") ; return(list(errmsg = errmsg, ok=FALSE)) }
    if(!all.equal2(dim(f),c(dimX,1))) {
      errmsg <- ("The dimension of the output from f is not c(dimX,1)") ; return(list(errmsg = errmsg, ok=FALSE)) }
    if(!all.equal2(dim(df),c(dimX,dimX))) {
      errmsg <- ("The dimension of the output from df is not c(dimX,dimX)") ; return(list(errmsg = errmsg, ok=FALSE)) }
    if(DataHasY) {
      dimY <- dim(Data$Y)[1]    
      if(!all.equal2(dim(g),c(dimY,1))) {
        errmsg <- ("The dimension of the output from g is not c(dimY,1)") ; return(list(errmsg = errmsg, ok=FALSE)) }
      if(!all.equal2(dim(dg),c(dimY,dimX))) {
        errmsg <- ("The dimension of the output from dg is not c(dimY,dimX)") ; return(list(errmsg = errmsg, ok=FALSE)) }
    } else {
      if(!dim(g)[1]==dim(dg)[1]){
        errmsg <- ("Number of rows in the output from g and dg are not equal") ; return(list(errmsg = errmsg, ok=FALSE)) }
      if(!dim(g)[2]==1){
        errmsg <- ("Number of columns in the output from g are not 1.") ; return(list(errmsg = errmsg, ok=FALSE)) }
      if(!dim(dg)[2]==dimX){
        errmsg <- ("Number of columns in the output from dg are not dimX.") ; return(list(errmsg = errmsg, ok=FALSE)) }
    }
    tmptest <- all.equal2(df,jacobian(func=Model$Functions$f,x=X0,u=Uk,time=Data$Time[1],phi=phi))
    if(!tmptest) {
        warning(paste("The output matrix from df does not appear to be equal to a numerical approximation of the jacobian of f.",
                      attributes(tmptest)$msg)) }
    tmptest <- all.equal2(dg,jacobian(func=Model$Functions$g,x=X0,u=Uk,time=Data$Time[1],phi=phi))
    if(!tmptest) {
      warning(paste("The output matrix from dg does not appear to be equal to a numerical approximation of the jacobian of g.",
                    attributes(tmptest)$msg)) }
  }


  # Test for positive semidefinit
  if(dim(SIG)[1]!=dim(SIG)[2] || dim(SIG)[1]==0) {
    errmsg <- ("Model$SIG is not a square matrix") ; return(list(errmsg = errmsg, ok=FALSE)) }
  if(dim(S)[1]!=dim(S)[2] || dim(S)[1]==0) {
    errmsg <- ("Model$S is not a square matrix") ; return(list(errmsg = errmsg, ok=FALSE)) }
      
  if( any( eigen(S)$values <0 ) ) {
    errmsg <- ("Model$S is not positiv semidefinit") ; return(list(errmsg = errmsg, ok=FALSE)) }

  
  # dimT
  dimT <- length(Data$Time)
  if(DataHasY) {
    if( dimT!=dim(Data$Y)[2]) {
      errmsg <- ("Data$Time and Data$Y does not have same length") ; return(list(errmsg = errmsg, ok=FALSE)) }
  }
  if(ModelHasInput) {
    if(dimT!=dim(Data$U)[2]) {
      errmsg <- ("Data$Time and Data$U does not have same length") ; return(list(errmsg = errmsg, ok=FALSE)) }
    } # ModelHasInput

  # dimX
  if(Linear) {
    dimX <- nrow(matA)
    if( dimX != ncol(matA) ) {
      errmsg <- ("A is not square [dimX dimX]") ; return(list(errmsg = errmsg, ok=FALSE)) }
    if( dimX!=ncol(matC) ) {
      errmsg <- ("A and C dimensions doesn't match") ; return(list(errmsg = errmsg, ok=FALSE)) }
    if( dimX!=length(X0) ) {
      errmsg <- ("X0 has incorrect size") ; return(list(errmsg = errmsg, ok=FALSE)) }
    if( dimX != nrow(SIG) | dimX != ncol(SIG) ) {
      errmsg <- ("SIG has incorrect size") ; return(list(errmsg = errmsg, ok=FALSE)) }
  }
  
  if(ModelHasInput && Linear) {
    if( dimX!=nrow(matB) ) {
      errmsg <- ("A and B dimensions doesn't match") ; return(list(errmsg = errmsg, ok=FALSE)) }
    } # ModelHasInput

  # dimY
  if(DataHasY) {
    dimY <- dim(Data$Y)[1]
    if (Linear && (dimY != nrow(matC)) ){
      errmsg <- ("C and Data$Y doesn't match") ; return(list(errmsg = errmsg, ok=FALSE)) }
    if( dimY != nrow(S) | dimY != ncol(S) ) {
      errmsg <- ("S has incorrect size") ; return(list(errmsg = errmsg, ok=FALSE)) }
  }

  if(Linear && ModelHasInput && DataHasY) {
    if (dimY != nrow(matD) ){
      errmsg <- ("Y and D doesn't match") ; return(list(errmsg = errmsg, ok=FALSE)) }
  }
  
  # dimU
  if( Linear && ModelHasInput) {
    dimU <- nrow(Data$U)
    if( dimU != ncol(matB)) {
      errmsg <- ("Data$U and B doesn't match") ; return(list(errmsg = errmsg, ok=FALSE)) }
    if( dimU != ncol(matD)) {
      errmsg <- ("Data$U and D doesn't match") ; return(list(errmsg = errmsg, ok=FALSE)) }
  }
  
  # Dose
  if( "Dose" %in% names(Model) ) {
      errmsg <- ("The dose specification $Dose in the model object has been moved to the data object (e.g. Data[[i]]$Dose) to allow for individual dosing.") ; return(list(errmsg = errmsg, ok=FALSE))
  }
  if( "Dose" %in% names(Data) ) {
    MD <- Data$Dose
    if( any(!(MD$Time %in% Data$Time))) {
      errmsg <- ("Dose times doesn't coincide with Data$Time") ; return(list(errmsg = errmsg, ok=FALSE)) }
      
    if( Linear && any( MD$State > dimX) ) {
      errmsg <- ("Dose states are larger than number of states") ; return(list(errmsg = errmsg, ok=FALSE)) }

    if( !all( (c(length(MD$State),length(MD$Amount))- length(MD$Time))==0  ) ) {
      errmsg <- ("Dose: Elements Time, State and Amount not of same length") ; return(list(errmsg = errmsg, ok=FALSE))}
  }

  return(list(errmsg = errmsg, ok=TRUE,Linear=Linear))
}
