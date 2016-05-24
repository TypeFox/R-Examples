## #=============================================================================
## FUN.Pformula calls FUN.with.trans
##
## Takes: fomula = a formula-object
##        effect = character-object ("none", "time", "individual", "twoways")
## 
## Gives: A List with P+1 components
##        (first for depend. variable second for indep. variables):
##        Each Component is a List with 6 Components:
##        1. Tr  (one of "none", "time", "individual", "twoways")
##        2. I   (TRUE, if intercept, FALSE, if not)                        
##        3. ODM (Original Data-Matrix)
##        4. TDM (Transformed Data-Matrix)
##        5. TDV (TDM as a vector)
##        6. TRm is a list with 3 Components:
##           1. OVc  (Overall constant)
##           2. InC  (Individual constants)
##           3. TiVC (Time variaing constants)
## #==============================================================================

FUN.Pformula <- function(formula, effect= c("none", "individual", "time", "twoways"))
  {
    data.fra <- model.frame(formula)
    dat.term <- attr(data.fra, "terms")	

    ## Construct data from formula
    ## dim(response) == TxN == dim(y.matrix) == TxN:
    y.matrix <- model.response(data.fra, "numeric")
    ## Check Data: Each Variable has to be a matrix
    if(!is.matrix(y.matrix)) stop("Each Variable has to be a TxN-matrix.")
    N  <- ncol(y.matrix)
    T  <- nrow(y.matrix)

    ## 1)Extract Regressors
    ## 2)Check the presence of a 'intercept' in the formula
    ## 3)And built the x.all.matrix (TxN*P)

    regressors.mat <- model.matrix(dat.term, data.fra)
    is.intercept   <- ifelse(colnames(regressors.mat)[1] == "(Intercept)", TRUE, FALSE)
    if(is.intercept){
      x.all.matrix       <- regressors.mat[,-1]
    }else{x.all.matrix   <- regressors.mat}

#    if(!is.intercept & effect=="twoways"){stop("Effects >> twoways << need an Intercept!")}
    
    ## Dimension parameters

    NT <- N*T
    P  <- as.integer(ncol(x.all.matrix)/N)
    if(P!= ncol(x.all.matrix)/N) stop("All the regressors must have the same dimension as the response variable Y")

    data.all.mat  <- cbind(y.matrix, x.all.matrix)
    
    ## Write the response variable, Y, and the 'p' regressors, X, in a list,
    ## where each component contains one of p+1 TxN-Matrices
    model.in.list <- lapply(1:(P+1), function(z, i) z[,seq((i-1)*N+1,i*N)], z = data.all.mat)
    
    ## Transform the response variable as well as the 'P' regressors and give them in a list where ech 
    ## componente contains also a list with:
    ##       1- "Tr" Name of the transformation
    ##       2- "I" Logical variable if ther is intercept or no
    ##       3- "ODM" Original Data matrix
    ##       4- "TDM" Transformed Data in a matrix
    ##       5- "TDV" Transformed Data in a NT x 1 Vector
    ##       6- "TRm" Sublist with
    ##           a- "OVc" Overall Constant
    ##           b- "InC" time constant individual effects
    ##           c- "TiVC" additive time varying effects


    data.in.list  <- sapply(model.in.list,
                            function(z) FUN.with.trans(z,
                                                       N            = N,
                                                       T            = T,
                                                       is.intercept = is.intercept,
                                                       effect       = effect),
                            simplify=FALSE)
    data.in.list
  }
