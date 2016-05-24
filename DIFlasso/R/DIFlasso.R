DIFlasso <-
function (Y, X, l.lambda = 20, trace = FALSE) 
{
  vardiffs <- abs(apply(X,2,var)-1)
  
  if(sum(vardiffs)>1e-6)
    stop("X has to be standardized")
  
  if(sum(!apply(Y,2,unique) %in% c(0,1))>0)
    stop("Y may only contain 0 or 1")
  
  if(!is.data.frame(X))
    stop("X has to be a data.frame")
  
  if(!is.data.frame(Y))
    stop("Y has to be a data.frame")
  
    # print trace?
    if (trace) {
        trace <- 1
    }else {
        trace <- 0
    }
    # number of persons
    P <- nrow(Y)
    # number of items
    I <- ncol(Y)
    # number of covariates
    l <- ncol(X)
    
    # Convert data.frames into matrices
    names.y <- names(Y)
    names.x <- names(X)
    
    Y <- matrix(as.double(as.matrix(Y)),ncol=I,byrow=FALSE)
    X <- matrix(as.double(as.matrix(X)),ncol=l,byrow=FALSE)
    
    # which persons have to be excluded? 
    exclu1 <- rowSums(Y) != I & rowSums(Y) != 0
    exclu <- rep(exclu1, each = I)
    
    # create the part of the design matrix containing the covariates
    xp <- matrix(0, ncol = l * I, nrow = P * I)
    suppressWarnings(for (q in 1:P) {
        xp[((q - 1) * I + 1):(q * I), ] <- matrix(as.double(c(X[q, 
            ], rep(0, l * I))), byrow = T, ncol = l * I, nrow = I)
    })
    xp <- xp[exclu, ]
    
    # create the response vector
    XP <- c()
    for (t in 1:P) {
        XP <- c(XP, Y[t, ])
    }
    XP <- XP[exclu]
    
    # new (reduced) number of persons
    ex <- sum(!exclu1)
    P <- P - ex
    
    # create the part of the design matrix containing the columns for person parameters theta
    help1 <- c(1, rep(0, P - 2))
    help2 <- c(rep(help1, I), 0)
    suppressWarnings(help3 <- matrix(help2, ncol = P - 1, nrow = P * 
        I, byrow = TRUE))
    help3[((P - 1) * I + 1):(P * I), ] <- 0
    
    # create the part of the design matrix containing the columns for item parameters beta
    help4 <- rep(diag(I), P)
    help5 <- matrix(help4, ncol = I, nrow = P * I, byrow = TRUE)
    
    # design matrix
    design.matrix <- cbind(help3, -help5, -xp)
    
    # index vector for group lasso penalization
    index = c(rep(NA, sum(exclu1) + I - 1), rep(1:I, each = l))
    
    # calculate maximal lambda
    lmax <- lambdamax(design.matrix, XP, index, penscale = sqrt, 
        model = LogReg(), center = FALSE, standardize = FALSE)
    
    # create sequence of lambdas
    lambda.seq <- seq(from = lmax * 1.05, to = 1e-06, length = l.lambda)
    
    # compute group lasso model
    suppressWarnings(m2 <- grplasso(design.matrix, XP, index, 
        lambda = lambda.seq, penscale = sqrt, model = LogReg(), 
        center = FALSE, standardize = FALSE, control = grpl.control(trace = trace)))
    
    # index for penalized parameters
    index2 <- index[!is.na(index)]
    # parameters for smallest lambda
    coef.unpen <- m2$coef[!is.na(index), l.lambda]
    # norm for parameters for smallest lambda
    norm.unpen <- c()
    for (o in 1:max(index2)) {
        norm.unpen[o] <- sqrt(sum(coef.unpen[index2 == o]^2))
    }
    # group sizes
    p.j <- table(index2)
    # matrix of coefficients
    coefs.unrest <- coefs <- m2$coef
    # matrix for gammas
    coefs.grp <- coefs[!is.na(index), ]
    # norms of gammas
    norm.gamma <- matrix(0, ncol = ncol(coefs), nrow = max(index2))
    for (j in 1:ncol(coefs)) {
        for (o in 1:max(index2)) {
            norm.gamma[o, j] <- sqrt(sum(coefs.grp[index2 == o, 
                j]^2))
        }
    }
    # degrees of freedom of all gammas
    df.grp <- colSums(norm.gamma > 0) + colSums((norm.gamma/norm.unpen) * 
        (c(p.j) - 1))
    # degrees of freedom for thetas and betas
    df.unpen <- sum(is.na(index))
    # total degrees of freedom
    df <- df.grp + df.unpen
    # predicted probabilities
    pi.i <- predict(m2, type = "response")
    # log likelihood
    loglik <- colSums(XP * log(pi.i) - XP * log(1 - pi.i) + log(1 - 
        pi.i))
    
    # BIC
    bic <- -2 * loglik + df * log(length(XP))
    # AIC
    aic <- -2 * loglik + df * 2
    
    ####################
    # which item will be the reference/anchor item
    ref.item <- NULL
    ind0 <- l.lambda
    while(is.null(ref.item)){
      cand <- which(colSums(matrix(coefs.grp[,ind0],nrow=l))==0)
      if(length(cand)>0){
        ref.item <- max(cand)
      }else{ind0 <- ind0-1}
    }
    
    restrict <- function(x){
      x <- x - rep(x[(l*ref.item-l+1):(l*ref.item)],I)
      return(x)}
    
    # center gamma around refernce item, only relevant when gamma_ref.item not zero
    coefs.grp <- apply(coefs.grp,2,restrict)
    
    # coefs2: add one row to coefs matrix because theta_P is not set zero anymore
    coefs2 <- matrix(0,nrow=nrow(coefs)+1,ncol=ncol(coefs))
    coefs2[(P+I+1):(nrow(coefs2)),] <- coefs.grp
    
    # estimated thetas, centered around reference item ref.item
    theta.hat <- predict(m2)[ref.item +seq(from=0, to= P*I-I, by=I),]
    
    # estimated betas, centered around reference item ref.item
    beta.hat <- coefs[(P):(P + I - 1), ]
    beta.hat <- t(apply(beta.hat,1,function(x){return(x - beta.hat[ref.item,])}))
    
    
    coefs2[1:P, ] <- theta.hat 
    # estimated betas
    coefs2[(P+1):(P + I), ] <- beta.hat 
    
    dif.mat <- matrix(coefs.grp[,which.min(bic)],nrow=l, dimnames=list(names.x,names.y))
    no.dif.cols <- which(coefs.grp[,which.min(bic)]==0)
    design.matrix <- insertCol(design.matrix,P,c(rep(0,(P-1)*I),rep(1,I)))
    design.matrix<- design.matrix[,-c(P+ref.item,P+I+no.dif.cols)]
    design.matrix <- cbind(design.matrix,XP)
    
    dif.items <- which(colSums(dif.mat)!=0)
    
    returns <- list(theta = theta.hat, beta = beta.hat, gamma = coefs.grp, P = P, I = I, m = l, 
                    logLik = logLik, BIC = bic, AIC = aic, df = df, refit.matrix = design.matrix, 
                    lambda = lambda.seq, ref.item = ref.item, dif.mat = dif.mat, dif.items = dif.items,
                    names.y = names.y, names.x = names.x, removed.persons = which(!exclu1))
                      
    class(returns) <- "DIFlasso"
                      
    return(returns)
}
