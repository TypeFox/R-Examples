#' enaControl --- control analyses
#' INPUT = network object
#' OUTPUT = list of control statistics
#' M. Lau | July 2011

#' ---------------------------------------------------

enaControl <- function(x, zero.na=TRUE,balance.override=FALSE){
                                        #Missing Data Check
  if (any(is.na(x%v%'storage'))){
    warning('This function requires quantified storage values.')
  }else{
                                        #Check for network class
    if (class(x) != 'network'){warning('x is not a network class object')}
                                        #Check for balancing
    if (balance.override){}else{
      if (any(list.network.attributes(x) == 'balanced') == FALSE){x%n%'balanced' <- ssCheck(x)}
      if (x%n%'balanced' == FALSE){warning('Model is not balanced'); stop}
    }
                                        #Unpack data
    Ti <- x%v%'input'
    exp <- x%v%'export';exp[is.na(exp)] <- 0
    res <- x%v%'respiration';res[is.na(res)] <- 0
    Tj <- exp + res
    flow <- t(as.matrix(x,attrname="flow"))
    Q <- enaStorage(x)$Q
    QP <- enaStorage(x)$QP

                                        #Input perspective
    T. <- Ti + apply(flow,1,sum)
    GP <- flow / T.
    if (zero.na){GP[is.na(GP)] <- 0}else{}
    I <- GP * 0; diag(I) <- 1
    NP <- ginv((I - GP))

                                        #Output perspective
    T. <- apply(flow,2,sum) + Tj
    G <- t(t(flow) / T.)
    if (zero.na){G[is.na(G)] <- 0}else{}
    I <- G * 0; diag(I) <- 1
    N <- ginv((I - G))

                                        #Calculate the control matrix
    CN <- N / NP
    if (zero.na){CN[is.na(CN)] <- 0}
    for (i in (1:nrow(CN))){
      for (j in (1:ncol(CN))){
        if (CN[i,j] < 1){CN[i,j] <- 1 - CN[i,j]}else{CN[i,j] <- 0}
      }
    }

                                        #Storage version
    CQ <- ginv(QP) %*% Q
                                     #SCHRAMSKY'S ALGORITHM
    #Algorithm Source: Schramsky et al Eco. Mod. I94(2006) 189-201########
    T_out <- N %*% Ti                               # Component Output Throughflow T_out
    T_in <- Tj %*% NP                               # Component Input Throughflow T_in given the Output y(=Tj)
    T_in_mat <- matrix(T_in,nrow=length(T_in), ncol=length(T_in),byrow='TRUE')# T_in_mat repeats T_in vector in rows
    T_out_mat <- matrix(T_out,nrow=length(T_out),ncol=length(T_out),byrow='FALSE') # T_out_mat repeats T_out in cols
    eta_1<-N/T_out_mat                              # eta matrix. Can be done by 2 different ways. Way1
    eta_2<-NP/T_in_mat                              # Way2 ::: Theoretically eta_1==eta_2 but Numerical error
    eta<-eta_1                                      # eta:=eta_1. Can be defined as eta_2 to observe variation
    eta2<-t(eta)                                    # Taking the transpose of eta
    CD = eta - eta2                                 # Calculating the Control Difference Matrix
    CR = CD/pmax(eta,eta2)                          # Control Ratio Matrix. Dividing terms by Max eta value
    sc = apply(CD,2,sum)                            # Calculating the System Control Vector
    names(sc) <- rownames(flow)                     # Assigning Names
    rownames(CR) <- colnames(CR) <- rownames(flow)
    rownames(CD) <- colnames(CD) <- rownames(flow)
                                        #Name nodes
    rownames(CN) <- colnames(CN) <- rownames(flow)
    rownames(CQ) <- colnames(CQ) <- rownames(flow)
                                        #Re-orient output
    orient <- get.orient()
    if (orient == 'rc'){
      CN <- t(CN)
      CQ <- t(CQ)
      CR <- t(CR)
      CD <- t(CD)
    }else{}
                                        #Package up output
    out <- list('CN'=CN,'CQ'=CQ,'CR'=CR,'CD'=CD,'sc'=sc)

    return(out)
  }
}
