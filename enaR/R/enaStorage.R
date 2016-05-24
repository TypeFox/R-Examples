#' enaStorage --- storage analysis
#' INPUT = network object
#' OUTPUT = list of storage statistics
#'
#' M. Lau | July 2011
#' ---------------------------------------------------

enaStorage <- function(x,balance.override=FALSE){
                                        #Missing Data Check
  if (any(is.na(x%v%'storage'))){
    warning('This function requires quantified storage values.')
  }else{
                                        #Check for network class
    if (class(x) != 'network'){warning('x is not a network class object')}

                                        #Check for balancing
    if (balance.override){}else{
      if (any(list.network.attributes(x) == 'balanced') == FALSE){x%n%'balanced' <- ssCheck(x)}
      if (x%n%'balanced'){}else{stop('Model is not balanced')}
    }
                                        #unpack data from x
    Flow <- t(as.matrix(x, attrname = 'flow'))  #flows
                                        #continue unpacking
    input <- x%v%'input' #inputs
    stor <- x%v%'storage' #storage values
    T. <- apply(Flow,1,sum) + input
    FD <- Flow - diag(T.) #flow matrix with negative throughflows on the diagonal
    I <- diag(1,nrow(Flow),ncol(Flow)) #create the identity matrix

                                        #Compute the Jacobian matrix
    C <- FD %*% ginv(diag(stor)) #output matrix
    CP <- ginv(diag(stor)) %*% FD #input matrix

                                        #smallest whole number to make diag(C) nonnegative
    dt <- -1 / floor(min(diag(C)))

                                        #calculating the storage-specific, output-oriented, intercompartmental flows (P)
    P <- I + C*dt
    PP <- I + CP*dt
                                        #calculating the dimensionalized integral output and input matrices
    S <- -ginv(C) #output
    SP <- -ginv(CP) #input

                                        #calculating the integral storage intensity matrix (Q)
    Q <- ginv(I - P) #output
    QP <- ginv(I - PP) #input
    dQ <- diag(Q) #diagonal of integral output storage matrix which is the same for input (i.e. diag(QP))

                                        #naming row and columns
    rownames(C) <- colnames(C) <- rownames(Flow)
    rownames(CP) <- colnames(CP) <- rownames(Flow)
    rownames(P) <- colnames(P) <- rownames(Flow)
    rownames(S) <- colnames(S) <- rownames(Flow)
    rownames(Q) <- colnames(Q) <- rownames(Flow)
    rownames(CP) <- colnames(CP) <- rownames(Flow)
    rownames(PP) <- colnames(PP) <- rownames(Flow)
    rownames(SP) <- colnames(SP) <- rownames(Flow)
    rownames(QP) <- colnames(QP) <- rownames(Flow)

    ##Storage Environ Properties
                                        #eigen analysis
    e <- eigen(P)$values
    lam1P <- e[1]
    rhoP <- e[1] / e[2]

    eP <- eigen(PP)$values
    lam1PP <- eP[1]
    rhoPP <- eP[1] / eP[2]

    TSS <- sum(stor) #total system storage
    TSScs <- sum(((dQ-1)/dQ)%*%stor) #cycled (mode 2) storage
    CIS <- TSScs / TSS #cucling index (storage)

                                        #Amplification parameter
    NAS <- length((Q-diag(diag(Q)))[(Q-diag(diag(Q))) > 1])
    NASP <- length((QP-diag(diag(QP)))[(QP-diag(diag(QP))) > 1])

                                        #Indirect effects parameter (srb fix 8.3.2011)
    ID.S.O <- sum(Q-I-P) /sum(P)
    ID.S.I <- sum(QP-I-PP) / sum(PP) #indirect to direct ratio (input matrix)

                                        #Indirect effects parameter (realized)  (srb fix 8.3.2011)
    ID.S <- sum(dt*(as.matrix((Q-I-P)) %*% input)) / sum(dt*as.matrix(P)%*% input ) #indirect to direct ratio (output matrix)

    BSI = sum( I %*% input *dt) / TSS
    DSI = sum( P %*% input *dt) / TSS
    ISI = sum( (Q-I-P) %*% input *dt) / TSS

                                        #Homogenization parameter
    CVP <- sd(as.numeric(P)) / mean(P) #Coefficient of variation for G
    CVQ <- sd(as.numeric(Q)) / mean(Q)  #Coefficient of variation for N
    HMG.S.O <- CVP / CVQ #homogenization parameter (output storage)

    CVPP <- sd(as.numeric(PP)) / mean(PP) #Coefficient of variation for GP
    CVQP <- sd(as.numeric(QP)) / mean(QP) #Coefficient of variation for NP
    HMG.S.I <- CVPP / CVQP #homogenization paraemeter (input storage)

                                        # Network Aggradation
    AGG.S <- TSS / sum(input) #network aggradation -- average amount of storage per system input

                                        # MODE Partition (Fath et al. 2001)
    z <- unpack(x)$z
    mode0.S = sum(z*dt)  # storage from boundary input flow
    mode1.S = sum( (ginv(diag(diag(Q))) %*% Q - I) %*% diag(z) * dt ) # storage from first-passage flow
    mode2.S = sum( diag(diag(Q) - 1) %*% (ginv(diag(diag(Q))) %*% Q) %*%  diag(z) * dt) # storage from compartment-wise dissipative flow
    mode3.S = mode1.S  # storage from first-passage flow (equal to mode 1 at Steady state)
    mode4.S = sum(unpack(x)$y * dt)  # storage from boundary output flow

    ##packing up network statistics for output
    ns <- cbind('TSS'=TSS,
                'CIS'=CIS,
                'BSI'= BSI,'DSI'= DSI,'ISI'= ISI,
                'ID.S'=ID.S, 'ID.S.I'=ID.S.I,'ID.S.O'=ID.S.O,
                'HMG.S.O'=HMG.S.O,'HMG.S.I'=HMG.S.I,
                'NAS'=NAS,'NASP'=NASP,
                mode0.S,mode1.S,mode2.S,mode3.S,mode4.S)

                                        #'lam1P'=abs(lam1P),'rhoP'=abs(rhoP),
                                        #'lam1PP'=abs(lam1PP),'rhoPP'=abs(rhoPP),'AGG.S'=AGG.S)
                                        #re-orientation
    orient <- get.orient()
    if (orient == 'rc'){
      C <- t(C)
      P <- t(P)
      S <- t(S)
      Q <- t(Q)
      CP <- t(CP)
      PP <- t(PP)
      SP <- t(SP)
      QP <- t(QP)
    }else{}

    out <- list('X'=stor,'C'=C,'P'=P,'S'=S,'Q'=Q,'CP'=CP,'PP'=PP,'SP'=SP,'QP'=QP,'dt'=dt,'ns'=ns)

    return(out)
  }
}
