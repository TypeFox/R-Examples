`%+%` <- function(a, b) paste0(a, b)

#------------------------------------
# sim_genk6dat_2.R
# author: Oleg Sofrygin 
# DATA SIMULATOR FOR NETWORKS WITH K=6
# Source the file to simulate network using these functions
#------------------------------------

#------------------------------------
# SEARCH TRIAL SIMULATION
#------------------------------------			
# Chance of getting infected P(Y=1) is logit linear in 
# no. of infected (W1=1) AND untreated (A=0) friends in the network
# For those who have no infected friends or
# no untreated infected friends P(Y=1)=0
# Thus W's have no direct effect on Y, 
# except for deterministically setting Y=1 if W1=1
# The source of confounding is W2=W_risk, 
# which is predictive of low chance of getting treated, P(A=1), as well as 
# acting on a lot of people (a lot of people have that person in their network)
# THUS, network N is constructed conditional on W2=W_risk
#------------------------------------	
# CONFOUNDING
#------------------------------------	
#1) W1[i] affects the total number of friends i has
#2) W1[i] affects the probability of having [i] selected as someone's friend

netvar <- function(varnm, fidx) { # OUTPUT format: Varnm_net.j
  cstr <- function(varnm, fidx) {
    slen <- length(fidx)
    rstr <- vector(mode = "character", length = slen)
    netidxstr <- !(fidx %in% 0L)
    rstr[netidxstr] <- stringr::str_c('_netF', fidx[netidxstr])  # vs. 1
    return(stringr::str_c(varnm, rstr))
  }
  if (length(varnm) > 1) {
    return(unlist(lapply(varnm, cstr, fidx)))
  } else {
    return(cstr(varnm, fidx))
  }
}

# get network A's & W's as a matrix
.f.redefineCov <- function(k, Var, VarNm, Net_vect, misval = 0L) {
  #get all friends Ws as a matrix of dim(n,k) filling unused cols with zeros
  .f.netCovar <-function(Covar, Net) sapply(Net, function(netwk)
    											c(Covar[netwk],
    												rep.int(misval, k-length(netwk))))
  netVar_full <- NULL
  netVar_names <- NULL
  sumnet_name <- NULL
  if (k>0) {
  	netVar_full <- .f.netCovar(Var, Net_vect)
  	if (k>1) netVar_full <- t(netVar_full)  
  	netVarNm <- paste("net", VarNm, "_", sep="")
  	netVar_names <- netvar(VarNm, c(1:k))
  }
  Var_names <- c(VarNm, netVar_names)	
  d <- cbind(Var, netVar_full)
  colnames(d) <- Var_names
  return(d)
}

#---------------------------------------------------------------------------------
# DEFINITION OF TREATMENT MECHANISM g(A|W) AND g*(A|W)
#---------------------------------------------------------------------------------
# Actual pop g(A|W) from NPSEM, logistic fcn with a coeff for network Ws and a coeff for individ W_i
# DEFINITION OF g(A_i|W)=P[A=1|cA]:  
# W1 - inf.RISK
# W2 - (0,1) infected at time t=0
# g_0 = b0 + b1*W1 + b2*(sum(netW1))
f.A <- function(k, data, ...) {
  .f_evalA_mtx <- function(netW1, netW2, netW3, n, ...) {
    sum_friendWs <- matrix(0, nrow = n, ncol = 1)
    for (k_ind in (1:dim(netW1)[2])) {
      sum_friendWs <- sum_friendWs + netW1[, k_ind] * coeff_A_W1friend +
                      netW2[, k_ind] * coeff_A_W2friend +
                      netW3[, k_ind] * coeff_A_W3friend
    }
    return(sum_friendWs)
  }
  
  # coefficients for true g(A|W)
  # Define coefficient for i individual's W
  Intercept_A <- 2
  coeff_A_W <- -0.5
  coeff_A_W1friend <- -0.1
  coeff_A_W2friend <- -0.4
  coeff_A_W3friend <- -0.7

  n <- nrow(data)
  W1 <- data[,"W1"]
  netW1 <- data[,netvar("W1", (1:k))]
  netW2 <- data[,netvar("W2", (1:k))]
  netW3 <- data[,netvar("W3", (1:k))]


  # P(A) component from individual infection risk (W1)
  indivW <- coeff_A_W * W1
  	sum_friendWs <- .f_evalA_mtx(netW1, netW2, netW3, n, ...)
  probA <- plogis(Intercept_A + indivW + sum_friendWs)
  	return(probA) # No deterministic treatments
}

# Set x% of community to A=1 (returns probability P(A=1))
f.A_x <- function(data, x, ...) rep(x, nrow(data))

# Set x% of community to A=1 only among W2=1
  # W1 - inf.RISK
  # W2 - (0,1) infected at time t=0
f.A_xlevelW2_1 <- function(data, x, ...) {
  # print("x"); print(x)
  n <- nrow(data)
  W2 <- data[, "W2"]
  pA <- rep(0,n)
  pA[which(W2==1)] <- x
  return(pA)
}

# Set x% of community to A=1 based on connectivity |N_i|= {low, high}, based on median
f.A_xlevelNi <- function(data, x_Ni_low, x_Ni_high, ...) {	
  n <- nrow(data)
  nFriends <- data[, "nFriends"]
  Ni_med = quantile(nFriends, 0.5)
  x <- rep(0,n)
  x[which(nFriends<= Ni_med)] <- x_Ni_low
  x[which(nFriends> Ni_med)] <- x_Ni_high
  return(x)  
}

# Deterministically set A=1 based on cutt-off value for Ws or if W_i = 1
f.A_cutt_offW <- function(k, data, cutt_offW, ...) {	
  n <- nrow(data)
  W1 <- data[,"W1"]
  netW1 <- data[,netvar("W1", (1:k))]
  sum_friendWs <- matrix(0, nrow=n, ncol=1)
  for (k_ind in (1:dim(netW1)[2])) {
    sum_friendWs <- sum_friendWs + netW1[, k_ind]
  }
  indivW <- W1
  A <- rep(0,n)
  A[which((indivW==1)|(sum_friendWs>cutt_offW))] <- 1 
  return(A)  	
}

#---------------------------------------------------------------------------------
# Generate the network population using structural equations model
#---------------------------------------------------------------------------------
EC <- 0.35
NONFIXED_N_FRIENDS <- TRUE 	# Simulate # of network friends from uniform (vs fixed)
SYMM_CONN <- FALSE 	# Set SYMM_CONN=T to generate symmetric connectivity mtx
Intercept_Y <- -1
coeff_Y_AW2 <- 2
coeff_Y_W3 <- -1.5

#Sample 1 community (C_j) with EC and f.g_A_W=A^C
gendata_Cj <- function(C_j = 1, n, k, EC, f.g_name=NULL, f.g_args=NULL) { 
  #n - number of individuals
  #k - max size of each individual's network
  #C_j - community # being sampled (j)
  #EC - community-level covariate, only influences W_i
  #f.g_name - fcn for community intervention on A^C, i.e. g(A|W)
  #f.g_args - additional args to be passed to f.g_name (as a list)

  #----------------------------------------------------------------
  # Defining structural equations 
  #----------------------------------------------------------------
  # W1 - categorical or continuous confounder (5 categories, 0-4)
  # W2 - baseline infection status (0/1)
  # W3 - binary confounder (0/1)
  nW1cat <- 6
  .f.W1 <- function(n) rbinom(n, 5, prob=c(0.4, 0.5, 0.7, 0.4))
  # W2 - binary infection status at t=0, positively correlated with W1
  .f.W2 <- function(Cj_prob, W1, n) {
  	# prob_W2 <- seq(0.25, 0.5, by=0.3/nW1cat)
  	prob_W2 <- seq(0.45, 0.8, by=0.3/nW1cat)
  	return(sapply(c(1:n), function(i) rbinom(1, 1, prob=prob_W2[W1[i]+1])))
  }
  .f.W3 <- function(Cj_prob2, W1, W2, n) {
  	# prob_W2 <- seq(0.25, 0.5, by=0.3/nW1cat)
  	prob_W3 <- 0.6
  	return(rbinom(n, 1, prob=prob_W3))
  }
	# Total number of friends for each i, influenced by W1 - inf. risk
  .f.Net_num <- function(W1, samp=FALSE, n, k) {
  # network confounding through W1
  k_arr <-c(1:k)
  pN_0 <- 0.02
  #1) W1[i] affects the total number of friends i has (prob goes up)
  # W1=0 probabilities of |F_i|
  prob_Ni_W1_0 <- c(pN_0, plogis(-3 - 0 - k_arr/2))
  # W1=1 probabilities of |F_i|
  prob_Ni_W1_1 <- c(pN_0, plogis(-1.5 - 0 - k_arr/3))
  # W1=2 probabilities of |F_i|
  prob_Ni_W1_2 <- c(pN_0, pnorm(-2*abs(2 - k_arr) / 5))
  # W1=3 probabilities of |F_i|
  prob_Ni_W1_3 <- c(pN_0, pnorm(-2*abs(3 - k_arr) / 5))
  # W1=4 probabilities of |F_i|
  prob_Ni_W1_4 <- c(pN_0, plogis(-4 + 2*(k_arr-2)))
  # W1=5 probabilities of |F_i|
  prob_Ni_W1_5 <- c(pN_0, plogis(-4 + 2*(k_arr-3)))

  prob_Ni <- list(prob_Ni_W1_0, prob_Ni_W1_1, prob_Ni_W1_2, prob_Ni_W1_3, prob_Ni_W1_4, prob_Ni_W1_5)
  prob_Ni <- lapply(prob_Ni, function(x) x/sum(x))
  Net_num <- sapply(c(1:n), function(i) sample(0:k, 1, replace = TRUE, prob = prob_Ni[[W1[i]+1]]))
  	return(Net_num)
  }
  # W1-based probs of i being selected as someone's friend
  W1cat_arr <- c(1:nW1cat)/2
  prob_f <- plogis(-4.5 + 2.5*W1cat_arr) / sum(plogis(-4.5 + 2.5*W1cat_arr))
  # Sample connectivity matrix (0,1), 
  # randomly selecting from the list of available individuals
  # so that we never exceed |N_i| for each i
  # may result in actual # of friends being < |N_i| for some
  .f.genConnectMatx_asym_biased <- function(W1, Net_num) {
    # I <- bigmemory::big.matrix(n,n, type="short", init=0, shared=FALSE)
    I <- matrix(0L, nrow = n, ncol = n)
    nFriendTot <- array(rep(0L,n))
    for (index in (1:n)) {
      I[index,index] <- 1
      FriendSampSet <- setdiff( c(1:n), index)  #set of possible friends to sample, anyone but itself
      nFriendSamp <- max(Net_num[index] -	nFriendTot[index], 0) #check i's network is not already filled to max
      if (nFriendSamp>0) {
        if (length(FriendSampSet)==1)  {  #To handle case where |FriendSampSet|=1
          friends_i <- FriendSampSet
        } else { #sample from the possible friend set, with prob based on W2=W_risk
          # W1[i] affects the probability of having [i] selected as someone's friend
          # W1-based probs of i being selected
          friends_i <- sort(sample(FriendSampSet, nFriendSamp, prob = prob_f[W1[FriendSampSet]+1]))
        }
        I[friends_i, index] <- 1
        nFriendTot[index] <- nFriendTot[index] + nFriendSamp
      }
    }
    return(I)
  }
  #Update the # of friends for each individual (given sampled connnectivity matrix)
  # .f.Net_num_update <- function(ConnectMatx) return(colsum(ConnectMatx, c(1:n)) - 1)
  .f.Net_num_update <- function(ConnectMatx) return(base::colSums(ConnectMatx) - 1)
  #Convert connectivity matx to a lists of vector of friend's ids
  .f.Net_vect <- function(ConnectMatx) {
    f.netwklist_i <- function(index)  {
      netwklist_i <- setdiff(which(ConnectMatx[, index]!=0), index)
      netwklist_i	}
    sapply(1:n, f.netwklist_i, simplify=F)
  }
  #Same, but using mwhich() instead of which() (faster for n>30,000)
  # .f.Net_vect_big <- function(ConnectMatx) {
  #   f.netwklist_i <- function(index) setdiff(mwhich(ConnectMatx, index, 0, 'neq'), index)
  #   sapply(1:n, f.netwklist_i, simplify=F)
  # }
  .f.Net_vect_ID <- function(ConnectMatx) {
    f.netwklist_i <- function(index)  {
      netwklist_i <- setdiff(which(ConnectMatx[, index]!=0), index)
      if (length(netwklist_i)>0) netwklist_i <- paste('I', netwklist_i, sep="")
      return(netwklist_i)}
    sapply(1:n, f.netwklist_i, simplify=F)
  }
  # .f.Net_vect_bigID <- function(ConnectMatx) {
  #   f.netwklist_i <- function(index)  {
  #     netwklist_i <- setdiff(mwhich(ConnectMatx, index, 0, 'neq'), index)
  #     if (length(netwklist_i)>0) netwklist_i <- paste('I', netwklist_i, sep="")
  #     return(netwklist_i)}
  #   sapply(1:n, f.netwklist_i, simplify=F)
  # }
  .f.mkstrNet <- function(Net) sapply(Net, function(Net_i) paste(Net_i, collapse=' '))
  .f.mkstrNetWi <- function(Net,i) sapply(Net, function(Net_i) paste(Net_i[[i]], collapse=" "))
    .f_g_wrapper <- function(k, W_netW, fcn_name, ...) {   #wrapper fcn for g(A|W)
      # assign(".f_gAW", get(fcn_name))
      args0 <- list(k=k, data=W_netW)
      args <- c(args0, ...)
      rbinom(n, 1, do.call(fcn_name, args))
    }
  # get all friends Ws as a list of vectors
  .f.cA_Pa <-function(W1, W2, W3, Net) sapply(Net, function(netwk) 
  										list(W1=W1[netwk], W2=W2[netwk], W3=W3[netwk]), simplify=F)
  # Calculate c^Y(Pa(Y_i)) - a function into R^(k+1), 
  # doesn't depend on n or i and is permutation invariant
  # for each i, get the W's and A's in the network => cY is a vector of (W,A),
  # not including individual (W_i,A_i)
  .f.cY_Pa <-function(W1, W2, W3, A, Net) {
    cY_Pa <- list(Pa_Ws=sapply(Net, function(netwk) list(W1=W1[netwk], W2=W2[netwk], W3=W3[netwk]), simplify=F),
                  Pa_As=sapply(Net, function(netwk) A[netwk], simplify=F))
  }
  #Define Y, using the network size |N|, netW's , netA, W_i's, A_i
  # Q0 = b0 + b1*I(N(inf & untrt)=0) + b2*N(inf & untrt) + b3*N_friends
  .f.Y <- function(W1, W2, W3, A, cY_Pa) {  
    # GET (total # of friends who are infected AND untreated):
    f.netY_AW2As <- function(i, coeff)  sum(coeff * cY_Pa$Pa_Ws[[i]]$W2 * (1-cY_Pa$Pa_As[[i]]))
    f.netY_W3 <- function(i, coeff)  sum(coeff * cY_Pa$Pa_Ws[[i]]$W3)

    sum_friendY_W2As <- sapply(c(1:n), f.netY_AW2As, coeff_Y_AW2)
    sum_friendY_W3s <- sapply(c(1:n), f.netY_W3, coeff_Y_W3)

    # N untreated in i's network
    f.netY_As <- function(i)  sum(1-cY_Pa$Pa_As[[i]])
    sum_friendA0 <- sapply(c(1:n), f.netY_As)
    # N infected in i's network
    f.netY_W2s <- function(i)  sum(cY_Pa$Pa_Ws[[i]]$W2)
    sum_friendW21 <- sapply(c(1:n), f.netY_W2s)

    p_YRisk <- plogis(Intercept_Y + sum_friendY_W2As + sum_friendY_W3s)
    Y <- rbinom(n, 1, p_YRisk)
    # No deterministic outcomes
    return(Y)
  }

  #-----------------------------------------------
  #Generating covars
  #-----------------------------------------------  
  print("---------------------------")
  # print("f.g:"); print(f.g_name); 
  print(paste("f.g_args:", f.g_args));
  #-----------------------------------------------    	
  W1 <- .f.W1(n)  	# categorical infection risk (confounder)
  W2 <- .f.W2(EC, W1, n) # infected at baseline
  W3 <- .f.W3(EC, W1, W2, n) # binary confounder
  # Generate # of friends, |N_i|
  nFriends <- .f.Net_num(W1, samp=NONFIXED_N_FRIENDS, n, k)
  #---------------------------------------------------------
  # Generate symmetric connectivity matrix of friends, by infection risk W1
  # if (SYMM_CONN) {
  #   ConnectMatx <- .f.genConnectMatx_sym(nFriends)
  # } else {
    ConnectMatx <-  .f.genConnectMatx_asym_biased(W1, nFriends)   
  # }
  # Generate asymmetric connectivity matrix of friends
 	#--------------------------------------------------------- 	
  #Update # of actual friends sampled, |N_i|
  nFriends <- .f.Net_num_update(ConnectMatx)
  #Get list of vectors of friend's ids for each i
  Net_vect <- .f.Net_vect(ConnectMatx)
  Net_vectID <- .f.Net_vect_ID(ConnectMatx)
  W1_netW <- data.frame(.f.redefineCov(k, W1, "W1", Net_vect))
  W2_netW <- data.frame(.f.redefineCov(k, W2, "W2", Net_vect))
  W3_netW <- data.frame(.f.redefineCov(k, W3, "W3", Net_vect))
  W_netW <- cbind(W1_netW, W2_netW, W3_netW)
  # SUM netW for each W separately
  .f_netWsum <- function(netW1, netW2, netW3) {
    sum_netWs <- matrix(0, nrow=n, ncol=3)
    for (k_ind in (2:dim(netW1)[2])) {
      sum_netWs[,1] <- sum_netWs[,1] + netW1[, k_ind]
      sum_netWs[,2] <- sum_netWs[,2] + netW2[, k_ind]
      sum_netWs[,3] <- sum_netWs[,3] + netW3[, k_ind]
    }
    colnames(sum_netWs) <- c("netW1_sum", "netW2_sum", "netW3_sum")
    return(sum_netWs)
  }
  netW123_sum <- .f_netWsum(W1_netW, W2_netW, W3_netW)

  cA_Pa <- .f.cA_Pa(W1, W2, W3, Net_vect)   #Get c^A - fcn for parents of A_i: (W)  	
  if (is.null(f.g_name)) {  
    A <- .f_g_wrapper(k, W_netW, "f.A") }
  else { 
    A <- .f_g_wrapper(k, W_netW, f.g_name, f.g_args)
  }

  #convert f.g_args to text format (for output)
  if (is.null(f.g_args)) {
    f.g_args_txt <- "NA" } 
  else {
    f.g_args_txt <- paste(f.g_args, collapse=" ") 
  } 
  #Get c^Y fcn for parents of Y_i: (A,W)   
  cY_Pa <- .f.cY_Pa(W1, W2, W3, A, Net_vect)
  A_netA <- data.frame(.f.redefineCov(k, A, "A", Net_vect))
  Y <- .f.Y(W1, W2, W3, A, cY_Pa)
  #Convert N_i, cA_Pa_i, cY_Pa_i to strings
  # Net_str <- .f.mkstrNet(Net_vect)
  Net_str <- .f.mkstrNet(Net_vectID)
  netW1_str <- .f.mkstrNetWi(cA_Pa, 1)
  netW2_str <- .f.mkstrNetWi(cA_Pa, 2)
  netW3_str <- .f.mkstrNetWi(cA_Pa, 3)
  netA_str <- .f.mkstrNet(cY_Pa$Pa_As)

  # Flag for no risk of infection (no infected untreated partners)
  NoRskF <- sapply(c(1:n), function(i) sum(cY_Pa$Pa_Ws[[i]]$W2 * (1-cY_Pa$Pa_As[[i]]))==0)
  IDs <- paste('I', seq(n), sep='')
  # Add N(infected & untreated) and net1(inf&untrt), ..., netk(inf&untrt) 	
  d <- data.frame(IDs=IDs, C_j, f.g_args_txt, EC, Y, nFriends, W_netW, netW123_sum,
                  A_netA, NoRskF=NoRskF, Net_str, netW1_str, netW2_str, netW3_str,
                  netA_str, stringsAsFactors = FALSE)
  #***************************   	
  # 1) Set A=0 when W2=0 (exclude W2=0 from modelling g_0)
  # 2) Set Y=1 when W2=1 (exclude W2=1 from modelling Q_0)   	
  #***************************
  # d$g_deterministic <- (d$W2==0)
  # d$Q_deterministic <- (d$W2==1)
  return(d) 
  }
  
#Sample the population (several communities) and combine in one dt.frame
gendata_pop <- function(nC = 1, n_arr, k_arr, EC_arr, f.g_list=NULL, f.g_args_list=NULL, rndseed = NULL) {
  #n_arr - number of individuals /per community
  #k_arr - max size of each individual's network /per community
  #EC_arr - community-level covariate, only influences W_i for i\in nC_j /per community  
  #nC - # communities to sample
  if (!is.null(rndseed)) set.seed(rndseed, kind = "Mersenne-Twister", normal.kind = "Inversion")
  .f.EC <- function(i, samp=FALSE, EC_arr) if (samp) sample(EC_arr, 1) else EC_arr[C_j]   
  EC_rand <- sapply(1:nC, .f.EC, T, EC_arr)  
  #sample a random assignment of community covars
  if ((nC>1) & (!(is.null(f.g_list)))) {
    fcn_ids <- sample(1:length(f.g_list), nC, replace=TRUE)
    f.g_rand <- f.g_list[fcn_ids]
    f.g_args_rand <-f.g_args_list[fcn_ids] }
  else {
    if (!(is.null(f.g_list))) {
      f.g_rand <- f.g_list
      f.g_args_rand <-f.g_args_list }
    else {
      f.g_rand <- list(NULL)
      f.g_args_rand <- list(NULL) }
    }
  if (is.null(f.g_args_list)) f.g_args_rand<-list(NULL)  
  d <- gendata_Cj(C_j=1, n=n_arr, k=k_arr, EC=EC_rand, f.g_name=f.g_rand, f.g_args=f.g_args_rand)
  # d <- d[, (names(d) %in% c('IDs','W1','W2', 'W3', "netW1_sum", "netW2_sum", "netW3_sum", 'A','nFriends','Y','Net_str'))]
  d <- d[, c('IDs','W1','W2', 'W3', 'A','Y','nFriends','Net_str')]
  return(d)
}

# #----------------------------------------------------------------------------------
# #DEFINE REGRESSION FORMS FORM Q AND g
# #---------------------------------------------------------------------------------- 
# # W1 - inf.RISK
# # W2 - (0,1) infected at time t=0                
# #------------------------------------------------------------------------------   
# # DEFINE Qform & gform
# # Qform Q0 = b0 + b1*I(N(inf & untrt)=0) + b2*N(inf & untrt) + b3*nFriends
# # Qform Q0 = b0 + b1*I(sum(netW2 & (1-netA))=0) + b2*N(netW2 & (1-netA)) + 
# # b3*nFriends
# .f.Qform <- function(k, qform_Vars=NULL, k_miss=1, miss=F) {
#   qform_NetVars <- NULL
#   if (miss) k <- k_miss
#   if (k > 0) { 
#     qform_NetVarsW2A <- stringr::str_c(netvar("W2", (1:k)), "*", "(1-",netvar("A", (1:k)), ")")
#     qform_NetVarsW2A <- paste("I(", paste(qform_NetVarsW2A, collapse="+"), ")", sep="")
#     qform_NetVarsW3 <- netvar("W3", (1:k))
#     qform_NetVarsW3 <- paste("I(", paste(qform_NetVarsW3, collapse="+"), ")", sep="")
#     qform_NetVars <- paste(c(qform_NetVarsW2A, qform_NetVarsW3), collapse="+")
#   }
#   Qform <- paste("Y ~ ", qform_NetVars, collapse = "")
#   return(Qform)
# }
# # gform  g_0 = b0 + b1*W1 + b2*(sum(netW1))
# .f.gform <- function(k, k_miss=1, gform_Vars, miss=F) {
#   gform_NetVars <- NULL
#   if (miss) k <- k_miss
#   if (k > 0) {
#     gform_NetVars1 <- netvar("W1", (1:k))
#     gform_NetVars2 <- netvar("W2", (1:k))
#     gform_NetVars3 <- netvar("W3", (1:k))
#     if (!miss) gform_NetVars <- c(gform_NetVars1, gform_NetVars2, gform_NetVars3, "nFriends")
#   }
#   gform <- paste("A ~ ", paste(c(gform_Vars, gform_NetVars), collapse = " + "),collapse = "")
#   return(gform)
# }
# .f.gformv2 <- function(k, k_miss=1, gform_Vars, miss=F) {
#   gform_NetVars <- NULL
#   if (miss) k <- k_miss
#   if (k > 0) {
#     if (!miss) gform_NetVars <- c("netW1_sum", "netW2_sum", "netW3_sum", "nFriends")
#   }
#   gform <- paste("A ~ ", paste(c(gform_Vars, gform_NetVars),collapse = " + "),collapse = "")
#   return(gform)
# }
