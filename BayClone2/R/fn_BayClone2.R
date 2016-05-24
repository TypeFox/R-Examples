#### function to generate N_tot and n_alt files
export_N_n <- function(vcf_in_file,N_tot_out_file,n_alt_out_file)
{
	##dyn.load("../src/BayClone2.so")
	.C("export_N_n",vcf_in_file,N_tot_out_file,n_alt_out_file);
}


BayClone2 <- function(min_C, max_C, SS, TT, Burn.in, N.sam, NN, nn, hpara, ave.B)
{
    ###########################
    #SPLITTING TRAINING AND TESTING
    ##FOR DETIALS, SEE THE REFERENCE LEE AT EL (2014)
    ##B: FRACTION USED TO SPLIT DATA (B IS THE PROPORTION FOR A TRAINING DATASET)
    ##FROM OUR EXPERIENCE, 2.5% WORKS FINE. HOWEVER!!!!! YOU MAY NEED TO CHANGE THIS ACCORDING TO YOUR DATA!!!!!!!!!!!!!!!!!!!!!
    aa <- 1000*ave.B
    bb <- 1000 - aa
    B <- array(rbeta(SS*TT, aa, bb), dim=c(SS, TT))  ##FRACTION USED TO SPLIT DATA (B IS THE PROPORTION FOR A TRAINING DATASET)
    #################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    n_Tr <- nn*B  ##n FOR TRAINING
    n_Te <- nn - n_Tr  ##n FOR TESTING
    
    N_Tr <- NN*B  ##N FOR TRAINING
    N_Te <- NN - N_Tr  ##N FOR TESTING
    
    ##PUT THE DATASETS IN A LONG SEQUENCE INSTEAD OF ARRAYS.
    n_Tr.tmp <- array(n_Tr, dim=c(1, SS*TT))[1,]
    N_Tr.tmp <- array(N_Tr, dim=c(1, SS*TT))[1,]
    N_Tr_sum <- apply(N_Tr, 2, sum)
    
    n_Te.tmp <- array(n_Te, dim=c(1, SS*TT))[1,]
    N_Te.tmp <- array(N_Te, dim=c(1, SS*TT))[1,]
    N_Te_sum <- apply(N_Te, 2, sum)
    
    B_Tr.tmp <- array(B, dim=c(1, SS*TT))[1,]
    B_Te.tmp <- array(1-B, dim=c(1, SS*TT))[1,]
    
    
    n.tmp <- array(nn, dim=c(1, SS*TT))[1,]
    N.tmp <- array(NN, dim=c(1, SS*TT))[1,]
    N_sum <- apply(NN, 2, sum)  ####COLUMN SUM OF N
    
    
    
    ###########################################################################
    ###DO MCMC FOR TRAINING DATASET
    ###FOR EACH VALUE OF C, SAMPLE THE OTHER PARAMETERS INCLUDING L, Z, W, PHI, PI, P0 USING TRAINING DATA FOR BURN-IN ITERATIONS
    print("Doing MCMC sampling for training dataset")
    sam.all.OR <- fn.Tr.CNV(hpara, SS, TT, N_Tr.tmp, n_Tr.tmp, N_Tr_sum, B_Tr.tmp, min_C, max_C, Burn.in)
    
    
    ###NOW WE DO TRANS-DIMENSIONAL MCMC (THAT IS, SAMPLE C AS WELL AS THE OTHER PARAMETERS) USING THE TEST DATASET.
    set.seed(115516131)
    print("Doing MCMC sampling using training and test datasets")
    print(date())
    MCMC.sam <- fn.MCMC.CNV(sam.all.OR, hpara, SS, TT, N_Tr.tmp, n_Tr.tmp, N_Tr_sum, B_Tr.tmp, N_Te.tmp, n_Te.tmp, N_Te_sum, B_Te.tmp, N.tmp, n.tmp, N_sum, min_C, max_C, N.sam)
    
    return(MCMC.sam)
}


###DO MCMC FOR EACH VALUE OF C BEWTEEN MIN_C AND MAX_C USING TRANING DATA SET ONLY.
###FOR EACH VALUE OF C SAVE THE POSTERIOR SAMPLES OF THE OTHER PARAMETERS FROM THE POSTERIOR BASED ON THE TRAINING ONLY.
fn.Tr.CNV <- function(hpara, SS, TT, NN_Tr, nn_Tr, NN_Tr_sum, BB_Tr, min_C, max_C, N_sam)
{
    #HYPERPARAMETERS
    a <- hpara$a; b <- hpara$b;  #PHI
    d0 <- hpara$d0; d <- hpara$d  #W
    
    #PRIOR FOR P_zO
    a_z0 <- hpara$a_z0; b_z0 <- hpara$b_z0 # FOR P_Z0
    
    #FOR L
    Q <- hpara$Q;
    alpha <- hpara$alpha
    beta <- hpara$beta
    gam <- hpara$gam
    
    #SAVE THE PROPOSALS FOR ALL C VALUES
    sam.all <- NULL
    
    sam.all$phi <- NULL
    sam.all$pi <- NULL
    
    sam.all$L <- NULL
    sam.all$Z <- NULL
    sam.all$A <- NULL
    
    sam.all$th <- NULL
    sam.all$p0_z <- rep(NA, max_C-min_C+1)
    sam.all$w <- NULL
    
    sam.all$M <- NULL
    sam.all$p <- NULL

    #INITIALIZATION FOR EACH C
    for(i.c in min_C:max_C)
    {
        print(date())
        print(paste("working on i.c=", i.c))
        set.seed(545754)
        
        ##L AND Z
        Z <- L <- array(0, dim=c(SS, i.c))
        
         #C--TO INITIALIZE
        p.ave <- apply(array(nn_Tr/NN_Tr, dim=c(SS, TT)), 1, mean)
        p.quant <- quantile(p.ave, probs=seq(0.10, 0.95, length.out=(i.c-1)))
        
        Z[,1] <- L[,1] <- 2 #CELL TYPE 0
        
        for(ii.c in 2:i.c) #EXCEPT CELL TYPE 0
        {
            Z[p.ave > p.quant[ii.c-1],ii.c] <- 2
        }

        for(ii.c in 2:i.c)
        for(i.s in 1:SS)
        {
            L[i.s, ii.c] <- sample((Z[i.s, ii.c]:Q), 1)
        }
        
        #COMPUTE A
        A <- array(NA, dim=c(i.c, (Q+1)))
        for(i.q in (0:Q))
        {
            A[,(i.q+1)] <- apply(L==i.q, 2, sum)
        }
        
        #PI
        ppi <- (A+5)/(SS + (Q+1)*5)
        
        #P0_Z
        p_z0 <- c(0.005, rep(1, i.c-1))

        #PHI
        phi <- rgamma(TT, a, b)
        
        
        #W_STAR AND W
        w_star <- array(NA, dim=c(TT, i.c))  #UN-NORMLAIZED WEIGHTS
        for(i_t in 1:TT) w_star[i_t,] <- rgamma(i.c, c(d0, rep(d, i.c-1)), 1)
        w <- w_star/t(array(rep(apply(w_star, 1, sum),each=i.c), dim=c(i.c, TT)))
        
        
        #COMPUTE M AND P
        M <- array(NA, dim=c(SS, TT))
        p <- array(NA, dim=c(SS, TT))
        for(i_s in 1:SS)
        {
            for(i_t in 1:TT)
            {
                M[i_s, i_t] <- sum(w[i_t,]*L[i_s,])
                p[i_s, i_t] <- sum(w[i_t,]*Z[i_s,]*p_z0)/M[i_s,i_t]
            }
        }
        
        #COPY CURRENT SAMPLE IN A DIFFERENT FORMAT
        A <- array(A, dim=c(1, i.c*(Q+1)))[1,]
        L <- array(L, dim=c(1, SS*i.c))[1,]
        Z <- array(Z, dim=c(1, SS*i.c))[1,]
        M <- array(M, dim=c(1, SS*TT))[1,]
        p <- array(p, dim=c(1, SS*TT))[1,]
        w <- array(w, dim=c(1, i.c*TT))[1,]
        w_star <- array(w_star, dim=c(1, TT*i.c))[1,]
        ppi <- array(ppi, dim=c(1, i.c*(Q+1)))[1,]

        ##CALL C FUNCTION TO DO MCMC FOR A GIVEN VALUE OF C
        output <- .C("fn_CNV_MCMC_1", alpha=as.double(alpha/(i.c-1)), beta=as.double(beta), gam=as.double(gam), QQ=as.integer(Q), a=as.double(a), b=as.double(b), dd0=as.double(d0), dd=as.double(d), aa0_z=as.double(a_z0), bb0_z=as.double(b_z0), CC=as.integer(i.c), ppi=as.double(ppi), phi=as.double(phi), L=as.double(L), Z=as.double(Z), A=as.double(A), p0_z=as.double(p_z0[1]), th=as.double(w_star), w=as.double(w), M=as.double(M), p=as.double(p), SS=as.integer(SS), TT=as.integer(TT), n=as.double(nn_Tr), N=as.double(NN_Tr), N_sum=as.double(NN_Tr_sum), BB=as.double(BB_Tr), NN_iter=as.integer(N_sam))
        
        sam.all$phi[[i.c-min_C+1]] <- output$phi
        sam.all$pi[[i.c-min_C+1]] <- output$ppi
        
        sam.all$L[[i.c-min_C+1]] <- output$L
        sam.all$Z[[i.c-min_C+1]] <- output$Z
        sam.all$A[[i.c-min_C+1]] <- output$A

        sam.all$th[[i.c-min_C+1]] <- output$th
        sam.all$p0_z[i.c-min_C+1] <- output$p0_z
        sam.all$w[[i.c-min_C+1]] <- output$w

        sam.all$M[[i.c-min_C+1]] <- output$M
        sam.all$p[[i.c-min_C+1]] <- output$p
    }

    return(sam.all)
}

###DO TRANS-DIMENSIONAL MCMC (MCMC INCLUDING UPDATE OF C)
####UPDATE PROPOSALS FROM THE TRAINING DATASET WHEN SAMPLING C
####FOR THE ACCEPTED C, UPDATE THE OTHER PARAMETERS USING THE FULL DATASET
fn.MCMC.CNV <- function(sam.all, hpara, SS, TT, NN_Tr, nn_Tr, NN_Tr_sum, BB_Tr, NN_Te, nn_Te, NN_Te_sum, BB_Te, NN, nn, NN_sum, min_C, max_C, N_sam)
{
    ############################################################################################################
    #HYPERPARAMETERS
    r <- hpara$r
    
    a <- hpara$a; b <- hpara$b;  #PHI
    d0 <- hpara$d0; d <- hpara$d  #W
    
    #PRIOR FOR P_zO
    a_l0 <- hpara$a_l0; b_l0 <- hpara$b_l0 # FOR P_L0
    a_z0 <- hpara$a_z0; b_z0 <- hpara$b_z0 # FOR P_Z0
    
    #FOR L
    Q <- hpara$Q;
    alpha <- hpara$alpha
    beta <- hpara$beta
    gam <- hpara$gam
    
    ############################################################################################################
    #SAVE MCMC SAMPLES
    sam <- NULL
    sam$C <- rep(0, N_sam)
    
    sam$L <- NULL
    sam$Z <- NULL
    
    sam$w <- NULL
    sam$th <- NULL
    
    sam$phi <- array(NA, dim=c(N_sam, TT))
    sam$pi  <- NULL
    
    sam$p0_z <- rep(0, N_sam)
    
    sam$M <- array(NA, dim=c(N_sam, SS, TT))
    sam$p <- array(NA, dim=c(N_sam, SS, TT))
    
    #sam$rec <- array(NA, dim=c(N_sam, 3))
    ############################################################################################################
    C_cur <- sample((min_C:max_C), 1)  ##CURRENT VALUE OF C
    for(i_sam in 1:N_sam)
    {
        if((i_sam%%100)==0)
        {
            print(paste(round(i_sam/N_sam*100,2), "% MCMC sampling has been done." ))
            print(date())
        }
        
        ############################################################################################################
        ##UPDATE (C, X) THROUGH M-H###############################################################
        ###CURRENT C_CUR
        output_cur <- .C("fn_CNV_MCMC_2", alpha=as.double(alpha/(C_cur-1)), beta=as.double(beta), gam=as.double(gam), QQ=as.integer(Q), a=as.double(a), b=as.double(b), dd0=as.double(d0), dd=as.double(d), aa0_z=as.double(a_z0), bb0_z=as.double(b_z0), CC=as.integer(C_cur), ppi=as.double(sam.all$pi[[C_cur-min_C+1]]), phi=as.double(sam.all$phi[[C_cur-min_C+1]]), L=as.double(sam.all$L[[C_cur-min_C+1]]), Z=as.double(sam.all$Z[[C_cur-min_C+1]]), A=as.double(sam.all$A[[C_cur-min_C+1]]), p0_z=as.double(sam.all$p0_z[[C_cur-min_C+1]]), th=as.double(sam.all$th[[C_cur-min_C+1]]), w=as.double(sam.all$w[[C_cur-min_C+1]]), M=as.double(sam.all$M[[C_cur-min_C+1]]), p=as.double(sam.all$p[[C_cur-min_C+1]]), SS=as.integer(SS), TT=as.integer(TT), n_tr=as.double(nn_Tr), N_tr=as.double(NN_Tr), N_tr_sum=as.double(NN_Tr_sum), BB_tr=as.double(BB_Tr), n_te=as.double(nn_Te), N_te=as.double(NN_Te), N_te_sum=as.double(NN_Te_sum),  BB_te=as.double(BB_Te), loglike=as.double(0), NN_iter=as.integer(0))
        alpha_cur <- output_cur$loglike  - (C_cur-2)*log(1-r)
        
        ###CURRENT C_PRO
        C_pro <- sample((min_C:max_C), 1)
        output_pro <- .C("fn_CNV_MCMC_2", alpha=as.double(alpha/(C_pro-1)), beta=as.double(beta), gam=as.double(gam), QQ=as.integer(Q), a=as.double(a), b=as.double(b), dd0=as.double(d0), dd=as.double(d), aa0_z=as.double(a_z0), bb0_z=as.double(b_z0), CC=as.integer(C_pro), ppi=as.double(sam.all$pi[[C_pro-min_C+1]]), phi=as.double(sam.all$phi[[C_pro-min_C+1]]), L=as.double(sam.all$L[[C_pro-min_C+1]]), Z=as.double(sam.all$Z[[C_pro-min_C+1]]), A=as.double(sam.all$A[[C_pro-min_C+1]]), p0_z=as.double(sam.all$p0_z[[C_pro-min_C+1]]), th=as.double(sam.all$th[[C_pro-min_C+1]]), w=as.double(sam.all$w[[C_pro-min_C+1]]), M=as.double(sam.all$M[[C_pro-min_C+1]]), p=as.double(sam.all$p[[C_pro-min_C+1]]), SS=as.integer(SS), TT=as.integer(TT), n_tr=as.double(nn_Tr), N_tr=as.double(NN_Tr), N_tr_sum=as.double(NN_Tr_sum), BB_tr=as.double(BB_Tr), n_te=as.double(nn_Te), N_te=as.double(NN_Te), N_te_sum=as.double(NN_Te_sum), BB_te=as.double(BB_Te), loglike=as.double(0), NN_iter=as.integer(1))
        alpha_pro <- output_pro$loglike  - (C_pro-2)*log(1-r)
        
        sam.all$phi[[C_pro-min_C+1]] <- output_pro$phi
        sam.all$pi[[C_pro-min_C+1]] <- output_pro$ppi
        
        sam.all$L[[C_pro-min_C+1]] <- output_pro$L
        sam.all$Z[[C_pro-min_C+1]] <- output_pro$Z
        sam.all$A[[C_pro-min_C+1]] <- output_pro$A
        
        sam.all$th[[C_pro-min_C+1]] <- output_pro$th
        sam.all$p0_z[C_pro-min_C+1] <- output_pro$p0_z
        sam.all$w[[C_pro-min_C+1]] <- output_pro$w
        
        sam.all$M[[C_pro-min_C+1]] <- output_pro$M
        sam.all$p[[C_pro-min_C+1]] <- output_pro$p
        
        ##ACCEPT??
        max_alpha <- max(alpha_cur, alpha_pro)
        alpha_cur <- exp(alpha_cur - max_alpha)
        alpha_pro <- exp(alpha_pro - max_alpha)
        
        #sam$rec[i_sam,] <- c(C_cur, C_pro, alpha_cur/(alpha_cur+alpha_pro))
        
        if(runif(1) < alpha_pro/(alpha_cur + alpha_pro))
        {
            output_cur <- output_pro
            C_cur <- C_pro
        }
        
        ############################################################################################################
        ##UPDATE X (THE OTHER PARAMETERS BUT C) USING ALL DATA
        output_cur <- .C("fn_CNV_MCMC_1", alpha=as.double(alpha/(C_cur-1)), beta=as.double(beta), gam=as.double(gam), QQ=as.integer(Q), a=as.double(a), b=as.double(b), dd0=as.double(d0), dd=as.double(d), aa0_z=as.double(a_z0), bb0_z=as.double(b_z0), CC=as.integer(C_cur), ppi=as.double(output_cur$ppi), phi=as.double(output_cur$phi), L=as.double(output_cur$L), Z=as.double(output_cur$Z), A=as.double(output_cur$A), p0_z=as.double(output_cur$p0_z), th=as.double(output_cur$th), w=as.double(output_cur$w), M=as.double(output_cur$M), p=as.double(output_cur$p), SS=as.integer(SS), TT=as.integer(TT), n=as.double(nn), N=as.double(NN), N_sum=as.double(NN_sum), BB=as.double(rep(1, SS*TT)), NN_iter=as.integer(10))
        
        
        ############################################################################################################
        ####SAVE SAMPLES
        sam$C[i_sam] <- C_cur
        sam$L[[i_sam]] <- array(output_cur$L, dim=c(SS, C_cur))
        sam$Z[[i_sam]] <- array(output_cur$Z, dim=c(SS, C_cur))
        
        sam$w[[i_sam]] <- array(output_cur$w, dim=c(TT, C_cur))
        sam$th[[i_sam]] <- array(output_cur$th, dim=c(TT, C_cur))
        
        sam$phi[i_sam,] <- output_cur$phi
        sam$pi[[i_sam]] <- array(output_cur$ppi, dim=c(C_cur, (Q+1)))
        
        sam$p0_z[i_sam] <- output_cur$p0_z
        
        sam$M[i_sam,,] <- array(output_cur$M, dim=c(SS, TT))
        sam$p[i_sam,,] <- array(output_cur$p, dim=c(SS, TT))
    }
    return(sam)
}


#COMPUTE THE POSTERIOR MARIGNAL DIST OF C
fn_post_C <- function(C_sam, min_C, max_C)
{
    if((min_C >= 2)&(max_C >= min_C))
    {
        ##NUMBER OF MCMC SAMPLES
        n_sam <- length(C_sam)
        post_dist <- array(0, dim=c(max_C, 2))
        post_dist[,1] <- (1:max_C)
        
        for(i.c in min_C:max_C) post_dist[i.c, 2] <- sum(C_sam==i.c)/n_sam
        
        colnames(post_dist) <- c("C", "Post_Prob")
        
        return(post_dist)
    }else{
        if(min_C < 2) print("Min of C should be greater than or equal to 2!")
        if(min_C > max_C) print("Min of C should be less than or equal to Max of C!")
    }
}


##FIND POSTERIOR POINT ESTIMATES OF L, Z, W, PHI, P0 FOR A GIVE C (HERE, C SHOULD BE LESS THAN OR EQUAL TO 10 DUE TO THE PERMUTATION)
fn_posterior_point <- function(CC, SS, TT, sam)
{
    #require(combinat)
    
    ##NUMBER OF MCMC SAMPLES
    n.sam <- length(sam$C)
    
    ###################################################################
    n.C <- sum(sam$C==CC)
    if((n.C > 2)&(CC <= 10))
    {
        ###COMPUTE DISTANCE AMONG L (SAMPLED THROUGH MCMC)
        dist.C <- NULL
        sam.C <- (1:n.sam)[sam$C==CC]
        D.C <- array(0, dim=c(n.C, n.C))
        
        ind_N <- factorial(CC-1)
        index_set <- array(NA, dim=c(ind_N, CC-1))
        tmp <- permn(CC-1)
        
        
        for(i.n in 1:ind_N)
        {
            index_set[i.n,] <- tmp[[i.n]]
        }
        
        index_set <- array(index_set-1, dim=c(1, ind_N*(CC-1)))[1,]
        
        print(paste("SAMPLE SIZE FOR THE GIVEN C =", n.C))
        
        for(i in 1:(n.C-1))
        {
            if((i%%100)==0)
            {
                print(paste(round(i/n.C*100, 2), "% has been done"))
            }
            
            for(j in (i+1):n.C)
            {
                i.sam <- sam.C[i]
                j.sam <- sam.C[j]
                
                L.i <- array(sam$L[[i.sam]], dim=c(1, SS*CC))[1,]
                L.j <- array(sam$L[[j.sam]], dim=c(1, SS*CC))[1,]
                
                output <- .C("fn_sum_Z_1", SS=as.integer(SS), CC=as.integer(CC-1), ind_NN=as.integer(ind_N), ind_set=as.integer(index_set), Z_1=as.integer(L.i), Z_2=as.integer(L.j), min_dist=as.double(0))
                
                D.C[i, j] <- D.C[j, i] <- output$min_dist
            }
        }
        
        dist.C[[CC]] <- D.C
        #################################################################################################################
        ###FIND THE SAMPLE OF L MINIMIZING THE DISTANCE FROM ALL THE OTHER SAMPLES OF L
        dist.tmp <- dist.C[[CC]]
        sum.dist <- apply(dist.tmp, 1, sum)
        ind.min <- which.min(sum.dist)
        ind.min <- sam.C[ind.min]
        
        #####################################################
        post.L <- sam$L[[ind.min]]  ### POSTERIOR POINT ESTIMATE OF L, L_HAT
        post.Z <- sam$Z[[ind.min]]  ### POSTERIOR POINT ESTIMATE OF Z, Z_HAT
        post.w <- sam$w[[ind.min]]  ### POSTERIOR POINT ESTIMATE OF W, W_HAT
        
        
        #####################################################
        sam.C <- (1:n.sam)[sam$C==CC]
        post.p0.z <- mean(sam$p0_z[sam.C])  ### POSTERIOR POINT ESTIMATE OF P0, P0_HAT
        post.phi <- apply(sam$phi[sam.C,], 2, mean)  ### POSTERIOR POINT ESTIMATE OF PHI, PHI_HAT
        
        #####################################################
        ####COMPUTE ESTIMATES OF M_ST AND P_ST CONDITIONAL ON C_HAT
        #####################################################
        p_st_post <- M_st_post <- array(0, dim=c(SS, TT))
        
        for(i.sam in sam.C)
        {
            phi_tmp <- sam$phi[i.sam,]
            M_tmp <- sam$M[i.sam,,]
            p_tmp <- sam$p[i.sam,,]
            
            M_st_post <- M_st_post + M_tmp
            p_st_post <- p_st_post + p_tmp
        }
        
        M_st_post <- M_st_post/length(sam.C)
        p_st_post <- p_st_post/length(sam.C)
        
        
        
        post.point.C <- NULL
        post.point.C$C <- CC
        post.point.C$L <- post.L
        post.point.C$Z <- post.Z
        post.point.C$w <- post.w
        post.point.C$p0 <- post.p0.z
        post.point.C$phi <- post.phi
        post.point.C$M <- M_st_post
        post.point.C$p <- p_st_post
        
        return(post.point.C)
    }else{
        if(n.C <= 2)print("ERROR: Not Enough Sample For Chosen C!")
        if(CC > 10) print("ERROR: The Chosen C Is Too Large!")
    }
    
}




