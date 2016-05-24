cooccur <-
function(mat, 
                    type = "spp_site", 
                    thresh = TRUE, 
                    spp_names = FALSE,
                    true_rand_classifier = 0.1,
                    prob = "hyper",
                    site_mask = NULL, ## NEW ADDITION
                    only_effects=FALSE,
                    eff_standard=TRUE,
                    eff_matrix=FALSE){
  
  # HANDLE ARGUEMENTS
    if (type == "spp_site"){spp_site_mat <- mat}
    if (type == "site_spp"){spp_site_mat <- t(mat)}
    if (spp_names == TRUE){
      spp_key <- data.frame(num=1:nrow(spp_site_mat),spp=row.names(spp_site_mat))  
    }
    if (!is.null(site_mask)){
      #if (nrow(site_mask)==nrow(spp_site_mat) & ncol(site_mask)==nrow(spp_site_mat)){
      #}else{
        if (nrow(site_mask)==nrow(spp_site_mat) & ncol(site_mask)==ncol(spp_site_mat)){
          N_matrix <- create.N.matrix(site_mask)
        }else{
          stop("Incorrect dimensions for site_mask, aborting.")
        }
      #}
    }else{
      site_mask <- matrix(data = 1,nrow = nrow(spp_site_mat),ncol = ncol(spp_site_mat)) 
      N_matrix <- matrix(data = ncol(spp_site_mat),nrow = nrow(spp_site_mat),ncol = nrow(spp_site_mat)) # needs doing, add nsite to all
    }
    
  # ORGANIZE & INITIALIZE FOR ANALYSIS 
    spp_site_mat[spp_site_mat>0] <- 1
    tsites <- ncol(spp_site_mat)
    nspp <- nrow(spp_site_mat)
    spp_pairs <- choose(nspp,2)
    #incidence <- prob_occur <- matrix(nrow=nspp,ncol=2) # needs expanding, pieces moved below
    incidence <- prob_occur <- obs_cooccur <- prob_cooccur <- exp_cooccur <- matrix(nrow=spp_pairs,ncol=3) # original obs_cooccur <- prob_cooccur <- exp_cooccur <- matrix(nrow=spp_pairs,ncol=3)
    # ORIGINAL: prob_share_site <- c(0:(nsite+1)) # NEW VERSION WILL BE LIKE sapply(nsite + 1, FUN = function(x){c(0:x)}, simplify = F)
    
  # CALCULATE PROBABILITIES
  
    incidence <- prob_occur <- matrix(nrow = nrow(N_matrix),ncol = ncol(N_matrix))
  
    # ORIGINAL incidence <- prob_occur <- cbind(c(1:nrow(spp_site_mat)),rowSums(spp_site_mat,na.rm=T)) # this must me expanded to allow pair specific incidence
    # prob_occur <- cbind(c(1:nrow(spp_site_mat)),rowSums(spp_site_mat,na.rm=T)/nsite) # ORIGINAL, now prob_occur moved up so nsite div later
    
    for (spp in 1:nspp){
      if (spp < nspp){
        for (spp_next in (spp + 1):nspp){
          incidence[spp,spp_next] <- sum(site_mask[spp,]*site_mask[spp_next,]*mat[spp,])
          incidence[spp_next,spp] <- sum(site_mask[spp,]*site_mask[spp_next,]*mat[spp_next,])
        }
      }
    }
  
    prob_occur <- incidence/N_matrix
  
    pb <- txtProgressBar(min = 0, max = (nspp + nrow(obs_cooccur)), style = 3)
    
    row <- 0
    for (spp in 1:nspp){
      if (spp < nspp){
        for (spp_next in (spp + 1):nspp){
          
          pairs <- sum(as.numeric(mat[spp,site_mask[spp,]*site_mask[spp_next,]==1]==1&mat[spp_next,site_mask[spp,]*site_mask[spp_next,]==1]==1))
          row <- row + 1
          #pairs <- 0
          #for (site in 1:tsites){
          #  if (spp_site_mat[spp,site] > 0 & spp_site_mat[spp_next,site] > 0){
          #    pairs <- pairs + 1
          #  }
          #}
          
          obs_cooccur[row,1] <- spp
          obs_cooccur[row,2] <- spp_next
          obs_cooccur[row,3] <- pairs
        
          prob_cooccur[row,1] <- spp
          prob_cooccur[row,2] <- spp_next
          prob_cooccur[row,3] <- prob_occur[spp,spp_next] * prob_occur[spp_next,spp] # ORIGINAL prob_occur[spp,2] * prob_occur[spp_next,2]
          
          exp_cooccur[row,1] <- spp
          exp_cooccur[row,2] <- spp_next
          exp_cooccur[row,3] <- prob_cooccur[row,3] * N_matrix[spp,spp_next] # original nsite
        
        }
      }
      
      setTxtProgressBar(pb, spp)
      
    }
    
  # SHOULD WE APPLY A THRESHOLD FOR EXPECTED COOCCURs
        
    if (thresh == TRUE){
      n_pairs <- nrow(prob_cooccur)
      prob_cooccur <- prob_cooccur[exp_cooccur[,3]>=1,]
      obs_cooccur <- obs_cooccur[exp_cooccur[,3]>=1,]
      exp_cooccur <- exp_cooccur[exp_cooccur[,3]>=1,]
      n_omitted <- n_pairs - nrow(prob_cooccur)
      pb <- txtProgressBar(min = 0, max = (nspp + nrow(obs_cooccur)), style = 3)
    }
    
  # PREPARE OUTPUT TABLE
    
    output <- data.frame(matrix(nrow=0,ncol=9))
    colnames(output) <- c("sp1","sp2","sp1_inc","sp2_inc","obs_cooccur","prob_cooccur","exp_cooccur","p_lt","p_gt")
    
  # COMBINATORICS
        
    for (row in 1:nrow(obs_cooccur)){
      
      sp1 <- obs_cooccur[row,1]
      sp2 <- obs_cooccur[row,2]
      
      sp1_inc <- incidence[sp1,sp2] # ORGIGINAL incidence[incidence[,1]==sp1,2]
      sp2_inc <- incidence[sp2,sp1] # ORGIGINAL incidence[incidence[,1]==sp2,2]
        
      max_inc <- max(sp1_inc,sp2_inc)
      min_inc <- min(sp1_inc,sp2_inc)
      
      nsite <- N_matrix[sp1,sp2] #ADDITION
      
      psite <- as.numeric(nsite+1)
      prob_share_site <- rep(x = 0,times = psite) 
      
  # CHOOSE TO CALCULATE PROBABILITIES USING THE
  # HYPERGEOMETRIC DISTRIBUTION OR VEECH 2013 COMBINATORICS
      
      ### Hypergeometric Implementation using phyper ###

          if (prob == "hyper"){
              if (only_effects == FALSE){
        
                all.probs <- phyper(0:min_inc, min_inc, nsite - min_inc, max_inc)
                prob_share_site[1]<-all.probs[1]
                for(j in 2:length(all.probs)){
                  prob_share_site[j] <- all.probs[j]-all.probs[j-1]
                }
              }else{
                for (j in 0:nsite){
                  if ((sp1_inc + sp2_inc) <= (nsite +j) ){
                    if (j <= min_inc){
                      
                      prob_share_site[(j+1)] <- 1
                      
                    }
                  }
                }
              }
                
          }
      
      ### Combinatorics Implementation from Veech 2013

          if (prob == "comb"){
              if (only_effects == FALSE){
                for (j in 0:nsite){
                  if ((sp1_inc + sp2_inc) <= (nsite +j) ){
                    if (j <= min_inc){
                      
                      prob_share_site[(j+1)] <- coprob(max_inc=max_inc,j=j,min_inc=min_inc,nsite=nsite) 
                      
                    }
                  }
                }
              }else{
                for (j in 0:nsite){
                  if ((sp1_inc + sp2_inc) <= (nsite +j) ){
                    if (j <= min_inc){
                      
                      prob_share_site[(j+1)] <- 1
                      
                    }
                  }
                }
              }    
          }
      
      ### END OF PROBABILITY CALCULATIONS
      
      p_lt <- 0
      p_gt <- 0
      
      for (j in 0:nsite){
        
        if (j <= obs_cooccur[row,3]){
          p_lt <- prob_share_site[(j+1)] + p_lt}
        if (j >= obs_cooccur[row,3]){
          p_gt <- prob_share_site[(j+1)] + p_gt}
        if (j == obs_cooccur[row,3]){
          p_exactly_obs <- prob_share_site[(j+1)]
        
        }
      }
      
      p_lt <- round(p_lt,5)
      p_gt <- round(p_gt,5)
      p_exactly_obs <- round(p_exactly_obs,5)
      prob_cooccur[row,3] <- round(prob_cooccur[row,3],3)
      exp_cooccur[row,3] <- round(exp_cooccur[row,3],1)
        
      output[row,] <- c(sp1,sp2,sp1_inc,sp2_inc,obs_cooccur[row,3],prob_cooccur[row,3],exp_cooccur[row,3],p_lt,p_gt)
     
      setTxtProgressBar(pb, nspp + row)
  
    }
        
    close(pb)
    
    if (spp_names == TRUE){
      
      sp1_name <- merge(x=data.frame(order=1:length(output$sp1),sp1=output$sp1),y=spp_key,by.x="sp1",by.y="num",all.x=T,sort=FALSE)
      sp2_name <- merge(x=data.frame(order=1:length(output$sp2),sp2=output$sp2),y=spp_key,by.x="sp2",by.y="num",all.x=T,sort=FALSE)
      
      output$sp1_name <- sp1_name[with(sp1_name,order(order)),"spp"]
      output$sp2_name <- sp2_name[with(sp2_name,order(order)),"spp"]
      
    } 
    
  # clasifying true randoms
    true_rand <- (nrow(output[(output$p_gt >= 0.05 & output$p_lt >= 0.05) & (abs(output$obs_cooccur - output$exp_cooccur) <= (tsites * true_rand_classifier)),]))
  
  # PREPARE AND DELIVER OUTPUT AS CLASS "cooccur"
    output_list <- list(call = match.call(),
                        results = output,
                        positive = nrow(output[output$p_gt < 0.05,]),
                        negative = nrow(output[output$p_lt < 0.05,]),
                        co_occurrences = (nrow(output[output$p_gt < 0.05 | output$p_lt < 0.05,])), #nrow(output[output$p_exactly_obs <= 0.05,]),
                        pairs = nrow(output),
                        random = true_rand,
                        unclassifiable = nrow(output) - (true_rand + nrow(output[output$p_gt < 0.05,]) + nrow(output[output$p_lt < 0.05,])),
                        sites = N_matrix,
                        species = nspp,
                        percent_sig = (((nrow(output[output$p_gt < 0.05 | output$p_lt < 0.05,]))) / (nrow(output))) * 100,
                        true_rand_classifier = true_rand_classifier
                        )
    
    if (spp_names == TRUE){
      output_list$spp_key <- spp_key
      output_list$spp.names = row.names(spp_site_mat)
    }else{
      output_list$spp.names = c(1:nrow(spp_site_mat))
    }
    if (thresh == TRUE){
      output_list$omitted <- n_omitted
      output_list$pot_pairs <- n_pairs
    }
    
    class(output_list) <- "cooccur"
  
    if (only_effects == F){
      output_list
    }else{
      effect.sizes(mod=output_list,standardized=eff_standard,matrix=eff_matrix)
    }
}
