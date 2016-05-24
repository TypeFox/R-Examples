# #**********************************************************************
# #**********************************************************************
# #*************           FCT3 EM_GR                    ****************
# #**********************************************************************
# #**********************************************************************
# #
# ###### Summary #######################
# #> Segmentation Function
# # fctP1 - sfMM : fMM segmentation function
# 
# #> Simulation Functions 
# # fctS1 - simul_fMM : simulation function for non spatial mixture model (including concomitant model)
# 
# #> Internal function for the fMM 
# # fctI1 - init1_fMM : pre-processing function of the fMMalgo arguments
# # fctI2 - test_fMM : test function of the fMMalgo arguments
# # fctI3 - init2_fMM : initialisation the EM algorithm
# 
# #*******************
# 
# #### 1- Segmentation function ####
# 
# sfMM <- function(M, G, data, Id = NULL, Wweight_SR = NULL, Wdist_LR = NULL, posterior_proba = NULL, 
#          formula_reg, offset_reg = NULL, family = gaussian(link = "identity"), prior_theta = "kmeans", prior_sigma = NULL, 
#          formula_group = NULL, prior_prevalence = FALSE, prior_proba = NULL, 
#          test.GR = FALSE, sigma_GR = "auto", proba_GR = NULL, proba_GRseed = NULL, seed = NULL, n.order = 3, 
#          test.ICM = FALSE, rho_ICM = "init", G_ICM = 1:G, prior_prevalenceICM = TRUE, rho_max = 10, update_rho = FALSE,                            
#          test.ICMregional = FALSE, coords = NULL, nbGroup_min = 100, threshold_potential = 0.1, distance_ref = NULL, multiV = FALSE, 
#          digit.pmin = 7, epsilon = 10^(-3), epsilon_corrSpat = epsilon * 10, iter_max = 100, 
#          verbose = TRUE, trace_iter = FALSE, trace_radius = FALSE, export.predicteur = FALSE)
# 		 
# { t0 <- date()
#   
#   initPackage("nnet", method = "sfMM")
#  
#   #### mise en forme des parametres ####
#   res.tempo <- init1_fMM(M = M, G = G, data = data, Id = Id, coords = coords, 
#                                   formula_reg = formula_reg, prior_theta = prior_theta, prior_sigma = prior_sigma, family = family, offset_reg = offset_reg, 
#                                   formula_group = formula_group, 
#                                   test.GR = test.GR, sigma_GR = sigma_GR, proba_GR = proba_GR, proba_GRseed = proba_GRseed, 
#                                   test.ICM = test.ICM, rho_ICM = rho_ICM, G_ICM = G_ICM, 
#                                   test.ICMregional = test.ICMregional)
#   
#   # chargement
#   n <- res.tempo$n
#   coords <- res.tempo$coords
#   formula_reg <- res.tempo$formula_reg
#   Var_reg <- res.tempo$Var_reg
#   character_reg <- res.tempo$character_reg
#   offset_reg <- res.tempo$offset_reg
#   Var_group <- res.tempo$Var_group
#   formula_group <- res.tempo$formula_group
#   if(!is.null(res.tempo$intercept_group)){intercept_group <- res.tempo$intercept_group}
#   Id <-  res.tempo$Id
#   if(!is.null(res.tempo$maj_paramGR)){maj_paramGR <- res.tempo$maj_paramGR}   
#   prior_theta <- res.tempo$prior_theta
#   prior_sigma <- res.tempo$prior_sigma
#   family <- res.tempo$family
#  
#   rm(res.tempo)
#   
#   
#   Lv_obs <- matrix(1, nrow = n, ncol = G) 
#   prior_proba <- matrix(1, nrow = n, ncol = G)
#   if(test.ICM == TRUE){Lv_spat <- matrix(NA, nrow = n, ncol = G)}
#   tempo_lvc <- matrix(0, nrow = n, ncol = G) 
#   
#   #### matrices de voisinage ####
#   # Wweight_SR
#   if(!is.null(Wweight_SR) && !is.list(Wweight_SR)){Wweight_SR <- list(Wweight_SR)}
#   
#   # W
#   if(!is.null(Wdist_LR) && !is.list(Wdist_LR)){Wdist_LR <- list(Wdist_LR)}
#   
#   if(test.ICMregional == TRUE && is.null(distance_ref)){
#     distance_ref <- c(-1, sort(unique(unlist(lapply(Wdist_LR, function(x){unique(x@x)})))) + 0.000001)    
#   }
#   
#   # test W normalise
#   W_x <- unique(unlist(lapply(Wdist_LR, function(x){unique(x@x)})))
#   if(any(W_x %in% seq(0, length(W_x) - 1) == FALSE)){
#     
#     if(verbose){cat("Wdist_LR normalisation \n")}
#     
#     Wdist_LR <- lapply(Wdist_LR, function(x){
#       test <- range(x@x)
#       if(test[1] <= min(distance_ref) || test[2] >= max(distance_ref)){
#         stop("sfMM[Fcts-fMM_segmentation.R] : \'distance_ref\' is not coherent with Wdist_LR \n", 
#              "observations out (or equal) of the range of distance_ref : ", min(distance_ref), " ", max(distance_ref), "\n", 
#              "range(Wdist_LR@x) : ", test[1], " ", test[2], "\n")
#       }       
#       x@x <- as.numeric(cut(x@x, distance_ref)) - 1
#       return(x)
#     })
#     
#     distance_ref <- sort(unique(unlist(lapply(Wdist_LR, function(x){unique(x@x)}))))       
#   }
#   
#   # Id 
#   if(test.GR + test.ICM > 0){
#     if(is.null(Id)){
#       Id <- sapply(1:length(Wweight_SR), function(x){rep(x, nrow(Wweight_SR[[x]]))})
#       Id <- as.character(unlist(Id))
#     }    
#     identifier <- unique(Id)
#     n.Id <- length(identifier)
#     index_pat <- tapply(X = 1:n, INDEX = Id, function(x){x})
#   }else{
#     index_pat <- NULL
#     n.Id <- NULL
#   }
#   
#   #### tests preliminaires #####
#   
#   res.tempo <- test_fMM(n = n, M = M, G = G, data = data, index_pat = index_pat, coords = coords, 
#                             Var_reg = Var_reg, family = family, offset_reg = offset_reg, 
#                             formula_group = formula_group, Var_group = Var_group, 
#                             prior_theta = prior_theta, prior_sigma = prior_sigma, prior_proba = prior_proba, posterior_proba = posterior_proba, 
#                             test.GR = test.GR, seed = seed, 
#                             test.ICM = test.ICM, test.ICMregional = test.ICMregional, Wdist_LR = Wdist_LR, Wweight_SR = Wweight_SR, G_ICM = G_ICM, rho_ICM = rho_ICM, 
#                             export.predicteur = export.predicteur)
#   
#   export.predicteur <- res.tempo$export.predicteur
#   rm(res.tempo)
#   
#   #### initialisation #####
#   
#   res.tempo <- init2_fMM(n = n, M = M, G = G, data = data, Id = Id, 
#                            formula_reg = formula_reg, Var_reg = Var_reg, offset_reg = offset_reg, family = family, 
#                            formula_group = formula_group, Var_group = Var_group, intercept_group = intercept_group, 
#                            prior_theta = prior_theta, prior_sigma = prior_sigma, prior_proba = prior_proba, posterior_proba = posterior_proba, 
#                            test.GR = test.GR, seed = seed, 
#                            test.ICM = test.ICM, Wweight_SR = Wweight_SR, prior_prevalenceICM = prior_prevalenceICM, rho_max = rho_max, rho_ICM = rho_ICM, 
#                            test.ICMregional = test.ICMregional, distance_ref = distance_ref, Wdist_LR = Wdist_LR, coords = coords, threshold = threshold_potential, nbGroup_min = nbGroup_min, multiV = multiV, 
#                            iter_max = iter_max, verbose = verbose)
#   
#   prior_theta <- res.tempo$prior_theta
#   prior_sigma <- res.tempo$prior_sigma
#   theta <- res.tempo$theta 
#   hist_theta <- res.tempo$hist_theta
#   sigma <- res.tempo$sigma
#   lm_adj <- res.tempo$lm_adj # a voir si c est necessaire
#   data_reg <- res.tempo$data_reg
#   if(!is.null(res.tempo$data_group)){data_group <- res.tempo$data_group}
#   if(!is.null(res.tempo$hist_beta)){hist_beta <- res.tempo$hist_beta}
#   if(!is.null(res.tempo$beta)){beta <- res.tempo$beta}
#   if(!is.null(res.tempo$nb.param_group)){nb.param_group <- res.tempo$nb.param_group}
#   if(!is.null(res.tempo$Lv_obs)){Lv_obs <- res.tempo$Lv_obs}
#   prior_proba <- res.tempo$prior_proba
#   prior_proba_spat <- res.tempo$prior_proba_spat
#   if(!is.null(res.tempo$init_corrSpat)){init_corrSpat <- res.tempo$init_corrSpat}     
#   if(!is.null(res.tempo$seed_GR)){seed_GR <- res.tempo$seed_GR}
#   if(!is.null(res.tempo$nbVois.max)){nbVois.max <- res.tempo$nbVois.max}
#   rho_ICM <- res.tempo$rho_ICM
#   
#   #### critere de convergence ####
#   critere_param <- rep(NA, iter_max)
#   critere_lv <- rep(NA, iter_max)
#   critere_lvbis <- rep(NA, iter_max)
#   critere_lv[1] <- -Inf
#   critere_lvbis[1] <- -Inf
#   critere_param[1] <- NA
#   critere_test <- NA
#   iter <- 1
#   
#   
#   if(verbose == TRUE && is.null(posterior_proba)){
#     if(is.null(posterior_proba)){
#       mat.out <- matrix(unlist(prior_theta), ncol = G, nrow = M, byrow = TRUE)
#       row.names(mat.out) <- paste("intercept", 1:M, ":")
#       colnames(mat.out) <- paste(" group", 1:G, "")
#       cat("# initialisation \n")
#       print(mat.out)
#     }        
#   }
#   
#   cat("cv criterion : ", epsilon, "\n")
#   if(test.ICM == TRUE){cat(" beginning regularization criterion : ", epsilon_corrSpat, "\n \n")}
#   
#   while(  ( is.na(critere_test) || critere_test > epsilon ) && (iter <= iter_max) ){
#     #### Etape E ####
#     if(iter > 1 || (iter == 1 && is.null(posterior_proba))){
#       
#       ## evaluation du posterior
#       if(!identical(init_corrSpat, FALSE) && test.ICM == TRUE){
#         posterior_proba <- prior_proba * Lv_spat * Lv_obs
#         posterior_proba <- posterior_proba / rowSums(posterior_proba)    
#       }else{
#         posterior_proba <- prior_proba * Lv_obs
#         posterior_proba <- posterior_proba /rowSums(posterior_proba)    
#       }
#       
#       
#       ## retraitement spatial par GR
#       if(test.GR == TRUE){
#         if(trace_iter == TRUE){cat("GR : ")}
#         dist_necrose <- rep(Inf, n)
#         
#         if(maj_paramGR){
#           index_GN <- which(apply(posterior_proba, 1, function(x){which.max(x) == G}) == 1)          
#           sigma_GR <- sd(posterior_proba[index_GN, G])
#         }
#         
#         if(trace_iter == TRUE){
#           cat("sigma = ", round(sigma_GR, 3), 
#               "/ proba_GR = ", if(!is.null(proba_GR)){round(proba_GR, 3)}, 
#               "/ proba_GRseed = ", if(!is.null(proba_GRseed)){round(proba_GRseed, 3)}, 
#               "/ n.order = ", n.order, "\n")
#         }
#         
#         for(iter_pat in 1:n.Id){
#           
#           # identification de la zone necrosee
#           # Region_necrose <- GRalgo(Intensity = posterior_proba[index_pat[[iter_pat]], G], 
#                                        # Wweight_SR = Wweight_SR[[iter_pat]], 
#                                        # seed = seed_GR[[iter_pat]], 
#                                        # sigma_max = sigma_GR, rescale = FALSE, seuil_min = proba_GR, seuil_min_seed = proba_GRseed)$GR
#           
# 		  Region_necrose <- GRalgo(contrast = posterior_proba[index_pat[[iter_pat]], G], 
# 		                           W = Wweight_SR[[iter_pat]], 
# 								   seed = seed_GR[[iter_pat]], 
# 								   sigma_max = sigma_GR, range = c(0, 1), breaks = seq(0, 1, 0.01), step = 0.01, operator = "stats::sd", 
# 								   iter_max = 100, keep.lower = FALSE, keep.upper = TRUE, 
#                                    history.sigma = FALSE, history.step = FALSE, history.front = FALSE)$GR
#           
#           if(is.null(Region_necrose)){
#             warning("sfMM[Fcts-fMM_segmentation.R] : lesion obtained by the GR step is NULL \n")
#           }else{ 
#             
#             for(iter_ordre in 1:n.order){
#               
#               #             subW_SR <- Wweight_SR[[iter_pat]][Region_necrose, -Region_necrose]
#               #             voisins <- which(colSums(subW_SR > 0) > 0)
#               Region_sain <- (1:length(index_pat[[iter_pat]]))[-Region_necrose]       
#               voisins <- Region_sain[unique(Wweight_SR[[iter_pat]][Region_sain,Region_necrose, drop = FALSE]@i + 1)]
#               
#               posterior_proba[index_pat[[iter_pat]][voisins], 4] <- posterior_proba[index_pat[[iter_pat]][voisins], 4]*EDK(x = iter_ordre, bandwidth = 1, power = 2)
#               
#               Region_necrose <- c(voisins, Region_necrose)
#             }
#             posterior_proba[index_pat[[iter_pat]][setdiff((1:length(index_pat[[iter_pat]])), Region_necrose)], 4] <-  0 
#           }
#         }
#         
#         posterior_proba <- posterior_proba / rowSums(posterior_proba)
#         
#       }
#           
#       # mise en ordre des groupes selon un Y[[1]] croissant
#       ordre <- order(apply(posterior_proba, 2, function(x){weighted.mean(data[, Var_reg[[1]][1]], w = x)}))
#       posterior_proba <- posterior_proba[, ordre]
#       
#       if(any(colSums(posterior_proba) < 2)){
#         warning("sfMM[Fcts-fMM_segmentation.R] : one group had insufficient number of observations \n", 
#                 "Sum of the proba by group : ", paste(round(colSums(posterior_proba), digits = 2), collapse = " "), "\n")
#         return(NA)
#       }
#       
#       if(any(rowSums(posterior_proba) < 10^(-digit.pmin))){
#         warning("sfMM[Fcts-fMM_segmentation.R] : all posterior were close to 0 for one observation  \n")
#         return(NA)
#       }
#     }
#    
#     #### Etape M ####
#     
#     # Lv_obs <- matrix(1, nrow = n, ncol = G) 
#     Lv_obs[] <- 1
#     
#     ### modele de reponse
#     for(iter_g in 1:G){
#       if(!is.null(digit.pmin)){        
#         index_pos <- which(posterior_proba[, iter_g] > 10^(-digit.pmin))
#       }else{
#         index_pos <- 1:n        
#       }
#       
#       for(iter_m in 1:M){
#         
#         ## WLS
#         offset <- offset_reg[[iter_m]][, iter_g]
#         if(family[[iter_m]][[iter_g]]$family == "gaussian"){
#           lm_adj[[iter_m]][[iter_g]] <- eval(parse(text = paste(
#             "stats:::lm(formula = ", character_reg[[iter_m]][[iter_g]], ", 
#             data = data_reg,                                                                       
#             weights = posterior_proba[, iter_g], 
#             subset = index_pos, 
#             offset = offset, model = FALSE)", 
#             sep = "")))
#         }
#         
#         if(family[[iter_m]][[iter_g]]$family %in% c("Gamma", "binomial", "quasibinomial")){
#           lm_adj[[iter_m]][[iter_g]] <- eval(parse(text = paste(
#             "stats:::glm(formula = ", character_reg[[iter_m]][[iter_g]], ", 
#             family = family[[iter_m]][[iter_g]], 
#             data = data_reg,                                                                       
#             weights = posterior_proba[, iter_g], 
#             subset = index_pos, 
#             offset = offset, 
#             model = FALSE)", 
#             sep = "")))
#         }
#         if(family[[iter_m]][[iter_g]]$family == "Beta"){
#           lm_adj[[iter_m]][[iter_g]] <- eval(parse(text = paste(
#             "betareg:::betareg(formula = ", character_reg[[iter_m]][[iter_g]], ", 
#             link = \"logit\", link.phi = \"identity\", type = \"ML\", 
#             data = data_reg,                                                                       
#             weights = posterior_proba[, iter_g], 
#             subset = index_pos, model = FALSE)", 
#             sep = "")))
#         }
#         
#         ## recuperation des parametres
#         if(length(coefficients(lm_adj[[iter_m]][[iter_g]])) > 0){
#           if(family[[iter_m]][[iter_g]]$family == "gaussian" )  
#           {theta[[iter_m]][names(coefficients(lm_adj[[iter_m]][[iter_g]])[1]), iter_g] <- coefficients(lm_adj[[iter_m]][[iter_g]])[1]}
#           if(family[[iter_m]][[iter_g]]$family == "Gamma")
#           {theta[[iter_m]][names(coefficients(lm_adj[[iter_m]][[iter_g]])[1]), iter_g] <- exp(coefficients(lm_adj[[iter_m]][[iter_g]])[1])}
#           if(family[[iter_m]][[iter_g]]$family %in% c("Beta", "binomial", "quasibinomial"))
#           {theta[[iter_m]][names(coefficients(lm_adj[[iter_m]][[iter_g]])[1]), iter_g] <- exp(coefficients(lm_adj[[iter_m]][[iter_g]])[1])}
#         }
#         
#         df <- sum(lm_adj[[iter_m]][[iter_g]]$weights) 
#         sigma[[iter_m]][iter_g] <- sqrt(sum(lm_adj[[iter_m]][[iter_g]]$weights * residuals(lm_adj[[iter_m]][[iter_g]], type = "response")^2) / df)
#         
#         ## calcul des predits
#         predit <- predict(lm_adj[[iter_m]][[iter_g]], 
#                           newdata = cbind(data_reg, offset = 0), 
#                           type = "response", se.fit = FALSE)
#         
#         if(!is.null(offset) && family[[iter_m]][[iter_g]]$family == "gaussian")
#         {predit <- predit + offset}
#         if(!is.null(offset) && family[[iter_m]][[iter_g]]$family == "Gamma")
#         {predit <- predit * exp(offset)}
#         if(!is.null(offset) && family[[iter_m]][[iter_g]]$family %in% c("Beta", "binomial", "quasibinomial"))
#         {predit <- logit(inv.logit(predit) + offset)}
#         
#         ## calcul de la vraisemblance
#         if(family[[iter_m]][[iter_g]]$family == "gaussian")
#         { Lv_obs[, iter_g] <- Lv_obs[, iter_g]  *  dnorm(x = data[, Var_reg[[iter_m]][1]], 
#                                                        mean = predit, 
#                                                        sd = sigma[[iter_m]][iter_g], 
#                                                        log = FALSE)}
#         
#         if(family[[iter_m]][[iter_g]]$family == "Gamma")
#         { Lv_obs[, iter_g] <- Lv_obs[, iter_g] *  dgamma(x = data[, Var_reg[[iter_m]][1]], 
#                                                        scale = sigma[[iter_m]][iter_g]^2 / predit, 
#                                                        shape = predit^2 / sigma[[iter_m]][iter_g]^2, 
#                                                        log = FALSE)}  
#         
#         if(family[[iter_m]][[iter_g]]$family %in% c("binomial", "quasibinomial"))
#         { Lv_obs[, iter_g] <- Lv_obs[, iter_g] * dbinom(x = data[, Var_reg[[iter_m]][1]], 
#                                                       size = 1, prob = predit, 
#                                                       log = FALSE)}
#         
#         if(family[[iter_m]][[iter_g]]$family == "Beta")
#         { Lv_obs[, iter_g] <- Lv_obs[, iter_g] * dbeta(x = data[, Var_reg[[iter_m]][1]], 
#                                                      shape1=predit * (predit * (1 - predit) / sigma[[iter_m]][iter_g]^2 - 1), 
#                                                      shape2=(1 - predit) * (predit * (1 - predit) / sigma[[iter_m]][iter_g]^2 - 1), 
#                                                      log = FALSE)}
#       }
#     }
#     
#     ### modele de structure
#     
#     if(!is.null(formula_group)){
#       if(!is.null(digit.pmin)){
#         index_pos <- which(as.vector(posterior_proba) > 10^(-digit.pmin))
#         index2_pos <- rep(1:n, G)[index_pos]
#         nb_pos <- cumsum(c(0, apply(posterior_proba, 2, function(x){sum(x > 10^{-digit.pmin})})))
#       }else{
#         index_pos <- 1:(G * n)  
#         index2_pos <- rep(1:n, G)
#         nb_pos <- cumsum(c(0, rep(n, G)))
#       }
#       
#       ## multinomiale     
#       multinom_beta <- eval(parse(text = paste(
#         "nnet::multinom(formula = ", formula_group, ", 
#                   weights = as.vector(posterior_proba), 
#                   subset = index_pos, 
#                   data = data_group, 
#                   wts = beta, 
#                   verbose = FALSE)", 
#         sep = "")))
#       
#       
#       # recuperation des parametres
#       beta <- coefficients(multinom_beta)
#       if(!is.matrix(beta)){beta <- as.matrix(beta, nrow = G - 1, ncol = nb.param_group)} ## beta des groupes 2:G
#       
#       # Calcul des priors
#       # prior_proba <- matrix(NA, nrow = n, ncol = G)
#       prior_proba[] <- NA 
#       
#       if(G == 2){
#         for(iter_g in 1:2){
#           index_g <- (1 + nb_pos[iter_g]):nb_pos[iter_g + 1]
#           prior_proba[index2_pos[index_g], 2] <- fitted.values(multinom_beta)[index_g]
#         }
#         prior_proba[,1] <- 1 - prior_proba[, 2]
#       }else{    
#         for(iter_g in 1:G){
#           index_g <- (1 + nb_pos[iter_g]):nb_pos[iter_g + 1]
#           prior_proba[index2_pos[index_g], ] <- fitted.values(multinom_beta)[index_g, ]
#         }
#       }
#       
#       if(prior_prevalence == FALSE){
#         prevalence <- colMeans(posterior_proba, na.rm = TRUE)
#         prior_proba <- prior_proba * matrix(1 / prevalence, byrow = TRUE, nrow = n, ncol = G)
#         prior_proba <- prior_proba / rowSums(prior_proba)
#       }
#     }
#     
#     ### ICM 
#     # test correction spatiale
#     if( test.ICM == TRUE && identical(init_corrSpat, FALSE) && iter > 1){
#       
#       ## critere de cv de la lv
#       critere_lv[iter + 1] <- sum(log(rowSums(prior_proba * Lv_obs) ))
#       
#       ## critere de cv des parametres sur les distributions
#       critere_param[iter + 1] <- max(max(-1, abs(unlist(theta) - unlist(hist_theta)) / unlist(hist_theta), na.rm = TRUE), 
#                                    max(abs(beta-hist_beta) / hist_beta))
#       
#       ## critere test
#       critere_test <- max(c(abs((critere_lv[iter + 1]-critere_lv[iter]) / critere_lv[iter]), critere_param[iter + 1]), na.rm = TRUE)
#       
#       if(!is.na(critere_test) && (critere_test <= epsilon_corrSpat) ){
#         init_corrSpat <- iter
#         cat("*** init. spatial regularization *** ")        
#         cat("\n")
#       }
#       
#     }
#     
#     # maj rho
#     if(update_rho == TRUE && test.ICM == TRUE && !identical(init_corrSpat, FALSE)){
#       
#       resICMp <- matrix(0, nrow = n, ncol = length(unique(G_ICM)))
#       resICMp[] <- 0
#       for(iter_g in 1:G){
#         resICMp[, G_ICM[iter_g]] <- resICMp[, G_ICM[iter_g]] + posterior_proba[, iter_g]
#       }    
#       
#       rho_ICM <- calcPotts(W_SR = Wweight_SR, W_LR = Wdist_LR, sample = NULL, Id = Id, 
#                            prior_prevalence = prior_prevalenceICM, test_Id = NULL, 
#                            rho = "auto", rho_max = rho_max, Y = resICMp, distance_ref = distance_ref, multiV = multiV, 
#                            test.ICMregional = test.ICMregional, coords = coords, nbGroup_min = nbGroup_min, distcat = TRUE, 
#                            trace_iter = FALSE, trace_radius = FALSE, trace_rho = FALSE, trace_pat = FALSE)$rho
#       
#     }
#     
#     # estimation du prior spatial
#     if(!identical(init_corrSpat, FALSE) && test.ICM == TRUE){
#       
#       if(trace_iter == TRUE){cat("ICM ")}
#       posterior_proba.prelim <- prior_proba * Lv_obs 
#       posterior_proba.prelim <- posterior_proba.prelim / rowSums(posterior_proba.prelim)
#       
#       resICMp <- matrix(0, nrow = n, ncol = length(unique(G_ICM)))
#       resICMp[] <- 0
#       for(iter_g in 1:G){
#         resICMp[, G_ICM[iter_g]] <- resICMp[, G_ICM[iter_g]] + posterior_proba.prelim[, iter_g]
#       }
#       
#       resICMp <- calcPotts(W_SR = Wweight_SR, W_LR = Wdist_LR, 
#                            sample = resICMp, 
#                            test.ICMregional = test.ICMregional, coords = coords, 
#                            rho = rho_ICM, test_Id = NULL, 
#                            Y = NULL, ls_index_pat = index_pat, distance_ref = distance_ref, distcat = TRUE, 
#                            Id = Id, nbGroup_min = nbGroup_min, multiV = multiV, 
#                            trace_iter = FALSE, trace_radius = trace_radius, trace_pat = trace_iter, iter_max = 150, critere = 0.005)   
#       gc()
#       
#       for(iter_g in 1:G){
#         Lv_spat[, iter_g] <- resICMp$pred_spatial[, G_ICM[iter_g]]
#       } 
#       
#       prior_proba_spat <- prior_proba * Lv_spat
#       prior_proba_spat <- prior_proba_spat / rowSums(prior_proba_spat)
#     }
#     
#     #### Critere de convergence ####
#     
#     ## log-vraisemblance
#     if(test.ICM == FALSE || identical(init_corrSpat, FALSE)){
#       critere_lv[iter + 1] <- sum(log(rowSums(prior_proba * Lv_obs) ))
#     }else{
#       critere_lv[iter + 1] <- sum(log(rowSums(prior_proba_spat  * Lv_obs) ))
#     }
#     
#     ## log-vraisemblance completee
#     #     tempo_lvc <- matrix(0, nrow = n, ncol = G) 
#     tempo_lvc[] <- 0
#     
#     if(!is.null(digit.pmin)){
#       index_pos <- which(posterior_proba > 10^(-digit.pmin))
#     }else{
#       index_pos <- 1:(G * n)
#     }
#     
#     if(test.ICM == FALSE || identical(init_corrSpat, FALSE)){
#       tempo_lvc[index_pos] <- posterior_proba[index_pos] * log(prior_proba[index_pos] * Lv_obs[index_pos] / posterior_proba[index_pos])  
#     }else{
#       tempo_lvc[index_pos] <- posterior_proba[index_pos] * log(prior_proba_spat[index_pos] * Lv_obs[index_pos] / posterior_proba[index_pos]) 
#     }
#     critere_lvbis[iter + 1] <- sum(tempo_lvc[!is.infinite(tempo_lvc)], na.rm = TRUE) 
#     
#     ## critere de cv des parametres sur les distributions
#     if(iter != 1)
#     {critere_param[iter + 1] <- max(max(-1, abs(unlist(theta) - unlist(hist_theta)) / unlist(hist_theta), na.rm = TRUE), 
#                                   max(abs(beta-hist_beta) / hist_beta))}else
#                                   {critere_param[iter + 1] <- NA}
#     
#     ## critere test
#     critere_test <- na.omit(c(abs((critere_lv[iter + 1]-critere_lv[iter]) / critere_lv[iter]), critere_param[iter + 1]))
#     if(length(critere_test) == 0)
#     {critere_test <- Inf}else{critere_test <- max(critere_test)}
#     if(!identical(init_corrSpat, FALSE) && iter == init_corrSpat){critere_test <- epsilon + 1}
#     
#     
#     
#     ## affichage iteration
#     if(trace_iter == TRUE){
#       cat("# iteration ", iter, " (ll = ", critere_lv[iter + 1], " | diff ll = ", (critere_lv[iter + 1] - critere_lv[iter]) / abs(critere_lv[iter]), " | diff param = ", critere_param[iter + 1], " | completed lv = ", critere_lvbis[iter + 1], ") \n")
# #       cat("# iteration ", iter, " (LL = ", critere_lv[iter + 1], " | diff param = ", critere_param[iter + 1], " | completed LL = ", critere_lvbis[iter + 1], ")\n")
#       
#       mat.out <- matrix(c(unlist(lapply(theta, function(x){x[1, ]})), unlist(sigma), colMeans(prior_proba), colMeans(posterior_proba)), nrow = 2 * M + 2, ncol = G, byrow = TRUE)
#       row.names(mat.out) <- c(paste("intercept", 1:M, ":"), paste("sigma", 1:M, ":"), "<prior > :", "<posterior > :")
#       colnames(mat.out) <- paste(" group", 1:G, "") 
#       print(mat.out)
#       if(test.ICM == TRUE && !identical(init_corrSpat, FALSE)){cat("ICM parameters : ", paste(round(rho_ICM, 2), collapse = " "), "\n")}
#       if(!is.null(formula_group)){print(beta)}
#     }
#     iter <- iter + 1
#     hist_theta <- theta
#     hist_beta <- beta
#   }
#   
#   if(verbose == TRUE){
#     
#     if(trace_iter == TRUE){cat("\n")}
#     #### affichage iteration finale ####
#     
#     cat("### FINAL iteration  ", iter - 1, " (ll = ", critere_lv[iter], " | diff ll = ", (critere_lv[iter] - critere_lv[iter - 1]) / abs(critere_lv[iter - 1]), " | diff param = ", critere_param[iter], " | completed lv = ", critere_lvbis[iter], ") \n")
#     mat.out <- matrix(c(unlist(lapply(theta, function(x){x[1,]})), unlist(sigma), colMeans(prior_proba), colMeans(posterior_proba)), nrow = 2 * M + 2, ncol = G, byrow = TRUE)
#     row.names(mat.out) <- c(paste("intercept", 1:M, ":"), paste("sigma", 1:M, ":"), "<prior > :", "<posterior > :")
#     colnames(mat.out) <- paste(" group", 1:G, "") 
#     print(mat.out)
#     if(test.ICM == TRUE && !identical(init_corrSpat, FALSE)){cat("ICM parameters : ", paste(round(rho_ICM, 2), collapse = " "), "\n")}
#     if(!is.null(formula_group)){print(beta)}
#     
#     if(is.na(critere_test) == FALSE && critere_test <= epsilon)
#     {cat("*** Convergence *** \n")}else
#     {cat("*** Convergence failure ***\n")}
#     
#   }
#   
#   #### Export ####
#   if(verbose == TRUE){cat("*** Export *** \n")}
#   
#   # attribution des groupes
#   group <- list()
#   for(iter_g in 1:G){
#       group[[iter_g]] <- which(max.col(posterior_proba) == iter_g)
#   }  
#   
#   # mise en forme du predicteur
#   if(export.predicteur == TRUE){
#     reg <- list(Y = lm_adj, 
#                 group = if(!is.null(formula_group)){multinom_beta}else{NULL})
#   }else{
#     reg <- list(Y = lapply(lm_adj, function(x){lapply(x, summary)}), 
#                 group = if(!is.null(formula_group)){coef(multinom_beta)}else{NULL})
#   }
#   
#   t1 <- date()
#   
#   param.export <- NULL  
#   if("index" %in% names(data)){param.export <- c(param.export, "index")}
#   
#   if(!is.null(Id)){
#     Id.export <- data.frame(Id = Id, data[,param.export, drop = FALSE])
#   }else{
#     Id.export <- data[,param.export, drop = FALSE]
#   }
#   
#   if(!is.null(coords)){Id.export <- data.frame(Id.export, coords)}
#   
#   
#   return(res <- list(Y = data[, unlist(lapply(Var_reg, function(x){x[1]})), drop = FALSE], 
#                      Id = Id.export, 
#                      theta = theta, 
#                      sigma = sigma, 
#                      reg = reg, 
#                      beta = if(!is.null(formula_group)){beta}else{NULL}, 
#                      rho_ICM = if(!identical(init_corrSpat, FALSE)){rho_ICM}else{NULL}, 
#                      posterior_proba = posterior_proba, 
#                      prior_proba = prior_proba, 
#                      proba_priorReg = if(!identical(init_corrSpat, FALSE)){prior_proba_spat}else{NULL}, 
#                      radius = if(!identical(init_corrSpat, FALSE)){resICMp$radius}else{NULL}, 
#                      group = group, 
#                      family = family, 
#                      critere_cv = list(lv = critere_lv, lv_completee = critere_lvbis, critere_param = critere_param), 
#                      cv = (is.na(critere_test) == FALSE && critere_test <= epsilon), 
#                      tps = list(t0, t1)))
# }
# 
# 
# 
# #### 2- Simulation Functions  ####
# simul_fMM <- function(n, G, 
#                       Ydist = "gaussian", Y_mu, Y_sd, nclass = "auto", 
#                       beta = FALSE, Xdist = NULL, X_mu = NULL, X_sd = NULL, 
#                       n.Id = NULL, Id_sd = NULL, 
#                       noise_sd = NULL, 
#                       display_Y = TRUE, display_X = TRUE, 
#                       window = FALSE, filename = "calc-echantillon", width = 1000, height = 700, path = NULL, unit = "px", res = NA)
# {
#   # tests preliminaires
#   if(length(n) == 1){n <- rep(n, G)}
#   if(length(noise_sd) == 1){noise_sd <- rep(noise_sd, G)}
#   if(length(Ydist) == 1){Ydist <- rep(Ydist, G)}
#   if(!is.null(Xdist) && length(Xdist) == 1){Xdist <- rep(Xdist, G)}
#   
#   if(length(n) != G)
#   {stop("simul_fMM[Fcts-fMM_segmentation.R] : \'n\' has incorrect length \n", 
#         "required length (\'G\') : ", G, "\n", 
#         "proposed length : ", length(n), "\n")}
#   
#   if(length(Ydist) != G)
#   {stop("simul_fMM[Fcts-fMM_segmentation.R] : \'Ydist\' has incorrect length \n", 
#         "required length (\'G\') : ", G, "\n", 
#         "proposed length : ", length(Ydist), "\n")}
#   
#   if(sum(Ydist %in% c("unif", "unif_segment", "gaussian", "Gamma", "Beta") == FALSE) > 0)
#   {stop("simul_fMM[Fcts-fMM_segmentation.R] : wrong specification of \'Ydist\' \n", 
#         "possible values : \"unif\" \"unif_segment\" \"gaussian\" \"Gamma\" \"Beta\" \n", 
#         "proposed value : ", paste(Ydist, collapse = " "), "\n")}
#   
#   if(length(Y_mu) != G)
#   {stop("simul_fMM[Fcts-fMM_segmentation.R] : \'Y_mu\' has incorrect length \n", 
#         "required length (\'G\') : ", G, "\n", 
#         "proposed length : ", length(Y_mu), "\n")}
#   
#   if(length(Y_sd) != G)
#   {stop("simul_fMM[Fcts-fMM_segmentation.R] : \'Y_sd\' has incorrect length \n", 
#         "required length (\'G\') : ", G, "\n", 
#         "proposed length : ", length(Y_sd), "\n")}
#   
#   if(length(beta) != 1 || beta %in% c(TRUE, FALSE) == FALSE)
#   {stop("simul_fMM[Fcts-fMM_segmentation.R] : wrong specification of \'beta\' \n", 
#         "possible \'beta\' : TRUE or FALSE\n", 
#         "proposed \'beta\' : ", beta, "\n")}
#   
#   if(beta == TRUE && (is.null(Xdist) || length(Xdist) != G))
#   {stop("simul_fMM[Fcts-fMM_segmentation.R] : \'Xdist\' has incorrect length \n", 
#         "required length (\'G\') : ", G, " \n", 
#         "proposed length : ", length(Xdist), "\n")}
#   
#   if(beta == TRUE && sum(Xdist %in% c("unif", "gaussian", "Gamma") == FALSE) > 0)
#   {stop("simul_fMM[Fcts-fMM_segmentation.R] : wrong specification of \'Ydist\'\n", 
#         "possible values : \"unif\" \"gaussian\" \"gamma\" \n", 
#         "proposed value : ", Ydist, "\n")}
#   
#   if(beta == TRUE && (is.null(X_mu) || length(X_mu) != G))
#   {stop("simul_fMM[Fcts-fMM_segmentation.R] : \'X_mu\' has incorrect length \n", 
#         "required length (\'G\') : ", G, "\n", 
#         "proposed length : ", length(X_mu), "\n")}
#   
#   if(beta == TRUE && (is.null(X_sd) || length(X_sd) != G))
#   {stop("simul_fMM[Fcts-fMM_segmentation.R] : \'X_sd\' has incorrect length \n", 
#         "required length (\'G\') : ", G, "\n", 
#         "proposed length : ", length(X_sd), "\n")}
#   
#   if(!is.null(n.Id) && (is.null(Id_sd) || length(Id_sd) != G))
#   {stop("simul_fMM[Fcts-fMM_segmentation.R] : \'XId_sd\'X has incorrect length \n", 
#         "required length (\'G\') : ", G, "\n", 
#         "proposed length : ", length(Id_sd), "\n")}
#   
#   if(!is.null(n.Id) && length(unique(n)) != 1)
#   {stop("simul_fMM[Fcts-fMM_segmentation.R] : \'n\' must be the same for each group \n", 
#         "proposed \'n\' : ", paste(n, collapse = " "), "\n")}
#   
#   if(window %in% c(TRUE, FALSE, "png", "eps", "svg") == FALSE)
#   {stop("simul_fMM[Fcts-fMM_segmentation.R] : \'window\' doit valoir TRUE FALSE \"eps\" \"png\" \"svg\" \n", 
#         "proposed \'window\' : ", window, "\n")}
#   
#   if(!is.null(noise_sd) && length(noise_sd) != G)
#   {stop("simul_fMM[Fcts-fMM_segmentation.R] : wrong specification of \'noise_sd\' \n", 
#         "required length (\'G\') : ", G, "\n", 
#         "proposed length : ", length(noise_sd), "\n")}
#   
#   if(unit %in% c("px", "in", "cm", "mm") == FALSE){
#     stop("simul_fMM[Fcts-fMM_segmentation.R] : wrong specification of \'unit\' \n", 
#          "possible \'unit\' : \"px\" \"in\" \"cm\" \"mm\" \n", 
#          "proposed \'unit\' : ", unit, "\n")
#   }
#   
#   # generation des priors
#   prior <- matrix(NA, nrow = sum(n), ncol = G)
#   
#   if(beta == FALSE){
#     X <- NULL
#     X_g <- NULL
#   }
#   
#   if(beta == TRUE){
#     
#     param1_loi <- rep(NA, G)
#     param2_loi <- rep(NA, G)
#     X_g <- list()
#     X <- numeric(0)
#     
#     for(iter_g in 1:G){
#       if(Xdist[iter_g] == "gaussian")
#       {param1_loi[iter_g] <- X_mu[iter_g] 
#        param2_loi[iter_g] <- X_sd[iter_g]
#        X_g[[iter_g]] <- rnorm(n = n[iter_g], mean = param1_loi[iter_g], sd = param2_loi[iter_g])
#       }
#       
#       if(Xdist[iter_g] == "Gamma")
#       {param1_loi[iter_g] <- X_mu[iter_g]^2 / X_sd[iter_g]^2
#        param2_loi[iter_g] <- 1 / (X_sd[iter_g]^2 / X_mu[iter_g])
#        X_g[[iter_g]] <- rgamma(n = n[iter_g], shape = param1_loi[iter_g], rate = param2_loi[iter_g])
#       }
#       
#       if(Xdist[iter_g] == "Beta"){
#         
#         if(X_mu[iter_g] >= 1 || X_mu[iter_g] <= 0){
#           stop("simul_fMM[Fcts-fMM_segmentation.R] : wrong specification of \'X_mu\'\n", 
#                "\'X_mu\' must be between 0 and 1 for a beta law \n", 
#                " proposed \'X_mu\'[", iter_g, "] : ", X_mu[iter_g], "\n")
#         }
#         
#         if(X_sd[iter_g] >= 0.25){
#           stop("simul_fMM[Fcts-fMM_segmentation.R] : wrong specification of \'X_sd\' \n", 
#                "\'X_sd\' must be between 0 and 1 / 4 for a beta law \n", 
#                "proposed \'X_sd\'[", iter_g, "] : ", X_sd[iter_g], "\n")
#         }
#         
#         param1_loi[iter_g] <- X_mu[iter_g] * ( X_mu[iter_g] * (1 - X_mu[iter_g]) / X_sd[iter_g]^2 -1)
#         param2_loi[iter_g] <- (1 - X_mu[iter_g]) * ( X_mu[iter_g] * (1 - X_mu[iter_g]) / X_sd[iter_g]^2 -1)
#         X_g[[iter_g]] <- rbeta(n = n[iter_g], shape1=param1_loi[iter_g], shape2=param2_loi[iter_g])
#       }
#       
#       if(Xdist[iter_g] == "unif")
#       {param1_loi[iter_g] <- X_mu[iter_g] - X_sd[iter_g]*sqrt(3)
#        param2_loi[iter_g] <- X_mu[iter_g] + X_sd[iter_g]*sqrt(3)
#        X_g[[iter_g]] <- runif(n = n[iter_g], min = param1_loi[iter_g], max = param2_loi[iter_g])
#       }
#       
#       if(Xdist[iter_g] == "unif_segment")
#       {  param1_loi[iter_g] <- X_mu[iter_g] 
#          param2_loi[iter_g] <- X_sd[iter_g]
#          X_g[[iter_g]] <- runif(n = n[iter_g], min = param1_loi[iter_g], max = param2_loi[iter_g])  }
#       
#       X_mu[iter_g] <- mean(X_g[[iter_g]]) 
#       X_sd[iter_g] <- sd(X_g[[iter_g]])
#       X <- c(X, X_g[[iter_g]])
#     }
#     
#     Y_tempo <- as.numeric(apply(rbind(1:G, n), 2, function(x){rep(x[1], x[2])}))
#     res_multinom  <- nnet::multinom(Y_tempo ~  X, verbose = FALSE)
#     beta <- coef(res_multinom)
#     prior[, 2:G] <- res_multinom$fitted.values
#     prior[, 1] <- 1 - rowSums(res_multinom$fitted.values)
#     
#   }
#   
#   
#   
#   # generation des lois
#   
#   famille <- list()
#   for(iter_g in 1:G){
#     if(Ydist[iter_g] == "gaussian"){famille[[iter_g]] <- gaussian(link = "identity") }
#     if(Ydist[iter_g] == "Gamma"){famille[[iter_g]] <- Gamma(link = "log") }
#   }
#   
#   param1_loi <- rep(NA, G)
#   param2_loi <- rep(NA, G)
#   Y_g <- list()
#   group <- matrix(0, nrow = sum(n), ncol = G)
#   cum_n <- c(0, cumsum(n))
#   
#   for(iter_g in 1:G){
#     
#     group[(cum_n[iter_g] + 1):cum_n[iter_g + 1], iter_g] <- 1
#     
#     if(Ydist[iter_g] == "gaussian")
#     {param1_loi[iter_g] <- Y_mu[iter_g] 
#      param2_loi[iter_g] <- Y_sd[iter_g]
#      Y_g[[iter_g]] <- rnorm(n[iter_g], param1_loi[iter_g], param2_loi[iter_g])
#     }
#     
#     if(Ydist[iter_g] == "Gamma")
#     {param1_loi[iter_g] <- Y_mu[iter_g]^2 / Y_sd[iter_g]^2
#      param2_loi[iter_g] <- 1 / (Y_sd[iter_g]^2 / Y_mu[iter_g])
#      Y_g[[iter_g]] <- rgamma(n[iter_g], param1_loi[iter_g], param2_loi[iter_g])
#     }
#     
#     if(Ydist[iter_g] == "Beta"){
#       
#       if(Y_mu[iter_g] >= 1 || Y_mu[iter_g] <= 0){
#         stop("simul_fMM[Fcts-fMM_segmentation.R] : wrong specification of \'Y_mu\' \n", 
#              "\'Y_mu\' must be between 0 and 1 for a beta law \n", 
#              "proposed \'Y_mu\'[", iter_g, "] : ", Y_mu[iter_g], "\n")
#       }
#       
#       if(Y_sd[iter_g] >= 0.25){
#         stop("simul_fMM[Fcts-fMM_segmentation.R] : wrong specification of \'Y_sd\' \n", 
#              "\'Y_sd\' must be between 0 and 1 / 4 for a beta law \n", 
#              "proposed \'Y_sd\'[", iter_g, "] : ", Y_sd[iter_g], "\n")
#       }
#       
#       param1_loi[iter_g] <- Y_mu[iter_g] * ( Y_mu[iter_g] * (1 - Y_mu[iter_g]) / Y_sd[iter_g]^2 -1)
#       param2_loi[iter_g] <- (1 - Y_mu[iter_g]) * ( Y_mu[iter_g] * (1 - Y_mu[iter_g]) / Y_sd[iter_g]^2 - 1)
#       Y_g[[iter_g]] <- rbeta(n = n[iter_g], shape1=param1_loi[iter_g], shape2=param2_loi[iter_g])
#     }
#     
#     if(Ydist[iter_g] == "unif")
#     {param1_loi[iter_g] <- Y_mu[iter_g] + Y_sd[iter_g]*sqrt(3)
#      param2_loi[iter_g] <- Y_mu[iter_g] - Y_sd[iter_g]*sqrt(3) 
#      Y_g[[iter_g]] <- runif(n[iter_g], param1_loi[iter_g], param2_loi[iter_g])
#     }
#     
#     if(Ydist[iter_g] == "unif_segment")
#     { param1_loi[iter_g] <- Y_mu[iter_g] 
#       param2_loi[iter_g] <- Y_sd[iter_g]
#       Y_g[[iter_g]] <- runif(n[iter_g], param1_loi[iter_g], param2_loi[iter_g])  }
#     
#     Y_mu[iter_g] <- mean(Y_g[[iter_g]]) 
#     Y_sd[iter_g] <- sd(Y_g[[iter_g]])
#     
#   }
#   
#   #### niveau hierarchique 
#   if(!is.null(n.Id)){
#     Id0 <- list()
#     n.cum <- c(0, cumsum(n))
#     u_G <- matrix(NA, nrow = n.Id, ncol = G)
#     
#     for(iter_g in 1:G){
#       Id0[[iter_g]] <- cut((n.cum[iter_g] + 1):n.cum[iter_g + 1], breaks = n.Id, labels = 1:n.Id)
#       u_G[, iter_g] <- rnorm(n.Id, mean = 0, sd = Id_sd[iter_g])
#       Y_g[[iter_g]] <- Y_g[[iter_g]] + u_G[Id0[[iter_g]], iter_g]
#     }
#     
#   }
#   
#   #### ajout de bruit
#   if(!is.null(noise_sd)){
#     for(iter_g in 1:G){
#       Y_g[[iter_g]] <- Y_g[[iter_g]] + rnorm(length(Y_g[[iter_g]]), 0, noise_sd[iter_g])
#     }
#   }
#   
#   #### assemblage
#   Y <- unlist(Y_g)  
#   
#   if(display_X || display_Y){
#     scale <- switch(unit, 
#                     "px" = 1, 
#                     "in" = 90, 
#                     "cm" = 35.43, 
#                     "mm" = 3.543)
#     
#     if(nclass == "auto")
#     {nclass = min(100, length(Y) / 100)}
#     col <- c("red", "blue", "green", "yellow", "purple", "brown")
#     
#     if(!identical(beta, FALSE) && display_X == TRUE){
#       
#       switch(window, 
#              T = dev.new(), 
#              eps = postscript(file = paste(path, filename, "_X.eps", sep = ""), width = width * scale / 90, height = height * scale / 90, horizontal = FALSE, onefile = FALSE, paper = "special"), 
#              svg = svg(filename = paste(path, filename, "_X.svg", sep = ""), width = width * scale / 90, height = height * scale / 90, onefile = FALSE), 
#              png = png(filename = paste(path, filename, "_X.png", sep = ""), width = width * scale, height = height * scale, res = res)
#       )
#       
#       layout(matrix(c(rep(1, G), 2:(G + 1)), byrow = FALSE, nrow = G, ncol = 2), c(2,1))
#       par(mar = c(3, 3, 1, 1), mgp = c(1.5, 0.5, 0), mai = c(0.5, 0.75, 0.25, 0.25))
#       
#       xlim <- c(NA, NA)
#       ylim <- c(NA, NA)
#       for(iter_g in 1:G){
#         xlim <- c(min(xlim[1], density(X_g[[iter_g]])$x, na.rm = TRUE), max(xlim[2], density(X_g[[iter_g]])$x, na.rm = TRUE))
#         ylim <- c(min(ylim[1], density(X_g[[iter_g]])$y, na.rm = TRUE), max(ylim[2], density(X_g[[iter_g]])$y, na.rm = TRUE))
#       }
#       
#       hist(1, xlim = xlim, ylim = ylim, col = "white", border = "white", nclass = nclass, freq = FALSE, 
#            main = "X_g distribution", xlab = "")
#       for(iter_g in 1:G)
#       {hist(X_g[[iter_g]], col = col[iter_g], nclass = nclass, freq = FALSE, add = TRUE, border = "grey")
#        points(density(X_g[[iter_g]])$x, density(X_g[[iter_g]])$y, 
#               type = "l", lwd = 3, col = paste(col[iter_g], "4", sep = ""))}
#       abline(v = lapply(X_g, "mean"), lwd = 3, col = paste(col, "4", sep = ""))
#       
#       for(iter_g in 1:G){
#         hist <- hist(prior[, iter_g], col = col[iter_g], freq = TRUE, min(100, length(Y) / 100), 
#                      main = paste("prior g", iter_g, sep = ""), xlab = "", xlim = c(-0.1, 1.1))
#         hist(c(rep(1.1, sum(prior[, iter_g]) * max(hist$count) / length(prior[, iter_g])), rep(-0.1, sum(1-prior[, iter_g]) * max(hist$count) / length(prior[, iter_g]))), 
#              freq = TRUE, col = "grey", add = TRUE, breaks = seq(-0.1, 1.1, 0.1))
#       }
#       
#       
#       switch(window, 
#              eps = dev.off(), 
#              svg = dev.off(), 
#              png = dev.off()
#       )
#     }
#     
#     if(display_Y == TRUE){
#       switch(window, 
#              T = dev.new(), 
#              eps = postscript(file = paste(path, filename, "_X.eps", sep = ""), width = width * scale / 90, height = height * scale / 90, horizontal = FALSE, onefile = FALSE, paper = "special"), 
#              svg = svg(filename = paste(path, filename, "_X.svg", sep = ""), width = width * scale / 90, height = height * scale / 90, onefile = FALSE), 
#              png = png(filename = paste(path, filename, "_X.png", sep = ""), width = width * scale, height = height * scale, res = res)
#       )
#       
#       par(mfrow = c(1, 2), mar = c(3, 3, 1, 1), mgp = c(1.5, 0.5, 0))
#       
#       xlim <- c(NA, NA)
#       ylim <- c(NA, NA)
#       for(iter_g in 1:G){
#         xlim <- c(min(xlim[1], density(Y_g[[iter_g]])$x, na.rm = TRUE), max(xlim[2], density(Y_g[[iter_g]])$x, na.rm = TRUE))
#         ylim <- c(min(ylim[1], density(Y_g[[iter_g]])$y, na.rm = TRUE), max(ylim[2], density(Y_g[[iter_g]])$y, na.rm = TRUE))
#       }
#       
#       hist(1, xlim = xlim, ylim = ylim, col = "white", border = "white", nclass = nclass, freq = FALSE, 
#            main = "Y_g distribution", xlab = "")
#       for(iter_g in 1:G)
#       {hist(Y_g[[iter_g]], col = col[iter_g], nclass = nclass, freq = FALSE, add = TRUE, border = "grey")
#        points(density(Y_g[[iter_g]])$x, density(Y_g[[iter_g]])$y, 
#               type = "l", lwd = 3, col = paste(col[iter_g], "4", sep = ""))}
#       abline(v = lapply(Y_g, "mean"), lwd = 3, col = paste(col, "4", sep = ""))
#       
#       hist(Y, col = "grey", border = "black", freq = FALSE, xlab = "", 
#            main = "Y distribution", nclass = nclass)
#       
#       switch(window, 
#              eps = dev.off(), 
#              svg = dev.off(), 
#              png = dev.off()
#       )
#     }
#   }
#   
#   return(list(Y = Y, 
#               X = X, 
#               moy = mean(Y), 
#               var = var(Y), 
#               Y_g = Y_g, 
#               X_g = X_g, 
#               moy_g = unlist(lapply(Y_g, "mean")), 
#               var_g = unlist(lapply(Y_g, "var")), 
#               Ydist = famille, 
#               Y_mu = Y_mu, 
#               Y_sd = Y_sd, 
#               beta = beta, 
#               Xdist = Xdist, 
#               X_mu = X_mu, 
#               X_sd = X_sd, 
#               group = group, 
#               u_Id = if(!is.null(n.Id)){u_G}else{NULL}, 
#               Id = if(!is.null(n.Id)){Id0}else{NULL}, 
#               prior = prior)
#   )
#   
# }
# 
# 
# #### 3- Internal function for the fMM algorithm ####
# init1_fMM <- function(M, G, data, Id, coords, 
#                       formula_reg, prior_theta, prior_sigma, family, offset_reg, 
#                       formula_group, 
#                       test.GR, sigma_GR, proba_GR, proba_GRseed, 
#                       test.ICM, rho_ICM, G_ICM, 
#                       test.ICMregional){
#   
#   #### mise en forme #####
#   n <- dim(data)[1]
#   
#   # modele de reponse
#   if(!is.list(formula_reg))
#   {formula_reg <- list(formula_reg)}
#   
#   if(length(formula_reg) != M){
#     stop("init1_fMM[Fcts-fMM_segmentation.R] : \'formula_reg\' is not coherent with M \n", 
#          "M : ", M, "\n", 
#          "length(formula_reg) : ", length(formula_reg), "\n")}
#   
#   Var_reg <- list()
#   character_reg <- list()
#   
#   if(M == 1 && class(family) == "family")
#   {family <- list(family)}
#   
#   
#   for(iter_m in 1:M){
#     
#     if(class(family[[iter_m]]) == "family"){
#       family[[iter_m]] <- rep(list(family[[iter_m]]), G)
#     }
#     
#     if(!is.list(formula_reg[[iter_m]]))
#     {formula_reg[[iter_m]] <- list(formula_reg[[iter_m]])}
#     
#     if(length(formula_reg[[iter_m]]) == 1)
#     {formula_reg[[iter_m]] <- rep(formula_reg[[iter_m]], G)}
#     
#     if(length(formula_reg[[iter_m]]) != G)
#     {stop("init1_fMM[Fcts-fMM_segmentation.R] : \'formula_reg\'[[", iter_m, "]] must contain a formula for each group \n", 
#           "required number of formula : ", G, "\n", 
#           "proposed number of formula : ", length(formula_reg[[iter_m]]), "\n")}
#     
#     if(any(unlist(lapply(formula_reg[[iter_m]], function(x){"formula" %in% is(x) == FALSE}))))
#     {stop("init1_fMM[Fcts-fMM_segmentation.R] : \'formula_reg\'[[", iter_m, "]] must be of type formula \n", 
#           "\'formula_reg\' type proposed : ", paste(is(formula_reg[[iter_m]]), collapse = " "))}
#     
#     
#     Var_reg[[iter_m]] <- unique(unlist(lapply(formula_reg[[iter_m]], function(x){strsplit(all.vars(x), split = ":")})))
#     character_reg[[iter_m]] <- lapply(formula_reg[[iter_m]], function(x){if("formula" %in% is(x)){as.character(list(x))}else{x}})
#   }
#   
#   if(!is.null(offset_reg) && is.list(offset_reg) == FALSE)
#   {offset_reg <- list(offset_reg)}
#   
#   if(!is.null(offset_reg)){
#     for(iter_m in 1:M)
#     { if(length(offset_reg[[iter_m]]) == 1){offset_reg[[iter_m]] <- rep(offset_reg[[iter_m]], G)}
#       if(length(offset_reg[[iter_m]]) == G){offset_reg[[iter_m]] <- data[, offset_reg[[iter_m]]]}
#     }
#   }
#   
#   
#   # modele de structure
#  if(!is.null(formula_group)){
#     if("formula" %in% is(formula_group) == FALSE)
#     {stop("init1_fMM[Fcts-fMM_segmentation.R] : \'formula_group\' must be of type formula \n", 
#           "type de \'formula_group\' proposed : ", paste(is(formula_group), collapse = " "))}
#     
#     X_group <- attributes(terms(formula_group))$term.labels
#     
#     if(attributes(terms(formula_group))$intercept == 1){
#       intercept_group <- 1
#       formula_group <- paste("Y_group", paste(X_group, collapse = "+"), sep = "~")
#       
#     }else{
#       intercept_group <- 0
#       formula_group <- paste("Y_group", paste(c("-1", X_group), collapse = "+"), sep = "~")
#       
#     }
#     
#     Var_group <- unlist(strsplit(c("Y_group", X_group), split = ":"))
#   }else{
#     Var_group <- NULL
#     intercept_group <- NULL
#   }
#   
#   # Coords 
#   if(all(is.character(coords)) && all(coords %in% names(data)))
#   {coords <- data[, coords]}
#   
#   if(!is.null(coords) && nrow(as.matrix(coords)) != nrow(data)){
#     stop("init1_fMM[Fcts-fMM_segmentation.R] : \'coords\' is incompatible with \'data\' \n", 
#          "nrow of \'coords\' : ", nrow(as.matrix(coords)), "\n", 
#          "nrow of \'data\' : ", nrow(data), "\n")
#   }
#   
#   # ICM
#   if(test.ICM == TRUE){
#     if(all(is.character(rho_ICM))  && !identical(rho_ICM, "init") && any(rho_ICM %in% names(data) == FALSE)){
#       stop("init1_fMM[Fcts-fMM_segmentation.R] : wrong specification of \'rho_ICM\' \n", 
#            "must be \"init\" or numeric or names in data : ", paste(names(data), collapse = " "), "\n", 
#            "missing names : ", paste(rho_ICM[rho_ICM %in% names(data) == FALSE], collapse = " "), "\n")    
#     }
#     
#   }else{
#     resICMp <- NULL
#   }
#   
#   # GR
#   if(test.GR == TRUE){
#     maj_paramGR <- sigma_GR == "auto"
#   }else{
#     maj_paramGR <- NULL
#   }
#   
#   # Id
#   if(length(Id) == 1 && Id %in% names(data))
#   {Id <- data[, Id]}
#   
#   if(!is.null(Id) && length(Id) != nrow(data)){
#     stop("init1_fMM[Fcts-fMM_segmentation.R] : \'Id\' is incompatible with \'data\' \n", 
#          "length of \'Id\' : ", length(Id), "\n", 
#          "nrow of \'data\' : ", nrow(data), "\n")
#   }
#   
#   # prior
#   if(!is.null(prior_theta) && !identical(prior_theta, "kmeans") && !is.list(prior_theta))
#   {prior_theta <- list(prior_theta)}
#   if(!is.null(prior_sigma) && !is.list(prior_sigma))
#   {prior_sigma <- list(prior_sigma)}
#   
#   #### Export ####
#   return(list(n = n, 
#               coords = coords, 
#               formula_reg = formula_reg, 
#               Var_reg = Var_reg, 
#               character_reg = character_reg,              
#               offset_reg = offset_reg, 
#               Var_group = Var_group, 
#               formula_group = formula_group, 
#               intercept_group = intercept_group, 
#               Id = Id,              
#               maj_paramGR = maj_paramGR, 
#               prior_theta = prior_theta, 
#               prior_sigma = prior_sigma, 
#               family = family))
# }
# 
# test_fMM <- function(n, M, G, data, index_pat, coords, posterior_proba, 
#                      Var_reg, family, offset_reg, 
#                      formula_group, Var_group, 
#                      prior_theta, prior_sigma, prior_proba, 
#                      test.GR, seed, 
#                      test.ICM, test.ICMregional, Wdist_LR, Wweight_SR, G_ICM, rho_ICM, 
#                      export.predicteur){
#   
#   #### general ####
#   
#   if(G < 2){
#     stop("test_fMM[Fcts-fMM_segmentation.R] : wrong specification of \'G\' \n", 
#          "the number of groups \'G\' must be greater or equal to 2 \n", 
#          "proposed G : ", G, "\n")}
#   
#    if(!is.null(formula_group) && any(Var_group[-1] %in% names(data) == FALSE))
#   {stop("test_fMM[Fcts-fMM_segmentation.R] : \'Var_group\' must be in data \n", 
#         "missing names : ", paste(Var_group[-1][Var_group[-1] %in% names(data) == FALSE], collapse = " "), "\n")}
#   
#   if(!is.null(formula_group) && length(Var_group[-1]) > 0 && any(is.na(data[, Var_group[-1]])))
#   {stop("test_fMM[Fcts-fMM_segmentation.R] : variables of \'formula_reg\' contains NA in data \n", 
#         "variables with NA : ", paste(Var_group[-1][which(colSums(is.na(data[, Var_group[-1]])) > 0)], collapse = " "), " \n")}
#   
#   if(is.null(prior_proba) == FALSE && (is.matrix(prior_proba) == FALSE || n != nrow(prior_proba) || G != ncol(prior_proba)) )
#   {stop("test_fMM[Fcts-fMM_segmentation.R] : \'prior_proba\' must be a matrix with correct size \n", 
#         "required size (n, G) : ", n, " ", G, "\n", 
#         "proposed size (type) : ", paste(dim(prior_proba), collapse = " "), "(", paste(is(prior_proba), collapse = " "), ") \n")}
#   
#   if(is.null(posterior_proba) == FALSE && (is.matrix(posterior_proba) == FALSE || n != nrow(posterior_proba) || G != ncol(posterior_proba)) )
#   {stop("test_fMM[Fcts-fMM_segmentation.R] : \'posterior_proba\' must be a matrix with correct size \n", 
#         "required size (n, G) : ", n, " ", G, "\n", 
#         "proposed size (type) : ", paste(dim(posterior_proba), collapse = " "), "(", paste(is(posterior_proba), collapse = " "), ") \n")}
#   
#   if(is.null(prior_theta) == FALSE && !identical(prior_theta, "kmeans") && (is.list(prior_theta) == FALSE || length(prior_theta) != M || any(unlist(lapply(prior_theta, length)) != G) ))
#   {stop("test_fMM[Fcts-fMM_segmentation.R] : \'prior_theta\' must be a list with correct length \n", 
#         "Can be \"kmeans\" or a list of length : ", M, " with vectors of size ", G, "\n", 
#         "length - type : ", length(prior_theta), "|", paste(unlist(lapply(prior_theta, length)), collapse = " "), " - ", paste(is(prior_theta), collapse = " "), "\n")}
#   
#   if(is.null(prior_sigma) == FALSE && (is.list(prior_sigma) == FALSE || length(prior_sigma) != M || any(unlist(lapply(prior_sigma, length)) != G) ))
#   {stop("test_fMM[Fcts-fMM_segmentation.R] : \'prior_sigma\' must be a list with correct length\n", 
#         "Can be NULL or a list of length : ", M, " with vectors of size ", G, "\n", 
#         "length - type : ", length(prior_sigma), "|", paste(unlist(lapply(prior_theta, length)), collapse = " "), " - ", paste(is(prior_sigma), collapse = " "), "\n")}
#   
#   #### modele de reponse ####
#   if(length(family) != M){
#     stop("test_fMM[Fcts-fMM_segmentation.R] : \'family\' must be of length M \n", 
#          "M : ", M, "\n", 
#          "length(family) : ", length(family), "\n")}
#   
#   for(iter_m in 1:M){
#     
#     # family 
#     
#     if(length(family[[iter_m]]) != G || is.list(family[[iter_m]]) != TRUE)
#     {stop("test_fMM[Fcts-fMM_segmentation.R] : \'family\'[[", iter_m, "]] must be a list with correct length \n", 
#           "required length : ", G, "\n", 
#           "proposed length (type) : ", length(family[[iter_m]]), " (", is(family[[iter_m]]), ")\n")}
#     
#     if(sum(unlist(lapply(family[[iter_m]], function(x){x$family %in% c("Gamma", "gaussian", "binomial", "quasibinomial", "Beta") == FALSE}))) > 0)
#     {stop("test_fMM[Fcts-fMM_segmentation.R] : wrong specification of \'family\'[[", iter_m, "]]  \n", 
#           "available families : gaussian Gamma(link = \"log\") binomial quasibinomial Beta \n", 
#           "proposed family : ", paste(lapply(family[[iter_m]], function(x){x$family}), collapse = " "), "\n")}
#     
#     if(any(unlist(lapply(family[[iter_m]], function(x){x$family == "Beta"})))){
#       
# 	    initPackage("betareg", method = "test_fMM")		
#   
#       if(!is.null(offset_reg))
#       {stop("test_fMM : \'offset_reg\' must be NULL for Beta regression \n")}
#       
#       if(export.predicteur == FALSE)
#       { export.predicteur <- TRUE
#         warning("test_fMM[Fcts-fMM_segmentation.R] : \'export.predicteur\' is TRUE and should be FALSE \n", 
#                 "(because of potential bug when applying summary to beta regression) \n")}
#     }
#     
#     # variables reg
#     if( any(Var_reg[[iter_m]] %in% names(data) == FALSE)){
#       stop("test_fMM[Fcts-fMM_segmentation.R] : \'formula_reg\'[[", iter_m, "]] contains variables missing in data \n", 
#            "missing names : ", Var_reg[[iter_m]][Var_reg[[iter_m]] %in% names(data) == FALSE], "\n")}
#     
#     if(any(is.na(data[, Var_reg[[iter_m]]]))){
#       stop("test_fMM[Fcts-fMM_segmentation.R] : \'data\'[, formula_reg[[", iter_m, "]]] contains NA \n", 
#            "variables with NA : ", paste(Var_reg[[iter_m]][which(colSums(is.na(data[,Var_reg[[iter_m]], drop = FALSE])) > 0)], collapse = " "), " \n")}
#     
#     if(any(unlist(lapply(family[[iter_m]], function(x){x$family == "Beta"}))) && (min(data[, Var_reg[[iter_m]][1]]) <= 0 || max(data[, Var_reg[[iter_m]][1]]) >= 1))
#     {stop("test_fMM[Fcts-fMM_segmentation.R] : \'Y\'[[", iter_m, "]] takes its values out of ]0;1[ which is incompatible with a Beta law \n", 
#           "range \'Y\' : ", paste(range(data[, Var_reg[[iter_m]][1]]), collapse = " "), "\n")}
#     
#     if(any(unlist(lapply(family[[iter_m]], function(x){x$family == "Gamma"}))) && min(data[, Var_reg[[iter_m]][1]]) <= 0)
#     {stop("test_fMM[Fcts-fMM_segmentation.R] : \'Y\'[[", iter_m, "]] takes negative values which is incompatible with a Gamma law \n", 
#           "min \'Y\' : ", min(data[, Var_reg[[iter_m]][1]]), "\n")}
#     
#     # offset
#     if(is.null(offset_reg) == FALSE && (is.matrix(offset_reg[[iter_m]]) == FALSE || n != nrow(offset_reg[[iter_m]]) || G != ncol(offset_reg[[iter_m]])) )
#     {stop("test_fMM[Fcts-fMM_segmentation.R] : \'offset_reg\'[[", iter_m, "]]  must be a matrix with correct length \n", 
#           "required size of \'offset_reg\' (n, G) : NULL or ", n, " ", G, "\n", 
#           "size of \'offset_reg\' (type) : ", paste(dim(offset_reg[[iter_m]]), collapse = " "), "(", paste(is(offset_reg[[iter_m]]), collapse = " "), ") \n")}
#     
#     
#   }
#   
#   #### spatial ####
#   
#   if(test.GR || test.ICM){
#     
#     if(is.null(Wweight_SR)){
#       stop("test_fMM[Fcts-fMM_segmentation.R] : \'Wweight_SR\' must be specified for GR \n")
#     }
#     
#     if(any(unlist(lapply(Wweight_SR, function(x){x@Dim})) != unlist(lapply(index_pat, function(x){rep(length(x), 2)}))))
#     { err <-  union(which(unlist(lapply(Wweight_SR, function(x){x@Dim[1]})) != unlist(lapply(index_pat, length))), 
#                     which(unlist(lapply(Wweight_SR, function(x){x@Dim[2]})) != unlist(lapply(index_pat, length))))
#       stop("test_fMM[Fcts-fMM_segmentation.R] : incorrect size of \'Wweight_SR\' \n", 
#            "patient : ", paste(err, collapse = " "), "\n", 
#            "required size : ", paste(unlist(lapply(index_pat[err], length)), collapse = " "), "\n", 
#            "proposed row size : ", paste(unlist(lapply(Wweight_SR, function(x){x@Dim[1]})), collapse = " "), "\n", 
#            "proposed column size : ", paste(unlist(lapply(Wweight_SR, function(x){x@Dim[2]})), collapse = " "), "\n")
#     }
#     
#   }
#   
#   if(test.GR && is.null(seed)){
#     stop("test_fMM[Fcts-fMM_segmentation.R] : \'seed\' must be specified for GR \n")
#   }
#   
#   if(test.ICM && (length(G_ICM) != G || any(G_ICM %in% 1:G == FALSE))){
#     stop("test_fMM[Fcts-fMM_segmentation.R] : incorrect specification of \'G_ICM\' that indicates combination of groups that should be used for ICM \n", 
#          "\'G_ICM\' must be of length ", G, " and contain numbers between 1 and ", G, "\n", 
#          "\'G_ICM\' proposed : ", paste(G_ICM, collapse = " "), "\n")
#   }
#   
#   if(test.ICM && test.ICMregional){
#     
#     if(is.null(Wdist_LR)){
#       stop("test_fMM[Fcts-fMM_segmentation.R] : \'Wdist_LR\' must be specified for regional ICM  \n")
#     }
#     
#     if(!is.null(Wdist_LR)){
#       if(any(unlist(lapply(Wdist_LR, function(x){x@Dim})) != unlist(lapply(index_pat, function(x){rep(length(x), 2)}))))
#       { err <-  union(which(unlist(lapply(Wdist_LR, function(x){x@Dim[1]})) != unlist(lapply(index_pat, length))), 
#                       which(unlist(lapply(Wdist_LR, function(x){x@Dim[2]})) != unlist(lapply(index_pat, length))))
#         stop("test_fMM[Fcts-fMM_segmentation.R] : incorrect size of \'Wdist_LR\' \n", 
#              "patient : ", paste(err, collapse = " "), "\n", 
#              "required size : ", paste(unlist(lapply(index_pat[err], length)), collapse = " "), "\n", 
#              "proposed row size : ", paste(unlist(lapply(Wdist_LR, function(x){x@Dim[1]})), collapse = " "), "\n", 
#              "proposed column size : ", paste(unlist(lapply(Wdist_LR, function(x){x@Dim[2]})), collapse = " "), "\n")
#       }
#     }
#     
#     if(is.null(coords)){
#       stop("test_fMM[Fcts-fMM_segmentation.R] : \'coords\' must be specified for regional ICM \n")
#     }
#     
#   }
#   
#   
#   #### export ####
#   
#   return(list(export.predicteur = export.predicteur))
#   
# }
# 
# init2_fMM <- function(n, M, G, data, Id, 
#                       formula_reg, Var_reg, offset_reg, family, 
#                       formula_group, Var_group, intercept_group, 
#                       prior_theta, prior_sigma, prior_proba, posterior_proba, 
#                       test.GR, seed, 
#                       test.ICM, Wweight_SR, prior_prevalenceICM, rho_max, rho_ICM, 
#                       test.ICMregional, distance_ref, Wdist_LR, coords, threshold, nbGroup_min, multiV, 
#                       iter_max, verbose){
#   
#   
#   #### preparation ####  
#   var <- unlist(lapply(Var_reg, function(x){var(data[, x[1]])}))
#   data <-data.frame("intercept" = 1, data)
#   
#   #### modele de reponse ####
#   nb.param_reg <- list()
#   theta <- list()
#   hist_theta <- list()
#   sigma <- list()
#   lm_adj <- list()
#   
#   for(iter_m in 1:M){  
#     nb.param_reg[[iter_m]] <- length(Var_reg[[iter_m]]) - 1 + max(unlist(lapply(formula_reg[[iter_m]], function(x){if("formula" %in% is(x)){attributes(terms(x))$intercept}})))
#     theta[[iter_m]] <- matrix(NA, nrow = nb.param_reg[[iter_m]], ncol = G)
#     row.names(theta[[iter_m]]) <- c(if(max(unlist(lapply(formula_reg[[1]], function(x){attributes(terms(x))$intercept}))) == 1){"(Intercept)"}else{NULL}, 
#                                     Var_reg[[iter_m]][-1])
#     hist_theta[[iter_m]] <- matrix(NA, nrow = nb.param_reg[[iter_m]], ncol = G)
#     sigma[[iter_m]] <- rep(NA, G)
#     lm_adj[[iter_m]] <- list() 
#     
#     # adaptation de l offset a la loi utilisee
#     if(!is.null(offset_reg[[iter_m]])){
#       for(iter_g in 1:G){
#         if(family[[iter_m]][[iter_g]]$family == "Gamma")
#         {offset_reg[[iter_m]][, iter_g] <- log(offset_reg[[iter_m]][, iter_g])}  
#       }
#     }
#   }
#   data_reg <- data[,unique(unlist(Var_reg)), drop = FALSE]
#   
#   #### modele de structure ####
#   if(!is.null(formula_group)){
#     nb.param_group <- length(Var_group) - 1 + intercept_group
#     beta <- matrix(NA, nrow = G - 1, ncol = nb.param_group)
#     hist_beta<- matrix(NA, nrow = G - 1, ncol = nb.param_group)
#     
#     Y_group <- matrix(0, nrow = G * n, ncol = G)
#     
#     data_group <- data.frame(matrix(0, nrow = dim(data)[1]*G, ncol = length(Var_group)))
#     
#     for(iter_g in 1:G){
#       Y_group[seq(1 + n * (iter_g - 1), n * iter_g, by = 1), iter_g] <- 1
#       data_group[seq(1 + n * (iter_g - 1), n * iter_g, by = 1), 1] <- iter_g
#       data_group[seq(1 + n * (iter_g - 1), n * iter_g, by = 1), 2:length(Var_group)] <- data[, Var_group[-1]]
#     }
#     names(data_group) <- c("Y_group", Var_group[-1])
#   }else{
#     beta <- -1000
#     hist_beta <- -1000
#     data_group <- NULL
#     nb.param_group <- NULL
#   }
#   
#   #### initialisation retraitement spatial ####
#   init_corrSpat <- FALSE
#   
#   # initialisation ICM
#   if(test.ICM == TRUE){
#     
#     if(all(is.character(rho_ICM)) && !identical(rho_ICM, "init")){
#       
#       rho_ICM <- calcPotts(W_LR = Wdist_LR, W_SR = Wweight_SR, 
#                            sample = NULL, test.ICMregional = test.ICMregional, threshold = threshold, 
#                            prior_prevalence = prior_prevalenceICM, 
#                            rho = "auto", rho_max = rho_max, 
#                            Y = data[, rho_ICM], distance_ref = distance_ref, coords = coords, 
#                            Id = Id, nbGroup_min = nbGroup_min, distcat = TRUE, multiV = multiV, 
#                            trace_iter = FALSE, trace_radius = TRUE, trace_rho = TRUE, iter_max = 150, critere = 0.005)$rho
#     }
#     
#   }else{
#     rho_ICM <- NULL
#   }
#   
#   
#   
#   # modele GR
#   if(test.GR == TRUE){
#     if(length(seed) == 1 && is.character(seed))
#     {seed_GR <- tapply(data[, seed], Id, function(x){which(x == 1)})}
#     
#     nbVois.max <- lapply(Wweight_SR, function(x){max(spam::rowSums(x > 0))})
#   }else{
#     seed_GR <- NULL
#     nbVois.max <- NULL
#   }
#   
#   
#   #### initialisation des valeurs ####
#   
#   # theta / sigma
#   if(!is.null(posterior_proba)){
#     cat("* initialisation using posterior proba *\n ")
#     prior_proba_spat <- NULL 
#     Lv_obs <- NULL
#   }else{
#     
#     if(is.null(prior_proba))
#     {prior_proba <- matrix(1 / G, nrow = n, ncol = G)}
#     
#     if(is.null(prior_theta)){
#       if(verbose){cat("* initialisation by random sampling *\n ")}
#       prior_theta <- list() 
#       
#       for(iter_m in 1:M)
#       {  iter_reinit <- 1
#          prior_theta[[iter_m]] <- rep(NA, G)
#          while(iter_reinit == 1 || (iter_reinit < 5 &&  crit_tempo ))
#          {prior_theta[[iter_m]] <- runif(G, min(data[, Var_reg[[iter_m]][1]]), max(data[, Var_reg[[iter_m]][1]]))
#           iter_reinit <- iter_reinit + 1
#           if(sum(!is.na(prior_theta[[iter_m]])) == 1)
#           {crit_tempo <- 0}else
#           {crit_tempo <- min(abs(dist(prior_theta[[iter_m]])), na.rm = TRUE) < sqrt(var[iter_m]) / (G - 1)}
#          }
#       }
#     }else{
#       if(identical(prior_theta, "kmeans")){    
#         if(verbose){cat("* initialisation by k means *\n ")}
#         res_kmeans <- kmeans(x = data[, unlist(lapply(Var_reg, function(x){x[1]}))], 
#                              centers = G)
#         ordre_kmeans <- order(res_kmeans$centers[, 1])
#         
#         prior_theta<- list()
#         for(iter_m in 1:M)
#         { prior_theta[[iter_m]] <- as.vector(res_kmeans$center[ordre_kmeans, iter_m]) }
#         
#         if(is.null(prior_sigma)){
#           prior_sigma <- list()
#           for(iter_m in 1:M)
#           { prior_sigma[[iter_m]] <- tapply(data[, Var_reg[[iter_m]][1]], res_kmeans$cluster, sd) }
#         }
#         
#       }else{
#         for(iter_m in 1:M){
#           if(is.null(offset_reg)){
#             if(min(data[, Var_reg[[iter_m]][1]], na.rm = TRUE) > min(prior_theta[[iter_m]], na.rm = TRUE) || max(data[, Var_reg[[iter_m]][1]], na.rm = TRUE) < max(prior_theta[[iter_m]], na.rm = TRUE) )
#             { warning("init2_fMM[Fcts-fMM_segmentation.R] : \'prior_theta\'[[", iter_m, "]] is not correctly specified \n", 
#                       "rang of \'Y\'[[", iter_m, "]] : ", paste(round(range(data[, Var_reg[[iter_m]][1]]), digits = 2), collapse = " "), " \n", 
#                       "proposed initialisation : ", paste(prior_theta[[iter_m]], collapse = " "), "\n") }
#           }
#         }
#         
#       }
#     }
#     
#     
#     if(is.null(prior_sigma)){
#       prior_sigma <- list()
#       for(iter_m in 1:M)
#       { prior_sigma[[iter_m]] <- rep(sqrt(var[iter_m]/G), G)}
#     }
#     
#     # Lv  
#     prior_proba_spat <- NULL 
#     Lv_obs <- matrix(1, nrow = n, ncol = G)  
#     
#     for(iter_g in 1:G){  
#       for(iter_m in 1:M){
#         
#         if(is.null(offset_reg[[iter_m]])){
#           lm_adj[[iter_m]][[iter_g]] <- data.frame(residuals = data[, Var_reg[[iter_m]][1]]-prior_theta[[iter_m]][iter_g], 
#                                                    fitted.values = prior_theta[[iter_m]][iter_g])
#         }else{
#           if(family[[iter_m]][[iter_g]]$family == "Gamma"){
#             lm_adj[[iter_m]][[iter_g]] <- data.frame(residuals = data[, Var_reg[[iter_m]][1]]-exp(offset_reg[[iter_m]][, iter_g]) - prior_theta[[iter_m]][iter_g], 
#                                                      fitted.values = exp(offset_reg[[iter_m]][, iter_g]) + prior_theta[[iter_m]][iter_g])
#           }else
#           {lm_adj[[iter_m]][[iter_g]] <- data.frame(residuals = data[, Var_reg[[iter_m]][1]]-offset_reg[[iter_m]][, iter_g]-prior_theta[[iter_m]][iter_g], 
#                                                     fitted.values = offset_reg[[iter_m]][, iter_g]+prior_theta[[iter_m]][iter_g])}
#         }
#         
#         
#         if(family[[iter_m]][[iter_g]]$family %in% c("binomial", "quasibinomial"))
#         {sigma[[iter_m]][iter_g] <- sqrt(prior_theta[[iter_m]][iter_g] * (1 - prior_theta[[iter_m]][iter_g]))}else
#         {sigma[[iter_m]][iter_g] <- prior_sigma[[iter_m]][iter_g]}
#         
#         if(family[[iter_m]][[iter_g]]$family == "gaussian")
#         { Lv_obs[, iter_g] <- Lv_obs[, iter_g]  *  dnorm(x = data[, Var_reg[[iter_m]][1]], 
#                                                        mean = fitted.values(lm_adj[[iter_m]][[iter_g]]), 
#                                                        sd = sigma[[iter_m]][iter_g], 
#                                                        log = FALSE)  }
#         
#         if(family[[iter_m]][[iter_g]]$family == "Gamma")
#         { Lv_obs[, iter_g] <- Lv_obs[, iter_g] *  dgamma(x = data[, Var_reg[[iter_m]][1]], 
#                                                        scale = sigma[[iter_m]][iter_g]^2 / fitted.values(lm_adj[[iter_m]][[iter_g]]), 
#                                                        shape = fitted.values(lm_adj[[iter_m]][[iter_g]])^2 / sigma[[iter_m]][iter_g]^2, 
#                                                        log = FALSE)}  
#         
#         if(family[[iter_m]][[iter_g]]$family == "Beta")
#         { Lv_obs[, iter_g] <- Lv_obs[, iter_g] *  dbeta(x = data[, Var_reg[[iter_m]][1]], 
#                                                       shape1=fitted.values(lm_adj[[iter_m]][[iter_g]]) * (fitted.values(lm_adj[[iter_m]][[iter_g]]) * (1 - fitted.values(lm_adj[[iter_m]][[iter_g]])) / sigma[[iter_m]][iter_g]^2 - 1), 
#                                                       shape2=(1-fitted.values(lm_adj[[iter_m]][[iter_g]])) * (fitted.values(lm_adj[[iter_m]][[iter_g]]) * (1 - fitted.values(lm_adj[[iter_m]][[iter_g]])) / sigma[[iter_m]][iter_g]^2 - 1), 
#                                                       log = FALSE)}  
#         
#         if(family[[iter_m]][[iter_g]]$family %in% c("binomial", "quasibinomial"))
#         { Lv_obs[, iter_g] <- Lv_obs[, iter_g] * dbinom(x = data[, Var_reg[[iter_m]][1]], 
#                                                       size = 1, prob = fitted.values(lm_adj[[iter_m]][[iter_g]]), 
#                                                       log = FALSE)}
#         
#       }
#     }
#   }
#   
#   #### Export ####
#   return(list(prior_theta = prior_theta, 
#               prior_sigma = prior_sigma, 
#               theta = theta, 
#               hist_theta = hist_theta, 
#               sigma = sigma, 
#               lm_adj = lm_adj, 
#               data_reg = data_reg, 
#               hist_beta = hist_beta, 
#               beta = beta, 
#               data_group = data_group, 
#               nb.param_group = nb.param_group, 
#               Lv_obs = Lv_obs, 
#               prior_proba = prior_proba, 
#               prior_proba_spat = prior_proba_spat, 
#               init_corrSpat = init_corrSpat, 
#               seed_GR = seed_GR, 
#               nbVois.max = nbVois.max, 
#               rho_ICM = rho_ICM))
#   
# }
# 
