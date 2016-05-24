# #### 1- User interface for the segmentation algorithm ####
# 
# launcher_fMM <- function(G, data, 
#                          Wweight_SR = NULL, distband.SR = NA, 
#                          Wdist_LR = NULL, distband.LR = NA, 
#                          var_reg, family = gaussian(link = "identity"), prior_theta = "kmeans", posterior_proba = NULL, 
#                          var_group = NULL, prior_prevalence = FALSE, 
#                          test.GR = FALSE, seed = NULL, 
#                          test.ICM = FALSE, rho_ICM = "auto", G_ICM = 1:G, update_rho = FALSE, 
#                          test.ICMregional = FALSE, coords = NULL, distance_ref = NULL, multiV = FALSE, 
#                          epsilon = 10^(-3), epsilon_corrSpat = epsilon * 10, iter_max = 100, 
#                          trace_iter = FALSE){
#   
#   #### gestion de Id et regresseurs ####
#   regresseurs <- "~ 1"
#   
#   #### modele de reponse ####
#   M <- length(var_reg)  
#   
#   if(is.list(var_reg) == FALSE){
#     formula_reg <- lapply(var_reg, 
#                           function(x){replicate(G, as.formula(paste(x, regresseurs, sep = "")), simplify = FALSE)}
#     )
#   }else{
#     formula_reg <- var_reg
#   }
#   
#   if(class(family) == "family"){
#     family <- replicate(M, family, simplify = FALSE)
#   }
#   if(length(family) != M){
#     stop("launcher_fMM[Fcts-fMM_segmentation.R] : \'family\' mispecified \n", 
#          "must have the same length as \'var_reg\' : ", M, " \n", 
#          "proposed length : ", length(family), "\n")
#   }
#   family  <- lapply(family, function(x){replicate(G, x, simplify = FALSE)})
#   
#   #### modele de structure ####
#   if(!is.null(var_group)){
#     formula_group <- as.formula(paste("~", paste(var_group, collapse = "+"), sep = ""))
#   }else{
#     formula_group <- NULL
#   }
#   
#   #### coords ####
#   if(all(is.character(coords))){
#     
#     if( any(coords %in% names(data) == FALSE)){
#       stop("launcher_fMM[Fcts-fMM_segmentation.R] : \'coords\' mispecified \n", 
#            "if character vector then must match the name of \'data\' : ", paste(names(data), collapse = " "), " \n", 
#            "proposed \'coords\' : ", paste(coords, collapse = " "), "\n")
#     }
#     
#     coords <- data[,coords]
#   }
#   
#   #### W ####
#   if(is.null(Wdist_LR)){
#     if(test.ICM && test.ICMregional){
#       
#       #Wdist_LR <- calcW_fMM(coords = coords, distband = distband.LR, 
#       #                     method = "euclidean", upper = NULL, format = "dgCMatrix", rowNorm = FALSE)
#       Wdist_LR <- calcW(coords, range == distband.LR, method = "euclidean", upper = NULL, format = "dgCMatrix", row.norm = FALSE)
# 	  
#     }else{
#       Wdist_LR <- NULL
#     }
#   }
#   
#   #### Wweight_SR ####  
#   if(is.null(Wweight_SR)){    
#     if(test.ICM || test.GR){
#       
# 	  Wweight_SR <- calcW(coords, range == distband.SR, method = "euclidean", upper = NULL, format = "dgCMatrix", row.norm = TRUE)
# 
#       #Wweight_SR <- calcW_fMM(coords = coords, distband = distband.SR, 
#       #                       method = "euclidean", upper = NULL, format = "dgCMatrix", rowNorm = TRUE)
#       
#     }else{
#       Wweight_SR <- NULL
#     }
#   }
#   
#   #### estimation ####
#   
#   res <- sfMM(M = M, G = G, data, Wweight_SR = Wweight_SR, Wdist_LR = Wdist_LR, posterior_proba = posterior_proba, 
#                          formula_reg = formula_reg, offset_reg = NULL, family = family, prior_theta = prior_theta, prior_sigma = NULL, 
#                          formula_group = formula_group, prior_prevalence = prior_prevalence, prior_proba = NULL, 
#                          test.GR = test.GR, sigma_GR = "auto", proba_GR = 0.1, proba_GRseed = 0.25, seed = seed, 
#                          test.ICM = test.ICM, rho_ICM = rho_ICM, G_ICM = G_ICM,                            
#                          test.ICMregional = test.ICMregional, coords = coords, nbGroup_min = 100, threshold_potential = 0.1, distance_ref = distance_ref, multiV = multiV, 
#                          prior_prevalenceICM = TRUE, rho_max = 10, update_rho = update_rho, 
#                          digit.pmin = 7, epsilon = epsilon, epsilon_corrSpat = epsilon_corrSpat, iter_max = iter_max, 
#                          verbose = TRUE, trace_iter = trace_iter, trace_radius = FALSE, export.predicteur = TRUE)
#   
#   #### export ####
#   return(res)
#   
# }
# 
# #### 2- postprocessing the results of the segmentation algorithm ####
# 
# plotProba_fMM <- function(res.EM, type = "posterior_proba", 
#                           main = "group", col = c("red", "blue", "green", "yellow", "purple", "brown"), border = "black", cex.main = 0.75, 
#                           window = FALSE, 
#                           filename = "EM.plotProb", width = 480, height = 480, path = NULL, unit = "px", res = NA)
# {  ## tests preliminaires
#   
#   if(type %in% c("prior_proba", "posterior_proba") == FALSE)
#   {stop("plotProba_fMM[Fcts-fMM_segmentation.R] : wrong specification of \'type\' \n", 
#         "available \'type\' :  \"prior_proba\" \"posterior_proba\" \n", 
#         "proposed \'type\' : ", type, "\n")}
#   
#   if(type %in% names(res.EM) == FALSE){
#     stop("plotProba_fMM[Fcts-fMM_segmentation.R] : wrong specification of \'res.EM\' \n", 
#          "\'res.EM\' must be a list containing the field : ", type , " \n", 
#          "proposed \'res.EM\' (names(res.EM) / is(res.EM)) : ", paste(names(res.EM), collapse = " "), " / ", paste(is(res.EM), collapse = " "), "\n")
#   }
#   
#   if(unit %in% c("px", "in", "cm", "mm") == FALSE){
#     stop("plotProba_fMM[Fcts-fMM_segmentation.R] : wrong specification of  \'unit\' \n", 
#          "available \'unit\' : \"px\" \"in\" \"cm\" \"mm\" \n", 
#          "proposed \'unit\' : ", unit, "\n")
#   }
#   
#   if(window %in% c(TRUE, FALSE, "eps", "svg", "png") == FALSE){
#     stop("plotProba_fMM[Fcts-fMM_segmentation.R] : wrong specification of \'window\' \n", 
#          "available \'window\' : TRUE FALSE \"eps\" \"svg\" \"png\" \n", 
#          "proposed \'window\' : ", window, "\n")
#   }
#   
#   prob <- res.EM[[type]]
#   G <- dim(prob)[2]
#   
#   scale <- switch(unit, 
#                   "px" = 1, 
#                   "in" = 90, 
#                   "cm" = 35.43, 
#                   "mm" = 3.543)
#   
#   switch(window, 
#          T = dev.new()(), 
#          eps = postscript(file = paste(path, filename, ".eps", sep = ""), width = width * scale / 90, height = height * scale / 90, horizontal = FALSE, onefile = FALSE, paper = "special"), 
#          svg = svg(filename = paste(path, filename, ".svg", sep = ""), width = width * scale / 90, height = height * scale / 90, onefile = FALSE), 
#          png = png(filename = paste(path, filename, ".png", sep = ""), width = width * scale, height = height * scale, res = res)
#   )
#   
#   if(length(col) < G){
#     stop("plotProba_fMM[Fcts-fMM_segmentation.R] : wrong specification of  \'col\' \n", 
#          "\'col\' must have length ", G, " \n", 
#          "proposed \'col\' : ", paste(col, collapse = " "), "\n")
#   }
#   
#   col <- col[1:G]
#   
#   if(length(border) == 1)
#   {border <- rep("black", G)}
#   
#   if(length(border) < G){
#     stop("plotProba_fMM[Fcts-fMM_segmentation.R] : wrong specification of \'border\' \n", 
#          "\'border\' must have length ", G, " \n", 
#          "proposed \'border\' : ", paste(border, collapse = " "), "\n")
#   }
#   
#   par(mfrow = c(G, 1), mai = rep(0.75, 4), mar = c(3, 3, 1, 1), mgp = c(2, 0.75, 0))
#   for(iter_g in 1:G){
#     hist(prob[,iter_g], 
#          main = paste(main, iter_g), cex.main = cex.main, 
#          col = col[iter_g], border = border[iter_g], 
#          xlab = "", xlim = c(0, 1))
#   }
#   
#   switch(window, 
#          eps = dev.off(), 
#          svg = dev.off(), 
#          png = dev.off()
#   )
#   
# }
# 
# plotCv_fMM <- function(res.EM, 
#                        window = FALSE, 
#                        filename = "EM.traceLv", width = 480, height = 480, path = NULL, unit = "px", res = NA){
#   
#   
#   #### tests ####
#   if("critere_cv" %in% names(res.EM)){
#     res.EM <- res.EM$critere_cv
#   }
#   
#   if(any(c("lv", "lv_completee", "critere_param") %in% names(res.EM) == FALSE)){
#     stop("plotCv_fMM[Fcts-fMM_segmentation.R] :  wrong specification of \'res.EM\' \n", 
#          "\'res.EM\' must be a list containing the fields : \"lv\" \"lv_completee\" \"critere_param\" \n", 
#          "proposed \'res.EM\' (names(res.EM) / is(res.EM)) : ", paste(names(res.EM), collapse = " "), " / ", paste(is(res.EM), collapse = " "), "\n")
#   }
#   
#   if(unit %in% c("px", "in", "cm", "mm") == FALSE){
#     stop("plotCv_fMM[Fcts-fMM_segmentation.R] : wrong specification of  \'unit\' \n", 
#          "available \'unit\' : \"px\" \"in\" \"cm\" \"mm\" \n", 
#          "proposed \'unit\' : ", unit, "\n")
#   }
#   
#   if(window %in% c(TRUE, FALSE, "eps", "svg", "png") == FALSE){
#     stop("plotCv_fMM[Fcts-fMM_segmentation.R] : wrong specification of  \'window\' \n", 
#          "available \'window\' : TRUE FALSE \"eps\" \"svg\" \"png\" \n", 
#          "proposed \'window\' : ", window, "\n")
#   }
#   
#   scale <- switch(unit, 
#                   "px" = 1, 
#                   "in" = 90, 
#                   "cm" = 35.43, 
#                   "mm" = 3.543)
#   
#   #### preparation ####
#   index_NA <- intersect(intersect(which(is.na(res.EM$lv)), which(is.na(res.EM$lv_completee))), 
#                         which(is.na(res.EM$critere_param)))
#   if(length(index_NA) > 0){
#     res.EM$lv <- res.EM$lv[-index_NA]
#     res.EM$lv_completee <- res.EM$lv_completee[-index_NA]
#     res.EM$critere_param <- res.EM$critere_param[-index_NA]
#   }
#   
#   iter <- length(na.omit(res.EM$lv))
#   
#   #### affichage
#   if(!is.null(window)){
#     switch(window, 
#            T = dev.new(), 
#            eps = postscript(file = paste(path, filename, ".eps", sep = ""), width = width * scale / 90, height = height * scale / 90, horizontal = FALSE, onefile = FALSE, paper = "special"), 
#            svg = svg(filename = paste(path, filename, ".svg", sep = ""), width = width * scale / 90, height = height * scale / 90, onefile = FALSE), 
#            png = png(filename = paste(path, filename, ".png", sep = ""), width = width * scale, height = height * scale, res = res)
#     )
#     
#     par(mfrow = c(2, 1), mar = c(0.5, 0.5, 0.5, 0.5), mai = c(0.75, 1, 0.5, 0.5), mgp = c(1.25, 0.4, 0))
#     
#     if(sum(is.na(res.EM$lv) == FALSE) > 0 && sum(is.infinite(res.EM$lv) == FALSE) > 0){
#       plot(1:iter, res.EM$lv[1:iter], xlab = "iteration number", ylab = expression(paste(L[v])), 
#            main = expression(paste("cv algo EM (", L[v], ")", sep = "")), 
#            cex.main = 0.9, cex.lab = 1, cex.axis = 0.6, xaxs = "r", type = "b", pch = 20, axes = FALSE)
#       points(1:iter, res.EM$lv_completee[1:iter], col = "grey", type = "b", lty = 3, cex = 0.75)
#       axis(1, at = 1:iter, labels = as.character(0:(iter - 1)))
#       axis(2)
#       legend("bottomright", legend = c(expression(L[v]), expression(paste(L[v], "completed", sep = ""))), col = c("black", "grey"), bty = "n", pch = c(20, 21))}else
#       {plot(1, 1, col = "white", xlab = "iteration number", ylab = expression(paste(Delta, L[v])), 
#             main = expression(paste("algo EM cv (crit = ", Delta, " ", L[v], ")", sep = "")))
#        legend("center", legend = "criterion lv unavailable")}
#     
#     if(sum(is.na(res.EM$critere_param) == FALSE) > 0){
#       plot(1:iter, res.EM$critere_param[1:iter], xlab = "iteration number", ylab = expression(paste(Delta, "param")), 
#            main = expression(paste("cv algo EM (crit = ", Delta, "param)", sep = "")), 
#            cex.main = 0.9, cex.lab = 1, cex.axis = 0.6, xaxs = "r", type = "b", pch = 20, axes = FALSE)
#       axis(1, at = 1:iter, labels = as.character(0:(iter - 1)))
#       axis(2)}else
#       {plot(1, 1, col = "white", xlab = "iteration number", ylab = expression(paste(Delta, L[v])), cex.main = 0.9, 
#             main = expression(paste("cv algo EM (crit = ", Delta, " ", L[v], ")", sep = "")))
#        legend("center", legend = "criterion param unavailable")}
#     
#     switch(window, 
#            eps = dev.off(), 
#            svg = dev.off(), 
#            png = dev.off()
#     )
#     
#   }
#   
#   #### test de croissance de la vraisemblance
#   diff_lv <- tail(res.EM$lv, -1) - head(res.EM$lv, -1)
#   test.croissance <- all(diff_lv >= 0)
#   
#   diff_lv_completee <- tail(res.EM$lv_completee, -1) - head(res.EM$lv_completee, -1)
#   test.croissance_completee <- all(diff_lv_completee >= 0)
#   
#   return(list(diff_lv = diff_lv, 
#               test.croissance = test.croissance, 
#               diff_lv_completee = diff_lv_completee, 
#               test.croissance_completee = test.croissance_completee))
# }
# 
# 
# plotY_fMM <- function(res.EM, prob = "posterior_proba", type = "density", 
#                       numY = 1, identifier = NULL, 
#                       nclass = NULL, x.legend = "topright", col = NULL, main = NULL, cex.main = 1, lwd = 1, lwd.mean = 2, lty.wmean = 2, 
#                       window = FALSE, 
#                       filename = "plotY_fMM", width = 480, height = 480, path = NULL, unit = "px", res = NA)
# { 
#   
#   EM.weighted.sd <- function(x, w){
#     
#     if(sum(w) <= 1){
#       stop("EM.weighted.sd[Fcts-fMM_segmentation.R] :  wrong specification of \'w\' \n", 
#            "sum(w) must be > 1\n", 
#            "proposed sum(w) : ", sum(w), "\n")
#     }
#     
#     mean <- weighted.mean(x = x, w = w, na.rm = T)
#     sd <- sqrt( sum(w * (x-mean)^2) / (sum(w) - 1) )
#     return(sd)
#   }
#   
#   ## tests preliminaires
#   if(prob %in% c("group", "prior_proba", "posterior_proba") == FALSE)
#   {stop("plotY_fMM[Fcts-fMM_segmentation.R] :  wrong specification of \'prob\' \n", 
#         "available \'prob\' :  \"group\" \"prior_proba\" \"posterior_proba\" \n", 
#         "proposed \'prob\' : ", prob, "\n")}
#   
#   if(type %in% c("hist", "density") == FALSE)
#   {stop("plotY_fMM[Fcts-fMM_segmentation.R] :  wrong specification of \'type\' \n", 
#         "available \'type\' :  \"hist\" \"density\" \n", 
#         "proposed \'type\' : ", type, "\n")}
#   
#   if(type == "hist" && prob != "group")
#   {stop("plotY_fMM[Fcts-fMM_segmentation.R] :  wrong specification of \'type\' \n", 
#         "\'hist\' is only available with \'prob\'=\"group\" \n", 
#         "proposed (\'type\'/\'prob\') : ", type, "/", prob, "\n")}
#   
#   if(any(c("Y", "family") %in% names(res.EM) == FALSE)){
#     stop("plotY_fMM[Fcts-fMM_segmentation.R] :  wrong specification of \'res.EM\' \n", 
#          "\'res.EM\' must be a list containing the fields : Y family \n", 
#          "proposed \'res.EM\' (names(res.EM) / is(res.EM)) : ", paste(names(res.EM), collapse = " "), " / ", paste(is(res.EM), collapse = " "), "\n")
#   }
#   
#   res.EM$Y <- as.matrix(res.EM$Y)
#   
#   if(numY > ncol(res.EM$Y)){
#     stop("plotY_fMM[Fcts-fMM_segmentation.R] :  wrong specification of \'numY\' \n", 
#          "there are ", ncol(res.EM$Y), " response variables \n", 
#          "proposed \'numY\' : ", numY, "\n")
#   }
#   
#   if(unit %in% c("px", "in", "cm", "mm") == FALSE){
#     stop("plotY_fMM[Fcts-fMM_segmentation.R] : wrong specification of  \'unit\' \n", 
#          "available \'unit\' : \"px\" \"in\" \"cm\" \"mm\" \n", 
#          "proposed \'unit\' : ", unit, "\n")
#   }
#   
#   if(window %in% c(TRUE, FALSE, "eps", "svg", "png") == FALSE){
#     stop("plotY_fMM[Fcts-fMM_segmentation.R] : wrong specification of  \'window\' \n", 
#          "available \'window\' : TRUE FALSE \"eps\" \"svg\" \"png\" \n", 
#          "proposed \'window\' : ", window, "\n")
#   }
#   
#   scale <- switch(unit, 
#                   "px" = 1, 
#                   "in" = 90, 
#                   "cm" = 35.43, 
#                   "mm" = 3.543)
#   
#   if(is.null(identifier) == FALSE && (length(identifier) > 1 || identifier %in% unique(res.EM$Id$Id) == FALSE))
#   {stop("\'identifier\' est mal specifie \n", 
#         "possible values  :  NULL", unique(res.EM$Id), "\n", 
#         "proposed value (nb) : ", identifier, " (", length(identifier), ")\n")}
#   
#   ## definition des variables d interet
#   if(is.null(identifier) == FALSE){
#     index_id <- which(res.EM$Id$Id == identifier)
#   }else{    
#     index_id <- 1:length(res.EM$Y[,1])
#   }
#   
#   G <-  switch(prob, 
#                "prior_proba" = ncol(res.EM$prior_proba), 
#                "posterior_proba" = ncol(res.EM$posterior_proba), 
#                "group" = length(res.EM$group)
#   )
#   
#   proba <-  switch(prob, 
#                    "prior_proba" = res.EM$prior_proba[index_id,], 
#                    "posterior_proba" = res.EM$posterior_proba[index_id,], 
#                    "group" = matrix(unlist(lapply(res.EM$group, 
#                                                 function(x){as.numeric(index_id %in% x)})), 
#                                   ncol = G, byrow = FALSE)
#   )
#   
#   Y <- res.EM$Y[index_id, numY]
#   
#   switch(window, 
#          T = dev.new(), 
#          eps = postscript(file = paste(path, filename, ".eps", sep = ""), width = width * scale / 90, height = height * scale / 90, horizontal = FALSE, onefile = FALSE, paper = "special"), 
#          svg = svg(filename = paste(path, filename, ".svg", sep = ""), width = width * scale / 90, height = height * scale / 90, onefile = FALSE), 
#          png = png(filename = paste(path, filename, ".png", sep = ""), width = width * scale, height = height * scale, res = res)
#   )
#   
#   if(is.null(col))
#   {col <- c("red", "blue", "green", "yellow", "purple", "brown")[1:G]}
#   if(is.null(nclass))
#   {nclass <- min(round(length(index_id) / 50), 100)}
#   
#   if(is.null(main)){
#     main <- switch(type, 
#                    "hist" = "histogram", 
#                    "density" = "density")
#   }
#   if(is.null(identifier) == FALSE)
#   {main <- paste(main, " Id = ", identifier, sep = "")}
#   
#   if(prob == "group" && type == "hist"){
#     
#     hist1 <- hist(Y, nclass = nclass, 
#                   main = main, xlab = colnames(res.EM$Y)[numY], 
#                   cex.main = cex.main)
#     
#     if(!is.null(x.legend)){
#       legend(x.legend, legend = paste("g", 1:G), col = col, pch = 15, bty = "n")
#     }
#     
#     for(iter_g in 1:G){
#       hist(Y[proba[,iter_g] == 1], 
#            col = col[iter_g], add = TRUE, breaks = hist1$breaks)
#     }
#     
#   }else{
#     
#     hist1 <- hist(Y, nclass = nclass, freq = FALSE, 
#                   main = main, xlab = colnames(res.EM$Y)[numY], 
#                   cex.main = cex.main)
#     
#     for(iter_g in 1:G){
#       seq_x <- seq(min(Y), max(Y), length.out = 1000)
#       mean_g <- weighted.mean(Y, w = proba[,iter_g])
#       sd_g <- EM.weighted.sd(Y, w = proba[,iter_g])
#       
#       if(res.EM$family[[numY]][[iter_g]]$family == "gaussian"){
#         points(seq_x, 
#                sum(proba[,iter_g]) / sum(proba) * dnorm(seq_x, mean = mean_g, sd = sd_g), 
#                type = "l", col = paste(col[iter_g], "4", sep = ""), lwd = lwd)        
#       }
#       
#       if(res.EM$family[[numY]][[iter_g]]$family == "Gamma"){
#         scale_g <- sd_g^2 / mean_g
#         shape_g <- mean_g^2 / sd_g^2
#         points(seq_x, 
#                sum(proba[,iter_g]) / sum(proba) * dgamma(seq_x, scale = scale_g, shape = shape_g), 
#                type = "l", col = paste(col[iter_g], "4", sep = ""))
#       }
#       
#       if(res.EM$family[[numY]][[iter_g]]$family == "Beta"){
#         shape1_g <- mean_g * (mean_g * (1 - mean_g) / sd_g^2 - 1)     
#         shape2_g <- (1 - mean_g) * (mean_g * (1 - mean_g) / sd_g^2 - 1)
#         points(seq_x, 
#                sum(proba[,iter_g]) / sum(proba) * dbeta(seq_x, shape1=shape1_g, shape2=shape2_g), 
#                type = "l", col = paste(col[iter_g], "4", sep = ""))
#       }
#       
#       abline(v = weighted.mean(Y, w = proba[,iter_g]), 
#              col = paste(col[iter_g], "4", sep = ""), lwd = lwd.mean, lty = lty.wmean)    
#       
#       if(!is.null(x.legend)){
#         legend(x.legend, legend = paste("g", 1:G), col = col, pch = 15, bty = "n")
#       }
#       
#     }
#   }
#   
#   switch(window, 
#          eps = dev.off(), 
#          svg = dev.off(), 
#          png = dev.off()
#   )
# }
