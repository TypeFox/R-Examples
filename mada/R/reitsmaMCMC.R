# reitsmaMCMC <- function(data, ...) UseMethod("reitsmaMCMC")
# 
# reitsmaMCMC.default <-
#   function(data = NULL, subset=NULL, formula = NULL,
#            TP="TP", FN="FN", FP="FP", TN="TN", 
#            alphasens = 1, alphafpr = 1, 
#            method = "MCMCglmm",  
#            n.samples = NULL, burnin = NULL, thin = NULL, verbose = FALSE,
#            ...){
#     call <- match.call()
#     mcall <- match.call(expand.dots = FALSE)
#     
#     if(is.null(burnin)){burnin <- 10000}
#     if(is.null(thin)){thin <- 100}
#     if(is.null(n.samples)){n.samples <- 1000}
# 
#     if(method == "MCMCglmm"){nitt <- burnin + thin*(n.samples)}
#     
#     stopifnot(0 <= alphasens, alphasens <= 2, 0 <= alphafpr, alphafpr <= 2,
#               method %in% c("MCMCglmm"),
#               is.numeric(TP) | (is.character(TP) & length(TP) == 1),
#               is.numeric(FP) | (is.character(FP) & length(FP) == 1),
#               is.numeric(TN) | (is.character(TN) & length(TN) == 1),
#               is.numeric(FN) | (is.character(FN) & length(FN) == 1))
#     if(!is.null(subset)){stopifnot(length(subset)>0, 
#                                    min(subset)>0, 
#                                    max(subset) < nrow(data)+1)}
#     if(is.null(subset)){subset <- 1:nrow(data)}
#         
#     if(!is.null(data) & is.character(c(TP,FP,TN,FN))){
#       X <- as.data.frame(data)
#       origdata <- data
#       TP <- getElement(X,TP)
#       FN <- getElement(X,FN)
#       FP <- getElement(X,FP)
#       TN <- getElement(X,TN)
#     }
#     
#     if(is.null(data) & !is.character(c(TP,FP,TN,FN))){
#       origdata <- data.frame(TP = TP, FN = FN, FP = FP, TN = TN)
#     }
#     
#     freqdata <- cbind(TP,FN,FP,TN)[subset,]
#     colnames(freqdata) <- c("TP", "FN", "FP", "TN")
#     checkdata(freqdata)
#     
#     N <- length(TP)  
#     
#   if(method == "MCMCglmm"){
#     if(is.null(formula)){formula <- cbind(TP,FN,FP,TN) ~ trait - 1}else{
#       cat("Make sure you understand the use of the reserved trait variable before specifying formulae.\n")
#     }
#     if(!class(formula) == "formula"){stop("formula must be of class formula")}
#     if(!formula[2] == (cbind(TP,FN,FP,TN)~1)[2]){
#       stop("The left hand side of formula must be cbind(TP,FN,FP,TN)")}
#     varnames <- all.vars(formula)
#     varnames <- varnames[! varnames %in% c("TP","FN","FP","TN", "trait")]  
#     if(!all(varnames %in% colnames(data))){
#       stop("Some variables in formula are not in data; supply them in a data frame")}
#     
#     fit <- MCMCglmm(formula, rcov =  ~ us(trait):mada_study_id, 
#                     data = cbind(data[subset, ], data.frame(freqdata),
#                                  data.frame(mada_study_id = 1:N)), 
#                     family = c("multinomial2","multinomial2"),
#                     burnin = burnin, thin = thin, nitt= nitt,
#                     verbose = verbose, ...)
#     } # end of method == "MCMCglmm"
#     class(fit) <- "reitsmaMCMC"
#     fit
#   }
# 
# print.reitsmaMCMC <- function(x, ...){print(str(x))}
# 
# summary.reitsmaMCMC <-  function(object, ...){summary(object, ...)}
# 
# print.summary.reitsmaMCMC <- function(x, ...){print(x, ...)}
