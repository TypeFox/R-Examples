 ##########################################################
 ## Example 1: Dataset containing nominal values on all attributes
 ##########################################################
 library(RoughSets)
 data(RoughSetData)
 decision.table <- RoughSetData$hiring.dt 

 ########## using fuzzy lower approximation ##############
 control <- list(t.implicator = "lukasiewicz", type.relation = c("tolerance", "eq.1"), type.aggregation = c("t.tnorm", "lukasiewicz"))
 reduct.1 <- FS.quickreduct.FRST(decision.table, type.method = "fuzzy.dependency", 
                             type.QR = "modified.QR", control = control)

 ########## using fuzzy boundary region ##############
 control <- list(t.implicator = "lukasiewicz", type.relation = c("tolerance", "eq.1"), type.aggregation = c("t.tnorm", "lukasiewicz"))
 reduct.2 <- FS.quickreduct.FRST(decision.table, type.method = "fuzzy.boundary.reg", 
                             type.QR = "modified.QR", control = control)
 
 ########## using vaquely quantified rough sets (VQRS) #########
 control <- list(alpha = 0.9, q.some = c(0.1, 0.6), q.most = c(0.2, 1), type.aggregation = c("t.tnorm", "lukasiewicz")) 
 reduct.3 <- FS.quickreduct.FRST(decision.table, type.method = "vqrs", 
                             type.QR = "modified.QR", control = control)

 ########## ordered weighted average (OWA) #########
 control <- list(t.implicator = "lukasiewicz", type.relation = c("tolerance", "eq.1"), m.owa = 3, type.aggregation = c("t.tnorm","lukasiewicz")) 
 reduct.4 <- FS.quickreduct.FRST(decision.table, type.method = "owa", 
                             type.QR = "modified.QR", control = control)

 ########## robust fuzzy rough sets (RFRS) #########
 control <- list(t.implicator = "lukasiewicz", type.relation = c("tolerance", "eq.1"), type.rfrs = "k.trimmed.min", 
                type.aggregation = c("t.tnorm", "lukasiewicz"), k.rfrs = 0) 
 reduct.5 <- FS.quickreduct.FRST(decision.table, type.method = "rfrs", 
                             type.QR = "modified.QR", control = control)

 ########## using min positive region (delta) ###########
 control <- list(alpha = 1, t.implicator = "lukasiewicz", type.relation = c("tolerance", "eq.1"), type.aggregation = c("t.tnorm", "lukasiewicz"))
 reduct.6 <- FS.quickreduct.FRST(decision.table, type.method = "min.positive.reg", 
                             type.QR = "modified.QR", control = control)

 ########## using FVPRS approximation ##############
 control <- list(alpha.precision = 0.05, t.implicator = "lukasiewicz", type.aggregation = c("t.tnorm", "lukasiewicz"), 
                type.relation = c("tolerance", "eq.1"))
 reduct.7 <- FS.quickreduct.FRST(decision.table, type.method = "fvprs", 
                             type.QR = "modified.QR", control = control)
 
 ########## using beta.PFRS approximation ##############
 control <- list(t.implicator = "lukasiewicz", type.relation = c("tolerance", "eq.1"), beta.quasi = 0.05, 
                 type.aggregation = c("t.tnorm", "lukasiewicz"))
 reduct.8 <- FS.quickreduct.FRST(decision.table, type.method = "beta.pfrs", 
                             type.QR = "modified.QR", control = control)

 ########## using fuzzy discernibility matrix ##############
 control <- list(alpha = 1)
 reduct.9 <- FS.quickreduct.FRST(decision.table, type.method = "fuzzy.discernibility", 
                             type.QR = "modified.QR", control = control)
							 
							 

 