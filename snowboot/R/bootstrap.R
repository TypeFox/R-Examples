# # This version allows growing patches without re-selecting the seeds for higher waves.
#
# # rm(list=ls()) SNOW BALL SAMPLING ######################################
#
#
# # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
# LSMI <- function(net, n.seeds = 10, n.neigh = 1, seeds = NULL) {
#   # this function returns the vertices samples by snowball up to wave n.neigh.
#   # net is object network (what is important is the component $edges and the length of $degree)
#   # this function randomly sample n.seeds and then select the neighbours up
#   # to wave n.neigh also give the index of those nodes last added
#   # and for which in the next stages of the sampling we assume we do not
#   # have their complete degree information. seed0 are the original seeds.
#   # sampleN are the possibly repeated elements in the sample.
#   # unodes are the no repeated elements in the sample.
#   # nodes.waves are the vertices added in each wave.
#   # Vertices may be present in more that one wave and more that once in a single wave.
#   # last.added are the vertices that are the most recently added into the set.
#
#
#   unodes <- nodes.waves <- as.list(rep(0, n.seeds))
#   # seeds selection: is without replacement and at random
#   if (is.null(seeds)) {
#     seed0 <- sort(sample(1:length(net$degree), n.seeds, rep = FALSE))
#   } else {
#     seed0 <- seeds
#   }
#   sampleN <- NULL
#   for (i in 1:n.seeds) {
#     res <- sample_about_one_seed(net, seed0[i], n.neigh)
#     sampleN <- c(sampleN, res$sampleN)
#     unodes[[i]] <- res$unodes
#     nodes.waves[[i]] <- res$nodes.waves
#   }
#   list(seeds = seed0, sampleN = sort(sampleN), unodes = unodes, nodes.waves = nodes.waves)
# }
# # Examples:
# # we are not really interested in running this function directly,
# # but within the next function called empdegree
# # distrib6 net<-local.network.MR.new5(n=100,distrib='pois',param=2)
# # a<-LSMI(net,n.seeds=3,n.neigh=3,seeds=NULL)
# # a<-sampleneigh(net,n.seeds=3,n.neigh=3,seeds=NULL)
#
# # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
#
# sample_about_one_seed <- function(net, seed0, n.neigh = 1) {
#   # this function returns the vertices samples by snowball up to wave n.neigh around a single seeds net is object network
#   # (what is important is the component $edges and the length of $degree) this function randomly sample n.seeds and then
#   # select the neighbours up to wave n.neigh seed0 the id of the seeds sampleN are the possibly repeated elements in the
#   # sample unodes are the no repeated elements in the sample nodes.waves are the vertices added in each wave. Vertices may
#   # be present in more that one wave and more that once in a single wave. last.added are the vertices that are the most
#   # recently added into the set.
#
#   sampleN <- nodes <- seed0
#   nodes.waves <- as.list(rep(0, n.neigh))
#   effEdges <- net$edges
#   more <- TRUE
#   nn <- n.neigh
#   new.nodes <- 0
#   # if(n.neigh==0) we only keep the seeds
#
#   wave <- 1
#   while (wave <= n.neigh & more) {
#     a <- is.element(effEdges, nodes)  #'nodes' will be accumulating all included vertices (non repeated)
#     if (any(a))
#       {
#         eedges <- which(matrix(a, dim(effEdges)[1], 2), arr.ind = TRUE)  #now it is the row number and column where they are in the edges matrix
#         nodes.waves[[wave]] <- arr.nodes <- sort(effEdges[cbind(eedges[, 1], sapply(eedges[, 2], FUN = switch, 2,
#           1))])  #the vertices we arrived to (duplicity is allowed)
#         # I need this specially to know which vertices were the last added:
#         if (!anyDuplicated(eedges[, 1])) {
#           new.nodes <- arr.nodes  #all the vertices we arrive to weren't already included in 'nodes'
#         } else {
#           new.nodes <- setdiff(arr.nodes, nodes)
#         }  #Then, already included vertices are not considered new because of inclusion of edge connecting them
#         ### subEdges<-effEdges[unique(a),] #the subset of edges. The repeated just have to be included once. maybe we are arriving
#         ### to the nodes more than once (due to small cycles) or we can get again to already included vertices (due to larger
#         ### cycles).  We want to include them as may times as they are neighbours of already included vertices. That is why I
#         ### consider arr.nodes.if a originally seeds vertex is included more than once, it is because it was selected also by
#         ### following one edge and then it also has the category of non seeds.
#
#         sampleN <- sort(c(sampleN, arr.nodes))
#         nodes <- unique(sampleN)
#         if (nn > 1)
#           effEdges <- effEdges[-unique(eedges[, 1]), ]  #I remove the 'used edges' to facilitate following searches within while,and assure
#         # we do not 'arrive' to a node more times than edges it has.
#         if (length(effEdges) > 0) {
#           if (is.vector(effEdges))
#           effEdges <- t(effEdges)  #I have to this very often in R to make sure effEdges is a matrix and not a vector
#         } else {
#           more <- FALSE
#         }  #when it reduces to become a matrix with one row.
#       }  #end if(any(a))
#     wave <- wave + 1
#   }  #end while
#   # browser() list(seeds=seed0,sampleN=sampleN,unodes=nodes,nodes.waves=nodes.waves,last.added=sort(new.nodes))
#   list(seeds = seed0, sampleN = sampleN, unodes = nodes, nodes.waves = nodes.waves)
# }
#
# # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
#
# ###### Simulataneously sampling ......interfererence among patches
# sampleneigh <- function(net, n.seeds = 10, n.neigh = 1, seeds = NULL) {
#   # this function returns the vertices samples by snowball up to wave n.neigh. net is object network (what is important is
#   # the component $edges and the length of $degree) this function randomly sample n.seeds and then select the neighbours up
#   # to wave n.neigh also give the index of those nodes last added and for which in the next stages of the sampling we
#   # assume we do not have their complete degree information. seed0 are the original seeds sampleN are the possibly repeated
#   # elements in the sample unodes are the no repeated elements in the sample nodes.waves are the vertices added in each
#   # wave. Vertices may be present in more that one wave and more that once in a single wave. last.added are the vertices
#   # that are the most recently added into the set.
#
#
#   # seeds selection: is without replacement and at random
#   if (is.null(seeds)) {
#     seed0 <- sort(sample(1:length(net$degree), n.seeds, rep = FALSE))
#   }
#   sampleN <- nodes <- seed0
#   nodes.waves <- as.list(rep(0, n.neigh))
#   effEdges <- net$edges
#   more <- TRUE
#   nn <- n.neigh
#   new.nodes <- 0
#   # if(n.neigh==0) we only keep the seeds
#   wave <- 1
#   while (wave <= n.neigh & more) {
#     a <- is.element(effEdges, nodes)  #'nodes' will be accumulating all included vertices (non repeated)
#     if (any(a))
#       {
#         eedges <- which(matrix(a, dim(effEdges)[1], 2), arr.ind = TRUE)  #now it is the row number and column where they are in the edges matrix
#         nodes.waves[[wave]] <- arr.nodes <- sort(effEdges[cbind(eedges[, 1], sapply(eedges[, 2], FUN = switch, 2,
#           1))])  #the vertices we arrived to (duplicity is allowed)
#         # I need this specially to know which vertices were the last added:
#         if (!anyDuplicated(eedges[, 1])) {
#           new.nodes <- arr.nodes  #all the vertices we arrive to weren't already included in 'nodes'
#         } else {
#           new.nodes <- setdiff(arr.nodes, nodes)
#         }  #Then, already included vertices are not considered new because of inclusion of edge connecting them
#         ### subEdges<-effEdges[unique(a),] #the subset of edges. The repeated just have to be included once. maybe we are arriving
#         ### to the nodes more than once (due to small cycles) or we can get again to already included vertices (due to larger
#         ### cycles).  We want to include them as may times as they are neighbours of already included vertices. That is why I
#         ### consider arr.nodes.if a originally seeds vertex is included more than once, it is because it was selected also by
#         ### following one edge and then it also has the category of non seeds.
#
#         sampleN <- sort(c(sampleN, arr.nodes))
#         nodes <- unique(sampleN)
#         if (nn > 1)
#           effEdges <- effEdges[-unique(eedges[, 1]), ]  #I remove the 'used edges' to facilitate following searches within while,and assure
#         # we do not 'arrive' to a node more times than edges it has.
#         if (length(effEdges) > 0) {
#           if (is.vector(effEdges))
#           effEdges <- t(effEdges)  #I have to this very often in R to make sure effEdges is a matrix and not a vector
#         } else {
#           more <- FALSE
#         }  #when it reduces to become a matrix with one row.
#       }  #end if(any(a))
#     wave <- wave + 1
#   }  #end while
#   # browser()
#   list(seeds = seed0, sampleN = sampleN, unodes = nodes, nodes.waves = nodes.waves, last.added = sort(new.nodes))
# }
# # Examples #we are not really interested in running this function directly but within the next function called empdegree
# # distrib6 net<-local.network.MR.new5(n=100,distrib='pois',param=2) a<-sampleneigh(net,n.seeds=3,n.neigh=1,seeds=NULL)
# # a<-sampleneigh(net,n.seeds=3,n.neigh=3,seeds=NULL)
#
# # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
#
# ###################################################################### SAMPLES ###########################
#
# Oempdegreedistrib <- function(net, n.seeds, n.neigh, num.sam, idname = "Temp", seeds) {
#   if (n.neigh == 0) {
#     # only information from the seeds
#     res <- Oempdegreedistrib0(net, n.seeds, n.neigh, num.sam, idname = "Temp", seeds)
#   } else {
#     res <- OempdegreedistribK(net, n.seeds, n.neigh, num.sam, idname = "Temp", seeds)
#   }
#   res
# }
#
# # ----------------------------------------------------------------------#
#
# Oempdegreedistrib0 <- function(net, n.seeds, n.neigh, num.sam, idname = "Temp", seeds) {
#   p0.seeds.array <- Oempd <- values.array <- val.seeds.array <- samples <- as.list(rep(NA, num.sam))
#   seeds1 <- matrix(NA, num.sam, n.seeds)
#   ## -------the 'real' parameters in the network:-------##
#   real <- summary.net(net)
#   realdd <- real$realdd
#   ## ---------------------------------------------------##
#   for (m in 1:num.sam) {
#     # if(m%%100==1)#cat('Obtaining empd of sample ',m,'\n')
#     neigh.seeds <- sort(sample(1:length(net$degree), n.seeds, rep = FALSE))  #n.neigh=0!!!!!!!
#     tab.seeds <- table(neigh.seeds)  #id seeds
#     seeds1[m, ] <- neigh.seeds
#     ###### degrees #####
#     deg.seeds <- realdd[rep(as.integer(names(tab.seeds)), tab.seeds)]  #duplicities allowed
#     # --------------------#
#     samples[[m]] <- list(freq.deg.seeds = freq.deg.seeds <- table(deg.seeds))
#     ##### resample and extract the degree of selected vertices######
#     values.array[[m]] <- values <- val.seeds.array[[m]] <- val.seeds <- sort(unique(deg.seeds))
#     p0.seeds.array[[m]] <- sum(deg.seeds == 0)/n.seeds
#     #### Frequency ##### (Not the relative frequency)
#     OFseed <- table.row(deg.seeds, values)
#     #### Empirical degree distribution ########
#     Oempd.seeds <- OFseed/n.seeds
#     Oempd[[m]] <- list(Oempd = Oempd.seeds)
#   }  # for(m in 1:num.sam)
#   # browser()
#   list(idname = idname, samples = samples, values = values.array, Oempd = Oempd, num.sam = num.sam, val.seeds = val.seeds.array,
#     n.seeds = n.seeds, n.neigh = n.neigh, p0.seeds = p0.seeds.array, seeds1 = seeds1)
# }
#
# # ----------------------------------------------------------------------#
#
# OempdegreedistribK <- function(net, n.seeds, n.neigh, num.sam, idname = "Temp", seeds) {
#   # This function obtains the empirical degree distribution from num.sam samples net is the network (only one) n.seeds is
#   # the number of seeds to set the neighbourhood sample n.neigh is the neighbouhood size around each seeds num.sam is the
#   # number of different samples taken from the same network idname is to identify from which nets we are sampling and
#   # resampling.
#   seeds1 <- matrix(0, num.sam, n.seeds)
#   p0.seeds.array <- Oempd <- ekseed.array <- values.array <- val.seeds.array <- val.nonseed.array <- samples <- as.list(rep(NA,
#     num.sam))
#
#   ## -------the 'real' parameters in the network:-------##
#   real <- summary.net(net)
#   realdd <- real$realdd
#   # rmeand<-real$mean(realdd) rquart<-real$rquart
#   rfreq <- real$rfreq
#   # rperc<-real$rperc -------------------------------------##
#
#   for (m in 1:num.sam) {
#     # if(m%%100==1)#cat('Obtaining empd of sample ',m,'\n') browser()
#     neigh <- LSMI(net, n.seeds = n.seeds, n.neigh = n.neigh, seeds = seeds[m, ])
#     seeds1[m, ] <- neigh$seeds
#     # nodes<-neigh$sampleN[!is.element(neigh$sampleN,neigh$last.added)] #vertices that are not included last (with their
#     # duplicities)
#     nodes <- neigh$sampleN  #vertices up to distance n.neigh(with their duplicities)!!!!!!!!!!!!!!
#     tab.nodes <- table(nodes)  #now it has the info of seeds and non-seeds
#     # Now we want to distinguish between seeds and non seeds. Remember that some seeds can also be non seeds if they were also
#     # selected by following one edge.
#     tab.seeds <- table(neigh$seeds)  #id seeds
#     a <- is.element(names(tab.nodes), names(tab.seeds))
#     if (all(names(tab.nodes[a]) == names(tab.seeds))) {
#       tab.nodes[a] <- tab.nodes[a] - tab.seeds
#     } else {
#       cat("no mismo orden")
#       browser()
#     }
#     tab.nodes <- tab.nodes[tab.nodes > 0]  #now I have only the id of non-seeds
#     ###################### degrees #####
#     deg.seeds <- realdd[rep(as.integer(names(tab.seeds)), tab.seeds)]  #in case seeds are present more than once as seeds
#     deg.nonseed <- realdd[rep(as.integer(names(tab.nodes)), tab.nodes)]  #to incorporate their duplicity
#     deg.nonseedU <- realdd[as.integer(names(tab.nodes))]  #it contains the non seeds only one time (I do not use in the rest of the code)
#     # --------------------#
#     samples[[m]] <- list(freq.deg.seeds = freq.deg.seeds <- table(deg.seeds), freq.deg.nonseed = freq.deg.nonseed <- table(deg.nonseed),
#       freq.deg.nonseedU = freq.deg.nonseedU <- table(deg.nonseedU))
#     ##### resample and extract the degree of selected vertices######
#     val.seeds.array[[m]] <- val.seeds <- sort(unique(deg.seeds))  #as.numeric(names(freq.deg.seeds))
#     val.nonseed.array[[m]] <- val.nonseed <- sort(unique(deg.nonseed))  #as.numeric(names(freq.deg.nonseed))
#     val.nonseedU <- sort(unique(deg.nonseedU))  #as.numeric(names(freq.deg.nonseedU))
#
#     p0.real <- rfreq[1]
#     p0.seeds <- 0
#     if (any(val.seeds == 0)) {
#       # if any seeds has degree zero
#       p0.seeds <- sum(deg.seeds == 0)/n.seeds
#     }
#     p0.seeds.array[[m]] <- p0.seeds
#     values <- sort(union(val.seeds, val.nonseed))  #all the possible degree values toresample
#     values.array[[m]] <- values
#
#     ###################### Frequency ##### (Not the relative frequency)
#     OFseed <- table.row(deg.seeds, values)
#     OFnonseed <- table.row(deg.nonseed, values)
#     #################################################################### combining information from seeds and nonseeds ###### mean degree computed from the original sampled seeds:
#     ekseed.array[[m]] <- ekseed <- sum(as.numeric(names(freq.deg.seeds)) * freq.deg.seeds)/n.seeds
#     colzero <- NULL
#     if (any(values == 0)) {
#       colzero <- which(values == 0)
#       vals <- values[-colzero]
#       Of.seeds <- OFseed[-colzero]
#       Of.nonseed.nw <- OFnonseed[-colzero]
#     } else {
#       vals <- values
#       Of.seeds <- OFseed
#       Of.nonseed.nw <- OFnonseed
#     }
#     ################################################ NWB # seeds and non weighted nonseeds #### p0 estimated from orginal sampled seeds#
#     Oempd.nw.p0sEks <- (Of.seeds + (1 - p0.seeds) * ekseed * Of.nonseed.nw/vals)/(n.seeds + ekseed * sum(Of.nonseed.nw/vals))
#     ######################################### p0 taken as known #
#     Oempd.nw.p0rEks <- (Of.seeds + (1 - p0.real) * ekseed * Of.nonseed.nw/vals)/(n.seeds + ekseed * sum(Of.nonseed.nw/vals))
#     if (any(values == 0)) {
#       Oempd.nw.p0sEks <- c(p0.seeds, Oempd.nw.p0sEks)  #5
#       Oempd.nw.p0rEks <- c(p0.real, Oempd.nw.p0rEks)  #8
#     }
#     Oempd[[m]] <- list(Oempd = Oempd.nw.p0sEks, Oempd.nw.p0rEks = Oempd.nw.p0rEks)
#   }  # for(m in 1:num.sam)
#   # browser()
#   list(idname = idname, samples = samples, values = values.array, Oempd = Oempd, num.sam = num.sam, val.seeds = val.seeds.array,
#     val.nonseed = val.nonseed.array, n.seeds = n.seeds, n.neigh = n.neigh, p0.real = p0.real, p0.seeds = p0.seeds.array,
#     ekseed = ekseed.array, seeds1 = seeds1)
# }
#
# ###################################################################### BOOTSTRAP SAMPLES ######################################
#
# bootdeg <- function(sam.out, num.sam, n.boot, idname = "Temp") {
#   if (sam.out$n.neigh == 0) {
#     # only information from the seeds
#     res <- bootdeg0(sam.out, num.sam, n.boot, idname = "Temp")
#   } else {
#     res <- bootdegK(sam.out, num.sam, n.boot, idname = "Temp")
#   }
#   res
# }
# # ----------------------------------------------------------------------#
# bootdeg0 <- function(sam.out, num.sam, n.boot, idname = "Temp") {
#   # This function obtains the bootstrap samples for each sample from a network sam.out is the output of Oempdegreedistrib
#   # num.sam is the number of different samples taken from the same network (Scalar or vector) n.boot is the number of
#   # bootstrap samples taken from each sample
#   n.seeds <- sam.out$n.seeds
#   n.neigh <- sam.out$n.neigh
#   if (length(num.sam) == 1)
#     num.sam <- 1:num.sam
#   empd <- as.list(rep(NA, length(num.sam)))
#   i <- 1
#   for (m in num.sam) {
#     # if(i%%100==1)#cat('Processing bootstrap samples of sample=',i,'\n') #print every 100
#     i <- i + 1
#     val.seeds <- sam.out$val.seeds[[m]]
#     freq.deg.seeds <- sam.out$samples[[m]]$freq.deg.seeds
#     bsam.seeds <- myBsample(val.seeds, n.seeds, n.boot, prob = freq.deg.seeds)
#     values <- sam.out$values[[m]]  #all the possible degree values toresample
#     #### Frequency ##### (Not the relative frequency)
#     Fseed <- t(apply(bsam.seeds, 1, table.row, vect = values))  #freq (sorted according to values)
#     # browser()
#     empd.seeds <- Fseed/n.seeds
#     empd[[m]] <- list(empd.seeds = empd.seeds)
#   }  # for(m in num.sam)
#   list(idname = idname, values = sam.out$values, empd = empd, num.sam = num.sam, n.boot = n.boot, n.neigh = n.neigh)
# }
#
# # ---------------------------------------------------------------------------------------#
#
# bootdegK <- function(sam.out, num.sam, n.boot, idname = "Temp") {
#   # This function obtains the bootstrap samples for each sample from a network sam.out is the output of Oempdegreedistrib
#   # num.sam is the number of different samples taken from the same network. Scalar o vector. n.boot is the number of
#   # bootstrap samples taken from each sample
#   n.seeds <- sam.out$n.seeds
#   n.neigh <- sam.out$n.neigh
#   if (length(num.sam) == 1)
#     num.sam <- 1:num.sam
#
#   empd <- as.list(rep(NA, length(num.sam)))
#   i <- 1
#   for (m in num.sam) {
#     # if(i%%100==1)#cat('Processing bootstrap samples of sample=',i,'\n') #print every 100
#     i <- i + 1
#     # Boostrap samples of seeds, nonseeds-noWeighted and nonseeds-Weighted:
#     val.seeds <- sam.out$val.seeds[[m]]
#     val.nonseed <- sam.out$val.nonseed[[m]]
#     freq.deg.seeds <- sam.out$samples[[m]]$freq.deg.seeds
#     freq.deg.nonseed <- sam.out$samples[[m]]$freq.deg.nonseed
#
#     bsam.seeds <- myBsample(val.seeds, n.seeds, n.boot, prob = freq.deg.seeds)  #matrix n.boot x n.seeds
#     bsam.nonseed.nw <- myBsample(val.nonseed, sum(freq.deg.nonseed), n.boot, prob = freq.deg.nonseed)  #matrix n.boot x sum(freq.deg.nonseed)
#     bsam.nonseed.w <- myBsample(val.nonseed, sum(freq.deg.nonseed), n.boot, prob = freq.deg.nonseed/val.nonseed)  #matrix
#
#     p0.B <- rep(0, n.boot)
#     if (any(val.seeds == 0)) {
#       # if any seeds has degree zero
#       p0.B <- rowSums(bsam.seeds == 0)/n.seeds  #the estimation from the bootstrap samples
#     }
#     p0.real <- sam.out$p0.real
#     p0.seeds <- sam.out$p0.seeds[[m]]
#
#     values <- sam.out$values[[m]]  #all the possible degree values to resample
#     # if (all(values==0)) {values<-0.1}
#
#     ###################### Frequency ##### (Not the relative frequency)
#     Fseed <- t(apply(bsam.seeds, 1, table.row, vect = values))  #frequency (sorted according to values)
#     # browser()
#     if (is.null(bsam.nonseed.nw)) {
#       # browser()
#       Fnonseed.nw <- 0
#     } else {
#       Fnonseed.nw <- t(apply(as.matrix(bsam.nonseed.nw), 1, table.row, vect = values))
#     }
#     if (is.null(bsam.nonseed.w)) {
#       Fnonseed.n <- 0
#     } else {
#       Fnonseed.w <- t(apply(as.matrix(bsam.nonseed.w), 1, table.row, vect = values))
#     }
#
#     # Fnonseed.nw<-t(apply(as.matrix(bsam.nonseed.nw),1,table.row,vect=values)) # '
#     # Fnonseed.w<-t(apply(as.matrix(bsam.nonseed.w),1,table.row,vect=values)) # '
#
#     #################################################################### combining information from seeds and nonseeds ######
#
#     # mean degree computed from the original sampled seeds:
#     ekseed <- sam.out$ekseed[[m]]
#
#     colzero <- NULL
#     if (any(values == 0)) {
#       colzero <- which(values == 0)
#       vals <- values[-colzero]
#       f.seeds <- Fseed[, -colzero]
#       f.nonseed.nw <- Fnonseed.nw[, -colzero]
#       f.nonseed.w <- Fnonseed.w[, -colzero]
#     } else {
#       vals <- values
#       f.seeds <- Fseed
#       f.nonseed.nw <- Fnonseed.nw
#       f.nonseed.w <- Fnonseed.w
#     }
#     # browser() WB # seeds and weighted nonseeds ###### consider the p0 fixed from the seeds information:
#     # empd.w.p0s<-(f.seeds+(1-p0.seeds)*f.nonseed.w)/(n.seeds+sum(freq.deg.nonseed))
#     empd.w.p0s <- (f.seeds + f.nonseed.w * (1 - p0.B))/(n.seeds + sum(freq.deg.nonseed))
#     #### consider the p0 fixed from the real information:
#     #### empd.w.p0r<-(f.seeds+(1-p0.real)*f.nonseed.w)/(n.seeds+sum(freq.deg.nonseed)) NWB # seeds and non weighted nonseeds ####
#     #### p0 estimated from orginal sampled seeds# E(K) estimated from bootstrap samples from the seeds
#     #### empd.nw.p0sEkb<-(f.seeds+(1-p0.seeds)*apply(bsam.seeds,1,FUN=mean)*t(t(f.nonseed.nw)/vals))/(n.seeds+
#     #### apply(bsam.seeds,1,FUN=mean)*rowSums(t(t(f.nonseed.nw)/vals)))
#     empd.nw.p0sEkb <- (f.seeds + t(t(f.nonseed.nw)/vals) * (1 - p0.B) * apply(bsam.seeds, 1, FUN = mean))/(n.seeds + rowSums(t(t(f.nonseed.nw)/vals)) *
#       apply(bsam.seeds, 1, FUN = mean))
#     # E(K) estimated from the original seeds sample
#     # empd.nw.p0sEks<-(f.seeds+(1-p0.seeds)*ekseed*t(t(f.nonseed.nw)/vals))/(n.seeds+ ekseed*rowSums(t(t(f.nonseed.nw)/vals)))
#     empd.nw.p0sEks <- (f.seeds + ekseed * t(t(f.nonseed.nw)/vals) * (1 - p0.B))/(n.seeds + ekseed * rowSums(t(t(f.nonseed.nw)/vals)))
#     ######################################### p0 taken as known # E(K) estimated from bootstrap samples from the seeds
#     ######################################### empd.nw.p0rEkb<-(f.seeds+(1-p0.real)*apply(bsam.seeds,1,FUN=mean)*t(t(f.nonseed.nw)/vals))/
#     ######################################### (n.seeds+apply(bsam.seeds,1,FUN=mean)*rowSums(t(t(f.nonseed.nw)/vals))) E(K) estimated from the original seeds sample
#     ######################################### empd.nw.p0rEks<-(f.seeds+(1-p0.real)*ekseed*t(t(f.nonseed.nw)/vals))/(n.seeds+ ekseed*rowSums(t(t(f.nonseed.nw)/vals)))
#
#     if (any(values == 0)) {
#       empd.w.p0s <- cbind(`0` = p0.B, empd.w.p0s)  #1
#       # empd.w.p0r<-cbind('0'=p0.real,empd.w.p0r) #2
#       empd.nw.p0sEkb <- cbind(`0` = p0.B, empd.nw.p0sEkb)  #3
#       empd.nw.p0sEks <- cbind(`0` = p0.B, empd.nw.p0sEks)  #4
#       # empd.nw.p0rEkb<-cbind('0'=p0.real,empd.nw.p0rEkb) #6 empd.nw.p0rEks<-cbind('0'=p0.real,empd.nw.p0rEks) #7
#     }
#     empd[[m]] <- list(empd.w.p0s = empd.w.p0s, empd.nw.p0sEkb = empd.nw.p0sEkb, empd.nw.p0sEks = empd.nw.p0sEks)
#   }  # for(m in num.sam)
#   list(idname = idname, values = sam.out$values, empd = empd, num.sam = num.sam, n.boot = n.boot, n.neigh = n.neigh)
# }
#
#
# ################################################################## Parameters Estimation from orignal samples #########
#
# OparametersEst <- function(outOempd) {
#   # This function obtains the parameters based on the empirical distritributions derived from the samples outOempd is
#   # output of Oempdegreedistrib
#   num.sam <- outOempd$num.sam
#   mean <- rep(NA, num.sam)
#   quartiles <- array(NA, c(num.sam, 3))  #1rst,2nd and 3rd quartile
#   rfreq <- array(NA, c(num.sam, 5))  #freq of 0,1,2,3,4
#   deciles <- array(NA, c(num.sam, 9))  #10,20,30,40,50,60,70,80,90
#   for (m in 1:outOempd$num.sam) {
#     distrib <- outOempd$Oempd[[m]]$Oempd
#     # vals<-as.numeric(names(distrib)) #checar diversos formatos
#     vals <- outOempd$values[[m]]
#     mean[m] <- sum(vals * distrib)
#     quartiles[m, ] <- vals[sapply(X = c(0.25, 0.5, 0.75), FUN = min.greater, v = cumsum(distrib), ge = TRUE)]
#     rfreq[m, ] <- c(distribvals(distrib, vals, 0), distribvals(distrib, vals, 1), distribvals(distrib, vals, 2), distribvals(distrib,
#       vals, 3), distribvals(distrib, vals, 4))
#     deciles[m, ] <- vals[sapply(X = seq(0.1, 0.9, by = 0.1), FUN = min.greater, v = cumsum(distrib), ge = TRUE)]
#   }
#   # browser()
#   list(mean = mean, quartiles = quartiles, rfreq = rfreq, deciles = deciles, m = m)
# }
# # ---------------------------------------------------------------------------------------#
#
# distribvals <- function(distrib, vals, x) {
#   if (any(vals == x)) {
#     res <- distrib[which(vals == x)]
#   } else {
#     res <- 0
#   }
#   res
# }
# # ---------------------------------------------------------------------------------------#
#
# table.row <- function(vecdat, vect) {
#   # vecdat es el vector que se quiere arreglar segun los valores de vect
#   table(c(vecdat, vect)) - 1
# }
# # example: m<-matrix(1:9,3,3) t(apply(m,1,table.row,1:9)) 1 2 3 4 5 6 7 8 9 [1,] 1 0 0 1 0 0 1 0 0 [2,] 0 1 0 0 1 0 0 1 0
# # [3,] 0 0 1 0 0 1 0 0 1
#
# # ---------------------------------------------------------------------------------------#
#
# min.greater <- function(v, x, ge = TRUE) {
#   # I construct this function to obtain percentiles from a empirical distribution.  This form is to allow obtaining it for
#   # a matrix of distributions. v is a vector x is a value
#   if (ge) {
#     res <- min(which(v >= x))
#   } else {
#     res <- min(which(v > x))
#   }
#   res
# }
#
# # ---------------------------------------------------------------------------------------#
#
# myBsample <- function(val, length, n.boot, prob = NULL) {
#   # val is a vector to resample the sample size is length*n.boot
#   if (length(val) == 0) {
#     res <- NULL
#   } else if (length(val) == 1) {
#     # problems with weighted sample if only one term
#     res <- matrix(rep(val, length * n.boot), n.boot, length, byrow = TRUE)
#   } else if (length(val) > 1) {
#     res <- matrix(sample(val, length * n.boot, replace = TRUE, prob = prob), n.boot, length, byrow = TRUE)
#   }
#   res
# }
#
#
# # ---------------------------------------------------------------------------------------#
#
# summary.net <- function(net) {
#   # this function obtains the real parameters in a network
#   realdd <- net$degree - net$degree.left
#   rmean <- mean(realdd)
#   rquart <- quantile(realdd, prob = c(0.25, 0.5, 0.7))
#   rfreq <- c(sum(realdd == 0), sum(realdd == 1), sum(realdd == 2), sum(realdd == 3), sum(realdd == 4))/length(net$degree)
#   rdeci <- quantile(realdd, prob = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))
#   list(realdd = realdd, rmean = rmean, rquart = rquart, rfreq = rfreq, rdeci = rdeci)
# }
#
# # ---------------------------------------------------------------------------------------#
#
# Bias <- function(opar.out, rpar.out) {
#   # opar.out the 'estimated parameters'. List with the names, mean, quetiles, rfreq, deciles rpar.out the 'real
#   # parameters'. Listi the the names rmean, rquart, rfreq, rdeci
#   meanBias <- mean(opar.out$mean - rpar.out$rmean)
#   quarBias <- rowMeans(t(opar.out$quartiles) - rpar.out$rquart)
#   rfreqBias <- rowMeans(t(opar.out$rfreq) - rpar.out$rfreq)
#   deciBias <- rowMeans(t(opar.out$deciles) - rpar.out$rdeci)
#   meanSMES <- sqrt(mean((opar.out$mean - rpar.out$rmean)^2))
#   quarSMES <- sqrt(rowMeans((t(opar.out$quartiles) - rpar.out$rquart)^2))
#   rfreqSMES <- sqrt(rowMeans((t(opar.out$rfreq) - rpar.out$rfreq)^2))
#   deciSMES <- sqrt(rowMeans((t(opar.out$deciles) - rpar.out$rdeci)^2))
#
#   rnam <- c("mean", "quart1", "quart2", "quart3", paste("rfreq", 0:4), paste("deci", 1:9))
#   BiasSMES <- data.frame(meanBias = c(meanBias, quarBias, rfreqBias, deciBias), sqrtSMES = c(meanSMES, quarSMES, rfreqSMES,
#     deciSMES))
#   rownames(BiasSMES) <- rnam
#   BiasSMES
# }
#
# # ---------------------------------------------------------------------------------------#
#
# Bias.SMSE.ij <- function(net, real.par, n.seeds, n.neigh, sam.size) {
#   All.biasSMSE <- 1
#   for (i in n.seeds) {
#     for (j in n.neigh) {
#       # cat('i= ',i,'\t','j= ',j,'\n')
#       Obs.distrib <- Oempdegreedistrib(net, n.seeds = i, n.neigh = j, num.sam = sam.size)
#       Oparam <- OparametersEst(Obs.distrib)
#       biasSMSE <- Bias(Oparam, real.par)
#       All.biasSMSE <- cbind(All.biasSMSE, rbind(biasSMSE, c(i, j)))
#     }
#   }
#   All.biasSMSE
# }
#
# # ---------------------------------------------------------------------------------------#
#
# BparametersEst <- function(outBempd) {
#   # outBempd is output of bootdeg
#   tn.sam <- length(outBempd$num.sam)
#   if (outBempd$n.neigh == 0) {
#     n.dist <- 1  #n.dist is the number of different emp distr.
#   } else {
#     n.dist <- 3
#   }
#
#   mean <- array(NA, c(tn.sam, outBempd$n.boot, n.dist))
#   quartiles <- array(NA, c(tn.sam, 3, outBempd$n.boot, n.dist))
#   rfreq <- array(NA, c(tn.sam, 5, outBempd$n.boot, n.dist))  #freq of 0,1,2,3,4
#   deciles <- array(NA, c(tn.sam, 9, outBempd$n.boot, n.dist))  #10,20,30,40,50,60,70,80,90
#
#   for (m in 1:tn.sam) {
#     w <- 1
#     in.while <- TRUE
#     while (in.while) {
#       if (outBempd$n.neigh == 0) {
#         empd <- outBempd$empd[[m]]$empd.seeds
#         in.while <- FALSE
#       } else {
#         if (w == 1) {
#           empd <- outBempd$empd[[m]]$empd.w.p0s
#         } else if (w == 2) {
#           empd <- outBempd$empd[[m]]$empd.nw.p0sEkb
#         } else if (w == 3) {
#           empd <- outBempd$empd[[m]]$empd.nw.p0sEks
#           in.while <- FALSE
#         }
#       }
#       vals <- outBempd$values[[m]]  #as.numeric(colnames(empd))  ###checar
#
#       if (dim(empd)[2] != length(as.vector(vals))) {
#         # browser()
#         mean[m, , w] <- t(empd) %*% as.vector(vals)
#       } else mean[m, , w] <- empd %*% as.vector(vals)
#
#       cempd <- t(apply(empd, 1, FUN = cumsum))
#       quartiles[m, , , w] <- t(sapply(X = c(0.25, 0.5, 0.75), FUN = cempdpercentile, cempd = cempd, vals = vals))
#       rfreq[m, , , w] <- t(sapply(X = 0:4, FUN = distribvalsmat, empd = empd, vals = vals))
#       deciles[m, , , w] <- t(sapply(X = seq(0.1, 0.9, by = 0.1), FUN = cempdpercentile, cempd = cempd, vals = vals))
#       w <- w + 1
#     }  #while(w<=3)
#   }  #for (m in 1:...)
#   # browser()
#   list(mean = mean, quartiles = quartiles, rfreq = rfreq, deciles = deciles, n.dist = n.dist, num.sam = outBempd$num.sam)
# }  #function
#
# # ---------------------------------------------------------------------------------------#
#
# distribvalsmat <- function(empd, vals, x) {
#   if (any(vals == x)) {
#     res <- empd[, which(vals == x)]
#   } else {
#     res <- rep(0, dim(empd)[1])
#   }
#   res
# }
#
# # empd 1 2 [1,] 0.5 0.5 [2,] 0.5 0.5 [3,] 0.5 0.5 [4,] 0.0 1.0 [5,] 0.5 0.5
# # sapply(X=0:4,FUN=distribvalsmat,empd=empd,vals=vals) [,1] [,2] [,3] [,4] [,5] [1,] 0 0.5 0.5 0 0 [2,] 0 0.5 0.5 0 0
# # [3,] 0 0.5 0.5 0 0 [4,] 0 0.0 1.0 0 0 [5,] 0 0.5 0.5 0 0
#
# # ---------------------------------------------------------------------------------------#
#
# cempdpercentile <- function(cempd, perc, vals) {
#   # cempd is a matrix perc is a scalar
#   res <- vals[apply(X = cempd, 1, FUN = min.greater, x = perc, ge = TRUE)]
#   res
# }
# # Ejemplo cempdpercentile(cempd,0.25,vals)
#
# # sapply(X=c(.25,.5,.75),FUN=cempdpercentile,cempd=cempd,vals=vals)
# # sapply(X=c(.25,.5,.75),FUN=cempdpercentile,cempd=cempd,vals=vals) [,1] [,2] [,3] [1,] 1 1 2 [2,] 1 1 2 [3,] 1 1 2 [4,]
# # 2 2 2 [5,] 1 1 2
#
# # ---------------------------------------------------------------------------------------#
#
# B.Bias <- function(bpar.out, opar.out) {
#   # bpar.out the bootstrap estimated parameters. List with the names mean, quartiles, rfreq, deciles. This elements are
#   # arrays opar.out the estimated parameters form the original samples. List with the names, mean, quetiles, rfreq, deciles
#
#   if (length(bpar.out$num.sam) == 1) {
#     sams <- 1:bpar.out$num.sam
#   } else {
#     sams <- bpar.out$num.sam
#   }  #the subset of samples
#   # apply(bpar.out$mean-opar.out$mean[num.sam],c(1,3),FUN=mean) #matrix of mean tn.sam X n.dist
#   meanBias <- apply(bpar.out$mean - opar.out$mean[sams], 3, FUN = mean)  #the mean error by emp distribution
#   if (bpar.out$n.dist > 1) {
#     q1Bias <- apply(bpar.out$quartiles[, 1, , ] - opar.out$quartiles[sams, 1], 3, FUN = mean)
#     q2Bias <- apply(bpar.out$quartiles[, 2, , ] - opar.out$quartiles[sams, 2], 3, FUN = mean)
#     q3Bias <- apply(bpar.out$quartiles[, 3, , ] - opar.out$quartiles[sams, 3], 3, FUN = mean)
#     ## quarBias<-apply(bpar.out$quartiles-opar.out$rquart[num.sam],c(2,4),FUN=mean)
#     quarBias <- cbind(q1Bias, q2Bias, q3Bias)
#
#     f0Bias <- apply(bpar.out$rfreq[, 1, , ] - opar.out$rfreq[sams, 1], 3, FUN = mean)
#     f1Bias <- apply(bpar.out$rfreq[, 2, , ] - opar.out$rfreq[sams, 2], 3, FUN = mean)
#     f2Bias <- apply(bpar.out$rfreq[, 3, , ] - opar.out$rfreq[sams, 3], 3, FUN = mean)
#     f3Bias <- apply(bpar.out$rfreq[, 4, , ] - opar.out$rfreq[sams, 4], 3, FUN = mean)
#     f4Bias <- apply(bpar.out$rfreq[, 5, , ] - opar.out$rfreq[sams, 5], 3, FUN = mean)
#     rfreqBias <- cbind(f0Bias, f1Bias, f2Bias, f3Bias, f4Bias)
#
#     d1Bias <- apply(bpar.out$deciles[, 1, , ] - opar.out$deciles[sams, 1], 3, FUN = mean)
#     d2Bias <- apply(bpar.out$deciles[, 2, , ] - opar.out$deciles[sams, 2], 3, FUN = mean)
#     d3Bias <- apply(bpar.out$deciles[, 3, , ] - opar.out$deciles[sams, 3], 3, FUN = mean)
#     d4Bias <- apply(bpar.out$deciles[, 4, , ] - opar.out$deciles[sams, 4], 3, FUN = mean)
#     d5Bias <- apply(bpar.out$deciles[, 5, , ] - opar.out$deciles[sams, 5], 3, FUN = mean)
#     d6Bias <- apply(bpar.out$deciles[, 6, , ] - opar.out$deciles[sams, 6], 3, FUN = mean)
#     d7Bias <- apply(bpar.out$deciles[, 7, , ] - opar.out$deciles[sams, 7], 3, FUN = mean)
#     d8Bias <- apply(bpar.out$deciles[, 8, , ] - opar.out$deciles[sams, 8], 3, FUN = mean)
#     d9Bias <- apply(bpar.out$deciles[, 9, , ] - opar.out$deciles[sams, 9], 3, FUN = mean)
#     deciBias <- cbind(d1Bias, d2Bias, d3Bias, d4Bias, d5Bias, d6Bias, d7Bias, d8Bias, d9Bias)
#
#     meanSMSE <- sqrt(apply((bpar.out$mean - opar.out$mean[sams])^2, 3, FUN = mean))
#
#     q1SMSE <- apply((bpar.out$quartiles[, 1, , ] - opar.out$quartiles[sams, 1])^2, 3, FUN = mean)
#     q2SMSE <- apply((bpar.out$quartiles[, 2, , ] - opar.out$quartiles[sams, 2])^2, 3, FUN = mean)
#     q3SMSE <- apply((bpar.out$quartiles[, 3, , ] - opar.out$quartiles[sams, 3])^2, 3, FUN = mean)
#     ## quarBias<-apply(bpar.out$quartiles-opar.out$rquart[num.sam],c(2,4),FUN=mean)
#     quarSMSE <- sqrt(cbind(q1SMSE, q2SMSE, q3SMSE))
#
#     f0SMSE <- apply((bpar.out$rfreq[, 1, , ] - opar.out$rfreq[sams, 1])^2, 3, FUN = mean)
#     f1SMSE <- apply((bpar.out$rfreq[, 2, , ] - opar.out$rfreq[sams, 2])^2, 3, FUN = mean)
#     f2SMSE <- apply((bpar.out$rfreq[, 3, , ] - opar.out$rfreq[sams, 3])^2, 3, FUN = mean)
#     f3SMSE <- apply((bpar.out$rfreq[, 4, , ] - opar.out$rfreq[sams, 4])^2, 3, FUN = mean)
#     f4SMSE <- apply((bpar.out$rfreq[, 5, , ] - opar.out$rfreq[sams, 5])^2, 3, FUN = mean)
#     rfreqSMSE <- cbind(f0SMSE, f1SMSE, f2SMSE, f3SMSE, f4SMSE)
#
#     d1SMSE <- apply((bpar.out$deciles[, 1, , ] - opar.out$deciles[sams, 1])^2, 3, FUN = mean)
#     d2SMSE <- apply((bpar.out$deciles[, 2, , ] - opar.out$deciles[sams, 2])^2, 3, FUN = mean)
#     d3SMSE <- apply((bpar.out$deciles[, 3, , ] - opar.out$deciles[sams, 3])^2, 3, FUN = mean)
#     d4SMSE <- apply((bpar.out$deciles[, 4, , ] - opar.out$deciles[sams, 4])^2, 3, FUN = mean)
#     d5SMSE <- apply((bpar.out$deciles[, 5, , ] - opar.out$deciles[sams, 5])^2, 3, FUN = mean)
#     d6SMSE <- apply((bpar.out$deciles[, 6, , ] - opar.out$deciles[sams, 6])^2, 3, FUN = mean)
#     d7SMSE <- apply((bpar.out$deciles[, 7, , ] - opar.out$deciles[sams, 7])^2, 3, FUN = mean)
#     d8SMSE <- apply((bpar.out$deciles[, 8, , ] - opar.out$deciles[sams, 8])^2, 3, FUN = mean)
#     d9SMSE <- apply((bpar.out$deciles[, 9, , ] - opar.out$deciles[sams, 9])^2, 3, FUN = mean)
#     deciSMSE <- sqrt(cbind(d1SMSE, d2SMSE, d3SMSE, d4SMSE, d5SMSE, d6SMSE, d7SMSE, d8SMSE, d9SMSE))
#   } else {
#     meanBias <- mean(bpar.out$mean - opar.out$mean[sams])
#
#     q1Bias <- mean(bpar.out$quartiles[, 1, , ] - opar.out$quartiles[sams, 1])
#     q2Bias <- mean(bpar.out$quartiles[, 2, , ] - opar.out$quartiles[sams, 2])
#     q3Bias <- mean(bpar.out$quartiles[, 3, , ] - opar.out$quartiles[sams, 3])
#     ## quarBias<-apply(bpar.out$quartiles-opar.out$rquart[num.sam],c(2,4),FUN=mean)
#     quarBias <- cbind(q1Bias, q2Bias, q3Bias)
#
#     f0Bias <- mean(bpar.out$rfreq[, 1, , ] - opar.out$rfreq[sams, 1])
#     f1Bias <- mean(bpar.out$rfreq[, 2, , ] - opar.out$rfreq[sams, 2])
#     f2Bias <- mean(bpar.out$rfreq[, 3, , ] - opar.out$rfreq[sams, 3])
#     f3Bias <- mean(bpar.out$rfreq[, 4, , ] - opar.out$rfreq[sams, 4])
#     f4Bias <- mean(bpar.out$rfreq[, 5, , ] - opar.out$rfreq[sams, 5])
#     rfreqBias <- cbind(f0Bias, f1Bias, f2Bias, f3Bias, f4Bias)
#
#     d1Bias <- mean(bpar.out$deciles[, 1, , ] - opar.out$deciles[sams, 1])
#     d2Bias <- mean(bpar.out$deciles[, 2, , ] - opar.out$deciles[sams, 2])
#     d3Bias <- mean(bpar.out$deciles[, 3, , ] - opar.out$deciles[sams, 3])
#     d4Bias <- mean(bpar.out$deciles[, 4, , ] - opar.out$deciles[sams, 4])
#     d5Bias <- mean(bpar.out$deciles[, 5, , ] - opar.out$deciles[sams, 5])
#     d6Bias <- mean(bpar.out$deciles[, 6, , ] - opar.out$deciles[sams, 6])
#     d7Bias <- mean(bpar.out$deciles[, 7, , ] - opar.out$deciles[sams, 7])
#     d8Bias <- mean(bpar.out$deciles[, 8, , ] - opar.out$deciles[sams, 8])
#     d9Bias <- mean(bpar.out$deciles[, 9, , ] - opar.out$deciles[sams, 9])
#     deciBias <- cbind(d1Bias, d2Bias, d3Bias, d4Bias, d5Bias, d6Bias, d7Bias, d8Bias, d9Bias)
#
#     meanSMSE <- sqrt(mean((bpar.out$mean - opar.out$mean[sams])^2))
#
#     q1SMSE <- mean((bpar.out$quartiles[, 1, , ] - opar.out$quartiles[sams, 1])^2)
#     q2SMSE <- mean((bpar.out$quartiles[, 2, , ] - opar.out$quartiles[sams, 2])^2)
#     q3SMSE <- mean((bpar.out$quartiles[, 3, , ] - opar.out$quartiles[sams, 3])^2)
#     ## quarBias<-apply(bpar.out$quartiles-opar.out$rquart[num.sam],c(2,4),FUN=mean)
#     quarSMSE <- sqrt(cbind(q1SMSE, q2SMSE, q3SMSE))
#
#     f0SMSE <- mean((bpar.out$rfreq[, 1, , ] - opar.out$rfreq[sams, 1])^2)
#     f1SMSE <- mean((bpar.out$rfreq[, 2, , ] - opar.out$rfreq[sams, 2])^2)
#     f2SMSE <- mean((bpar.out$rfreq[, 3, , ] - opar.out$rfreq[sams, 3])^2)
#     f3SMSE <- mean((bpar.out$rfreq[, 4, , ] - opar.out$rfreq[sams, 4])^2)
#     f4SMSE <- mean((bpar.out$rfreq[, 5, , ] - opar.out$rfreq[sams, 5])^2)
#     rfreqSMSE <- cbind(f0SMSE, f1SMSE, f2SMSE, f3SMSE, f4SMSE)
#
#     d1SMSE <- mean((bpar.out$deciles[, 1, , ] - opar.out$deciles[sams, 1])^2)
#     d2SMSE <- mean((bpar.out$deciles[, 2, , ] - opar.out$deciles[sams, 2])^2)
#     d3SMSE <- mean((bpar.out$deciles[, 3, , ] - opar.out$deciles[sams, 3])^2)
#     d4SMSE <- mean((bpar.out$deciles[, 4, , ] - opar.out$deciles[sams, 4])^2)
#     d5SMSE <- mean((bpar.out$deciles[, 5, , ] - opar.out$deciles[sams, 5])^2)
#     d6SMSE <- mean((bpar.out$deciles[, 6, , ] - opar.out$deciles[sams, 6])^2)
#     d7SMSE <- mean((bpar.out$deciles[, 7, , ] - opar.out$deciles[sams, 7])^2)
#     d8SMSE <- mean((bpar.out$deciles[, 8, , ] - opar.out$deciles[sams, 8])^2)
#     d9SMSE <- mean((bpar.out$deciles[, 9, , ] - opar.out$deciles[sams, 9])^2)
#     deciSMSE <- sqrt(cbind(d1SMSE, d2SMSE, d3SMSE, d4SMSE, d5SMSE, d6SMSE, d7SMSE, d8SMSE, d9SMSE))
#   }
#   # browser()
#   BiasSMSE <- cbind(meanBias, quarBias, rfreqBias, deciBias, meanSMSE, quarSMSE, rfreqSMSE, deciSMSE)
#   rnam <- colnames(BiasSMSE)
#   BiasSMSE <- t(BiasSMSE)
#   if (bpar.out$n.dist == 1) {
#     colnames(BiasSMSE) <- "only seeds"
#   } else {
#     colnames(BiasSMSE) <- c("empd.w.p0s", "empd.nw.p0sEkb", "empd.nw.p0sEks")
#   }
#   BiasSMSE
# }
#
# # ---------------------------------------------------------------------------------------#
#
# B.Bias.RMSE.ij.S <- function(net, n.seeds, n.neigh, sam.size, n.boot, otherNetParameters = FALSE) {
#   # sam.size is the number of different samples taken from the network for each i and j otherNetParameters is true if
#   # intervals and fallins for the rest of the parmeters (other than mean) are required.
#   seeds2 <- array(0, dim = c(length(as.vector(n.neigh)), length(as.vector(n.seeds)), sam.size, max(n.seeds)))
#   All.biasRMSE <- 1
#   Mean.intervals.list <- Mean.fallins <- B.mean <- as.list(rep(NA, length(n.seeds) * length(n.neigh)))
#   realparam <- summary.net(net)
#   OtherPar.intervals.list <- NULL
#   if (otherNetParameters) {
#     OtherPar.intervals.list <- OtherPar.fallins <- as.list(rep(NA, length(n.seeds) * length(n.neigh)))
#   }
#   b <- 1
#   for (i in n.seeds) {
#     for (j in n.neigh) {
#
#       if (j == 0) {
#         n.dist <- 1  #n.dist is the number of different emp distr.
#       } else {
#         n.dist <- 3
#       }
#
#       # cat('i= ',i,'\t','j= ',j,'\n') browser()
#       if (j == 0) {
#         Obs.distrib <- Oempdegreedistrib(net, n.seeds = i, n.neigh = j, num.sam = sam.size)
#         TMP <- Obs.distrib$seeds1
#       } else {
#         Obs.distrib <- Oempdegreedistrib(net, n.seeds = i, n.neigh = j, num.sam = sam.size, seeds = TMP)
#       }
#
#       Oparam <- OparametersEst(Obs.distrib)
#       # browser()
#       seeds2[which(n.neigh == j), which(n.seeds == i), , 1:dim(Obs.distrib$seeds1)[2]] <- Obs.distrib$seeds1
#
#       B.distrib <- bootdeg(Obs.distrib, num.sam = sam.size, n.boot = n.boot)
#       Bparam <- BparametersEst(B.distrib)
#       # browser()
#       B.mean[[b]] <- Bparam$mean
#
#
#       biasRMSE <- B.Bias(Bparam, Oparam)
#       All.biasRMSE <- cbind(All.biasRMSE, biasRMSE)
#
#       # mean intervals and fallins
#
#       Mean.intervals.list[[b]] <- intervals <- BMean.intervals(Opar = Oparam$mean, Bpar = Bparam$mean, n.dist)
#       Mean.fallins[[b]] <- list(ij = c(i, j), fallins = fallInsMean(ints = intervals, rpar = realparam$rmean, n.dist))
#       if (otherNetParameters)
#         {
#           OtherPar.intervals.list[[b]] <- OtherPar.intervals <- list(quartiles = Bintervalsi(Opari = Oparam$quart,
#           Bpari = Bparam$quart, n.dist, num.par = 3), rfreq = Bintervalsi(Opari = Oparam$rfreq, Bpari = Bparam$rfreq,
#           n.dist, num.par = 5), deciles = Bintervalsi(Opari = Oparam$deciles, Bpari = Bparam$deciles, n.dist, num.par = 9))
#           Other.fallins[[b]] <- list(ij = c(i, j), fi.quart = NA, fi.rfreq = NA, fi.deciles = NA)
#         }  #if(OtherNetParameters)
#
#       b <- b + 1
#     }
#   }
#   list(All.biasRMSE = All.biasRMSE, Mean.intervals = Mean.intervals.list, Mean.fallins = Mean.fallins, B.mean = B.mean,
#     seeds2 = seeds2)
# }
#
# # ---------------------------------------------------------------------------------------#
#
# BMean.intervals <- function(Opar, Bpar, n.dist) {
#
#   # Opar is output of OparametersEst Opar$mean Bpar is output of BparametersEst Bpar$mean n.dist is the number of empirical
#   # distributions considered browser()
#   if (n.dist > 1) {
#     mean.int <- apply(Bpar, c(1, 3), FUN = quantile, prob = c(0.025, 0.975))  #debe ser una matriz 2x tn.sam x n.dist
#     mean.int1 <- cbind(2 * Opar - mean.int[2, , ], 2 * Opar - mean.int[1, , ])[, c(1, 4, 2, 5, 3, 6)]  #firts to columns dist1, then dist2 and finally dist3
#     mean.int2 <- cbind(Opar + apply(Bpar, c(1, 3), FUN = mean) - mean.int[2, , ], Opar + apply(Bpar, c(1, 3), FUN = mean) -
#       mean.int[1, , ])[, c(1, 4, 2, 5, 3, 6)]
#     mean.int3 <- cbind(Opar^2/mean.int[2, , ], 2 * Opar^2/mean.int[1, , ])[, c(1, 4, 2, 5, 3, 6)]
#     mean.int <- cbind(t(mean.int[, , 1]), t(mean.int[, , 2]), t(mean.int[, , 3]))
#   } else if (n.dist == 1)
#     {
#       mean.int <- t(apply(Bpar, c(1), FUN = quantile, prob = c(0.025, 0.975)))  #debe ser una matriz tn.sam x2
#       mean.int1 <- cbind(2 * Opar - mean.int[, 2], 2 * Opar - mean.int[, 1])
#       mean.int2 <- cbind(Opar + apply(Bpar, c(1), FUN = mean) - mean.int[, 2], Opar + apply(Bpar, c(1), FUN = mean) -
#         mean.int[, 1])
#       mean.int3 <- cbind(Opar^2/mean.int[, 2], 2 * Opar^2/mean.int[, 1])
#     }  #if (n.dist>1)
#   list(mean.int = mean.int, mean.int1 = mean.int1, mean.int2 = mean.int2, mean.int3 = mean.int3)
# }
#
# # ---------------------------------------------------------------------------------------#
#
# fallInsMean <- function(ints, rpar, n.dist) {
#   # ints in output of Bintervals (mean.int, mean.int1,mean.int2,mean.int3) rpar is output of summary.net rpar$rmean
#   # n.dist is the number of empirical distributions browser()
#   if (n.dist > 1) {
#     fi.mean <- fi1.mean <- fi2.mean <- fi3.mean <- rep(0, 3)
#     col <- c(1, 3, 5)
#     for (k in 1:3) {
#       fi.mean[k] <- sum(ints$mean.int[, col[k]] <= rpar & rpar <= ints$mean.int[, col[k] + 1])/dim(ints$mean.int)[1]
#       fi1.mean[k] <- sum(ints$mean.int1[, col[k]] <= rpar & rpar <= ints$mean.int1[, col[k] + 1])/dim(ints$mean.int)[1]
#       fi2.mean[k] <- sum(ints$mean.int2[, col[k]] <= rpar & rpar <= ints$mean.int2[, col[k] + 1])/dim(ints$mean.int)[1]
#       fi3.mean[k] <- sum(ints$mean.int3[, col[k]] <= rpar & rpar <= ints$mean.int3[, col[k] + 1])/dim(ints$mean.int)[1]
#     }
#   } else {
#     fi.mean <- sum(ints$mean.int[, 1] <= rpar & rpar <= ints$mean.int[, 2])/dim(ints$mean.int)[1]
#     fi1.mean <- sum(ints$mean.int1[, 1] <= rpar & rpar <= ints$mean.int1[, 2])/dim(ints$mean.int)[1]
#     fi2.mean <- sum(ints$mean.int2[, 1] <= rpar & rpar <= ints$mean.int2[, 2])/dim(ints$mean.int)[1]
#     fi3.mean <- sum(ints$mean.int3[, 1] <= rpar & rpar <= ints$mean.int3[, 2])/dim(ints$mean.int)[1]
#   }
#   fi <- rbind(fi.mean = fi.mean, fi1.mean = fi1.mean, fi2.mean = fi2.mean, fi3.mean = fi3.mean)
#   fi
# }
# # ---------------------------------------------------------------------#
#
#
# Bintervalsi <- function(Opari, Bpari, n.dist, num.par) {
#   # Opari is the parameter(s) of interesout from OparametersEst (either: mean,quartiles,rfreq,deciles) Bpari is the
#   # parameter(s) of interes from BparametersEst (either: mean,quartiles,rfreq,deciles) n.dist is the number of empirical
#   # distributions considered
#
#   reorder <- as.vector(matrix(rep(1:(num.par), each = 2), 2, num.par) + c(0, num.par))
#   # browser()
#   if (n.dist > 1) {
#     par.int <- apply(Bpari, c(1, 2, 4), FUN = quantile, prob = c(0.025, 0.975))  #debe ser una matriz: 2 x num.sam x num param x n.dist
#
#     par.int1.a <- apply(-par.int[2, , , ], 3, FUN = "+", 2 * Opari)  #matrix (num.sam)(num.par) x ndist #the first col is the vect of a matrix num.sam x num.par con dist1
#     par.int1.b <- apply(-par.int[1, , , ], 3, FUN = "+", 2 * Opari)  # '
#     par.int1 <- cbind(par.int1.a, par.int1.b)[, reorder]  #the first two columns are the lower and upper limit according to n.dist1...
#
#     ##### test for the quartiles plot(1:60,rep(5,60),t='n',ylim=c(-3,6)) segments(1:60,par.int[,1],1:60,par.int[,2])
#
#     par.int2.a <- apply(apply(Bpari, c(1, 2, 4), FUN = mean) - par.int[2, , , ], 3, FUN = "+", Opari)
#     par.int2.b <- apply(apply(Bpari, c(1, 2, 4), FUN = mean) - par.int[1, , , ], 3, FUN = "+", Opari)
#     par.int2 <- cbind(par.int2.a, par.int2.b)[, reorder]
#
#     par.int3.a <- apply(1/par.int[2, , , ], 3, FUN = "*", Opari)
#     par.int3.b <- apply(1/par.int[1, , , ], 3, FUN = "*", Opari)
#     par.int3 <- cbind(par.int3.a, par.int3.b)[, reorder]
#
#     ##### segments(1:60,par.int3[,1],1:60,par.int3[,2])
#   } else if (n.dist == 1)
#     {
#       par.int <- apply(Bpari, c(1, 2), FUN = quantile, prob = c(0.025, 0.975))  #debe ser una matriz 2 x num.sam x num.par
#       par.int1 <- cbind(2 * Opari - par.int[2, , ], 2 * Opari - par.int[1, , ])[, reorder]  #debe ser una matriz num.sam x 2(num.par)
#       par.int2 <- cbind(Opari + apply(Bpari, c(1, 2), FUN = mean) - par.int[2, , ], Opari + apply(Bpari, c(1, 2), FUN = mean) -
#         par.int[1, , ])[, reorder]
#       par.int3 <- cbind(Opari^2/par.int[2, , ], 2 * Opari^2/par.int[1, , ])[, reorder]
#       mpar.int <- NULL
#       for (w in 1:num.par) mpar.int <- cbind(mpar.int, t(par.int[, , w]))
#     }  #if (n.dist>1)
#   list(par.int = mpar.int, par.int1 = par.int1, part.int2 = par.int2, par.int3 = par.int3, num.par = num.par)
# }
#
# # ---------------------------------------------------------------------#
# fallInsi <- function(intsi, rpari, n.boot, n.dist) {
#   # intsi is output of Bintervalsi (for quarts, rfreq or deciles) rpari is the parameter(s) of interes from summary.net
#   # (either:rquart,rfreq,rdeci) n.dist is the number of empirical distributions considered
#
#   if (n.dist > 1) {
#     fi.mean <- fi1.mean <- fi2.mean <- fi3.mean <- rep(0, 3)
#     col <- c(1, 3, 5)
#     for (k in 1:3) {
#       fi.mean[k] <- sum(ints$mean.int[, col[k]] <= rpar & rpar <= ints$mean.int[, col[k] + 1])/dim(ints$mean.int)[1]
#       fi1.mean[k] <- sum(ints$mean.int1[, col[k]] <= rpar & rpar <= ints$mean.int1[, col[k] + 1])/dim(ints$mean.int)[1]
#       fi2.mean[k] <- sum(ints$mean.int2[, col[k]] <= rpar & rpar <= ints$mean.int2[, col[k] + 1])/dim(ints$mean.int)[1]
#       fi3.mean[k] <- sum(ints$mean.int3[, col[k]] <= rpar & rpar <= ints$mean.int3[, col[k] + 1])/dim(ints$mean.int)[1]
#     }
#   } else {
#     fi.mean <- sum(ints$mean.int[, 1] <= rpar & rpar <= ints$mean.int[, 2])/dim(ints$mean.int)[1]
#     fi1.mean <- sum(ints$mean.int1[, 1] <= rpar & rpar <= ints$mean.int1[, 2])/dim(ints$mean.int)[1]
#     fi2.mean <- sum(ints$mean.int2[, 1] <= rpar & rpar <= ints$mean.int2[, 2])/dim(ints$mean.int)[1]
#     fi3.mean <- sum(ints$mean.int3[, 1] <= rpar & rpar <= ints$mean.int3[, 2])/dim(ints$mean.int)[1]
#   }
#   fi <- rbind(fi.mean = fi.mean, fi1.mean = fi1.mean, fi2.mean = fi2.mean, fi3.mean = fi3.mean)
#   fi
# }
#
# ###############################################################################
#
# #
# # load(paste(distr, "10000x", order, ".RData", sep = ""))
# #
# # MC <- 200  #length(networks)
# # nsam <- 2
# # n.boot <- 500
# # # n.seeds<-c(3,5,7,10,20,50,100)#
# # n.neigh <- c(0, 1, 2, 3, 4, 5)  #
# # # OBSME <- MC.BBRMSE <- 0
# # # MC.int.mean<-matrix(0,(length(n.seeds)*length(n.neigh)*4),6,dimnames=list(c(rep.int(rep.int(c('int','int1','int2','int3'),length(n.neigh)),length(n.seeds))),c(rep.int(c('2.5%','97.5%'),3))))
# # # MC.FI<-array(0,dim=c((length(n.seeds)*length(n.neigh)*4),3,MC),
# # # dimnames=list(c(rep.int(rep.int(c('int','int1','int2','int3'),length(n.neigh)),length(n.seeds))),c('dist1','dist2','dist3')))
# # Bmean <- array(NA, dim = c(length(n.neigh), length(n.seeds), n.boot, nsam, MC))
# # # sample.means<-array(NA, dim=c(length(n.neigh), length(n.seeds), nsam, MC))
# # seeds <- array(0, dim = c(MC, length(n.neigh), length(n.seeds), nsam, max(n.seeds)))
# # for (mc in (n * MC + 1):(n * MC + MC)) {
# #   # mc=1
# #   net <- networks[[mc]]
# #   # rp<-summary.net(net) OBSME<-OBSME+Bias.SMSE.ij.S(net,rp,n.seeds,n.neigh,nsam)/MC
# #
# #   # Bootstrap
# #   BBRMSE <- B.Bias.RMSE.ij.S(net, n.seeds, n.neigh, nsam, n.boot)
# #   # MC.BBRMSE<-BBRMSE$All.biasRMSE/MC+MC.BBRMSE
# #
# #   seeds[(mc - n * MC), , , , ] <- BBRMSE$seeds2
# #
# #   # INTERVALS for(i in 1:(length(n.seeds)*length(n.neigh))){
# #   # MC.int.mean[(i*4-3),1:ncol(BBRMSE$Mean.intervals[[i]]$mean.int)]<-MC.int.mean[(i*4-3),1:ncol(BBRMSE$Mean.intervals[[i]]$mean.int)]+apply(BBRMSE$Mean.intervals[[i]]$mean.int,2,mean)/MC
# #   # MC.int.mean[(i*4-2),1:ncol(BBRMSE$Mean.intervals[[i]]$mean.int1)]<-MC.int.mean[(i*4-2),1:ncol(BBRMSE$Mean.intervals[[i]]$mean.int1)]+apply(BBRMSE$Mean.intervals[[i]]$mean.int1,2,mean)/MC
# #   # MC.int.mean[(i*4-1),1:ncol(BBRMSE$Mean.intervals[[i]]$mean.int2)]<-MC.int.mean[(i*4-1),1:ncol(BBRMSE$Mean.intervals[[i]]$mean.int2)]+apply(BBRMSE$Mean.intervals[[i]]$mean.int2,2,mean)/MC
# #   # MC.int.mean[(i*4),1:ncol(BBRMSE$Mean.intervals[[i]]$mean.int3)]<-MC.int.mean[(i*4),1:ncol(BBRMSE$Mean.intervals[[i]]$mean.int3)]+apply(BBRMSE$Mean.intervals[[i]]$mean.int3,2,mean)/MC
# #   # }
# #
# #   # #FALLINS for(i in 1:(length(n.seeds)*length(n.neigh)))
# #   # MC.FI[(i*4-3):(i*4),1:ncol(BBRMSE$Mean.fallins[[i]]$fallins),(mc-n*MC)]=BBRMSE$Mean.fallins[[i]]$fallins
# #
# #   # #Intermediate output print(c('Monte Carlo', mc)); print(apply(MC.FI[,,1:(mc-n*MC)], c(1,2), sum)/(mc-n*MC))
# #   print(mc)
# #
# #   # B.means[[(mc-n*MC)]]<-BBRMSE$B.mean
# #   b <- 1
# #   for (j in 1:length(n.seeds)) {
# #     for (i in 1:length(n.neigh)) {
# #       # B.mean is a list of lists: 1-5mc <= 1-24combinations of n.neigh and n.seeds <= <=nsam*n.boot*3 matrix, where 1-3 are
# #       # different distributions, and only 1 distribution exists for n.neigh=0 (i.e. nsam*n.boot matrix). So, we take 1st
# #       # distribution for each of 24 combinations of seeds and neighbours, apply mean over bootstraps and get 1*nsam vector of
# #       # means, which we sort and embed in the sample.means: sample.means[i,j,,(mc-n*MC)]<-sort(apply(BBRMSE$B.mean[[b]][,,1],
# #       # 1, mean))
# #       Bmean[i, j, , , (mc - n * MC)] <- t(BBRMSE$B.mean[[b]][, , 1])
# #       b <- b + 1
# #     }
# #   }
# # }
# # save(seeds, file = paste(distr, order, "_seeds_", n, ".RData", sep = ""))
# # save(Bmean, file = paste(distr, order, "_BMean_", n, ".RData", sep = ""))
# # # save(MC.FI, file = paste(distr, order, '_MCFI_', n, '.RData', sep=''))
# #
# # # write.table(MC.BBRMSE,file = paste(distr, order, '_MCBBRMSE_', n,'.txt', sep='')) write.table(MC.int.mean,file =
# # # paste(distr, order,'_MCintmean_', n,'.txt', sep=''))
