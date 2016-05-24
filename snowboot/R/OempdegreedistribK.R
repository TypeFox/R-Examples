OempdegreedistribK <- function(net, n.seeds, n.neigh, num.sam, seeds) {
      # This function obtains the empirical degree distribution from num.sam samples net is the network (only one) n.seeds is
      # the number of seeds to set the neighbourhood sample n.neigh is the neighbouhood size around each seeds num.sam is the
      # number of different samples taken from the same network idname is to identify from which nets we are sampling and
      # resampling.

      seeds1 <- matrix(0, num.sam, n.seeds)
      nodes_of_LSMI <- list()
      p0.seeds.array <- Oempd <- ekseed.array <- values.array <- val.seeds.array <- val.nonseed.array <- samples <- as.list(rep(NA,
                                                                                                                              num.sam))

      ## -------the 'real' parameters in the network:-------##
      real <- summary_net(net)
      realdd <- real$realdd
      # rmeand<-real$mean(realdd) rquart<-real$rquart
      rfreq <- real$rfreq
      # rperc<-real$rperc -------------------------------------##

      for (m in 1:num.sam) {
            # if(m%%100==1)#cat('Obtaining empd of sample ',m,'\n') browser()
            neigh <- LSMI(net, n.seeds = n.seeds, n.neigh = n.neigh, seeds = seeds[m, ])
            seeds1[m, ] <- neigh$seeds
            # nodes<-neigh$sampleN[!is.element(neigh$sampleN,neigh$last.added)] #vertices that are not included last (with their
            # duplicities)
            nodes <- neigh$sampleN  #vertices up to distance n.neigh(with their duplicities)!!!!!!!!!!!!!!
            nodes_of_LSMI <- c(nodes_of_LSMI, list(nodes))
            tab.nodes <- table(nodes)  #now it has the info of seeds and non-seeds
            # Now we want to distinguish between seeds and non seeds. Remember that some seeds can also be non seeds if they were also
            # selected by following one edge.
            tab.seeds <- table(neigh$seeds)  #id seeds
            a <- is.element(names(tab.nodes), names(tab.seeds))
            if (all(names(tab.nodes[a]) == names(tab.seeds))) {
                  tab.nodes[a] <- tab.nodes[a] - tab.seeds
            } else {
                  cat("no mismo orden")
                  browser()
            }
            tab.nodes <- tab.nodes[tab.nodes > 0]  #now I have only the id of non-seeds
            ###################### degrees #####
            deg.seeds <- realdd[rep(as.integer(names(tab.seeds)), tab.seeds)]  #in case seeds are present more than once as seeds
            deg.nonseed <- realdd[rep(as.integer(names(tab.nodes)), tab.nodes)]  #to incorporate their duplicity
            deg.nonseedU <- realdd[as.integer(names(tab.nodes))]  #it contains the non seeds only one time (I do not use in the rest of the code)
            # --------------------#
            samples[[m]] <- list(freq.deg.seeds = freq.deg.seeds <- table(deg.seeds), freq.deg.nonseed = freq.deg.nonseed <- table(deg.nonseed),
                                 freq.deg.nonseedU = freq.deg.nonseedU <- table(deg.nonseedU))
            ##### resample and extract the degree of selected vertices######
            val.seeds.array[[m]] <- val.seeds <- sort(unique(deg.seeds))  #as.numeric(names(freq.deg.seeds))
            val.nonseed.array[[m]] <- val.nonseed <- sort(unique(deg.nonseed))  #as.numeric(names(freq.deg.nonseed))
            val.nonseedU <- sort(unique(deg.nonseedU))  #as.numeric(names(freq.deg.nonseedU))

            p0.real <- rfreq[1]
            p0.seeds <- 0
            if (any(val.seeds == 0)) {
                  # if any seeds has degree zero
                  p0.seeds <- sum(deg.seeds == 0)/n.seeds
            }
            p0.seeds.array[[m]] <- p0.seeds
            values <- sort(union(val.seeds, val.nonseed))  #all the possible degree values toresample
            values.array[[m]] <- values

            ###################### Frequency ##### (Not the relative frequency)
            OFseed <- table.row(deg.seeds, values)
            OFnonseed <- table.row(deg.nonseed, values)
            #################################################################### combining information from seeds and nonseeds ###### mean degree computed from the original sampled seeds:
            ekseed.array[[m]] <- ekseed <- sum(as.numeric(names(freq.deg.seeds)) * freq.deg.seeds)/n.seeds
            colzero <- NULL
            if (any(values == 0)) {
                  colzero <- which(values == 0)
                  vals <- values[-colzero]
                  Of.seeds <- OFseed[-colzero]
                  Of.nonseed.nw <- OFnonseed[-colzero]
            } else {
                  vals <- values
                  Of.seeds <- OFseed
                  Of.nonseed.nw <- OFnonseed
            }
            ################################################ NWB # seeds and non weighted nonseeds #### p0 estimated from orginal sampled seeds#
            Oempd.nw.p0sEks <- (Of.seeds + (1 - p0.seeds) * ekseed * Of.nonseed.nw/vals)/(n.seeds + ekseed * sum(Of.nonseed.nw/vals))
            ######################################### p0 taken as known #
            Oempd.nw.p0rEks <- (Of.seeds + (1 - p0.real) * ekseed * Of.nonseed.nw/vals)/(n.seeds + ekseed * sum(Of.nonseed.nw/vals))
            if (any(values == 0)) {
                  Oempd.nw.p0sEks <- c(p0.seeds, Oempd.nw.p0sEks)  #5
                  Oempd.nw.p0rEks <- c(p0.real, Oempd.nw.p0rEks)  #8
            }
            Oempd[[m]] <- list(Oempd = Oempd.nw.p0sEks, Oempd.nw.p0rEks = Oempd.nw.p0rEks)
      }  # for(m in 1:num.sam)
      # browser()
      list(samples = samples, values = values.array, Oempd = Oempd, num.sam = num.sam, val.seeds = val.seeds.array,
           val.nonseed = val.nonseed.array, n.seeds = n.seeds, n.neigh = n.neigh, p0.real = p0.real, p0.seeds = p0.seeds.array,
           ekseed = ekseed.array, seeds1 = seeds1, nodes_of_LSMI=nodes_of_LSMI)
}
