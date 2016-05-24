Oempdegreedistrib0 <- function(net, n.seeds, n.neigh, num.sam, seeds) {
      p0.seeds.array <- Oempd <- values.array <- val.seeds.array <- samples <- as.list(rep(NA, num.sam))
      seeds1 <- matrix(NA, num.sam, n.seeds)
      ## -------the 'real' parameters in the network:-------##
      real <- summary_net(net)
      realdd <- real$realdd
      ## ---------------------------------------------------##
      for (m in 1:num.sam) {
            # if(m%%100==1)#cat('Obtaining empd of sample ',m,'\n')
            neigh.seeds <- sort(sample(1:length(net$degree), n.seeds, replace = FALSE))  #n.neigh=0!!!!!!!
            tab.seeds <- table(neigh.seeds)  #id seeds
            seeds1[m, ] <- neigh.seeds
            ###### degrees #####
            deg.seeds <- realdd[rep(as.integer(names(tab.seeds)), tab.seeds)]  #duplicities allowed
            # --------------------#
            samples[[m]] <- list(freq.deg.seeds = freq.deg.seeds <- table(deg.seeds))
            ##### resample and extract the degree of selected vertices######
            values.array[[m]] <- values <- val.seeds.array[[m]] <- val.seeds <- sort(unique(deg.seeds))
            p0.seeds.array[[m]] <- sum(deg.seeds == 0)/n.seeds
            #### Frequency ##### (Not the relative frequency)
            OFseed <- table.row(deg.seeds, values)
            #### Empirical degree distribution ########
            Oempd.seeds <- OFseed/n.seeds
            Oempd[[m]] <- list(Oempd = Oempd.seeds)
      }  # for(m in 1:num.sam)
      # browser()
      list(samples = samples, values = values.array, Oempd = Oempd, num.sam = num.sam, val.seeds = val.seeds.array,
           n.seeds = n.seeds, n.neigh = n.neigh, p0.seeds = p0.seeds.array, seeds1 = seeds1)
}
