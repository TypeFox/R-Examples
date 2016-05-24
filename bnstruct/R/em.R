#' @rdname em
#' @aliases em,InferenceEngine,BNDataset
setMethod("em",
          c("InferenceEngine","BNDataset"),
          function(x, dataset, threshold = 0.001, max.em.iterations = 10, ess = 1)
          {
            # We assume:
            # 1) there is a BN with learnt or known parameters
            # 2) there is a raw dataset with missing values (raw complete data case is not treated)
            #
            # We do:
            # 1) belief propagation
            # 2) re-guess imputed data, substituting NAs with most probable values obtained with BP step above
            # 3) re-learn parameters for BN from imputed data from step 2
            # 4) if convergence condition is met (parameters vary by less than threshold) stop, otherwise go back to step 1
            
            # for steps 1,2:
            # - do bp
            # - get computed network
            # - get cpts of computed network
            # - for NA values:
            #   - select corresponding cpts
            #   - set parents to observed values
            #     - in case of parents with missing values: take variables in topological order
            #   - compute values according to distribution

            bnstruct.start.log("starting EM algorithm ...")
            
            rawdata  <- raw.data(dataset)
            if (test.updated.bn(x))
              orig.bn <- updated.bn(x)
            else
              orig.bn <- bn(x)
            
            bak.bn     <- orig.bn
            bn         <- orig.bn
            eng        <- x
            num.items  <- num.items(dataset)
            num.nodes  <- num.nodes(orig.bn)
            node.sizes <- node.sizes(orig.bn)
            var.names  <- variables(orig.bn)
            
            cliques    <- jt.cliques(x)
            
            # ndataset <- dataset
            no.iterations <- 1
            difference    <- threshold + 1
            prev.log.lik  <- 0
            while((difference > threshold && no.iterations <= max.em.iterations) || no.iterations < 2)
            {
              observations(eng) <- list(NULL, NULL) # clean observations, if needed
              eng      <- belief.propagation(eng)            
              imp.data <- matrix(data=rep(0, prod(dim(rawdata))),
                                 nrow=nrow(rawdata),
                                 ncol=ncol(rawdata))
              
              still.has.NAs <- c()
              for (row in 1:num.items)
              {
                y             <- rawdata[row,]
                # mpv           <- rep(0, num.nodes)
                obsd.vars     <- which(!is.na(y))
                obsd.vals     <- y[obsd.vars]
                non.obsd.vars <- setdiff(1:num.nodes, obsd.vars)
                mpv           <- c(y)
                
                if (length(non.obsd.vars) == 0)
                {
                  imp.data[row,] <- mpv
                  next
                }
  
                j                 <- jpts(eng)
                to.evaluate       <- non.obsd.vars
                overall.obsd.vars <- obsd.vars
                overall.obsd.vals <- obsd.vals
                
                
                while(length(to.evaluate) > 0)
                {
                  to.evaluate.next <- c()
                  for (i in to.evaluate)
                  {
                    # TODO: PROFILING: which is really better?
                    #target.cliques <- which(!is.na(sapply(cliques, function(cl) {match(i, c(cl))})))
                    target.cliques <- which(sapply(cliques, function(cl) {i %in% cl}))
                    
                    tc <- 1
                    while (tc  <= length(target.cliques))
                    {
                      target.clique  <- target.cliques[tc]
                      jpt            <- j[[target.clique]]
                      # TODO: PROFILING: can this be improved?
                      d              <- c(match(names(dimnames(jpt)),var.names))
                      
                      # TODO: PROFILING: is there anything better than intersect?
                      if (length(intersect(d, overall.obsd.vars)) == length(d)-1)
                      {
                        dd <- d
                        for (v in intersect(d, overall.obsd.vars))
                        {
                          #dmnms <- dimnames(jpt)
                          #nms   <- names(dimnames(jpt))
                          cpt         <- rep(0, node.sizes[v])
                          cpt[mpv[v]] <- 1
                          out         <- mult(jpt, dd, cpt, c(v), node.sizes)
                          jpt         <- out$potential
                          dd          <- out$vars
                          jpt         <- jpt / sum(jpt)
                          #dimnames(jpt) <- dmnms
                          #names(dimnames(jpt)) <- nms
                          # out         <- marginalize(jpt, dd, v)
                          # jpt         <- out$potential
                          # dd          <- out$vars
                          remaining <- (1:length(dd))[-which(dd == v)]
                          dd <- dd[remaining]
                          jpt <- apply(jpt, remaining, sum)
                        }
                        
                        if (length(dd) == 1 && !is.element(NaN,jpt) && !is.element(NA,jpt))
                        {
                          tc <- length(target.cliques) + 100
                          
                          # TODO: PROFILING: can this be done more efficiently?
                          wm <- which(!is.na(match(c(jpt),max(jpt))))
                          if (length(wm) == 1)
                          {
                            mpv[i] <- wm # jpt[wm]
                          }
                          else
                          {
                            mpv[i] <- sample(wm,1) #,replace=TRUE
                          }
                          overall.obsd.vars <- sort(c(overall.obsd.vars,i))
                        }
                      }
                      tc <- tc + 1
                    }
  
                    if(tc <= length(target.cliques)+2 || is.element(NaN,jpt) || is.element(NA, jpt))
                    {
                      to.evaluate.next <- c(to.evaluate.next, i)
                    }
                  }

                  if (length(to.evaluate.next) < length(to.evaluate))
                  {
                    to.evaluate <- to.evaluate.next
                  }
                  else
                  {
                    # ok, for now skip, or will loop.
                    # find out how to deal with this.
                    to.evaluate   <- c()
                    still.has.NAs <- c(still.has.NAs, row)
                  }
                }
  
                imp.data[row,] <- mpv
              }
              
              # Fill in holes left by rows with too many NAs for being all identified.
              # Use the network
              if (length(still.has.NAs) > 0)
              {
                bis.data <- imp.data[-still.has.NAs,]
                bis.dataset <- dataset
                imputed.data(bis.dataset) <- bis.data
                num.items(bis.dataset) <- num.items - length(still.has.NAs)
                bis.net <- learn.params(bn, bis.dataset, ess=ess, use.imputed.data=T)

                for (bis.row in still.has.NAs)
                {
                  bis.ie <- InferenceEngine(bis.net)
                  ov     <- which(!is.na(imp.data[bis.row]))
                  bis.ie.1 <- belief.propagation(bis.ie, list("observed.vars" = ov,
                                                              "observed.vals" = (imp.data[bis.row])[ov]))
                  imp.data[bis.row,] <- get.most.probable.values(bis.ie.1)
                }
              }
              
              storage.mode(imp.data) <- "integer"
              imputed.data(dataset)  <- imp.data
              
              bn <- learn.params(bn, dataset, ess=ess, use.imputed.data=T)

              no.iterations <- no.iterations + 1
              curr.log.lik  <- log.likelihood(dataset, bn, use.imputed.data = T)
              difference    <- prev.log.lik - curr.log.lik
              
              orig.bn      <- bn
              prev.log.lik <- curr.log.lik
            }
            
            updated.bn(x) <- bn
            
            bnstruct.end.log("EM algorithm completed.")
            
            return(list("InferenceEngine" = x, "BNDataset" = dataset))
          })


# Compute log likelihood of the dataset given the network.
# Takes a BNDataset and a BN as input.
log.likelihood <- function(dataset, net, use.imputed.data = FALSE)
{
  dag        <- dag(net)
  node.sizes <- node.sizes(net)
  n.nodes    <- num.nodes(net)
  variables  <- variables(dataset)
  
  parents    <- lapply(1:n.nodes, function(x) which(dag[,x] != 0))
  jpts       <- lapply(cpts(net), function(x) x / sum(x))
  dim.vars   <- lapply(1:n.nodes,
                       function(x)
                         match(
                           c(unlist(
                             names(dimnames(jpts[[x]])), F, F
                           )),
                           variables
                         )
  )
  
  if (use.imputed.data)
    data <- imputed.data(dataset)
  else
    data <- raw.data(dataset)
  
  sorted.nodes <- topological.sort(dag)
  vals         <- rep(0, n.nodes)
  ll           <- 0
  
  for (row in 1:nrow(data))
    for (node in sorted.nodes)
    {
      if (length(parents[[node]]) == 0) {
        vals[node] <- data[row,node]
        ll <- ll + log(jpts[[node]][data[row,node]])
      } else {
        jpt  <- jpts[[node]]
        vars <- c(unlist(dim.vars[[node]]))
        for (p in parents[[node]]) {
          sumout          <- rep(0, node.sizes[p])
          sumout[vals[p]] <- 1
          out  <- mult(jpt, vars, sumout, c(p), node.sizes)
          jpt  <- out$potential
          vars <- out$vars
        }
        jpt <- c(jpt[which(jpt != 0)])
        vals[node] <- data[row,node]
        ll <- ll + log(jpt[data[row,node]])
      }
    }
  
  return(ll)
}

