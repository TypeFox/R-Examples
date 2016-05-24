# ' @rdname sem
# ' @aliases sem,InferenceEngine,BNDataset
setMethod("sem",
          c("BN","BNDataset"),
          function(x, dataset, struct.threshold = 0, param.threshold = 0, 
                   max.sem.iterations = 25, max.em.iterations = 10, scoring.func = "BDeu",
                   initial.network = NULL, alpha = 0.05, ess = 1, bootstrap = FALSE,
                   layering = c(), max.fanin.layers = NULL,
                   max.fanin = num.variables(dataset), cont.nodes = c(), use.imputed.data = FALSE,
                   use.cpc = TRUE, ...)
          {
            net <- x
            
            num.nodes <- num.nodes(net)

            if (is.character(scoring.func))
              scoring.func <- match(tolower(scoring.func), c("bdeu", "aic", "bic"))

            if (is.na(scoring.func))
            {
              message("scoring function not recognized, using BDeu")
              scoring.func <- 0
            }
            else
              scoring.func <- scoring.func - 1
            # scoring.func(bn) <- c("BDeu", "AIC", "BIC")[scoring.func + 1]

            # starting from an empty network: learn a starting point using MMHC
            if (is.null(initial.network) || is.na(sum(dag(net))) || sum(dag(net)) == 0)
            {
              #w.net <- net
              #w.net <- learn.network(w.net, dataset, "mmhc", c("bdeu", "aic", "bic")[scoring.func+1],
              #                       alpha=alpha, ess = ess, bootstrap = bootstrap,
              #                       layering = layering, max.fanin.layers = max.fanin.layers,
              #                       max.fanin = max.fanin, cont.nodes = cont.nodes,
              #                       use.imputed.data = F, use.cpc = use.cpc, ...)
              w.net <- sample.chain(dataset)
              
            }
            else
            {
              # start from an already learnt network
              w.net     <- initial.network #net
            }
            
            w.dataset <- dataset
            w.eng     <- InferenceEngine(w.net)
            
            sem.iterations <- 0

            repeat
            {
              out <- em(w.eng, dataset, param.threshold, max.em.iterations, ess)
              
              new.eng     <- out$InferenceEngine
              new.dataset <- out$BNDataset
              
              new.net <- learn.network(new.dataset, "mmhc", c("bdeu", "aic", "bic")[scoring.func+1],
                                       initial.network = w.net, # NULL,
                                       alpha = alpha, ess = ess, bootstrap = bootstrap,
                                       layering = layering, max.fanin.layers = max.fanin.layers,
                                       max.fanin = max.fanin, cont.nodes = cont.nodes,
                                       use.imputed.data = T, use.cpc = use.cpc, ...)
              
              difference <- shd(dag(w.net), dag(new.net))
              
              w.net     <- new.net
              w.dataset <- new.dataset
              
              sem.iterations <- sem.iterations + 1
              
              if (difference <= struct.threshold || sem.iterations >= max.sem.iterations)
                break
              else
                w.eng     <- InferenceEngine(w.net)
            }

            
            return(w.net)
          })
