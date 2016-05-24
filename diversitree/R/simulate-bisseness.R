## 12. tree.bisseness
tree.bisseness <- function(pars, max.taxa=Inf, max.t=Inf,
                           include.extinct=FALSE, x0=NA) {
  check.pars.bisseness(pars)
  if ( is.na(x0) )
    x0 <- as.integer(runif(1) > stationary.freq.bisseness(pars)[1])
  else if (length(x0) != 1 || !(x0 == 0 || x0 == 1)) 
    stop("Invalid root state")

  info <- make.tree.bisseness(pars, max.taxa, max.t, x0)
  phy <- me.to.ape.bisse(info[-1, ], info$state[1])
  if ( include.extinct || is.null(phy) )
    phy
  else
    prune(phy)
}

## 13. make.tree.bisseness
make.tree.bisseness <- function(pars, max.taxa=Inf, max.t=Inf, x0,
                                single.lineage=TRUE)  {
  p.pars<- matrix(c(1-pars[c(7,9)],
                    pars[c(7,9)]*pars[c(8,10)],
                    pars[c(7,9)]*(1-pars[c(8,10)])), 2, 3,
                  byrow=FALSE)
  pars <- matrix(pars[1:6], 2, 3)
  
  extinct <- FALSE
  split <- FALSE
  parent <- 0
  n.i <- c(0, 0)
  r.i <- rowSums(pars)
  len <- 0
  t <- 0
  hist <- list()
  
  if ( single.lineage ) {
    states <- x0
    n.taxa <- lineages <- n.i[x0 + 1] <- 1
    start <- 0
  } else {
    stop("Nope.")
  }
  
  while ( n.taxa <= max.taxa && n.taxa > 0 ) {
    ## When does an event happen?
    r.n <- r.i * n.i
    r.tot <- sum(r.n)
    dt <- rexp(1, r.tot)
    t <- t + dt
    if ( t > max.t ) {
      dt <- dt - (t - max.t)
      len[lineages] <- len[lineages] + dt
      t <- max.t
      break
    }

    len[lineages] <- len[lineages] + dt

    ## Proceed.  What state does an event happen to?
    state <- as.integer(runif(1) > r.n[1]/r.tot)
    state.i <- state + 1

    ## Pick a lineage for that state:
    j <- sample(n.i[state.i], 1)
    lineage <- lineages[states[lineages] == state][j]

    ## Pick an event: 1= speciation, 2= extinction, 3= state change
    type <- sample(3, 1, FALSE, pars[state.i, ])

    if ( type == 1 ) {
      ## Speciating:      
      if (n.taxa == max.taxa) 
        break
      new.i <- length(extinct) + 1:2
      split[lineage] <- TRUE
      extinct[new.i] <- split[new.i] <- FALSE
      
      if ( p.pars[state.i,1]==1 ) {
        states[new.i] <- state
        n.i[state.i] <- n.i[state.i] + 1
      } else {
        ## The daughter states be inherited from the parent state? 
        ## 1= complete inheritance (i.e. traditional bisse), 
        ## 2= partial inheritance, 
        ## 3= no inheritance.
        inherit.type<- sample(3, 1, FALSE, p.pars[state.i,])

        ## RGF: I suspect that moving the statements from this switch
        ## into the three-way if below will be more efficient.
        states[new.i] <- switch(inherit.type,
                                state * c(1,1),
                                {new1<-sample(0:1, 1); c(new1, 1-new1)},
                                rep(1-state, 2))
        
        if ( inherit.type == 1 )
          n.i[state.i] <- n.i[state.i] + 1
        else if ( inherit.type == 2 )
          n.i[(1-state)+1] <- n.i[(1-state)+1] + 1
        else if ( inherit.type == 3 ) {
          n.i[(1-state)+1] <- n.i[(1-state)+1] + 2
          n.i[state.i] <- n.i[state.i] - 1
        }
      }
      
      parent[new.i] <- lineage
      start[new.i] <- t
      len[new.i] <- 0
      n.taxa <- n.taxa + 1
      lineages <- which(!split & !extinct)
    } else if ( type == 2 ) {
      extinct[lineage] <- TRUE
      lineages <- which(!split & !extinct)
      n.i[state.i] <- n.i[state.i] - 1
      n.taxa <- n.taxa - 1
    } else {
      n.i <- n.i + if (state == 0) c(-1, 1) else c(1, -1)
      states[lineage] <- 1 - state
      hist[[length(hist) + 1]] <- c(lineage, t, state, 1 - state)
    }
  }

  info <- data.frame(idx=seq_along(extinct), len=len, parent=parent, 
                     start=start, state=states, extinct=extinct,
                     split=split)
  hist <- as.data.frame(do.call(rbind, hist))
  if (nrow(hist) == 0)
    hist <- as.data.frame(matrix(NA, 0, 4))
  names(hist) <- c("idx", "t", "from", "to")
  hist$x0 <- info$start[match(hist$idx, info$idx)]
  hist$tc <- hist$t - hist$x0

  attr(info, "t") <- t
  attr(info, "hist") <- hist
  info
}
