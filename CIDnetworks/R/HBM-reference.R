
#library(Rcpp); library(mvtnorm); library(msm); sourceCpp ("../src/cid.cpp"); source("CID-basefunctions.R");

# Hierarchical Block Model: Reference Class
# This adapts and extends the "Hierarchical Random Graph" model from Clauset, Moore and Newman 2008 in Nature.


closest.ancestor.from.parents <- function (parents.vector) {
  #parents.vector=c(0,1,1,2,2)

  nodes <- length(parents.vector)
  history <- function(kk) {
    out <- kk
    while (parents.vector[kk] > 0) {kk <- parents.vector[kk]; out <- c(kk,out)}
    return(out)
  }
  history.all <- lapply(1:length(parents.vector), history)
  length.set <- sapply(history.all, length)

  out.table <- array(0, c(nodes,nodes))
  for (ii in 1:nodes)
    for (jj in ii:nodes) {
      ll <- min(length.set[ii], length.set[jj])
      pick <- history.all[[ii]][max(which(history.all[[ii]][1:ll] ==
                                          history.all[[jj]][1:ll]))]
      out.table[ii,jj] <- out.table[jj,ii] <- pick
    }

  return(out.table)
}



HBMcid <-
  setRefClass(
    "HBMcid",
    fields = list(
      n.groups="numeric",

      block.value="numeric",
      block.value.m="numeric",
      block.value.v="numeric",

      membership="numeric",
      tree.parent="numeric",
      common.anc="matrix",

      #membership.edge="matrix",   #looks like edge list: has the membership number of each participant.
      #membership.node="matrix",   #each column is the Dirichlet distribution for each membership.
      #membership.alpha0="numeric",  #prior strength for the Dirichlet -- alpha0*vector(1)

      shift="numeric",
      restrict.and.shift="logical",
      group.pairs="matrix",

      #single.membership="logical",

      #inherited from main.
      node.names="character",
      n.nodes="numeric",
      outcome="numeric",
      edge.list="matrix",
      residual.variance="numeric",
      edge.list.rows="list"    #,
      ),

    methods=list(

      initialize = function (

        n.groups=2,

        n.nodes=10,
        edge.list=make.edge.list(n.nodes),
        edge.list.rows=row.list.maker(edge.list),
        residual.variance=1,
        outcome=numeric(0),

        block.value=rep(0, n.groups), #*(n.groups+1)/2),
        block.value.m=rep(0, n.groups), #*(n.groups+1)/2),
        block.value.v=rep(10000, n.groups), #*(n.groups+1)/2),

        membership=sample(n.groups, n.nodes, replace=TRUE),
        tree.parent=c(0, sapply(1:(n.groups-1), function(gg) sample(gg, 1))),
        #membership.alpha0=0.1,
        #membership.node=rdirichlet.block (matrix(membership.alpha0, nrow=n.groups, ncol=n.nodes)),
        #membership.edge=draw.MMSB.from.nodes(edge.list, membership.node),

        shift=0,

        restrict.and.shift=FALSE,
        generate=FALSE

        ) {

        .self$n.nodes <<- n.nodes
        .self$edge.list <<- edge.list
        .self$edge.list.rows <<- edge.list.rows
        .self$node.names <<- as.character(1:.self$n.nodes)

        .self$n.groups <<- n.groups

        .self$block.value <<- block.value
        .self$block.value.m <<- block.value.m
        .self$block.value.v <<- block.value.v

        if (n.groups != length(block.value)) stop(paste("block.value",paste(block.value,collapse=","),"does not have length specified by n.groups,",n.groups))

        .self$membership <<- membership
        .self$tree.parent <<- tree.parent

        if (n.groups != length(tree.parent)) stop(paste("tree.parent",paste(tree.parent,collapse=","),"does not have length specified by n.groups,",n.groups))

        .self$common.anc <<- closest.ancestor.from.parents(tree.parent)

        .self$residual.variance <<- residual.variance
        .self$restrict.and.shift <<- restrict.and.shift
        #.self$single.membership <<- FALSE

        .self$group.pairs <<- makeEdgeListSelfies(n.groups)

        .self$shift <<- shift
        if (generate) .self$generate() else .self$outcome <<- outcome
      },
      center.me = function () if (restrict.and.shift) {
        shift <<- mean(block.value)
        block.value <<- block.value - shift
      },

      reinitialize = function (n.nodes=NULL,
        edge.list=NULL, node.names=NULL) {
        if (!is.null(n.nodes)) n.nodes <<- n.nodes
        if (!is.null(edge.list)) {
          edge.list <<- edge.list
          edge.list.rows <<- row.list.maker(edge.list)
        }


        if (n.groups > n.nodes) {
          warning ("HBM: Resetting number of groups to one less than the number of nodes.")
          n.groups <<- n.nodes - 1
          membership <<- sample(n.groups, n.nodes, replace=TRUE)
          block.value <<- block.value[1:n.groups]
        }


        if (length(membership) != n.nodes) {
          message ("Reinitializing HBM Membership Vector")
          membership <<- sample(n.groups, n.nodes, replace=TRUE)
        }


        if (!is.null(node.names)) {
          if (length(node.names) == .self$n.nodes) node.names <<- node.names
        } else node.names <<- as.character(1:.self$n.nodes)

      },

      pieces = function (include.name=FALSE) {
        out <- list (block.value=block.value,
                     membership=membership,
                     tree.parent=tree.parent)
        class(out) <- "HBMout"
        out
      },

      show = function () {
        message("block.value:"); print(block.value)
        message("membership:"); print(membership)
        message("tree.parent:"); print(tree.parent)
      },
      plot = function (memb=membership, tree.par=tree.parent, blockval=block.value) {
        circular.dendrogram (memb, tree.par, blockval, node.labels=node.names)
      },
      plot.network = function (color=outcome, ...) {
        image.netplot (edge.list, color, node.labels=node.names, ...)
      },



      value = function () {
        #sbm.matrix <- symBlock(block.value)
        block.value[common.anc[membership[edge.list[,1]] + n.groups*(membership[edge.list[,2]]-1)]]
      },
      value.ext = function (parameters=pieces(), edges=1:nrow(edge.list)) {   #slightly slower.
        common.anc.temp <- closest.ancestor.from.parents(parameters[[3]])
        parameters[[1]][common.anc.temp[parameters[[2]][edge.list[edges,1]] + n.groups*(parameters[[2]][edge.list[edges,2]]-1)]]

      },



      generate = function () {outcome <<- rnorm(nrow(edge.list), value(), sqrt(residual.variance))},

      log.likelihood = function(parameters=pieces(), edges=1:nrow(edge.list)) {
        meanpart <- value.ext (parameters, edges)
        sum(dnorm(outcome[edges], meanpart, sqrt(residual.variance), log=TRUE))
      },


      random.start = function () {
        tree.parent <<- c(0, sapply(1:(n.groups-1), function(gg) sample(gg, 1)))
        membership <<- sample(n.groups, n.nodes, replace=TRUE)
        block.value <<- rnorm(n.groups, 0, 0.5)
      },


      rotate = function () {
        rotation <- SBM.ID.rotation(membership, n.groups)
        membership <<- rotation[membership]
        block.value <<- block.value[rotation]

        inv.rotation <- sapply(1:length(rotation), function(rr) which(rotation==rr))
        tp.temp <- tree.parent[inv.rotation]
        tp.temp[tp.temp>0] <- rotation[tp.temp[tp.temp>0]]
        tree.parent <<- tp.temp
      },

      draw = function (verbose=0) {

        if (length(outcome) != nrow(edge.list)) stop ("HBM: outcome and edge.list have different lengths.")

        #Hold me!
        b.memb <- membership
        common.anc.temp <- closest.ancestor.from.parents (tree.parent)

        # draw node memberships in random order, so we don't favor escapes of lower-order nodes.
        for (ii in sample(1:n.nodes)) if (sum(b.memb==b.memb[ii]) + sum(tree.parent==b.memb[ii]) > 2) {  #all internal nodes must have two children; can't let it leave if it would leave 1 remaining.

          log.pp.vec <- sapply(1:n.groups, function(gg) {
            b.memb[ii] <- gg
            piece <- block.value[common.anc.temp[b.memb[edge.list[edge.list.rows[[ii]],1]] +
                                                 n.groups*(b.memb[edge.list[edge.list.rows[[ii]],2]]-1)]]

            sum(dnorm(outcome[edge.list.rows[[ii]]], piece, sqrt(residual.variance), log=TRUE))
          })
          log.pp.vec <- log.pp.vec - max(log.pp.vec)
          b.memb[ii] <- sample (1:n.groups, 1, prob=exp(log.pp.vec))
        } else if (sum(b.memb != b.memb[ii]) > 0) { #if it would cause a problem, and it's an option, find a leaf on another internal node and propose a swap.
          other.node <- sample((1:n.nodes)[b.memb != b.memb[ii]], 1)
          b.memb.temp <- b.memb; b.memb.temp[ii] <- b.memb[other.node]; b.memb.temp[other.node] <- b.memb[ii]
          rowset <- intersect(edge.list.rows[[ii]], edge.list.rows[[other.node]])

          piece0 <- block.value[common.anc.temp[b.memb[edge.list[rowset,1]] +
                                                n.groups*(b.memb[edge.list[rowset,2]]-1)]]
          log.pp.0 <- sum(dnorm(outcome[rowset], piece0, sqrt(residual.variance), log=TRUE))

          piece1 <- block.value[common.anc.temp[b.memb.temp[edge.list[rowset,1]] +
                                                n.groups*(b.memb.temp[edge.list[rowset,2]]-1)]]
          log.pp.1 <- sum(dnorm(outcome[rowset], piece1, sqrt(residual.variance), log=TRUE))

          if (log.pp.1 - log.pp.0 > -rexp(1)) b.memb <- b.memb.temp  #Standard Metropolis step.
        }
        membership <<- b.memb


        #update block value.
        edge.group <- common.anc.temp[membership[edge.list[,1]] + n.groups*(membership[edge.list[,2]]-1)]
        block.value <<- sapply(1:length(block.value), function(bb) {
          #common ancestor values
          picks <- which(edge.group == bb)
          if (length(picks) > 0) {
            var.b <- 1/(length(picks)/residual.variance + 1/block.value.v[bb])
            mean.b <- var.b*(sum(outcome[picks])/residual.variance +
                             block.value.m[bb]/block.value.v[bb])
            output <- rnorm(1, mean.b, sqrt(var.b))
          } else output <- rnorm(1, 0, 0.5)
        })

        #update tree structure. How do we do this? Find another internal node that is not its own descendant.
        b.tree <- tree.parent
        common.anc.temp <- closest.ancestor.from.parents (tree.parent)
        for (gg in sample(1:n.groups)) {
          possibles <- which(common.anc.temp[gg,] != gg)
          if (length(possibles) > 1 & sum(membership==b.tree[gg]) + sum(b.tree==b.tree[gg]) > 2 ) {
            log.pp.vec <- sapply(1:length(possibles), function(pp) {
              b.tree[gg] <- possibles[pp]
              common.anc.temp <- closest.ancestor.from.parents (b.tree)
              piece <- block.value[common.anc.temp[membership[edge.list[,1]] +
                                                   n.groups*(membership[edge.list[,2]]-1)]]

              sum(dnorm(outcome, piece, sqrt(residual.variance), log=TRUE))
            })

            log.pp.vec <- log.pp.vec - max(log.pp.vec)
            b.tree[gg] <- sample (possibles, 1, prob=exp(log.pp.vec))
            common.anc.temp <- closest.ancestor.from.parents (b.tree)
          }
        }
        tree.parent <<- b.tree
        common.anc <<- closest.ancestor.from.parents(tree.parent)

        if (restrict.and.shift) {center.me()}
        rotate()

      },


      gibbs.full = function (report.interval=0, draws=100, burnin=0, thin=1,
        make.random.start=FALSE) {
        out <- list()

        if (make.random.start) random.start()
        for (kk in 1:(draws*thin+burnin)) {
          draw();
          index <- (kk-burnin)/thin
          if (kk > burnin & round(index)==index) {
            out[[index]] <- c(pieces(), list(log.likelihood=log.likelihood()))
            if (report.interval > 0) if (index %% report.interval == 0) message("HBM ",index)
          }
        }
        return(out)
      },

      gibbs.value = function (gibbs.out) sapply(gibbs.out, function(gg) {
        value.ext (gg)
      }),

      gibbs.summary = function (gibbs.out) {
        #message("Note: HBM is not built for summaries over the posterior, only maximal draws.")
        return(gibbs.out[[length(gibbs.out)]])
      },
      print.gibbs.summary = function (gibbs.out) {
        get.sum <- gibbs.summary(gibbs.out)
        message ("Last state of the Markov Chain: Membership")
        print(get.sum$membership)
        message ("Internal tree parents:")
        print(get.sum$tree.parent)
        message ("Block values:")
        print(get.sum$block.value)
      },

      gibbs.mean = function(gibbs.out){
        get.sum <- gibbs.summary(gibbs.out)

        return(HBM(n.groups=n.groups,n.nodes=n.nodes,
                   edge.list=edge.list,
                   edge.list.rows=edge.list.rows,
                   residual.variance=residual.variance,
                   outcome=outcome,
                   block.value=get.sum$block.value,
                   block.value.m=block.value.m,
                   block.value.v=block.value.v,
                   membership=get.sum$membership,
                   tree.parent=get.sum$tree.parent,
                   shift=shift,
                   restrict.and.shift=restrict.and.shift))
      },

      gibbs.plot = function (gibbs.out) {
        #message("Note: HBM is not built for summaries over the posterior, only maximal draws.")
        #print(gibbs.out[[length(gibbs.out)]])
        get.sum <- gibbs.summary(gibbs.out)
        plot (get.sum$membership,
              get.sum$tree.parent,
              get.sum$block.value)
      },

      gibbs.node.colors = function (gibbs.out, colors=(1:n.groups) + 1) {
        get.sum <- gibbs.summary(gibbs.out)
        return(colors[get.sum$membership])
      }


      )
    )

