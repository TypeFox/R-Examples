dag.init <-
function (outcome = NULL, exposure = NULL, covs = c(), arcs = c(), 
    assocs = c(), xgap = 0.04, ygap = 0.05, len = 0.1, y.name = NULL, 
    x.name = NULL, cov.names = c(), symbols = NULL, ...) 
{ # covs: 1 for a covariable, 2 for an unknown;
  # arcs: the numbering refers to the covs vector, i.e. it
  #       differs from the later numbering in the DAG objects;
  #       exposure X is 0, outcome Y is -1;
    dag.out <- c()
  # creating first coordinates for covariates
    i1 <- 0
    cl <- length(covs)
    node.x <- c(0)
    node.y <- c(0)
    while (i1 < cl) {
        i1 <- i1 + 1
        i1.deg <- i1 * pi/(cl + 1)
        node.x <- c(node.x, 0.5 - cos(i1.deg) * 0.7)
        node.y <- c(node.y, sin(i1.deg) * 0.7)
    }
    node.x <- c(node.x, 1)
    node.y <- c(node.y, 0)
    arcs[arcs < 0] <- cl + 1
    arcs <- arcs + 1
    curve.x <- rep(NA, length(arcs)/2)
    curve.y <- rep(NA, length(arcs)/2)

  # set arc.type to default 0, and to 1 for associations
    arc.type <- rep(0, length(arcs)/2)
    i2 <- 0
    while (i2 < length(assocs)) {
        i2 <- i2 + 1
        arc.type[assocs[i2]] <- 1
    }

  # use standard names for X and Y if none provided;
    if (is.null(y.name)) 
        y.name <- "outcome"
    if (is.null(x.name)) 
        x.name <- "exposure"

  # if not enough labels for the covs, use standard
    if (length(cov.names) < cl) {
        standard.names <- c("covariable", "unknown")
        while (length(cov.names) < cl) {
            cov.names <- c(cov.names, standard.names[covs[length(cov.names) + 
                1]])
        }
    }
    node.names <- c(x.name, cov.names, y.name)

  # create custom symbols array if necessary;
    if(is.null(symbols)) symbols<-rep(NA, length(node.x));

    dag.out$cov.types <- c(0, covs, -1)
    dag.out$x <- node.x
    dag.out$y <- node.y
    dag.out$arc <- matrix(arcs, ncol = 2, byrow = TRUE)
    dag.out$arc.type <- arc.type
    dag.out$curve.x <- curve.x
    dag.out$curve.y <- curve.y
    dag.out$xgap <- xgap
    dag.out$ygap <- ygap
    dag.out$len <- len
    dag.out$names <- node.names
    dag.out$adj <- c()
    dag.out$symbols <- symbols;
    dag.out$version<-installed.packages()["dagR","Version"];
    class(dag.out)<-"dagRdag";
    return(dag.out)
}
