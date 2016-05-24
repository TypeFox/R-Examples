lambda <- function(tab, rp=rep(1/nrow(tab), nrow(tab)), cp=rep(1/ncol(tab), ncol(tab))) {
  logp <- log(prop.table(tab))

  lambda <- sweep(logp, 1, rowSums(sweep(logp, 2, cp, "*")), "-")
  lambda <- sweep(lambda, 2, colSums(sweep(logp, 1, rp, "*")), "-")
  lambda + sum(logp * rp %o% cp)
}

maor <- function(tab, phi=FALSE, cell=FALSE,
                 weighting=c("marginal", "uniform", "none"), norm=2,
                 component=c("total", "symmetric", "antisymmetric"),
                 row.weights=NULL, col.weights=NULL) {
  weighting <- match.arg(weighting)
  component <- match.arg(component)

  if(!length(dim(tab)) %in% 2:3) {
      stop("Only two- and three-way tables are supported.")
  }

  if(!norm %in% 1:2)
     stop("'norm' must be 1 or 2")

  if(norm != 2 && !phi)
      stop("Only norm=2 is currently supported when phi=FALSE")

  if(!(is.null(row.weights) || !is.null(col.weights)) && !missing(weighting))
      warning("Argument 'weighting' is ignored when custom row/column weights are specified.")

  if(any(is.na(tab))) {
      warning("NA cells are currently not supported, returning NA.")
      return(NA)
  }

  if(length(dim(tab)) == 3) {
      rp <- margin.table(tab, 1)
      cp <- margin.table(tab, 2)

      if(!is.null(row.weights))
          rp <- row.weights

      if(!is.null(col.weights))
          cp <- col.weights

      if(any(tab == 0)) {
          tab <- tab + 0.5
          warning("Cells with zero counts found: adding 0.5 to each cell of the table.")
      }

      if(weighting == "marginal")
          return(apply(tab, 3, maor,
                       phi=phi, cell=cell, norm=norm,
                       component=component,
                       row.weights=rp, col.weights=cp))
      else
          return(apply(tab, 3, maor,
                       phi=phi, cell=cell,
                       weighting=weighting, norm=norm,
                       component=component,
                       row.weights=row.weights, col.weights=col.weights))
  }

  if(any(tab == 0)) {
      if(all(tab == 0)) {
          warning("Table contains only empty cells, returning NA.")
          return(NA)
      }
      else if(sum(tab == 0)/length(tab) > .25) {
          warning("More than 25% of cells are empty, the value of the index may not be reliable.")
      }

      tab <- tab + 0.5
      warning("Cells with zero counts found: adding 0.5 to each cell of the table.")
  }

  if(weighting == "marginal") {
      p <- prop.table(tab)
      rp <- margin.table(p, 1)
      cp <- margin.table(p, 2)
  }
  else if(weighting == "uniform") {
      rp <- rep(1/nrow(tab), nrow(tab))
      cp <- rep(1/ncol(tab), ncol(tab))
  }
  else {
      rp <- rep(1, nrow(tab))
      cp <- rep(1, ncol(tab))
  }

  if(!is.null(row.weights))
      rp <- prop.table(row.weights)

  if(!is.null(col.weights))
      cp <- prop.table(col.weights)

  if(component %in% c("symmetric", "antisymmetric"))
      rp <- cp <- (rp + cp)/2

  rp1 <- prop.table(rp)
  cp1 <- prop.table(cp)

  l <- lambda(tab, rp1, cp1)

  if(component == "symmetric")
      l <- (l + t(l))/2
  else if(component == "antisymmetric")
      l <- (l - t(l))/2

  lambda.norm <- abs(l^norm * rp %o% cp)

  if(phi) {
      if(cell)
          lambda.norm
      else
          sum(lambda.norm)^(1/norm)
  }
  else {
      if(weighting == "none") {
          if(cell)
              4 * nrow(tab) * ncol(tab) * lambda.norm
          else
              exp((4 * nrow(tab) * ncol(tab) * sum(lambda.norm))^(1/norm))
      }
      else {
          if(cell)
              4/sum((rp1 * (1 - rp1)) %o% (cp1 * (1 - cp1))) * lambda.norm
          else
              exp((4/sum((rp1 * (1 - rp1)) %o% (cp1 * (1 - cp1))) * sum(lambda.norm))^(1/norm))
      }
  }
}
