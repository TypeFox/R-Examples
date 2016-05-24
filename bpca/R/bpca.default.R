bpca.default <- function(x,
                         d=1:2,
                         center=2,
                         scale=TRUE,
                         method=c('hj', 'sqrt', 'jk', 'gh'),
                         iec=FALSE,
                         var.rb=FALSE,
                         var.rd=FALSE,
                         limit=10, ...)
{
  stopifnot(is.matrix(x) || is.data.frame(x))

  li <- d[1]
  le <- d[length(d)]
  n.lambda <- (le - li + 1)

  if(n.lambda < 2 || n.lambda > 3)
    stop('Please, check the parameter d:\n',
         'The (d[1] - d[length(d)] + 1) must equal to 2 (for bpca.2d) or 3 (for bpca.3d).\n\n')

  x <- as.matrix(x)

  x.cent <- x                                             # 0: no centering
  switch(center,                                          # of course, if center=1:3
         x.cent <- sweep(x, 1, mean(x)),                  # 1: global-centered
         x.cent <- sweep(x, 2, apply(x, 2, mean)),        # 2: column-centered
         x.cent <- sweep(sweep(x, 1, apply(x, 1, mean)),  # 3: double-centered
                         2, apply(x, 2, mean)) + mean(x))

  if(scale)
    x.scal <- sweep(x.cent, 2, apply(x.cent, 2, sd), '/')
  else
    x.scal <- x.cent

  svdx.scal <- svd(x.scal)

  rownames(svdx.scal$u) <- rownames(x) # obj
  rownames(svdx.scal$v) <- colnames(x) # var
  colnames(svdx.scal$v) <- paste('PC',
                                 1:length(svdx.scal$d),
                                 sep='')

  s2.scal <- diag(svdx.scal$d)

  switch(match.arg(method),
         hj = {
           g.scal  <- svdx.scal$u %*% s2.scal
           h.scal  <- s2.scal %*% t(svdx.scal$v)
           hl.scal <- t(h.scal)
         },
         sqrt = {
           g.scal  <- svdx.scal$u %*% sqrt(s2.scal)
           h.scal  <- sqrt(s2.scal) %*% t(svdx.scal$v)
           hl.scal <- t(h.scal)
         },
         jk = {
           g.scal  <- svdx.scal$u %*% s2.scal
           hl.scal <- svdx.scal$v

           ## The below is the GGEBiplot aproach
           #d1 <- (max(hl.scal[,li]) - min(hl.scal[,li])) /
           #(max(g.scal[,li]) - min(g.scal[,li]))
           #d2 <- (max(hl.scal[,le]) - min(hl.scal[,le])) /
           #(max(g.scal[,le]) - min(g.scal[,le]))
           #d <- max(d1,d2)

           #hl.scal / d
         },
         gh = {
           g.scal  <- sqrt(nrow(x)-1) * svdx.scal$u
           h.scal  <- 1/sqrt(nrow(x)-1) * s2.scal %*% t(svdx.scal$v)
           hl.scal <- t(h.scal)

           ## The below is the GGEBiplot aproach
           #g.scal  <- svdx.scal$u
           #hl.scal  <- svdx.scal$v %*% s2.scal

           #d1 <- (max(g.scal[,li]) - min(g.scal[,li])) / 
           #(max(hl.scal[,li]) - min(hl.scal[,li]))
           #d2 <- (max(g.scal[,le]) - min(g.scal[,le])) / 
           #(max(hl.scal[,le]) - min(hl.scal[,le]))
           #d <- max(d1,d2)
           #g.scal <- g.scal / d
         })

  pc.names <- paste('PC',
                    1:ncol(hl.scal),
                    sep='')

  if(is.null(rownames(x.scal)))
    rownames(g.scal) <- 1:nrow(x.scal)
  else
    rownames(g.scal) <- rownames(x.scal)
  colnames(g.scal) <- pc.names
  if(is.null(colnames(x.scal)))
    rownames(hl.scal) <- paste('V',
                               1:ncol(x),
                               sep='')
  else
    rownames(hl.scal) <-colnames(x.scal)
  colnames(hl.scal) <- pc.names

  # variables
  if(var.rb)
    var.rb.res <- var.rbf(hl.scal[,d[1]:d[length(d)]])
  else
    var.rb.res <- NA

  if(var.rb & var.rd)
    var.rd.res <- var.rdf(x.scal, var.rb.res, limit)
  else
    var.rd.res <- NA

  if(iec){
    svdx.scal$v <- (-1) * svdx.scal$v
    g.scal      <- (-1) * g.scal     
    hl.scal     <- (-1) * hl.scal    
  }   

  res <- list(call=match.call(),
              eigenvalues=svdx.scal$d,
              eigenvectors=svdx.scal$v,
              number=seq(li, le, 1),
              importance=rbind(general=round(sum(svdx.scal$d[li:le]^2) /
                                             sum(svdx.scal$d^2), 3),
                               partial=round(sum(svdx.scal$d[li:le]^2) /
                                             sum(svdx.scal$d[li:length(svdx.scal$d)]^2), 3)),
              coord=list(objects=g.scal,
                         variables=hl.scal),
              var.rb=var.rb.res,
              var.rd=var.rd.res)

  colnames(res$importance) <- 'explained'

  if(n.lambda == 2)
    class(res) <- c('bpca.2d', 'bpca', 'list')
  else if(n.lambda == 3)
    class(res) <- c('bpca.3d', 'bpca', 'list')

  invisible(res)
}
