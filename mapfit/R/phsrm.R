## PH package for R

# setAs("cf1", "cf1srm", function(from, to){
#   new("cf1srm", size=from@size, alpha=from@alpha, Q=from@Q, xi=from@xi, rate=from@rate, omega=c(0))
# })

# setMethod("phfit.init", signature(model = "cf1srm"),
#   function(model, data, ...) {
#     cf1srm <- as(cf1.param(model@size, mean(data), ...), "cf1srm")
#     cf1srm@omega <- sum(data@data$counts)
#     cf1srm
#   }
# )

# setMethod("phfit.estep", signature(model = "cf1srm", data = "phdata.group"),
#   function(model, data, ufact = 1.01, eps = 1.0e-8, ...) {
#   nnz <- length(model@Q@x)
#   p <- model@Q@p + 1L
#   j <- model@Q@j + 1L
#   qtrans <- .Fortran("em_sparse_phase_qtrans", PACKAGE="mapfit",
#     n=as.integer(model@size),
#     nnz=as.integer(nnz),
#     Q=as.double(model@Q@x),
#     rowptr=as.integer(p),
#     colind=as.integer(j),
#     P=as.double(numeric(nnz)),
#     qv=as.double(numeric(1)),
#     ufactor=as.double(ufact))

#   ba <- -solve(t(as.matrix(model@Q)), model@alpha)
#   data@data$instant[is.na(data@data$counts)] <- 0
#   data@data$counts[is.na(data@data$counts)] <- -1
#   l <- data@size
#   if (is.infinite(data@data$time[l])) {
#     gdatlast <- data@data$counts[l]
#     data@data <- data@data[-l,]
#     data@size <- data@size - 1
#   } else {
#     gdatlast <- 0
#   }

#   res <- .Fortran("em_sparse_phase_estep_trunc_group_poi", PACKAGE="mapfit",
#     n=as.integer(model@size),
#     alpha=as.double(model@alpha),
#     baralpha=as.double(ba),
#     xi=as.double(model@xi),
#     vone=as.double(rep(1,model@size)),
#     nnz=as.integer(nnz),
#     Q=as.double(model@Q@x),
#     P=as.double(qtrans$P),
#     rowptr=as.integer(p),
#     colind=as.integer(j),
#     qv=as.double(qtrans$qv),
#     omega=as.double(model@omega),
#     eps=as.double(eps),
#     m=as.integer(data@size),
#     tdat=as.double(data@data$time),
#     gdat=as.integer(data@data$counts),
#     idat=as.integer(data@data$instant),
#     gdatlast=as.integer(gdatlast),
#     eb=as.double(numeric(model@size)),
#     en=as.double(numeric(nnz)),
#     ey=as.double(numeric(model@size)),
#     llf=as.double(numeric(1)))
#   list(eres=list(eb=res$eb, ey=res$ey, en=res$en, eomega=res$omega), llf=res$llf)
#   })

# #### mstep

# setMethod("phfit.mstep", signature(model = "cf1srm"),
#   function(model, eres, data, ...) {
#   nnz <- length(model@Q@x)
#   p <- model@Q@p + 1L
#   j <- model@Q@j + 1L
#   res <- .Fortran("em_sparse_phase_mstep", PACKAGE="mapfit",
#     n=as.integer(model@size),
#     alpha=as.double(model@alpha),
#     xi=as.double(model@xi),
#     nnz=as.integer(nnz),
#     Q=as.double(model@Q@x),
#     rowptr=as.integer(p),
#     colind=as.integer(j),
#     eb=as.double(eres$eb),
#     en=as.double(eres$en),
#     ey=as.double(eres$ey))
#   rate <- -res$Q[res$Q < 0]
#   res2 <- .Fortran("em_cf_phase_sort", PACKAGE="mapfit",
#     n=as.integer(model@size),
#     alpha=as.double(res$alpha),
#     rate=as.double(rate))
#   cf1srm <- as(cf1(alpha=res2$alpha, rate=res2$rate), "cf1srm")
#   cf1srm@omega <- eres$eomega
#   cf1srm
#   })

# ########
