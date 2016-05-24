## sigma estimation for RBF kernels
## author: alexandros

setGeneric("sigest", function(x, ...) standardGeneric("sigest"))
setMethod("sigest",signature(x="formula"),
function (x, data=NULL, frac = 0.5, na.action = na.omit, scaled = TRUE){
  call <- match.call()
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
 ## m$... <- NULL
  m$formula <- m$x
  m$x <- NULL
  m$scaled <- NULL
  m$frac <- NULL
  m[[1L]] <- quote(stats::model.frame)
  m <- eval(m, parent.frame())
  Terms <- attr(m, "terms")
  attr(Terms, "intercept") <- 0
  x <- model.matrix(Terms, m)
  if (length(scaled) == 1)
    scaled <- rep(scaled, ncol(x))
  if (any(scaled)) {
    remove <- unique(c(which(labels(Terms) %in% names(attr(x, "contrasts"))),
                       which(!scaled)
                       )
                     )
    scaled <- !attr(x, "assign") %in% remove
  }
   ret <- sigest(x, scaled = scaled, frac = frac, na.action = na.action)
   return (ret)
})
setMethod("sigest",signature(x="matrix"),
function (x,
          frac = 0.5,
          scaled    = TRUE,
          na.action = na.omit)
          {
            x <- na.action(x)
            
            if (length(scaled) == 1)
              scaled <- rep(scaled, ncol(x))
            if (any(scaled)) {
              co <- !apply(x[,scaled, drop = FALSE], 2, var)
              if (any(co)) {
                scaled <- rep(FALSE, ncol(x))
                warning(paste("Variable(s)",
                              paste("`",colnames(x[,scaled, drop = FALSE])[co],
                                    "'", sep="", collapse=" and "),
                              "constant. Cannot scale data.")
                        )
              } else {
                xtmp <- scale(x[,scaled])
                x[,scaled] <- xtmp
              }
            }
  
            m <- dim(x)[1]
            n <- floor(frac*m)
            index <- sample(1:m, n, replace = TRUE)
            index2 <- sample(1:m, n, replace = TRUE)
            temp <- x[index,, drop=FALSE] - x[index2,,drop=FALSE]
            dist <- rowSums(temp^2)
            srange <- 1/quantile(dist[dist!=0],probs=c(0.9,0.5,0.1))

            ## ds <- sort(dist[dist!=0])
            ## sl <- ds[ceiling(0.2*length(ds))]
            ## su <- ds[ceiling(0.8*length(ds))]
            ## srange <- c(1/su,1/median(ds), 1/sl)
            ##   names(srange) <- NULL

            return(srange)
          })
