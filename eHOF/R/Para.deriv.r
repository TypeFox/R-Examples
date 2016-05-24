Para_deriv <- function (
		resp, # response data object of HOF function
		newdata=NULL, # specify new x value vector
		model=NULL, 
		p, 
		type=c('slope','inflection'), 
		...)
{
  type <- match.arg(type, c('slope','inflection'))
  if(is.null(newdata)) x <- scale01(seq(resp$range[1] - diff(resp$range), resp$range[2] + diff(resp$range), length.out=10000)) else x <- scale01(newdata, ...)  
  if(is.null(model)) model <- pick.model(resp, gam=FALSE, ...)
  if(missing(p)) p <- resp$models[[model]]$par
  a = p[1]; b = p[2]; c = p[3]; d = p[4]; f = p[5]
  M <- resp$M
  
  fv <- switch(as.character(model),
     'I'  = expression(M/(1 + exp(a))),
     'II' = expression(M/(1 + exp(a + b * x))),
     'III'= expression(M/(1 + exp(a + b * x))/(1 + exp(c))),
     'IV' = expression(M/(1 + exp(a + b * x))/(1 + exp(c - b * x))), 
     'V'  = expression(M/(1 + exp(a + b * x))/(1 + exp(c - d * x))),
     'VI' = expression(M/(1 + exp(a + b * x)) * 1/(1 + exp(c - b * x)) + M/(1 + exp(a + b * (x-d))) * 1/(1 + exp(c - b * (x-d)))),
    'VII' = expression(M/(1 + exp(a + b * x)) * 1/(1 + exp(c - b * x)) + M/(1 + exp(a + b * (x-d))) * 1/(1 + exp(c - f * (x-d))))
    )

  if(type == 'slope') {
    if(model=='I') out <- rep(0, length(x)) else {
    	fv1 <- D(fv, 'x')
    	out <- eval(fv1)
    }
  }

  if(type == 'inflection') {
    optima <- suppressMessages(Para_opt(resp, model)$opt)
    optima <- scale01(optima, resp$range)
    pessima <- suppressMessages(Para_opt(resp, model)$pess)
    pessima <- scale01(pessima, resp$range)
    fv1 <- D(fv,'x')
    fv2 <- deriv(fv1, 'x')
    fv3 <- as.function(list(fv2[[1]]))
    formals(fv3) <- alist(x=,a=a,b=b,c=c,d=d,f=f)
    out <- switch(as.character(model),
      I  = NA,
      II = if(a>0) optimize(fv3, c(0,1), a,b,c,d, maximum=TRUE)$maximum else optimize(fv3, c(0,1), a,b,c,d, maximum=FALSE)$minimum,
      III= {# range <- if(a<0) c(max(optima[2]), min(pessima)) else c(max(pessima),min(optima))
            if(a<0) optimize(fv3, c(0,1), a,b,c,d, maximum=FALSE)$minimum else optimize(fv3, c(0,1), a,b,c,d, maximum=TRUE)$maximum
            },
      IV = c(optimize(fv3, c(0, optima), a,b,c,d, maximum=TRUE)$maximum, optimize(fv3, c(optima, 1), a,b,c,d, maximum=FALSE)$minimum),
      V  = c(optimize(fv3, c(0, optima), a,b,c,d, maximum=TRUE)$maximum, optimize(fv3, c(optima, 1), a,b,c,d, maximum=FALSE)$minimum),
      VI = c(optimize(fv3, c(0, optima[1]), a,b,c,d, maximum=TRUE)$maximum, 
             optimize(fv3, c(optima[1], pessima), a,b,c,d,maximum=FALSE)$minimum, 
             optimize(fv3, c(pessima, optima[2]), a,b,c,d, maximum=TRUE)$maximum, 
             optimize(fv3, c(optima[2], max(x)), a,b,c,d, maximum=FALSE)$minimum),

      VII = c(optimize(fv3, c(0, optima[1]), a,b,c,d,f, maximum=TRUE)$maximum, 
              optimize(fv3, c(optima[1], pessima), a,b,c,d,f, maximum=FALSE)$minimum,
              optimize(fv3, c(pessima, optima[2]), a,b,c,d,f, maximum=TRUE)$maximum, 
              optimize(fv3, c(optima[2], max(x)), a,b,c,d,f, maximum=FALSE)$minimum)
   )
   out <- rescale01(out, resp$range)
## Filter true inflection points
   if(model == 'II') {
      if(a < 0 & out < (resp$range[1] + diff(resp$range)/100)) out <- NA
      if(a > 0 & out > (resp$range[2] - diff(resp$range)/100)) out <- NA
   }
   if(model == 'IV') {
      if(out[1] < (resp$range[1] + diff(resp$range)/100)) out <- NA
      if(length(out) > 1) if(out[2] > (resp$range[2] - diff(resp$range)/100) ) out[2] <- NA
   }
   if(model == 'V') {
      if(out[1] < (resp$range[1] + diff(resp$range)/100)) out[1] <- NA
      if(length(out)>1) if(out[2] > (resp$range[2] - diff(resp$range)/100)) out[2] <- NA
   }
}

return(out) 

}

