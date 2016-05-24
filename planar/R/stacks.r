##' invert the description of a multilayer to simulate the opposite direction of incidence
##'
##' inverts list of epsilon and thickness of layers
##' @title invert_stack
##' @param p list
##' @return list
##' @export
##' @family helping_functions
##' @author Baptiste Auguie
invert_stack <- function(p){
  p[["epsilon"]] <- rev(p[["epsilon"]])
  p[["thickness"]] <- rev(p[["thickness"]])
  p
}


check_stack <- function(s){
  inherits(s, "stack") && length(s[["thickness"]] == length(s[["epsilon"]]))
}

##' invert the description of a multilayer to simulate the opposite direction of incidence
##'
##' inverts list of epsilon and thickness of layers
##' @title rev.stack
##' @param x stack
##' @return stack
##' @family helping_functions user_level stack
##' @export
##' @author Baptiste Auguie
rev.stack <- function(x) {
  x[["epsilon"]] <- rev(x[["epsilon"]])
  x[["thickness"]] <- rev(x[["thickness"]])
  x
}


##' @export
c.stack <- function(..., recursive = FALSE){
  sl <- list(...)
  
  tl <- lapply(sl, "[[", "thickness")
  el <- lapply(sl, "[[", "epsilon")
  
  ll <- list(epsilon=do.call(c, el), 
             thickness=do.call(c, tl))
  structure(ll, class="stack")
}

##' @importFrom utils str
##' @export
print.stack <- function(x, ...){
  str(x)
}

##' @importFrom graphics plot rect
##' @export
plot.stack <- function(x, ...){
  xx <- c(0, cumsum(x[['thickness']]))
  material <- epsilon_label(x[['epsilon']])
  plot(xx, 0*xx, t="n", ylim=c(0,1), yaxt="n", ylab="", bty="n",
       xlab=expression(x/nm), ...)
  rect(xx[-length(xx)], 0, xx[-1], 1, col=material)
}

##' @importFrom ggplot2 fortify
##' @export
fortify.stack <- function(model, data, ...){
  
  xx <- c(0, cumsum(model[['thickness']]))
  if(!missing(data)) xx <- xx - data # shift origin
  material <- epsilon_label(model[['epsilon']])
  N <- length(xx)
  data.frame(xmin=xx[-N],
             xmax=xx[-1],
             material = material)
  
}

# silly visible binding 
xmin <- xmax <- material <- NULL

##' @importFrom ggplot2 autoplot aes fortify expand_limits theme scale_x_continuous scale_y_continuous ggplot geom_rect
##' @export
autoplot.stack <- function(object, ...){
  d <- fortify(object, ...)
  ggplot(d) + 
    expand_limits(y=c(0,1))+
    geom_rect(aes(xmin=xmin, xmax=xmax, 
                  ymin=-Inf, ymax=Inf, fill=material))+
    scale_x_continuous(expression(x/nm),expand=c(0,0))+
    scale_y_continuous("", expand=c(0,0), breaks=NULL)+
    theme()
  
}

##' Single-layer stack structure
##'
##' returns a stack describing a single layer
##' @title layer_stack
##' @export
##' @param epsilon dielectric function (numeric, character, or complex)
##' @param thickness layer thickness in nm
##' @param ... ignored
##' @return list of class 'stack'
##' @author baptiste Auguie
##' @family stack user_level
layer_stack <- function(epsilon="epsAu", thickness=50, ...){
  ll <- list(epsilon=list(epsilon), thickness=thickness)
  structure(ll, class="stack")
}

##' Embed stack structure
##'
##' embeds a stack in semi-infinite media
##' @title embed_stack
##' @export
##' @param s stack (finite structure)
##' @param nleft real refractive index on the left side
##' @param nright real refractive index on the right side
##' @param dleft dummy layer thickness in nm
##' @param dright dummy layer thickness in nm
##' @param ... ignored
##' @return list of class 'stack'
##' @author baptiste Auguie
##' @family stack user_level
embed_stack <- function(s, nleft=1.0, nright=1.0,
                        dleft= 200, dright=200, ...){
  
  Nlay <- length(s[["thickness"]])
  if(s[["thickness"]][1] == 0L && s[["thickness"]][Nlay] == 0L)
    warning("already embedded?")
  
  s[["thickness"]] <- c(0, dleft, s[["thickness"]], dright, 0)
  s[["epsilon"]] <- c(nleft^2, nleft^2, s[["epsilon"]], 
                      nright^2, nright^2)
  
  s
}

##' DBR stack structure
##'
##' periodic structure of dielectric layers
##' @title dbr_stack
##' @export
##' @param lambda0 central wavelength of the stopband
##' @param n1 real refractive index for odd layers
##' @param n2 real refractive index for even layers
##' @param d1 odd layer thickness in nm
##' @param d2 even layer thickness in nm
##' @param N number of layers, overwrites pairs
##' @param pairs number of pairs
##' @param ... ignored
##' @return list of class 'stack'
##' @author baptiste Auguie
##' @family stack user_level
dbr_stack <- function(lambda0=630, 
                      n1=1.28, n2=1.72, 
                      d1=lambda0/4/n1, d2=lambda0/4/n2,
                      N=2*pairs, pairs=4, ...){
  epsilon <- as.list(rep(c(n1^2, n2^2), length.out=N))
  thickness <- rep(c(d1,d2), length.out=N)
  
  ll <- list(epsilon=epsilon, thickness=thickness)
  structure(ll, class="stack")
}

##' DBR-metal stack structure
##'
##' periodic structure of dielectric layers against metal film
##' @title tamm_stack
##' @export
##' @param lambda0 central wavelength of the stopband
##' @param n1 real refractive index for odd layers
##' @param n2 real refractive index for even layers
##' @param d1 odd layer thickness in nm
##' @param d2 even layer thickness in nm
##' @param N number of layers, overwrites pairs
##' @param pairs number of pairs
##' @param dx1 variation of last odd layer thickness in nm
##' @param dx2 variation of last even layer thickness in nm
##' @param dm thickness of metal layer
##' @param metal character name of dielectric function
##' @param nleft refractive index of entering medium
##' @param nright refractive index of outer medium
##' @param dleft distance from the left side for visualisation
##' @param dright distance from the right side for visualisation
##' @param position metal position relative to DBR
##' @param incidence direction of incidence
##' @param ... ignored
##' @return list of class 'stack'
##' @author baptiste Auguie
##' @family stack user_level
tamm_stack <- function(lambda0=630, 
                       n1=1.28, n2=1.72, 
                       d1=lambda0/4/n1, d2=lambda0/4/n2,
                       N=2*pairs, pairs=4, 
                       dx1 = 0, dx2 = 0,
                       dm=50, metal="epsAu",
                       position=c("after", "before"),
                       incidence = c("left", "right"),
                       nleft = 1.5, nright=1.0,
                       dleft= 200, dright=200,
                       ...){
  position <- match.arg(position)
  incidence <- match.arg(incidence)
  
  dbr <- dbr_stack(lambda0=lambda0, 
                   n1=n1, n2=n2, 
                   d1=d1, d2=d2,
                   N=N-2, pairs=pairs-1)
  
  variable <- dbr_stack(lambda0=lambda0, 
                        n1=n1, n2=n2, 
                        d1=d1 + dx1, d2=d2 + dx2,
                        N=2)
  
  met <- layer_stack(epsilon=metal, thickness=dm)
  
  ## case odd number of layers
  ## N-2 is odd, so end with same as start
  ## fine if variable is before DBR
  ## wrong if after, needs to be flipped
  if(N %% 2 && position == "after") variable <- rev(variable)
  
  struct <- switch(position,
                   before = c(met, variable, dbr),
                   after = c(dbr, variable, met))
  
  s <- embed_stack(struct, nleft=nleft, nright=nright, 
                   dleft= dleft, dright=dright) 
  
  switch(incidence,
         left = s,
         right = rev(s))
  
}

##' DBR-metal stack structure
##'
##' periodic structure of dielectric layers against metal film
##' @title tamm_stack_ir
##' @export
##' @param lambda0 central wavelength of the stopband
##' @param n1 real refractive index for odd layers
##' @param n2 real refractive index for even layers
##' @param d1 odd layer thickness in nm
##' @param d2 even layer thickness in nm
##' @param N number of layers, overwrites pairs
##' @param pairs number of pairs
##' @param dx1 variation of last odd layer thickness in nm
##' @param dx2 variation of last even layer thickness in nm
##' @param dm thickness of metal layer
##' @param metal character name of dielectric function
##' @param nleft refractive index of entering medium
##' @param nright refractive index of outer medium
##' @param position metal position relative to DBR
##' @param incidence direction of incidence
##' @param ... ignored
##' @return list of class 'stack'
##' @author baptiste Auguie
##' @family stack user_level
tamm_stack_ir <- function(lambda0=950, 
                          n1=3, n2=3.7, 
                          d1=lambda0/4/n1, d2=lambda0/4/n2,
                          N=2*pairs, pairs=4, 
                          dx1 = 0, dx2 = 0,
                          dm=50, metal="epsAu",
                          position="after",
                          incidence = "left",
                          nleft = n2, nright=1.0,
                          ...){
  tamm_stack(lambda0=lambda0, 
             n1=n1, n2=n2,
             d1=d1, d2=d2, 
             N=N, pairs=pairs, 
             dx1=dx1, dx2=dx2, 
             dm=dm, metal=metal,
             position=position,
             incidence=incidence,
             nleft=nleft, nright=nright,
             ...)
}

##' DBR-metal stack structure
##'
##' periodic structure of dielectric layers against metal film
##' @title tamm_stack_porous
##' @export
##' @param lambda0 central wavelength of the stopband
##' @param n1 real refractive index for odd layers
##' @param n2 real refractive index for even layers
##' @param d1 odd layer thickness in nm
##' @param d2 even layer thickness in nm
##' @param N number of layers, overwrites pairs
##' @param pairs number of pairs
##' @param dx1 variation of last odd layer thickness in nm
##' @param dx2 variation of last even layer thickness in nm
##' @param dm thickness of metal layer
##' @param metal character name of dielectric function
##' @param nleft refractive index of entering medium
##' @param nright refractive index of outer medium
##' @param position metal position relative to DBR
##' @param incidence direction of incidence
##' @param ... ignored
##' @return list of class 'stack'
##' @author baptiste Auguie
##' @family stack user_level
tamm_stack_porous <- function(lambda0=600, 
                          n1=1.72, n2=1.28, 
                          d1=lambda0/4/n1, d2=lambda0/4/n2,
                          N=2*pairs, pairs=4, 
                          dx1 = 0, dx2 = 0,
                          dm=20, metal="epsAu",
                          position="before",
                          incidence = "right",
                          nleft = 1.5, nright=1.0,
                          ...){
  tamm_stack(lambda0=lambda0, 
             n1=n1, n2=n2,
             d1=d1, d2=d2, 
             N=N, pairs=pairs, 
             dx1=dx1, dx2=dx2, 
             dm=dm, metal=metal,
             position=position,
             incidence=incidence,
             nleft=nleft, nright=nright,
             ...)
}

