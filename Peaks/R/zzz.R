.First.lib <- function(lib, pkg) {
       library.dynam("Peaks", pkg, lib)
}

SpectrumBackground <- function(y,iterations=100,decreasing=FALSE,order=c("2","4","6","8"),smoothing=FALSE,window=c("3","5","7","9","11","13","15"),compton=FALSE){
  p <- .Call("R_SpectrumBackground",
     as.vector(y),
     as.integer(iterations),
     as.integer(decreasing),
     as.integer(as.integer(match.arg(order))/2-1),
     as.integer(smoothing),
     as.integer(as.integer(match.arg(window))),
     as.integer(compton),
     PACKAGE="Peaks"
     )
  return(p)
}

SpectrumSmoothMarkov <- function(y,window=3){
  p <- .Call("R_SpectrumSmoothMarkov",
     as.vector(y),
     as.integer(window),
     PACKAGE="Peaks"
     )
  return(p)
}
  
SpectrumDeconvolution <- function(y,response,iterations=10,repetitions=1,boost=1.0,method=c("Gold","RL")){
  method <- match.arg(method)
  if (length(as.vector(response))<length(as.vector(y))){
    response <- c(response,rep(0,length(y)-length(response)))
  }
  if (length(as.vector(response))>length(as.vector(y))){
    stop("response length should be shorter or equal y length")
  }
  switch(method,
         Gold={
           p1 <- "R_SpectrumDeconvolution"
         },
         RL={
           p1 <- "R_SpectrumDeconvolutionRL"
         })
  p <- .Call(p1,
             as.vector(y),
             as.vector(response),
             as.integer(iterations),
             as.integer(repetitions),
             as.numeric(boost),
             PACKAGE="Peaks")

  return(p)
}

SpectrumSearch <-  function(y,sigma=3.0,threshold=10.0,background=FALSE,iterations=13,markov=FALSE,window=3){
  p <- .Call("R_SpectrumSearchHighRes",
     as.vector(y),
     as.numeric(sigma),
     as.numeric(threshold),
     as.integer(background),
     as.integer(iterations),
     as.integer(markov),
     as.integer(window),
     PACKAGE="Peaks"
     )
  return(p)
}
