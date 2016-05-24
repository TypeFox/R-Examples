# Methods for classes *StepParam defined in AllClass.R 

GLMStepParam <- function(cuts){
      ncuts <- as.integer(length(cuts))
      nbands <- as.integer(length(cuts)-1)
      
  new("GLMStepParam",
      cuts=cuts,
      min=min(cuts),
      max=max(cuts),
      ncuts=ncuts,
      nbands=nbands,
      steps=cuts[2:ncuts]- cuts[1:nbands],
      points=(cuts[2:ncuts] + cuts[1:nbands])/2
      )
}

NCStepParam <- function(step, min=0, max=+Inf){
  new("NCStepParam",
      min=min,
      max=max,
      step=step
      )
}


NCAdaptedStepParam <- function(step, fromT=0, toT, mult=2, min=0, max=+Inf){
  cutt <- cutTfromto(fromT, toT, step, mult, min, max)
  new("NCAdaptedStepParam",
      Nstep=cutt$Nstep,
      theSteps=cutt$stepT,
      min=min,
      max=max,
      step=step,
      from=fromT,
      to=toT
      )
}


