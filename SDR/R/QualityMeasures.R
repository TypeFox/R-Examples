

#------------------------------------------------------------------------
# 
#        Devuelve el .coverage ( N(Cond) / Ns )  de la regla
#        
#        Donde:
#        ? N(Cond) es el n?mero de ejemplos que cumplen el antecedente
#        ? Ns es el n?mero de ejemplos del dataset
#        
#------------------------------------------------------------------------

.coverage <- function(x){
  c <- x[[1]] / x[[4]]
  if(is.nan(c)){
    0
  } else {
    c
  }
}


#------------------------------------------------------------------------
# 
#        Devuelve el soporte ( N(TargetValue ? Cond) / Ns )  de la regla
#        
#        Donde:
#        ? N(TargetValue ? Cond) es el n?mero de ejemplos que cumplen el antecedente y tienen como consecuente nuestra clase objetivo
#        ? Ns es el n?mero de ejemplos del dataset
#        
#------------------------------------------------------------------------


.Csupport <- function(x){
  c <- x[[2]] / x[[4]]
  if(is.nan(c)){
    0
  } else {
    c
  }
}

.Fsupport <- function(x){
  c <- x[[11]] / x[[4]]
  if(is.nan(c)){
    0  
  } else {
    c
  }
}

.FLocalSupport <- function(x){
  c <- x[[12]] / x[[9]]
  if(is.nan(c)){
    0  
  } else {
    c
  }
}

#------------------------------------------------------------------------
# 
#        Devuelve la confianza ( N(TargetValue ? Cond) / N(Cond) )  de la regla
#        
#        Donde:
#        ? N(TargetValue ? Cond) es el n?mero de ejemplos que cumplen el antecedente y tienen como consecuente nuestra clase objetivo
#        ? N(Cond) es el n?mero de ejemplos que cumplen el antecedente
#        
#------------------------------------------------------------------------

.confianza <- function(x){
  
  c <- x[[2]] / x[[1]]
  if(is.nan(c)){
    0
  } else {
    c
  }
  
}

.confianzaDifusa <- function(x){
  c <- x[[11]] / x[[10]]
  if(is.nan(c)){
    0
  } else {
    c
  }
}



.unusualness <- function(x){
  .coverage(x) * ( .confianza(x) - (x[[5]] / x[[4]]) )
}

.norm_unusualness <- function(x){
  .coverage(x) * .confianza(x)
  
}


.sensitivity <- function(x){
  c <- x[[2]] / x[[5]]
  if(is.nan(c)){
    0
  } else {
    c
  }
}


.accuracy <- function(x){
  c <-  (x[[2]] + 1) / (x[[1]] + NROW(x[[7]]))
  if(is.nan(c)){
    0
  } else {
    c
  }
}



.significance <- function(x){
  cove <- .coverage(x = x)
  if(cove > 0){
    SumsingClase <- 0
    for(i in 1:length(x[[6]])){
      if(x[[6]][i] > 0 && x[[7]][i] > 0){
        SumsingClase <- SumsingClase + x[[6]][i] * log10(x = (x[[6]][i] / (x[[7]][[i]] * cove)  ))
      }
    }
    
    2 * SumsingClase
  
  } else {
    0
  }
}



.LocalSupport <- function(x){
  
  c <- x[[8]] / x[[9]]
  if(is.infinite(c)){
    0
  } else {
    c
  }
}





