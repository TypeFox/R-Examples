is.installed <- function(mypkg){
  is.element(mypkg, installed.packages()[,1])
} 
library(rkafkajars)
producerProperties <- .jnew("com/musigma/producer/ProducerProperties")
condition=rJava::.jinstanceof(producerProperties,"com/musigma/producer/ProducerProperties")
RUnit::checkTrue(condition,"check1")
