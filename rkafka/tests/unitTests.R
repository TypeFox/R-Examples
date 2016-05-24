is.installed <- function(mypkg){
  is.element(mypkg, installed.packages()[,1])
} 
library(rkafka)
  
producer1=rkafka.createProducer("127.0.0.1:9092")
condition=rJava::.jinstanceof(producer1,"com/musigma/producer/MuProducer")
RUnit::checkTrue(condition,"check1")

simpleConsumer1=rkafka.createSimpleConsumer("127.0.0.1","9092","10000","10000","test")
condition2=rJava::.jinstanceof(simpleConsumer1,"com/musigma/consumer/MuSimpleConsumer")
RUnit::checkTrue(condition2,"check2")

