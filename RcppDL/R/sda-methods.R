setMethod("summary", 
          signature=c("sda"),
          function(object) return(object@backend$summary()))

setMethod("pretrain", 
          signature=c("sda"),
          function(object) return(object@backend$pretrain()))

setMethod("finetune", 
          signature=c("sda"),
          function(object) return(object@backend$finetune()))

setMethod("PretrainLearningRate",
          signature=c("sda"),
          function(object) {
            info <- summary(object)           
            return(info$PretrainLearningRate)
          })
          
setMethod("CorruptionLevel",
          signature=c("sda"),
          function(object) {
            info <- summary(object)           
            return(info$CorruptionLevel)
          })

setMethod("PretrainingEpochs",
          signature=c("sda"),
          function(object) {
            info <- summary(object)           
            return(info$PretrainingEpochs)
          })

setMethod("FinetuneLearningRate",
          signature=c("sda"),
          function(object) {
            info <- summary(object)           
            return(info$FinetuneLearningRate)
          })

setMethod("FinetuneEpochs",
          signature=c("sda"),
          function(object) {
            info <- summary(object)           
            return(info$FinetuneEpochs)
          })
          
setMethod("predict",
          signature=c("sda","matrix"),
          function(object,test) {
            return(object@backend$predict(test))
          })
          
setMethod("setPretrainLearningRate",
          signature=c("sda","numeric"),
          function(object, x) {
            return(object@backend$setPretrainLearningRate(x))
          })

setMethod("setPretrainEpochs",
          signature=c("sda","numeric"),
          function(object, x) {
            return(object@backend$setPretrainEpochs(x))
          })
          
setMethod("setFinetuneLearningRate",
          signature=c("sda","numeric"),
          function(object, x) {
            return(object@backend$setFinetuneLearningRate(x))
          })     
          
setMethod("setFinetuneEpochs",
          signature=c("sda","numeric"),
          function(object, x) {
            return(object@backend$setFinetuneEpochs(x))
          })                

setMethod("setCorruptionLevel",
          signature=c("sda","numeric"),
          function(object, x) {
            return(object@backend$setCorruptionLevel(x))
          })  

setMethod("LearningRate",
          signature=c("sda"),
          function(object) {
            info <- summary(object)   
            ll <- list()
            ll$'PretrainLearningRate' <- info$PretrainLearningRate
            ll$'FinetuneLearningRate' <- info$FinetuneLearningRate
            return(ll)
          })

Rsda <- function(x, y, hidden){
    sdaModule <- new(Sda)
    sdaModule$init(x, y, hidden)
    return(new("sda",
               backend=sdaModule))

}
