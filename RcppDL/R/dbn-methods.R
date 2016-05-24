setMethod("summary", 
          signature=c("dbn"),
          function(object) return(object@backend$summary()))

setMethod("pretrain", 
          signature=c("dbn"),
          function(object) return(object@backend$pretrain()))

setMethod("finetune", 
          signature=c("dbn"),
          function(object) return(object@backend$finetune()))

setMethod("PretrainLearningRate",
          signature=c("dbn"),
          function(object) {
            info <- summary(object)           
            return(info$PretrainLearningRate)
          })
          
setMethod("PretrainingEpochs",
          signature=c("dbn"),
          function(object) {
            info <- summary(object)           
            return(info$PretrainingEpochs)
          })

setMethod("FinetuneLearningRate",
          signature=c("dbn"),
          function(object) {
            info <- summary(object)           
            return(info$FinetuneLearningRate)
          })

setMethod("FinetuneEpochs",
          signature=c("dbn"),
          function(object) {
            info <- summary(object)           
            return(info$FinetuneEpochs)
          })
          
setMethod("predict",
          signature=c("dbn","matrix"),
          function(object,test) {
            return(object@backend$predict(test))
          })
          
setMethod("setPretrainLearningRate",
          signature=c("dbn","numeric"),
          function(object, x) {
            return(object@backend$setPretrainLearningRate(x))
          })

setMethod("setPretrainEpochs",
          signature=c("dbn","numeric"),
          function(object, x) {
            return(object@backend$setPretrainEpochs(x))
          })
          
setMethod("setFinetuneLearningRate",
          signature=c("dbn","numeric"),
          function(object, x) {
            return(object@backend$setFinetuneLearningRate(x))
          })     
          
setMethod("setFinetuneEpochs",
          signature=c("dbn","numeric"),
          function(object, x) {
            return(object@backend$setFinetuneEpochs(x))
          })                

setMethod("setStep",
          signature=c("dbn","numeric"),
          function(object, x) {
            return(object@backend$setStep(x))
          })

setMethod("Step",
          signature=c("dbn"),
          function(object) {
            info <- summary(object)           
            return(info$ContrastiveDivergenceStep)
          })

setMethod("LearningRate",
          signature=c("dbn"),
          function(object) {
            info <- summary(object)   
            ll <- list()
            ll$'PretrainLearningRate' <- info$PretrainLearningRate
            ll$'FinetuneLearningRate' <- info$FinetuneLearningRate
            return(ll)
          })

Rdbn <- function(x, y, hidden){
    dbnModule <- new(Dbn)
    dbnModule$init(x, y, hidden)
    return(new("dbn", backend=dbnModule))

}
