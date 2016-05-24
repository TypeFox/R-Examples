setMethod("summary", 
          signature=c("rbm"),
          function(object) return(object@backend$summary()))

setMethod("train", 
          signature=c("rbm"),
          function(object) return(object@backend$train()))

setMethod("reconstruct",
          signature=c("rbm","matrix"),
          function(object,test) {
            return(object@backend$reconstruct(test))
          })

setMethod("setLearningRate",
          signature=c("rbm","numeric"),
          function(object, x) {
            return(object@backend$setLearningRate(x))
          })

setMethod("LearningRate",
          signature=c("rbm"),
          function(object) {
            info <- summary(object)           
            return(info$LearningRate)
          })

setMethod("setTrainingEpochs",
          signature=c("rbm","numeric"),
          function(object, x) {
            return(object@backend$setTrainingEpochs(x))
          })

setMethod("TrainingEpochs",
          signature=c("rbm"),
          function(object) {
            info <- summary(object)           
            return(info$TrainingEpochs)
          })

setMethod("setHiddenRepresentation",
          signature=c("rbm","numeric"),
          function(object, x) {
            return(object@backend$setHiddenRepresentation(x))
          })

setMethod("HiddenRepresentation",
          signature=c("rbm"),
          function(object) {
            info <- summary(object)           
            return(info$HiddenRepresentation)
          })

setMethod("setStep",
          signature=c("rbm","numeric"),
          function(object, x) {
            return(object@backend$setStep(x))
          })

setMethod("Step",
          signature=c("rbm"),
          function(object) {
            info <- summary(object)           
            return(info$ContrastiveDivergenceStep)
          })

Rrbm <- function(x){
    rbmModule <- new(Rbm)
    rbmModule$init(x)
    return(new("rbm", backend=rbmModule))
}
