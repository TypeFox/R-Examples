setMethod("summary", 
          signature=c("da"),
          function(object) return(object@backend$summary()))

setMethod("train", 
          signature=c("da"),
          function(object) return(object@backend$train()))

setMethod("reconstruct",
          signature=c("da","matrix"),
          function(object,test) {
            return(object@backend$reconstruct(test))
          })

setMethod("setLearningRate",
          signature=c("da","numeric"),
          function(object, x) {
            return(object@backend$setLearningRate(x))
          })

setMethod("LearningRate",
          signature=c("da"),
          function(object) {
            info <- summary(object)           
            return(info$LearningRate)
          })

setMethod("setTrainingEpochs",
          signature=c("da","numeric"),
          function(object, x) {
            return(object@backend$setTrainingEpochs(x))
          })

setMethod("TrainingEpochs",
          signature=c("da"),
          function(object) {
            info <- summary(object)           
            return(info$TrainingEpochs)
          })

setMethod("setHiddenRepresentation",
          signature=c("da","numeric"),
          function(object, x) {
            return(object@backend$setHiddenRepresentation(x))
          })

setMethod("HiddenRepresentation",
          signature=c("da"),
          function(object) {
            info <- summary(object)           
            return(info$HiddenRepresentation)
          })

setMethod("setCorruptionLevel",
          signature=c("da","numeric"),
          function(object, x) {
            return(object@backend$setCorruptionLevel(x))
          })  

setMethod("CorruptionLevel",
          signature=c("da"),
          function(object) {
            info <- summary(object)           
            return(info$CorruptionLevel)
          })

Rda <- function(x){
    daModule <- new(dA)
    daModule$init(x)
    return(new("da", backend=daModule))
}
