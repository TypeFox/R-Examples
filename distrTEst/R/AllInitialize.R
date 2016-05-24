## Initialization method
setMethod("initialize", "Evaluation",
          function(.Object, name = NULL, filename = NULL,
                    call.ev = NULL, result = NULL,
                    estimator, Data) {
             if(missing(Data)) Data <- Simulation();
             if(missing(estimator)) estimator <- mean;
             if(is.null(name))
                {if(!is.null(name(Data)))
                    name <- name(Data)
                 else
                    name <- "Default-name"}
             if(is.null(filename))
                {if(!is.null(filename(Data)))
                     filename <- filename(Data)
                 else
                     filename <- name}
            .Object@name <- name
            .Object@filename <- filename
            .Object@call.ev <- call.ev
            .Object@result <- result
            .Object@estimator <- estimator
            .Object@Data <- Data
            .Object
          })


