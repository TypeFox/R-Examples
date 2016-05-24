ImputateGaps <-
function(dataset, impt.times, responses, method) {
  dataset <- dataset[dataset[, 2] <= max(impt.times), ]
  if(method == "locf") {
    LastObsCarriedForward <- function(data.id) {
      for(k in 1:length(responses)) {
        meassure <- data.id[, responses[k]]
        last.obs <- max(which(!is.na(meassure)))
        meassure[1:last.obs] <- na.locf(meassure, na.rm=FALSE, xout=1:last.obs)
        data.id[, responses[k]] <- meassure 
      }
      return(data.id)
    }
    dataset.split.id <- split(dataset, dataset[, 1])
    dataset.split.id.impt <- lapply(dataset.split.id, LastObsCarriedForward)
    dataset <- do.call("rbind", dataset.split.id.impt)
  }
  if(method == "approx") {
    imptDfApprox <- function(data.id) {
      for(k in 1:length(responses)) {
        if(length(data.id[!is.na(data.id[, responses[k]]), responses[k]]) > 1) {
          data.id[, responses[k]] <- na.approx(data.id[, responses[k]],
                                          na.rm=FALSE)
        }
      }
      return(data.id)   
    }
    dataset.split.id <- split(dataset, dataset[, 1])
    dataset.split.id.impt <- lapply(dataset.split.id, imptDfApprox)
    dataset <- do.call("rbind", dataset.split.id.impt)
  }
  return(dataset)
}
