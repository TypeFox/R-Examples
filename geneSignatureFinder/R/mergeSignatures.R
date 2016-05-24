mergeSignatures <-
function(signature1, signature2, cpuCluster = NULL, stopCpuCluster = TRUE) {
  if(is.null(cpuCluster)) {
    print("!!!")
    return(NULL)
  }
  message("... merging: ", signature1$signatureName, " and ",
signature2$signatureName)

  startingSignatureIDs <- intersect(signature1$signatureIDs, signature2$signatureIDs)
  if(length(startingSignatureIDs) == 0) {
    message("No common genes have benn found in the two signatures.\nUsing both original starting genes")
    startingSignatureIDs <- c(signature1$signatureIDs[signature1$startingSignature],
                              signature2$signatureIDs[signature2$startingSignature])
  }
  
  message(paste("starting gene/s: ", paste(colnames(geData)[startingSignatureIDs], collapse = ", "), sep = ""))
  featuresSpace <- union(signature1$signatureIDs, signature2$signatureIDs)
  result <- signatureFinder(startingSignatureIDs,
                            subsetToUse = featuresSpace,
                            cpuCluster = cpuCluster,
                            stopCpuCluster = stopCpuCluster)
  result$featureSpace <- colnames(geData)[featuresSpace]
  
  result$signatureName <-  signature1$signatureName
  # c'è un problema di riniminazione della signature!!!!!!!!!!!!!!!!!! al momento la funzione non è visibile
  if(is.null(result$mergedBy)) 
      result$mergedBy <- paste(signature1$signatureName, "*", signature2$signatureName, sep = "") else {
        result$mergedBy <- paste("(", result$mergedBy, ")", sep = "")
        result$mergedBy <- paste(result$mergedBy, "*", signature1$signatureName, "*", signature2$signatureName, sep = "")
        }
  #result$mergedBy <- paste(signature1$startingSignature, "*", signature2$startingSignature, sep = "")
  result$firstSignature <- signature1$signature
  result$firstTValue <- signature1$tValue
  result$firstClassification <- signature1$classification
  result$secondSignature <- signature2$signature
  result$secondTValue <- signature2$tValue
  result$secondClassification <- signature2$classification
  result$dissimilarityBetweenSignatures <- table(result$firstClassification, result$secondClassification)
  result$dissimilarityBetweenSignatures <- 1 - sum(diag(result$dissimilarityBetweenSignatures))/
    sum(result$dissimilarityBetweenSignatures)

#  if(sum(result$signature == signature1$startingSignature) > 0)
#    result$startingSignature <- signature1$startingSignature else
#  result$startingSignature <- signature2$startingSignature
#    message("signature renamed as: ", result$startingSignature, "\n")
  return(result)
}
