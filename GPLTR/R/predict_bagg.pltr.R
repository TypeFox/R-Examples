predict_bagg.pltr <- function(bag_pltr, Y.name, newdata, type = "response", thresshold = seq(0, 1, by = 0.1))
{
  
    predict_glm <- lapply(bag_pltr$Glm_BAG, function(uw)
    {
      pred <- predict.glm(uw, newdata = newdata, type = type)
      return(sapply(thresshold, function(wz) as.numeric(pred > wz)))
    })
    
    Bag <- length(bag_pltr$Glm_BAG)
    PRED_IND <- list()
    for(jj in seq(length(thresshold))){
      PRED_IND[[jj]] <- sapply(1:Bag, function(ww) predict_glm[[ww]][,jj])
    }
    
    FINAL_PRED_IND1 <- lapply(PRED_IND, function(www) apply(www, 1, function(zzz) as.numeric(mean(zzz) > 0.5)))
    names(FINAL_PRED_IND1) <- paste('CUT', 1: length(thresshold), sep = '')
    PRED_ERROR1 <- sapply(FINAL_PRED_IND1, function(uuu) mean( uuu != newdata[, Y.name]))    
  
  
  PRED_ERRORS_PBP <- lapply(1: length(thresshold),function(vvv)
  { 
    return(sapply(1: Bag, function(ww) mean(newdata[, Y.name] != predict_glm[[ww]][, vvv])))
  })
  names(PRED_ERRORS_PBP) <- paste('CUT', 1: length(thresshold), sep = '')
  PRED_ERROR_PBP <- sapply(PRED_ERRORS_PBP, function(uu) mean(uu))
  
  
  PROB_LIST <- lapply(bag_pltr$Glm_BAG, function(uu){
  pred <- predict.glm(uu, newdata = newdata, type = type)
  return(pred)  
  })
  PROB_MAT <- matrix(unlist(PROB_LIST), ncol = Bag, byrow = FALSE)
  PROB_VECT <- apply(PROB_MAT, 1, mean)
  FINAL_PRED_IND2 <- sapply(thresshold, function(ttt) as.numeric(PROB_VECT > ttt))
  FINAL_PRED_IND2 <- as.list(as.data.frame( FINAL_PRED_IND2))
  names(FINAL_PRED_IND2) <- paste('CUT', 1: length(thresshold), sep = '')
  confusion2 <-lapply(FINAL_PRED_IND2, function(cc) table(cc,newdata[, Y.name], dnn = c("Predicted Class", "Observed Class")))
  PRED_ERROR2 <- sapply(FINAL_PRED_IND2, function(uuu) mean( uuu != newdata[, Y.name]))    
  confusion1 <-lapply(FINAL_PRED_IND1, function(cc) table(cc,newdata[, Y.name], dnn = c("Predicted Class", "Observed Class")))
  
  return(list(FINAL_PRED_IND1 = FINAL_PRED_IND1, FINAL_PRED_IND2 = FINAL_PRED_IND2,  PRED_ERROR1 = PRED_ERROR1, 
              PRED_ERROR2 = PRED_ERROR2, CONF1 = confusion1, CONF2 = confusion2, PRED_ERRORS_PBP = PRED_ERRORS_PBP,
              PRED_ERROR_PBP = PRED_ERROR_PBP))
  
}