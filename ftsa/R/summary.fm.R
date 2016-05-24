summary.fm = function (object, ...) 
{
    if(class(object)[1] == "fm"|class(object)[1] == "ftsm"){
       print(object)
       junk <- MISE(object$y, object$fitted)
       junk1 <- cbind(junk$ME, junk$MSE, junk$MPE, junk$MAPE)
       colnames(junk1) <- c("ME", "MSE", "MPE", "MAPE")
       junk2 <- cbind(junk$MIE, junk$MISE, junk$MIPE, junk$MIAPE)
       colnames(junk2) = c("IE", "ISE", "IPE", "IAPE")
       cat(paste("\nAverages across x:\n"))
       print(round(colMeans(junk1), 5))
       cat(paste("\nAverages across time:\n"))
       print(round(colMeans(junk2), 5))
       cat("\n")
    }
    else{
       warning("object is neither a functional time series model nor a functional model.") 
    }
}
