summary.agreement <-
function(object,dec,...){
    print <- cbind(round(object$CCC,4), round(object$Precision,4),round(object$Accuracy,4),round(object$TDI,dec),round(object$CP,4),round(object$RBS,2));
    colnames(print) <- c("CCC","Precision","Accuracy","TDI","CP","RBS");
    rownames(print) <- "";
    print(print,quote=FALSE, na.print=".",...);
}

