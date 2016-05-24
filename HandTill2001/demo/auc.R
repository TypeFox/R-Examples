library(HandTill2001)
data(ht01.twoclass)
data(ht01.multipleclass)
message("AUC for a two class response")
auc(bincap(
	    response = as.factor(ht01.twoclass$observed)
	    , predicted = ht01.twoclass$predicted
	    , true = "1"
	    )
)
message("AUC for a multiple class response")
auc(multcap(
	     response = ht01.multipleclass$observed
	     , predicted = as.matrix(ht01.multipleclass[, levels(ht01.multipleclass$observed)])
	     )
)

