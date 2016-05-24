evalucc <-
function(CC, EC)
{
index<-list()
cat("\n","Evaluation of core collection","\n","\n")

index[["DD"]] <- DD(CC, EC)
cat("1) DD%      ",index[["DD"]],"\n","\n")

index[["MD"]] <- MD(CC, EC)
cat("2) MD%      ",index[["MD"]],"\n","\n")


index[["VD"]] <- VD(CC, EC)
cat("3) VD%      ",index[["VD"]],"\n","\n")

index[["RR"]] <- RR(CC, EC)
cat("4) RR%      ",index[["RR"]],"\n","\n")

invisible(index)
}
