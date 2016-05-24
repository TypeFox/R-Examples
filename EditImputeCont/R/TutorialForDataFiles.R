makeDemoTxt <- function(NestedEx){
	
	write.table(NestedEx$D.obs, file = "Demo_D_obs.txt", row.names = FALSE, col.names = TRUE)
	write.table(NestedEx$Ratio.edit, file = "Demo_Ratio_edit.txt", row.names = FALSE, col.names = FALSE)
	write.table(NestedEx$Balance.edit, file = "Demo_Balance_edit.txt", row.names = FALSE, col.names = FALSE)
	write.table(NestedEx$Range.edit, file = "Demo_Range_edit.txt", row.names = FALSE, col.names = FALSE)
	
} 
