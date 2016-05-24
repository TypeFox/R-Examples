# To create geometric sequences:
#	gseq(.4,4E10,10,addzero=TRUE)
#	0 .4 4 4E1 4E1 4E3 4E4 4E5 4E6 4E7 4E8 4E9 4E10
 
gseq <- function(from, to, by, addrem = FALSE, addzero = FALSE) {
	retvec = c()
	addval = from
	
	if (addzero) {
		retvec[1] = 0
		count = 2
	} else count = 1
	
	while (addval <= to) {
		retvec[count] = addval
		addval = by * addval
		count = count + 1
	}
	
	if (addval != to && addrem) retvec[count] = to
	
	retvec
}