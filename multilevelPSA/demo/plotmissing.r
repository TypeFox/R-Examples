require(multilevelPSA)
if(require(pisa, quietly=TRUE)) {
	data(pisa.student)
	data(pisa.psa.cols)
	pisa.student = pisa.student[,c('CNT', pisa.psa.cols)]
	pisa.student$CNT = as.character(pisa.student$CNT)
	missing.plot(pisa.student, pisa.student$CNT)
} else {
	message("pisa package not installed. Try\nrequire(devtools)\ninstall_github('pisa','jbryer')")
}
