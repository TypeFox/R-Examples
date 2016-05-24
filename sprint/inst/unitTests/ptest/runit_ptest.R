test.mpiIsWorkingOnSeveralProcesses <- function()
{	
	result <- ptest()
	print(paste(result))
	one_result <- paste(result, collapse = '')
	checkTrue(grepl("1", one_result), "MPI is not running across multiple processes. Did you run this test from the command line using 'mpiexec -n'?")
}