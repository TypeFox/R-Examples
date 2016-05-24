library(testthat)
library(GetoptLong)

is.solaris = function()
	grepl('SunOS',Sys.info()['sysname'])


if(!is.solaris()) {
	if(Sys.which("perl") != "") {
		test_check("GetoptLong", filter = "qq|GetoptLong")
	} else {
		test_check("GetoptLong", filter = "qq")
	}
}
