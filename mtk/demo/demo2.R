
##############
# testExpWorkflow,  start a session of sensitivity analysis
# with  ishigami_morris.xml, ishigami_fast.xml, WWDM_morris.xml as example.

	xmlFile <- "ishigami_morris.xml"

	xmlFilePath <- paste(path.package("mtk", quiet = TRUE),
		"/extdata/",xmlFile,sep = "") ## find where the examples are held. 
	
	workflow <- mtkExpWorkflow(xmlFilePath=xmlFilePath)
	
	cat("\n\n **** Begin of the workflow: *****\n\n")
	run(workflow)
	cat("\n\n **** End of the workflow: *****\n\n")
	
	cat("\n\n **** Begin of the reporting: *****\n\n")
	summary(workflow,maximum=1, digits=2)
	
	cat("\n\n **** End of the reporting: *****\n\n")


