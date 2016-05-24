
.onAttach = function(libname, pkgname) {
	.jpackage(pkgname);
	if(canUseJava() && .Platform$OS!="unix"){
		lf<-.jcall("javax.swing.UIManager","S","getSystemLookAndFeelClassName");
		.jcall("javax.swing.UIManager",,"setLookAndFeel",lf);
		.jcall("javax.swing.JDialog",,"setDefaultLookAndFeelDecorated",TRUE);
	}

	desc <- packageDescription(pkgname)
	DQdate <-  desc$Date
	DQVersion =  desc$Version
	packageStartupMessage("This is ", pkgname, " ", desc$Version, " ", desc$Date)
	return(invisible(NULL))
}

