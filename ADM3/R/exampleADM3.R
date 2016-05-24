`exampleADM3` <- function(outfile="example_ouput_ADM3.txt") {
	exFile= paste(system.file("extdata", package = "ADM3"), "/", "example.bed", sep="");
	plotADM3(ADM3(file=exFile, outfile=outfile));
}