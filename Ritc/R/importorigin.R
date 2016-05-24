# import data from origin ITC data sheet exported to .csv files
# input: file name as string for x
# output: list of vectors for each column in the .csv file, with names maintained. 
importorigin=function(x){
	options(warn=-1); # suppress warning messages
	itcdata1=read.csv(x, colClasses="character");
	itcdata2=lapply(itcdata1, as.numeric);
	# generates warning messages as filling uneven columns with NAs
	itcdata3=lapply(itcdata2, rmlastna);
}
