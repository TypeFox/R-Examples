#################################
#'  Read QuantiNemo extended format for genotype files
#'  
#'  Read QuantiNemo (\url{http://www2.unil.ch/popgen/softwares/quantinemo/}) genotype files extended format (option 2)
#'  
#'  @usage qn2.read.fstat(fname, na.s = c("NA","NaN"))
#'  @param fname quantinemo file name
#'  @param na.s  na string used
#'  @return dat  a data frame with nloc+1 columns, the first being the population
#'  to which the individual belongs and the next being the genotypes, one column per locus; 
#'  and ninds rows
#'  @return sex  the sex of the individuals
#'  @author Jerome Goudet \email{jerome.goudet@@unil.ch}
#'  @seealso \code{\link{read.fstat}}
#'  @references \href{http://www2.unil.ch/popgen/softwares/quantinemo/2008_Neuenschwander_et_al_BioInf_quantiNEMO.pdf}{Neuenschwander S, Hospital F, Guillaume F, Goudet J (2008)} 
#'  quantiNEMO: an individual-based program to simulate quantitative traits with explicit 
#'  genetic architecture in a dynamic metapopulation Bioinformatics 24, 1552-1553.
#'  @examples 
#'   dat<-qn2.read.fstat(system.file("extdata","qn2_sex.dat",package="hierfstat"))
#'   sexbias.test(dat[[1]],sex=dat[[2]])
#'  @export
########################################################################################  
qn2.read.fstat<-function (fname, na.s = c("NA","NaN")) {
  #written to allow direct reading of quantinemo genotype files extended format (option 2) in R
  #eliminates juveniles and ignores the columns after age (ind and parents id)
  #split the data set in correct format (dat) and the vector of sexes (1:M and 2:F)
  x <- scan(fname, n = 4)
  nloc <- x[2]
  lnames <- scan(fname, what = character(), skip = 1, nlines = nloc)
  lnames <- c("Pop", lnames)
  dat <- scan(fname, skip = nloc + 1, na.strings = na.s,comment.char="_")
  dat <- data.frame(matrix(dat, ncol = nloc + 4, byrow = TRUE))
  age<-dat[,nloc+2]
  sex<-dat[age==2,nloc+3]
  asex<-character(length(sex))
  asex[sex==0]<-"M"
  asex[sex==1]<-"F"
  dat<-dat[age==2,1:(nloc+1)]
  names(dat) <- lnames
  return(list(dat=dat,sex=asex))
}

