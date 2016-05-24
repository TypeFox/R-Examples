# Hinda Haned, December 2008
# haned@biomserv.univ-lyon1.fr

##################################
# forensim  Classes Definitions 
# & Auxilary Functions
##################################




# Specifying attributes types : unions of classes
setClassUnion("vectorOrNULL", c("vector","NULL"))
setClassUnion("vectorOrdataframe", c("vector","data.frame"))
setClassUnion("matrixOrdataframe", c("matrix","data.frame"))
setClassUnion("factorOrNULL", c("factor","NULL"))
setClassUnion("listOrdataframe", c("list","data.frame"))
setClassUnion("characterOrNULL", c("character","NULL"))



#simugeno: simulated profiles object
setClass("simugeno",representation(tab.freq="list",nind="numeric",pop.names="factorOrNULL",popind="factorOrNULL",which.loc="characterOrNULL",
tab.geno="matrix",indID="character"))

# simumix: simulated mixture objects
setClass("simumix",representation(ncontri="numeric",mix.prof="matrix",mix.all="list",which.loc="characterOrNULL", popinfo="factorOrNULL"))

#tabfreq:  simulated allele frequencies 
setClass("tabfreq",representation(tab="listOrdataframe",which.loc="characterOrNULL",pop.names="factorOrNULL"))


#####################
# Auxilary Functions
#####################

# Function names
setMethod("names", signature(x = "simugeno"), function(x){
    return(slotNames(x))
})

setMethod("names", signature(x = "simumix"), function(x){
    return(slotNames(x))
})


setMethod("names", signature(x = "tabfreq"), function(x){
    return(slotNames(x))
})



#Method  show 
# Method show  for the simugeno class
setMethod("show", "simugeno", function(object){
  x <- object
  cat("\n")
  cat("   # Simugeno object: simulated genotypes # \n")
  cat("\n")
  cat("\n@which.loc: vector of ", length(x@which.loc), " locus names")
  cat("\n@nind:", ifelse(length(x@nind)==0,0, x@nind))
  cat("\n@indID:  vector of the individuals ID")
  cat("\n@tab.geno: ", nrow(x@tab.geno),"x", ncol(x@tab.geno), "data frame of genotypes")
  cat("\n@tab.freq:  allele frequencies for the", length(x@which.loc), " loci")
  cat("\n\n")
  cat("Population-related information: ")
  cat("\n@pop.names:", ifelse(is.null(x@pop.names), "- empty -", "population names"))
  cat("\n@popind:", ifelse(is.null(x@popind), "- empty -", "factor giving the population of each individual"))
  cat("\n")
} 
)


# Method show for  the simumix class
setMethod("show", "simumix", function(object){
  x <- object
  cat("\n")
  cat("   # Simumix object: simulated mixture # \n")
  cat("\n")
  cat("\n@which.loc: vector of ", length(x@which.loc), " locus names")
  cat("\n@ncontri:", ifelse(length(x@ncontri)==0, 0, x@ncontri))
  cat("\n@mix.prof: ", nrow(x@mix.prof),"x",ncol(x@mix.prof),"data frame of the contributors genotypes")
  cat("\n@mix.all: list of the alleles found in the mixture")
  cat("\n@popinfo:",  ifelse(is.null(x@popinfo), "- empty -", "populations of the contributors"))
  cat("\n")
} 
)

# Method show for the tabfreq class
setMethod("show","tabfreq", function(object) {
	x <- object
	cat("\n")
	cat("   # Tabfreq object: allele frequencies  # \n")
	cat("\n")
	#cat("\n@inputab: input data ")
	cat("\n@tab: list of allele frequencies")
	#cat("\n\nOptionnal : ")
	#code copie confrome de celui de Tibo, a changer
	cat("\n@which.loc: vector of ", length(x@which.loc), " locus names")
	cat("\n@pop.names: ", ifelse(is.null(x@pop.names), "- empty -", "populations names"))
	#cat("\n@theta: ", ifelse(is.null(x@theta), "- empty -", "coancestry coefficient(s)"))
	cat("\n")
}
)


##################
# Basic Methods
##################

is.simugeno <- 
function(x)
{
  res <- ( is(x, "simugeno") & validObject(x))
  return(res)
}

is.simumix <- 
function(x)
{
  res <- ( is(x, "simumix") & validObject(x))
  return(res)
}

is.tabfreq <- 
function(x)
{
	res <- ( is(x, "tabfreq") & validObject(x))
	return(res)
}





