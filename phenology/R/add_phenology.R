#' add_phenology creates a new dataset.
#' @title Create a new dataset or add a timeserie to a previous dataset.
#' @author Marc Girondot
#' @return Return a list of formated data that can be used ith fit_phenology()
#' @param previous Previous data formated with add_phenology or NULL [default] if no previous data exist
#' @param add The data to be added. It can be a set of several entities that uses the same reference and date format
#' @param name The name of the monitored site
#' @param reference as.Date('2001-12-31') The date used as 1st date
#' @param sep.dates Separator used to separate dates when incertitude is included
#' @param month_ref If no reference date is given, use this month as a reference
#' @param header If the data is read from a file, can be used to force header or not
#' @param format The format of the date in the file. Several format can be set and the last one that give compatible result is used
#' @param silent Does information about added timeseries is shown
#' @description To create a new dataset, the syntaxe is :\cr
#' data <- add_phenology(add=newdata, name="Site", reference=as.Date('2001-12-31'), 
#' format='\%d/\%m/\%y')\cr\cr
#' To add a dataset to a previous one, the syntaxe is :\cr
#' data <- add_phenology(previous=previousdata, add=newdata, name='Site', \cr
#' reference=as.Date('2001-12-31'), adjust_ref=TRUE, format="\%Y-\%m-\%d") \cr\cr
#' To add several timeseries at the same time with '\%d/\%m/\%y' or '\%d/\%m/\%Y' date format:\cr
#' data<-add_phenology(add=list(newdata1, newdata2), name=c('Site1', 'Site2'),\cr 
#' reference=as.Date('2001-12-31'), format=c('\%d/\%m/\%y', '\%d/\%m/\%Y'))\cr\cr
#' The dataset to be added must include 2 or 3 columns.\cr
#' The first one is the date in the format specified by 
#' the parameter format=. If the number of nests is known  
#' for an exact data, then only one date must be indicated.\cr  
#' If the number of nests is known for a range of date, the 
#' first and last dates must be separated but a - (dash).\cr
#' For example: 1/2/2000-10/2/2000\cr\cr
#' The second column is the number of nests observed for 
#' this date or this range of dates.\cr
#' The third column is optional and is the name of the rookery.\cr\cr
#' If only two columns are indicated, the name can be indicated as  
#' a parameter of the function with name=. If no name is indicated,  
#' the default name Site will be used, but take care, only one 
#' rookery of this name can be used.\cr\cr
#' Several rookeries can be included in the same file but in this case 
#' the rookery name is obligatory at the third column.\cr\cr
#' The simplest use of this function is just: \cr
#' phen <- add_phenology()\cr\cr
#' Some problems that can occur:\cr
#' If a name is defined as a third column of a data.frame and a name is defined also with name=, the third column has priority.\cr
#' Two different timeseries MUST have different name and character _ is forbiden in timeseries names.
#' @examples
#' \dontrun{
#' # Get the lastest version at:
#' # install.packages("http://www.ese.u-psud.fr/epc/conservation/CRAN/phenology.tar.gz", 
#'      repos=NULL, type="source")
#' library(phenology)
#' # Read a file with data
#' Gratiot<-read.delim("http://max2.ese.u-psud.fr/epc/conservation/BI/Complete.txt", header=FALSE)
#' data(Gratiot)
#' # Generate a formatted list nammed data_Gratiot 
#' refdate <- as.Date("2001-01-01")
#' data_Gratiot <- add_phenology(Gratiot, name="Complete", 
#' 	reference=refdate, format="%d/%m/%Y")
#' # Generate initial points for the optimisation
#' parg <- par_init(data_Gratiot, parametersfixed=NULL)
#' # Run the optimisation
#' result_Gratiot <- fit_phenology(data=data_Gratiot, parametersfit=parg, 
#' 	parametersfixed=NULL, trace=1)
#' data(result_Gratiot)
#' # Plot the phenology and get some stats
#' output <- plot(result_Gratiot)
#' }
#' @export


add_phenology <-
function(add=file.choose(), name=NULL, reference=NULL, month_ref= NULL, sep.dates="-", 
         header=NULL, format=NULL, previous=NULL, silent=FALSE) {


# name=NULL; reference=NULL; month_ref= NULL; header=NULL; sep.dates="-"; format=NULL; previous=NULL; silent=FALSE
# add=file.choose()
  
if (class(previous)!="phenologydata" && !is.null(previous)) {
  stop("The previous dataset must be already formated!")
}
  
  if (!is.null(format)) {
    if (any(grepl(sep.dates, format))) {
    stop("Separator between day, month and year cannot be the same as the separator between two dates")
    }
  }
  

if(class(add)=="try-error") {
	stop("No file has been chosen!")
}	

nm <- name

if (class(add)=="character") {
# j'ai utilisé le file.choose
	nm <- ifelse(is.null(name), add, name)
	add <- lapply(add,readLines, warn=FALSE)
  addec <- add[[1]]
	add[[1]] <- addec[addec!=""]
}
	
## if (is.null(add) || !exists(as.character(substitute(add)))) {
if (is.null(add)) {
	stop("Data to be added does not exist!")
}


# rp <- phenology:::.read_phenology(add, header, reference, month_ref, format, nm, sep.dates)

rp <- .read_phenology(add, header, reference, month_ref, format, nm, sep.dates)

if (any(grepl(sep.dates, rp$format))) {
  stop("The date separator is used also within a date. It is not possible.")
}

add <- rp$DATA
reference <- rp$reference
format <- rp$format


if (dim(add)[2]==3) {
  if (all(add[1,3]==add[,3])) name <- add[1,3]
}
# je ne comprends rien 14/3/2014
# là je crée une liste
	nbdatasets <- 1
	addlist <- list(add)
if (is.null(name)) {
	if (is.null(nm)) {
    if (length(deparse(substitute(add)))==1) {
		names(addlist) <- deparse(substitute(add))
	    } else {
	  names(addlist) <- paste("Rookery", runif(1, 1, 10000))
	  }
	} else {
		names(addlist) <- basename(nm)
	}
} else {
		names(addlist) <- name
}



# print(paste("2: ", name))


for (kk in 1:nbdatasets) {

	add <- addlist[[kk]]
	name<-names(addlist)[kk]

# print(paste("3: ", name))

# fichier pnd
# si 5 colonnes, j'en retire les 2 dernières
	if (dim(add)[2]==5) {
		add<-add[,-c(3:5)]
	}
#Si j'ai 3 colonnes et rien sur la 3ème, je n'en garde que deux
	if (dim(add)[2]==3) {
		if (all(is.na(add[,3]))) {
			add<-add[,-3]
		}
	}

  # 23122002
  fin  <- FALSE
	if (dim(add)[2]==2) {fin <- TRUE} else {
    if (length(levels(factor(add[,3])))==1) fin <- TRUE
	}
  
  
	if (fin) {


# J'ai deux colonnes et le nom des séries dans name
	if (!silent) message(name)
# Je n'ai pas de nom de site. Il n'y a donc qu'une seule série
		colnames(add)=c("Date", "Nombre")	
		add$Date<-as.character(add$Date)
		add[, 2]<-as.character(add[, 2])

# je retire toutes les lignes pour lesquelles je n'ai pas d'observation de ponte - 19012012
		add<-add[(add[,1]!="") & (!is.na(add[,1])) & (gsub("[0-9 ]", "", add[,2])=="") & (!is.na(add[,2])) & (add[,2]!=""),1:2]
		
		if (dim(add)[1]==0) {
		  if (!silent) warning(paste("The timeseries", name, "is empty; check it"))
			return(invisible())
		} 

		
		# Je le mets en numérique: 21/3/2012
		add[, 2]<-as.numeric(add[, 2])
		
		addT<-data.frame(Date=rep(as.Date("1900-01-01"), dim(add)[1]), Date2=rep(as.Date("1900-01-01"), dim(add)[1]), nombre=rep(NA, dim(add)[1]), ordinal=rep(NA, dim(add)[1]), ordinal2=rep(NA, dim(add)[1]), Modeled=rep(NA, dim(add)[1]), LnL=rep(NA, dim(add)[1]))
		

		if (is.null(reference)) {		
			warning("No refence date can be calculated")
			return(invisible())
		}
		
		
    if (!silent) message(paste("Reference: ", reference))
		
		# dans i la ligne en cours
		for(i in 1:dim(add)[1]) {
		
			essai<-unlist(strsplit(add$Date[i], sep.dates))
			
			dtcorrect <- NULL
			for (fd in 1:length(format)) {
				dtencours <- as.Date(essai[1], format[fd])
				ref <- as.numeric(dtencours-reference+1)
				if (ref>=0 & ref<=366) {
					dtcorrect <- dtencours
				}
			}
			
			addT$Date[i] <- dtcorrect
			if (is.na(addT$Date[i])) {
				warning(paste("Error date ", essai[1], sep=""))
			} else {
				addT$ordinal[i]<-as.numeric(addT$Date[i]-reference+1)
				if (length(essai)==2) {
					dtcorrect <- NULL
					for (fd in 1:length(format)) {
						dtencours <- as.Date(essai[2], format[fd])
						ref <- as.numeric(dtencours-reference+1)
						if (ref>=0 & ref<=366) {
							dtcorrect <- dtencours
						}
					}
							
					addT$Date2[i]<-dtcorrect
					addT$ordinal2[i]<-as.numeric(addT$Date2[i]-reference+1)
				} else {
					addT$Date2[i]<-NA
					addT$ordinal2[i]<-NA
				}
			}
		}
		addT$nombre<-add$Nombre
		nb<-length(previous)
		previous$kyYI876Uu<-addT
		names(previous)[nb+1]<-name

		

	} else {
# J'ai plus d'un nom de site. Il y a donc plusieurs séries - 21012012
		colnames(add)=c("Date", "Nombre", "Site")	
		add$Date<-as.character(add$Date)

# Je le mets en numérique: 21/3/2012
		add[, 2]<-as.numeric(add[, 2])
		
# Maintenant j'ai tous les sites remplis
# mais je dois avoir la liste des sites
		fac <- factor(add$Site)
# je génère un nouveau data.frame avec facteur par facteur
		dtaorigin=previous
		for(i in 1:length(levels(fac))) {
			siteencours<-levels(fac)[i]
			addencours<-add[add[,3]==siteencours, 1:2]
			dtaorigin<-add_phenology(previous=dtaorigin, add=addencours, name=siteencours, reference=NULL, month_ref= month_ref, format=format)
		}
		previous<-dtaorigin

	}
	problem<-FALSE
	for(i in 1:length(previous)) {
		problem<-(problem) || (any(previous[[i]]$ordinal[!is.na(previous[[i]]$ordinal)]>366)) || (any(previous[[i]]$ordinal2[!is.na(previous[[i]]$ordinal2)]>366)) || (any(previous[[i]]$ordinal[!is.na(previous[[i]]$ordinal)]<0)) || (any(previous[[i]]$ordinal2[!is.na(previous[[i]]$ordinal2)]<0))

	}
#	print(problem)
	if (problem) {
		warning(paste("Problem in at least one date; check them. Take care about format used.\n", 
                  "Within a file, all dates must conform to the same format.\n",
                  "Data should not be longer than one year for a site at the same time."))
	}	
	

}

class(previous) <- "phenologydata"


if (any(grepl("_", names(previous)))) {
  print("The character _ is forbinden in names of timesseries. It has been changed to -.")
  names(previous) <- gsub("_", "-", names(previous))
}

if (any(duplicated(names(previous)))) {
  stop("The names of timesseries must be unique.")
}


return(previous)
}
