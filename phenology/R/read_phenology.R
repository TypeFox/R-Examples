#' .read_phenology translate file format
#' @title The function .read_phenology
#' @author Marc Girondot \email{marc.girondot@@u-psud.fr}
#' @return Return a dataframe
#' @param obj_list A list with different datasets
#' @param header Does the timeseries has header.
#' @param reference Date used as reference. Is the day 1.
#' @param month_ref Reference month. Generally will be 1 (January) or 7 (July).
#' @param format Format for dates.
#' @param nm Path for file or name
#' @description Function for package Phenology
#' @export


.read_phenology <-
function(obj_list=NULL, header=NULL, reference=NULL, month_ref= NULL, format=NULL, 
         nm=NULL, sep.dates=NULL) {

# obj_list=NULL; header=NULL; reference=NULL; month_ref= NULL; format=NULL; nm=NULL
# rp <- .read_phenology(add, header, reference, month_ref, format, nm, sep.dates)
# obj_list=add; header=header;reference=reference;month_ref=month_ref; format=format; nm=nm;sep.dates=sep.dates

if (!is.null(format)) {
dtaspl <- substr(gsub("[%dmYy]", "", format), 1, 1)
  } else {
dtaspl <- "/"
warning("Separator between day, month and year is supposed to be the / character")
}
  
if (class(obj_list)=="list") {

				analyse <- unlist(obj_list)[1]
				if (length(grep("\t", analyse))!=0) {
					spl <- "\t"				
				
				} else {
					if (length(grep(";", analyse))!=0) {
						spl <- ";"
										
					} else {
						spl <- ","
					}
				}
				lsp <- lapply(obj_list, strsplit, split=spl)[[1]]
				comble <- max(unlist(lapply(lsp, length)))
				for (i in 1:length(lsp)) {
					lsp[[i]] <- c(lsp[[i]], rep(NA, comble-length(lsp[[i]])))
				}
				
				# il faut que je vérifie si un header				
				obj_list_df <- as.data.frame(t(as.data.frame(lsp, stringsAsFactors=FALSE)), stringsAsFactors=FALSE)
				row.names(obj_list_df) <- NULL
				
				if (any(obj_list_df[,1]=="#")) 
            obj_list_df <- obj_list_df[1:(which(obj_list_df[,1]=="#")-1), ]
				
				obj_list_df <- obj_list_df[complete.cases(obj_list_df[,1]),]
				obj_list_df <- obj_list_df[complete.cases(obj_list_df[,2]),]
				obj_list_df <- obj_list_df[obj_list_df[,2]!="",]
				obj_list_df <- obj_list_df[obj_list_df[,1]!="",]
								
				if (is.null(header)) {
				
				dta <- obj_list_df[2,1]
				dta <- obj_list_df[1,1]
				if (length(grep(dtaspl, dta))!=0) {header <- FALSE} else {header <- TRUE}				
				}
				
				if (header) {
					colnames(obj_list_df) <- obj_list_df[1,]
					obj_list_df <- obj_list_df[-1,]
				}
				
} else {
  obj_list_df <- obj_list
	if (any(obj_list_df[,1]=="#")) 
	  obj_list_df <- obj_list_df[1:(which(obj_list_df[,1]=="#")-1), ]
  
	obj_list_df <- obj_list_df[complete.cases(obj_list_df[,1]),]
	obj_list_df <- obj_list_df[complete.cases(obj_list_df[,2]),]
	obj_list_df <- obj_list_df[obj_list_df[,2]!="",]
	obj_list_df <- obj_list_df[obj_list_df[,1]!="",]
}

#        if (analyse1 =="") {
#          print("No data in first column !")
#          return(invisible())
#        }
      
  DATA <- obj_list_df

  col2 <- FALSE
  if (!is.null(nm)) {
    if (toupper(substring(nm, nchar(nm)-2))=="PND") col2 <- TRUE
  }

  if (is.factor(DATA[,1])) DATA[,1] <- as.character(DATA[,1])
  if (dim(DATA)[2]==3) {
    if (is.factor(DATA[,3])) DATA[,3] <- as.character(DATA[,3])
  }
  
	      if (col2) {
				  DATA  <- DATA[,1:2]
				} else {
          DATA  <- DATA[,1:min(3, dim(DATA)[2])]
          if (dim(DATA)[2]==3) {
            if (all(is.na(DATA[,3]))) DATA <- DATA[,1:2]
          }
				}

  
      # DATA[,1] j'ai les dates
			
        essai <- unlist(as.data.frame(strsplit(DATA[,1], sep.dates), stringsAsFactors=FALSE))

				es <- as.data.frame(strsplit(essai, dtaspl), stringsAsFactors=FALSE)

				nf <- length(levels(as.factor(as.numeric(es[1,]))))
				nf <- c(nf, length(levels(as.factor(as.numeric(es[2,])))))
				nf <- c(nf, length(levels(as.factor(as.numeric(es[3,])))))
				
				
				
				if (is.null(format)) {
				  
				  y <- "%y"
				  an <- max(c(max(nchar(es[1,])), max(nchar(es[2,])), max(nchar(es[3,]))))
				  if (an==4) y <- "%Y"
				  
				  maxdt <- c(max(as.numeric(es[1,])), 
				             max(as.numeric(es[2,])), 
				             max(as.numeric(es[3,])))
				  
				  ypos <- which(maxdt>31)
				  if (identical(ypos, integer(0))) ypos <- 3
				  
				  dpos <- which(maxdt>12)
				  if (length(dpos)==2) dpos <- dpos[dpos!=ypos]
				  
				  mpos <- c(1:3)[-c(ypos, dpos)]
				
				if (ypos == 1) {format <- y} else {
					if (dpos==1) {format <- "%d"} else {format <- "%m"}}
				format <- paste(format, dtaspl, sep="")
				if (ypos == 2) {format <- paste(format, y, sep="") } else {
					if (dpos==2) {format <- paste(format, "%d", sep="")} else {format <- paste(format, "%m", sep="")}}
				format <- paste(format, dtaspl, sep="")
				if (ypos == 3) {format <- paste(format, y, sep="") } else {
					if (dpos==3) {format <- paste(format, "%d", sep="")} else {format <- paste(format, "%m", sep="")}}				
				} else {
				  # j'ai format. Je dois prendre mpos, dpos et ypos de format
				  dpos <- which(strsplit(format, dtaspl)[[1]]=="%d")
				  mpos <- which(strsplit(format, dtaspl)[[1]]=="%m")
				if (identical(which(strsplit(format, dtaspl)[[1]]=="%y"), integer(0))) {
				  ypos <- which(strsplit(format, dtaspl)[[1]]=="%Y")
				} else {
				  ypos <- which(strsplit(format, dtaspl)[[1]]=="%y")
				}
				}
# j'ai le fichier mais peut-etre 3 colonnes

				if (dim(DATA)[2]==3) {
					DATA[,3] <- ifelse(DATA[,3]=="", NA, DATA[,3])
					DATA[,3] <- na.locf(DATA[,3])
					plage <- levels(as.factor(DATA[,3]))
				} else {
 					if (is.null(nm)) {
						plage <- paste("Roockery", runif(1, 1, 10000))
					} else {
						plage <- basename(nm)
					}
					DATA <- cbind(DATA, Site=rep(plage, dim(DATA)[1]))
				}
				names(DATA) <- c("Date", "Number", "Site")
				
				
				if (is.null(reference)) {
				year <- min(as.numeric(es[ypos,]))
				if (year == 0 & max(as.numeric(es[ypos,])) == 99) year=1999
				if (year<100) {
					if (year<50) {year <- year + 2000} else {year <- year + 1900}
				}
						
				year <- as.character(year)
				
				if (is.null(month_ref)) {
				
				if (nf[ypos] ==2) {
					# 2 annees donc juin
					reference <- as.Date(paste(year, "-07-01", sep=""))
				} else {
					reference <- as.Date(paste(year, "-01-01", sep=""))
				}
				} else {
# si month_ref > premier mois de la série, c'est l'année n-1
				mois <- as.numeric(es[mpos, 1])
				if (mois<month_ref) year <- as.character(as.numeric(year)-1)
				reference <- as.Date(paste(year, "-", month_ref, "-01", sep=""))
				
				}
				
				}
					
				
				return(list(DATA=DATA, reference=reference, format=format))

}
