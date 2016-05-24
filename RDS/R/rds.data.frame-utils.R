


#was get.recruiter.id
#' Determines the recruiter.id from recruitment coupon information
#' @param data a data.frame
#' @param subject.coupon The variable representing the coupon returned by subject
#' @param coupon.variables The variable representing the coupon ids given to the subject
#' @param subject.id The variable representing the subject's id
#' @param seed.id The recruiter.id to assign to seed subjects.
#' @export
#' @examples
#' fpath <- system.file("extdata", "nyjazz.csv", package="RDS")
#' dat <- read.csv(fpath)
#' dat$recruiter.id <- rid.from.coupons(dat,"own.coupon",
#'                       paste0("coupon.",1:7),"id")
#' 
#' #create and rds.data.frame
#' rds <- as.rds.data.frame(dat,network.size="network.size")
rid.from.coupons <- function(data, subject.coupon=NULL, coupon.variables, 
		subject.id=NULL,seed.id="seed"){
	max.coupons <- length(coupon.variables)
	column.names <- colnames(data)
	named.columns <- column.names[column.names != "" && !is.na(column.names)]
	if(is.null(subject.id))
		subject.id <- as.char(1:nrow(data))
	else
		subject.id <- as.char(data[[subject.id]])
	if(any(duplicated(subject.id)))
		stop("subject.id is not a unique identifier")
	rcv <- data[,coupon.variables, drop=FALSE]
	for(i in 1:ncol(rcv)){
		rcv[[coupon.variables[i]]] <- as.char(rcv[[coupon.variables[i]]])
	}
	subject.coupon <- as.char(data[[subject.coupon]])
	seed.code <- seed.id
	i <- 0
	while(seed.code %in% subject.id){
		seed.code <- paste0("seed",i)
		i <- i+1
	}
	
	recruiter.id <- rep(as.character(NA),nrow(data))
	for (row in 1:nrow(data)) {
		coupons.given <- rcv[row, ]
		recruitees <- which(subject.coupon %in% coupons.given[!is.na(coupons.given)])
		recruiter.id[recruitees] <- subject.id[row]
	}
	recruiter.id[is.na(recruiter.id)] <- seed.code
	recruiter.id
}


#' Calculates the root seed id for each node of the recruitement tree.
#' @param data An rds.data.frame
#' @export
#' @examples
#' data(fauxmadrona)
#' seeds <- get.seed.id(fauxmadrona)
#' #number recruited by each seed
#' barplot(table(seeds))
get.seed.id <- function(data){
	if(!is.rds.data.frame(data))
		stop("data must be an rds.data.frame")
	id <- get.id(data)
	recruiter.id <- get.rid(data)
	sid <- get.seed.rid(data)
	get.seed <- function(i, history) {
		row <- match(i, id)
		rec.id <- recruiter.id[row]
		if(rec.id==i){
			stop(sprintf("Yikes! The data says that the person with id %s recruited themselves! Please check that the coupon information in the data for that person is correct :-)",i),call.=FALSE)}
		if(rec.id %in% history){
			stop("Loop found in recruitment tree.")
		}
		#print(sprintf("i %s rec.id %s ",i,rec.id))
		if (rec.id == sid) {
			return(i)
		}
		else {
			get.seed(rec.id,history=c(history,i))
		}
	}
	seed <- sapply(id, get.seed,history=c())
	seed
}


#' Calculates the depth of the recruitment tree (i.e. the recruitment wave) 
#' at each node.
#' @param data An rds.data.frame
#' @export
#' @examples
#' data(fauxmadrona)
#' #number subjects in each wave
#' w <- get.wave(fauxmadrona)
#' #number recruited in each wave
#' barplot(table(w))
get.wave <- function(data){
	if(!is.rds.data.frame(data))
		stop("data must be an rds.data.frame")
	id <- get.id(data)
	recruiter.id <- get.rid(data)
	sid <- get.seed.rid(data)
	get.w <- function(i, history) {
		row <- match(i, id)
		rec.id <- recruiter.id[row]
		if(rec.id==i){
			stop(sprintf("Yikes! The data says that the person with id %s recruited themselves! Please check that the coupon information in the data for that person is correct :-)",i),call.=FALSE)}
		if(rec.id %in% history){
			stop("Loop found in recruitment tree.")
		}
		if (rec.id == sid) {
			return(length(history))
		}
		else {
			get.w(rec.id,history=c(history,i))
		}
	}
	seed <- sapply(id, get.w,history=c())
	seed
}



#' Calculates the number of (direct) recuits for each respondent.
#' @param data An rds.data.frame
#' @export
#' @examples
#' data(fauxmadrona)
#' nr <- get.number.of.recruits(fauxmadrona)
#' #frequency of number recruited by each id
#' barplot(table(nr))
get.number.of.recruits <- function(data)
    {
    if(!is(data,"rds.data.frame"))
              stop("data must be of type rds.data.frame")

    # Make sure that we can compute the number of recruits for each respondent's respondents.
  
    number.of.recruits <- sapply(as.character(data[[attr(data,"id")]]),function(i){length(which(as.character(data[[attr(data,"recruiter.id")]]) == i))})
    as.integer(number.of.recruits)
    }


