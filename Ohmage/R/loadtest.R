enquote <- function(string){
  return(paste('"', string, '"', sep=""));
}

setXmlValues <- function(xmlfile, n.users, prefix="loadtest"){
	xmldoc <- XML::xmlParse(xmlfile);
	urnnode <- XML::getNodeSet(xmldoc, "/campaign/campaignUrn")[[1]];
	namenode <- XML::getNodeSet(xmldoc, "/campaign/campaignName")[[1]];
	oldname <- XML::xmlValue(namenode);
	XML::xmlValue(urnnode) <- paste("urn:campaign", prefix, oldname, n.users, sep=":");
	XML::xmlValue(namenode) <- paste(prefix, oldname, n.users, sep=".");
	newfile <- tempfile();
	XML::saveXML(xmldoc, newfile);
	attr(newfile, "oldname") <- oldname;
	return(newfile);
}

stripspace <-function(x) {
	sub("^\\s*([^ ]*.*[^ ])\\s*$", "\\1",x);
}

stripuuid <- function(datastring){
	newstring <- "12345678-1234-1234-1234-1234567890ab"
	myexpr <- "[a-f0-9]{8}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{12}"
	hasphoto <- length(grep(myexpr, datastring)>0);
	if(hasphoto){
		newdatastring <- gsub(myexpr, newstring, datastring);
	} else {
		newdatastring <- datastring;
	}
	attr(newdatastring, "hasphoto") <- hasphoto;
	return(newdatastring);
}

getuuids <- function(datastring, prefix){
	myexpr1 <- paste("\"", prefix, "\":\"[a-f0-9]{8}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{12}\"", sep="");
	myexpr2 <- "[a-f0-9]{8}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{12}"
	photomatches <- (regmatches(datastring, gregexpr(myexpr1, datastring))[[1]]);
	return(unlist(regmatches(photomatches, gregexpr(myexpr2, photomatches))))
}

makeuuidparams <- function(datastring, testimage, prefix){
	uuids <- getuuids(datastring, prefix)
	outlist <- list()
	for(i in uuids){
		outlist[[i]] <- testimage;
	}
	return(outlist);
}


numtolet <- function(number, numlength=6){
	num2let <- function(number) return(letters[number])
	number <- as.integer(number);
	number <- formatC(number, digits=0, width=numlength, flag=0, format="f");
	allindexes <- lapply(strsplit(number,""), as.numeric);
	allindexes <- lapply(allindexes, "+", 1)
	strings <- sapply(lapply(allindexes, num2let), paste, collapse="");
	return(strings);
}

#' Data generations for loadtesting
#' @param n.users number of users to generate
#' @param n.days number of days per user
#' @param n.responses number of responses per day
#' @param xmlfile xml file or string
#' @param recycle if https connection should be kept alive (recommended)
#' @param verbose verbose output
#' @param shareall should data be shared? T/F
#' @param user.prefix prefix used for campaign and user names.
#' @export
loadtest <- function(n.users = 10, n.days=5, n.responses=2, xmlfile = system.file(package="Ohmage", "files/jeroen.xml"), recycle=TRUE, verbose=FALSE, shareall=TRUE, user.prefix="loadtest"){

	#load xml package
	library(XML);

	#check for java
	if(system('java', ignore.stdout=TRUE, ignore.stderr=TRUE) != 0){
    stop("java was not found");
  }

	#statics
	class_urn <- "urn:class:loadtest";
	password <- "Test.123";
	datagenjar <- system.file(package="Ohmage", "files/andwellness-survey-generator-2.11.jar");
	ohmage_username <- getOption("ohmage_username")

	#jsonfile <- gsub("Rtmp.*", "datagen.json", tempdir()); # "/tmp/datagen.json"
	jsonfile <- path.expand("~/.datagen.json");

	#remove old file
	if(file.exists(jsonfile)){
		if(!file.remove(jsonfile)){
			stop("Could not write edit tempfile. Please manually delete: ", jsonfile);
		}
	}

	#new xml file
	myxmlfile <- tempfile();
	if(isTRUE(regexpr("<?xml", xmlfile, fixed=T))){
		#write to file.
		write(xmlfile, file=myxmlfile);
	} else if(file.exists(xmlfile)) {
		#copy existing xml file
		file.copy(xmlfile, myxmlfile);
	} else {
		stop("xmlfile argument needs to be xml string or a filename.")
	}

	myxmlfile <- setXmlValues(myxmlfile, n.users, user.prefix);

	#parse xml
	oldname <- attr(myxmlfile, "oldname")
	xmldoc <- XML::xmlParse(myxmlfile);
	campaignUrn <- stripspace(XML::xmlValue(XML::getNodeSet(xmldoc, "/campaign/campaignUrn")[[1]]));

	#try to create class (might already exist);
	mytry <- try(oh.class.create(class_urn, "Loadtestclass", recycle=recycle), silent=T);
	if(class(mytry) == "try-error") message(stripspace(paste(strsplit(mytry[1],":")[[1]][-1], collapse=":")));
	oh.class.update(class_urn, user_role_list_add=paste(ohmage_username,";privileged", sep=""));

	#create campaign
	oh.campaign.create(xml=paste(readLines(myxmlfile), collapse="\n"), class_urn_list=class_urn);
	creationtime <- as.character(oh.campaign.read(output_format="short")[campaignUrn,"creation_timestamp"]);

	#generate usernames.
	usernames <- paste(user.prefix,".", substring(oldname, 1,5), "." ,n.users, ".",numtolet(0:(n.users-1), ceiling(log(n.users,10))), sep="");

	#create users.
	for(thisuser in usernames){
		oh.user.create(thisuser, password, recycle=recycle)
	}

	#add them to a class
	user.role.list <- paste(usernames, ";restricted", collapse=",", sep="");
	message("Adding users to class...")
	oh.class.update(class_urn, user_role_list_add=user.role.list);

	#build the system command
	command <- paste("java -jar", enquote(datagenjar), enquote(myxmlfile), n.days, n.responses, enquote(jsonfile), "upload");
	testimage <- fileUpload(system.file("files/lolcat.jpg", package="Ohmage"), contentType="image/jpg")

	#for all users...
	for(thisuser in usernames){

		#generate some data
		unlink(jsonfile);
		system(command, intern=TRUE);
		jsondata <- readChar(jsonfile, file.info(jsonfile)$size);
		uuidparams <- makeuuidparams(jsondata, testimage, prefix="value");

		#get hashed passwd
		hashedpass <- oh.user.auth(thisuser, "Test.123", recycle=recycle);

		#call upload function
		surveyargs <- list(campaign_urn=campaignUrn, user=thisuser, password=hashedpass, campaign_creation_timestamp=creationtime,
				surveys=jsondata, recycle=recycle, verbose=verbose);
		do.call("oh.survey.upload", c(surveyargs, uuidparams));

		#share all data
		if(shareall==TRUE){
			surveykeys <- getuuids(jsondata, prefix="survey_key");
			for(thiskey in surveykeys){
				oh.survey_response.update(campaignUrn, thiskey);
			}
		}

	}

	message("Successfully generated data for ", n.users, " users, ", n.days, " days, and ", n.responses, " responses.\n\n");
}


#' Delete campaign generated by loadtest function
#' @param n.users number of users that was generated
#' @param xmlfile path or string of orriginal xml
#' @param recycle to recycle the connection
#' @param user.prefix a character string used in the usernames and campaign name. Make it short.
#' @export
loadtest.wipe <- function(n.users = 10, xmlfile = system.file(package="Ohmage", "files/jeroen.xml"), recycle=TRUE, user.prefix = "loadtest"){

	#load xml package
	library(XML);

	#new xml file
	myxmlfile <- tempfile();
	if(isTRUE(regexpr("<?xml", xmlfile, fixed=T))){
		#write to file.
		write(xmlfile, file=myxmlfile);
	} else if(file.exists(xmlfile)) {
		#copy existing xml file
		file.copy(xmlfile, myxmlfile);
	} else {
		stop("xmlfile argument needs to be xml string or a filename.")
	}

	myxmlfile <- setXmlValues(myxmlfile, n.users, user.prefix);

	#parse xml
	oldname <- attr(myxmlfile, "oldname")
	xmldoc <- XML::xmlParse(myxmlfile);
	campaignUrn <- stripspace(XML::xmlValue(XML::getNodeSet(xmldoc, "/campaign/campaignUrn")[[1]]));

	#usernames
	usernames <- paste(user.prefix,".", substring(oldname, 1,5), "." ,n.users, ".",numtolet(0:(n.users-1), ceiling(log(n.users,10))), sep="");

	#remove all
	message("trying to wipe: ", campaignUrn);
	try(oh.campaign.delete(campaignUrn));

	for(thisuser in usernames){
		message("Deleting user: ", thisuser);
		try(oh.user.delete(thisuser, recycle=recycle));
	}
}
