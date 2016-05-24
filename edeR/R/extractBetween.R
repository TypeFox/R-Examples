extractBetween <-
function(username,password,folder,startDate,endDate,nmail=-1)
{
 #...Initializing "imaplib" from python
 imap<-rJython(modules="imaplib")
 imap$exec("import imaplib")

 #...Authorized login using username and passowrd
 imap$exec(paste("user_email = imaplib.IMAP4_SSL('","imap.gmail.com", "')", sep = ""))
 imap$exec(paste("user_email.login(\'",username, "\',\'",password, "\')", sep = ""))

 #...Extracting unique email message id from the specified date range
 imap$exec(paste("nmsg = user_email.select(\'",folder,"\')",sep=""))
 imap$exec(paste("stts,msgid = user_email.uid('search', None, '(SINCE ",startDate, " BEFORE ",endDate,")')",sep=""))
 imap$exec(paste("msgid=msgid[0].split()",sep=""))
 msg_uid <- .jstrVal(imap$get("msgid"))
 msg_uid <- gsub(" ","",strsplit(gsub("'","",gsub("\\]","",gsub("\\[","",msg_uid))),",")[[1]])

 #...Extracting number of available email in the specified "folder"
 nmsg <- length(msg_uid)


 if(nmail<0 | nmail>nmsg)
  nmail<-nmsg
 #...Initialize email data frame
 email_data <- NULL

 #...Setting counter so that we can control how many message to extract
 counter<-1
 for(uid in msg_uid)
 {
  if(counter<=nmail)
  {
   #...Extracting the email address where the email come from
   imap$exec(paste("stts,addr_frm=user_email.uid('fetch',\'",uid,"\','(BODY[HEADER.FIELDS (FROM)])')",sep=""))
   imap$exec("addr_frm = addr_frm[0][1]")
   from <- .jstrVal(imap$get("addr_frm"))
   from <- unlist(strsplit(from, "[<>\r\n, \"]"))
   from <- sub("from: ", "", from, ignore.case = TRUE)
   from <- grep("@", from, value = TRUE)

   #...Extracting address where the email send to (there could be more than one receipients)
   imap$exec(paste("stts,addr_to=user_email.uid('fetch',\'",uid,"\','(BODY[HEADER.FIELDS (TO)])')",sep=""))
   imap$exec("addr_to = addr_to[0][1]")
   addrs_to <- .jstrVal(imap$get("addr_to"))
   addrs_to <- unlist(strsplit(addrs_to, "[<>\r\n, \"]"))
   addrs_to <- sub("addrs_to: ", "", addrs_to, ignore.case = TRUE)
   addrs_to <- grep("@", addrs_to, value = TRUE)

   #...Extracting email subject line
   imap$exec(paste("stts,subj=user_email.uid('fetch',\'",uid,"\','(BODY[HEADER.FIELDS (SUBJECT)])')",sep=""))
   imap$exec("subj = subj[0][1]")
   subj <- .jstrVal(imap$get("subj"))
   subj <- unlist(strsplit(subj,"\r"))[1]

   #...Extracting timestamp
   imap$exec(paste("stts,timestamp=user_email.uid('fetch',\'",uid,"\','(BODY[HEADER.FIELDS (DATE)])')",sep=""))
   imap$exec("timestamp = timestamp[0][1]")
   timestamp <- .jstrVal(imap$get("timestamp"))
   timestamp <- unlist(strsplit(unlist(strsplit(timestamp,"Date:"))[2],"\r"))[1]

   #...Producing dataframe using extracted information
   id <- rep(uid,length(addrs_to))
   from <- rep(from,length(addrs_to))
   subj <- rep(subj,length(addrs_to))
   timestamp <- rep(timestamp,length(addrs_to))
   email_data_temp <- data.frame(email_uid=id,email_from=from,email_to=addrs_to,email_subj=subj,email_date=timestamp,stringsAsFactors=FALSE)
   email_data <- rbind(email_data,email_data_temp)

   rm(email_data_temp); rm(from); rm(addrs_to); rm(timestamp); rm(subj);rm(id)
   counter<-counter+1
  }
  else
   break
 }
 #...Logging out from the account
 imap$exec("user_email.shutdown()")
 outlist <- list(n_message=nmsg,data=email_data)
 return(outlist)
}
