sendmail = function(recipient, subject="Notification from R",
		message="Calculation finished!", password="rmail")
{
	x = sendmail_(recipient, subject, message, password)
	return(x)
}

sendmail_ = function(recipient, subject, message, password)
{
	# Leerteichen kodiert als <
	# Zeilenumbruch "\n" oder "\r" kodiert als >
	pattern = "[^0-9,a-z,A-Z,.,:,;_!,?,%,<,>,=,@,*]"
	recipient = gsub(" "    ,"" , recipient)
	recipient = gsub(pattern,"_", recipient)
	password  = gsub(pattern,"_", password)
	subject   = gsub("[<]"  , "_", subject)
	subject   = gsub(" "    , "<", subject)
	subject   = gsub(pattern, "_", subject)
	message   = gsub("[<>]", "_", message)
	message   = gsub(" "   , "<", message)
	message   = gsub("\n"  , ">", message)
	message   = gsub("\r"  , ">", message)
	message   = gsub(pattern,"_", message)
	path = paste("http://rmail.linhi.de/rmail.php?Email=",recipient,
			"&Passwort=",password,
			"&Betreff=",subject,
			"&Nachricht=",message,
			sep="")
	x = readLines(path)
	return(x)
}
