message<-function(messageType,names=""){
	messageList<-list(
			INFO_WELCOME='',
			FILE_CREDENTIALS_NO_EXISTS='Function util_setCredentials() should be called before. To retrieve your account credentials, please log in to your account or register for a free account at https://www.portfolioeffect.com/registration.',
			WRONG_SETTINGS_ARGUMENTS='',
			NOT_ESTIMATOR_CLASS='The object must have the class "estimator". Use the function estimator_create to create it.',
			OBJECT_NOT_CHARACTER_CLASS=paste('The',names,'must have the class "character".'),
			OBJECT_NOT_NUMERIC_CLASS=paste('The',names,'must have the class "numeric".'),
			OBJECT_NOT_SINGLE_NUMBER=paste('The',names,'must be single number.'),
			OBJECT_NOT_POSITIVE_NUMBER=paste('The',names,'must be a positive number.')
	)
	return(messageList[[messageType]])
}

stopMessage<-function(messageType,names=""){
	stop(message(messageType,names))
}
infoMessage<-function(messageType,names=""){
	cat(message(messageType,names))
}