############
## Utility functions for sending and receiving real time messages. 
############

aqInitMessaging <- function(host = "localhost", port = 61618){
  .C("aqInit", host, port, PACKAGE="aqr")
}

aqEnableDebugMessages <- function(){
  .C("aqEnableDebugMessages", PACKAGE="aqr")
}

aqDisableDebugMessages <- function(){
  .C("aqDisableDebugMessages", PACKAGE="aqr")
}

aqSubscribeChannel <- function(channel){
  .C("aqSubscribe", paste("/topic/", channel, sep=""), PACKAGE="aqr")
}

aqUnsubscribeChannel <- function(channel){
  .C("aqUnsubscribe", paste("/topic/", channel, sep=""), PACKAGE="aqr")
}

aqPoll <- function(){
  return(.C("aqPollAll", PACKAGE="aqr"))
}

# waits for data and returns a list of channels for which data is available. 
# this is a synchronous call and thus blocks. 
aqDataReady <- function(){
  return(.C("aqDataReady", PACKAGE="aqr"))
}

aqWaitForData <- function(){
  return(.C("aqWaitForData", PACKAGE="aqr"))
}

aqSend <- function(channel, message){
  .C("aqSend", paste("/topic/", channel, sep=""), message, PACKAGE="aqr")
}

aqTestCallToDynLib<- function(testMessage){
  return(.C("testCall", testMessage, PACKAGE="aqr"))
}
