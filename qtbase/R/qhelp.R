
getAssistant <- local({
  assistants <- list()
  function(package) {
    assistant <- assistants[[package]]
    if (is.null(assistant) || !isOpen(assistant)) {
      cmd <- paste("assistant", "-collectionFile",
                   system.file("qhelp", paste(package, "qhc", sep = "."),
                               package),
                   "-enableRemoteControl")
      assistants[[package]] <<- pipe(cmd, "wb")
    }
    assistants[[package]]
  }
})

activateAssistantKeyword <- function(assistant, keyword) {
  writeBin(paste("activateKeyword", keyword), assistant)
}

qhelp <- function(topic, package = "qtbase") {
  if (is.name(substitute(topic)))
    topic <- as.character(substitute(topic))
  assistant <- getAssistant(package)
  activateAssistantKeyword(assistant, topic)
}
