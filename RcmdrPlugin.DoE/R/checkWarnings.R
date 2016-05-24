checkWarnings <- function (messages)
{
  ### copied from Rcmdr version 2.0-0
    if (getRcmdr("suppress.X11.warnings")) {
        X11.warning <- grep("X11 protocol error|Warning in structure",
            messages)
        if (length(X11.warning) > 0) {
            messages <- messages[-X11.warning]
        }
        if (length(messages) == 0)
            Message()
        else if (length(messages) > 10) {
            messages <- c(paste(length(messages), "warnings."),
                gettextRcmdr("First and last 5 warnings:"), head(messages,
                  5), ". . .", tail(messages, 5))
            Message(message = paste(messages, collapse = "\n"),
                type = "warning")
        }
        else {
            if (length(grep("warning", messages, ignore.case = TRUE)) >
                0)
                Message(message = paste(messages, collapse = "\n"),
                  type = "warning")
            else Message(message = paste(messages, collapse = "\n"),
                type = "note")
        }
    }
    else {
        if (length(messages) == 0)
            Message()
        else if (length(messages) > 10) {
            messages <- c(paste(length(messages), "warnings."),
                gettextRcmdr("First and last 5 warnings:"), head(messages,
                  5), ". . .", tail(messages, 5))
            Message(message = paste(messages, collapse = "\n"),
                type = "warning")
        }
        else {
            if (length(grep("warning", messages, ignore.case = TRUE)) >
                0)
                Message(message = paste(messages, collapse = "\n"),
                  type = "warning")
            else Message(message = paste(messages, collapse = "\n"),
                type = "note")
        }
    }
    tkfocus(CommanderWindow())
}
