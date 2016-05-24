## Author: Ingo Feinerer

## E-mail document
MailDocument <-
function(x = character(0),
         author = character(0),
         datetimestamp = as.POSIXlt(Sys.time(), tz = "GMT"),
         description = character(0),
         header = character(0),
         heading = character(0),
         id = character(0),
         language = character(0),
         origin = character(0),
         ...,
         meta = NULL)
    PlainTextDocument(x, author, datetimestamp, description, heading, id,
                      language, origin, header = header, ..., meta,
                      class = "MailDocument")
