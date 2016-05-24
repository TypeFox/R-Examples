## Author: Ingo Feinerer

readMail <- structure(
function(DateFormat = "%d %B %Y %H:%M:%S")
{
    stopifnot(is.character(DateFormat))

    format <- DateFormat
    function(elem, language, id) {
        mail <- elem$content

        ## The header is separated from the body by a blank line.
        ## http://en.wikipedia.org/wiki/E-mail#Internet_e-mail_format
        for (index in seq_along(mail))
            if (identical(mail[index], "")) break

        header <- mail[1:index]
        content <- mail[(index + 1):length(mail)]

        author <- gsub("From: ", "",
                       grep("^From:", header, value = TRUE, useBytes = TRUE))
        datetimestamp <-
            strptime(gsub("Date: ", "",
                          grep("^Date:", header, value = TRUE, useBytes = TRUE)),
                     format = format,
                     tz = "GMT")
        mid <- gsub("Message-ID: ", "",
                    grep("^Message-ID:", header, value = TRUE, useBytes = TRUE))
        origin <- gsub("Newsgroups: ", "",
                       grep("^Newsgroups:", header, value = TRUE, useBytes = TRUE))
        heading <- gsub("Subject: ", "",
                        grep("^Subject:", header, value = TRUE, useBytes = TRUE))

        MailDocument(content, author, datetimestamp, character(0), header, heading,
                     if (length(mid)) mid[1] else id, language, origin)
    }
}, class = c("FunctionGenerator", "function"))
