convert_mbox_eml <- function(mbox, dir) {
    dir.create(dir, recursive = TRUE)
    content <- readLines(mbox)
    counter <- start <- end <- 1
    nMessages <- length(grep("^From ", content, useBytes = TRUE))
    fmt <- paste("%0", nchar(nMessages), "d", sep = "")
    needWrite <- FALSE
    for (i in seq_along(content)) {
        if (length(grep("^From ", content[i], useBytes = TRUE)) > 0) {
            end <- i - 1
            if (needWrite && start <= end) {
                con <- file(file.path(dir, sprintf(fmt, counter)))
                writeLines(content[start:end], con)
                close(con)
                needWrite <- FALSE
                counter <- counter + 1
            }
            start <- i
            needWrite <- TRUE
        }
    }
    if (needWrite) {
        con <- file(file.path(dir, sprintf(fmt, counter)))
        writeLines(content[start:length(content)], con)
        close(con)
    }
    invisible(TRUE)
}

# Remove e-mail citations beginning with >
removeCitation <- function(x, removeQuoteHeader) UseMethod("removeCitation", x)
removeCitation.character <- function(x, removeQuoteHeader = FALSE) {
    citations <- grep("^[[:blank:]]*>", x, useBytes = TRUE)
    if (removeQuoteHeader) {
      headers <- grep("wrote:$|writes:$", x)
      # A quotation header must immediately preceed a quoted part,
      # possibly with one empty line in between
      headers <- union(headers[(headers + 1) %in% citations],
                       headers[(headers + 2) %in% citations])
      citations <- union(headers, citations)
    }
    if (length(citations)) x[-citations] else x
}
removeCitation.MailDocument <-
    content_transformer(removeCitation.character)

# Remove non-text parts from multipart e-mail messages
removeMultipart <- function(x) UseMethod("removeMultipart", x)
removeMultipart.character <- function(x) {
    # http://en.wikipedia.org/wiki/Multipart_message#Multipart_Messages
    # We are only interested in text/plain parts
    i <- grep("^Content-Type: text/plain", x, useBytes = TRUE)
    r <- character(0)
    k <- 2
    for (j in i) {
        end <- if (k <= length(i)) i[k]-1 else length(x)
        content <- x[j:end]
        ## Find boundary (starting with "--")
        # In most cases the boundary is just one line before the Content-Type header
        start <- j - 1
        while (j > 0) {
            if (substr(x[j], 1, 2) == "--") {
                start <- j
                break
            }
            else
                j <- j - 1
        }
        index <- grep(x[start], content, useBytes = TRUE, fixed = TRUE)
        index <- if (length(index) == 0) length(content) else (index[1] - 1)
        content <- content[1:index]
        # Now remove remaining headers
        index <- grep("^$", content, useBytes = TRUE)
        index <- if (length(index) == 0) 1 else (index[1] + 1)
        r <- c(r, content[index:length(content)])
        k <- k + 1
    }

    if (length(r) == 0) x else r
}
removeMultipart.MailDocument <-
    content_transformer(removeMultipart.character)

# Remove e-mail signatures
removeSignature <- function(x, marks) UseMethod("removeSignature", x)
removeSignature.character <- function(x, marks = character(0)) {
    # "---" is often added to Sourceforge mails
    # "___" and "***" are also common, i.e.,
    # marks <- c("^_{10,}", "^-{10,}", "^[*]{10,}")

    # "-- " is the official correct signature start mark
    marks <- c("^-- $", marks)

    signatureStart <- length(x) + 1
    for (m in marks)
        signatureStart <- min(grep(m, x, useBytes = TRUE), signatureStart)

    if (signatureStart <= length(x)) x[-(signatureStart:length(x))]
    else x
}
removeSignature.MailDocument <-
    content_transformer(removeSignature.character)

get.thread.id <- function(parentID, ht) {
    threadID <- NA
    threadLevel <- 2

    if (!identical(parentID, "") && length(parentID) != 0 && is.numeric(ht[[parentID]][1]))
        threadID <- as.integer(ht[[parentID]][1])

    if (!identical(parentID, "") && length(parentID) != 0 && is.numeric(ht[[parentID]][2]))
        threadLevel <- as.integer(ht[[parentID]][2] + 1)

    list(threadID = threadID, threadLevel = threadLevel)
}

# Compute thread IDs and thread depth (based on heuristics)
# See http://www.jwz.org/doc/threading.html for a more complete and correct approach
threads <- function(x)
{
    # Hash table storing (thread ID, thread level) for a given message ID
    ht <- new.env()
    tid <- 1
    threadIDs <- threadLevels <- integer(length(x))
    for (i in seq_along(x)) {
        messageID <- meta(x[[i]], "id")
        parentID <- gsub("In-Reply-To: ", "", grep("^In-Reply-To:", attr(x[[i]], "Header"), value = TRUE, useBytes = TRUE))
        refID <- gsub("References: ", "", grep("^References:", attr(x[[i]], "Header"), value = TRUE, useBytes = TRUE))
        refID <- sub(",$", "", refID)

        # Generate new thread
        if (!length(parentID) & !length(refID)) {
            ht[[messageID]] <- c(tid, 1)
            threadIDs[i] <- tid
            threadLevels[i] <- 1
            tid <- tid + 1
        }
        # Use existing thread
        else {
            info <- get.thread.id(refID, ht)
            if (!is.na(info$threadID)) {
                ht[[messageID]] <- c(info$threadID, info$threadLevel)
                threadIDs[i] <- info$threadID
                threadLevels[i] <- info$threadLevel
                next
            }

            parentID <- unique(parentID)
            if (length(parentID) > 1) {
                for (i in 1:length(parentID)) {
                    info <- get.thread.id(parentID[[i]], ht)
                    if (!is.na(info$threadID))
                        next
                }
            } else
                info <- get.thread.id(parentID, ht)

            if (is.na(info$threadID)) {
                # The message is a reply to some other message, but the parent
                # e-mail is not available --- for instance, it could be a follow-up
                # to a discussion from a separate mailing list. Create a new thread.
                ht[[messageID]] <- c(tid, 1)
                threadIDs[i] <- tid
                threadLevels[i] <- 1
                tid <- tid + 1
            }
            else {
                ht[[messageID]] <- c(info$threadID, info$threadLevel)
                threadIDs[i] <- info$threadID
                threadLevels[i] <- info$threadLevel
            }
        }
    }
    list(ThreadID = threadIDs, ThreadDepth = threadLevels)
}
