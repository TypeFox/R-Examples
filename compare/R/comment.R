
# TODO:
# - make commentQuestions() work with comparisionList
#   as well as comparisonListList

# Convert a comparison result into a set of comments

transformComment <- function(pattern, comment) {
    tr <- list(pattern=pattern, comment=comment)
    class(tr) <- "transformComment"
    tr
}

as.character.transformComment <- function(x, ...) {
    paste("if '", x$pattern, "' then '", x$comment, "'", sep="")
}

print.transformComment <- function(x, ...) {
    cat(as.character(x), "\n")
}
    
comments <- function(answerName, ...) {
    comment <- list(answer=answerName, 
                 transforms=list(...))
    class(comment) <- "markingComments"
    comment
}

print.markingComment <- function(x, ...) {
    cat(paste("Answer ", x$answer, ": ",
              paste(sapply(x$transforms, as.character), collapse=", "),
              sep=""), "\n")
}

# Given a "comparison" and a "tranformComment"
# determine a comment
comparisonComment <- function(tr, comp) {
    if (length(grep(tr$pattern, comp$transform)) > 0)
        tr$comment
    else
        ""
}

# Given a "comparison" and a "markingComment"
# determine comments
comparisonComments <- function(comment, comp) {
    if (length(comment$transforms) == 0) {
        ""
    } else {
        comments <- sapply(comment$transforms, comparisonComment, comp)
        nonblanks <- which(nchar(comments) > 0)
        paste(comments[nonblanks], collapse=", ")
    }
}

questionComments <- function(answerNames,
                             ...) {
    answerNames <- as.character(answerNames)
    comments <- list(...)
    if (!all(sapply(comments, inherits, "markingComments")))
        stop("all '...' arguments must be 'markingComments' objects")
    quest <- list(names=answerNames,
                  comments=comments)
    class(quest) <- "questionComments"
    quest
}
                      
# Given a "comparisonList" and the names of answer objects
# and a list of comments ...
# for each answer, determine all relevant comments
commentQuestion <- function(result, question) {
    UseMethod("commentQuestion")
}

growComments <- function(comments, newComments) {
    if (length(comments) == 0) {
        newComments
    } else {
        if (nchar(comments) == 0) {
            newComments
        } else {
            paste(comments, newComments, sep=", ")
        }
    }
}

commentQuestion.comparisonList <- function(result, question) {
    comments <- character()
    for (i in question$names) {
        comp <- result[[i]]
        # If there is a comment for this question ...
        comment <- grep(i, question$comments)
        if (length(comment) > 0) {
            newcomments <- comparisonComments(question$comments[[comment]],
                                              comp)
            comments <- growComments(comments, newcomments)
        }
    }
    comments
}

commentQuestion.comparisonListList <- function(result, question) {
    sapply(result,
           function(x) {
               commentQuestion(x, question)
           })
}

commentQuestions <- function(result, ...) {
    questions <- list(...)
    if (!all(sapply(questions, inherits, "questionComments")))
        stop("all '...' arguments must be 'questionComments' objects")
    comments <- do.call("cbind",
                        lapply(questions,
                               function(q) { commentQuestion(result, q) }))
    colnames(comments) <- sapply(questions,
                                 function(q) { paste(q$names, collapse="-") })
    comments
}

