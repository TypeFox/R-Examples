
# TODO:
# - make markQuestions() work with comparisionList
#   as well as comparisonListList

# Convert a comparison result into a mark

transformRule <- function(pattern, mark) {
    tr <- list(pattern=pattern, mark=mark)
    class(tr) <- "transformRule"
    tr
}

as.character.transformRule <- function(x, ...) {
    paste("'", x$pattern, "' deducts ", x$mark, sep="")
}

print.transformRule <- function(x, ...) {
    cat(as.character(x), "\n")
}
    
rule <- function(answerName, falseMark, ...) {
    rule <- list(answer=answerName, false=falseMark,
                 transforms=list(...))
    class(rule) <- "markingRule"
    rule
}

print.markingRule <- function(x, ...) {
    cat(paste("Answer ", x$answer, ": FALSE deducts ",
              x$false, ", ",
              paste(sapply(x$transforms, as.character), collapse=", "),
              sep=""), "\n")
}

# Given a "comparison" and a "tranformRule"
# determine a deduction
transformDeduction <- function(tr, comp) {
    if (length(grep(tr$pattern, comp$transform)) > 0)
        -tr$mark
    else
        0
}

# Given a "comparison" and a "markingRule"
# determine a deduction
ruleDeduction <- function(rule, comp, name) {
    if (rule$answer == name) {
        if (comp$result) {
            if (length(rule$transforms) == 0) {
                0
            } else {
                sum(sapply(rule$transforms, transformDeduction, comp))
            }
        } else {
            -rule$false
        }
    } else {
        0
    }
}

questionMarks <- function(answerNames,
                          maxMark,
                          ...) {
    answerNames <- as.character(answerNames)
    rules <- list(...)
    if (!all(sapply(rules, inherits, "markingRule")))
        stop("all '...' arguments must be 'markingRule' objects")
    mark <- as.numeric(maxMark)
    quest <- list(names=answerNames,
                  mark=maxMark,
                  rules=rules)
    class(quest) <- "questionMarks"
    quest
}
                      
# Given a "comparisonList" and the names of answer objects
# and the maximum mark, and a list of rules ...
# for each answer, apply all relevant rules, deducting
# appropriate mark (until get to 0)
markQuestion <- function(result, question) {
    UseMethod("markQuestion")
}

markQuestion.comparisonList <- function(result, question) {
    mark <- as.numeric(question$mark)
    for (i in question$names) {
        comp <- result[[i]]
        deductions <- sapply(question$rules, ruleDeduction,
                             comp, i)
        mark <- mark + sum(deductions)
    }
    if (mark <= 0)
        mark <- 0
    mark
}

markQuestion.comparisonListList <- function(result, question) {
    sapply(result,
           function(x) {
               if (is.null(x))
                   0
               else
                   markQuestion(x, question)
           })
}

markQuestions <- function(result, ...) {
    questions <- list(...)
    if (!all(sapply(questions, inherits, "questionMarks")))
        stop("all '...' arguments must be 'questionMarks' objects")
    marks <- do.call("cbind",
                     lapply(questions,
                            function(q) { markQuestion(result, q) }))
    colnames(marks) <- sapply(questions,
                              function(q) { paste(q$names, collapse="-") })
    marks
}

# Tests:
function() {
    transformRule("coerced", 1)
    rule("ageSummary", 2,
         transformRule("coerced", 1),
         transformRule("rounded", 1))
    rule("ageSummary", 2,
         transformRule("coerced|rounded", 1))
    comp <- compare(letters, factor(letters))
    transformDeduction(transformRule("coerced", 1), comp)
    ruleDeduction(rule("IndianMothers", 1),
                  comp)
    ruleDeduction(rule("ageSummary", 2,
                       transformRule("coerced|rounded", 1)),
                  comp)
    comp <- compare(letters, factor(letters), allowAll=TRUE)
    transformDeduction(transformRule("coerced", 1), comp)
    ruleDeduction(rule("ageSummary", 2,
                       transformRule("coerced|rounded", 1)),
                  comp)
    ruleDeduction(rule("ageSummary", 2,
                       transformRule("coerced", 1),
                       transformRule("rounded", 1)),
                  comp)
    questionMarks(answerNames=c("id", "age", "edu", "class"),
                  maxMark=1,
                  rule("id", 1),
                  rule("age", 1),
                  rule("edu", 1),
                  rule("class", 1))
}
