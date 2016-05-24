#' Grade exams built using rte.grade.exams
#'
#' This function will take as input the information from the exam and grade it using the framework of \pkg{RndTexExams}
#'
#' @param exam.names A character vector with the names of the students, obtained from the test
#' @param exam.version A numeric vector with the version of the exam for each student, obtained from the exam
#' @param exam.answer.matrix A matrix with the answers of the students where the rows represent each student and the columns are the answers to each question
#' @param list.build.rdn.exam A list with several information of the random exams (output from rte.build.rdn.text)
#' @param question.points A numeric vector with the score for each question (default = 1/n.question)
#'
#' @return A list with the following items: \describe{
#' \item{df.grade}{A dataframe with the partial resuts from grading}
#'   \item{df.final.score}{A dataframe with the final results for each student} }
#' @examples
#' # define some options
#' latex.dir.out = 'latexOut' # Name of folder where latex files are going (will create if not exists)
#' pdf.dir.out = 'PdfOut'     # Name of folder where resulting pdf files are going
#' f.out <- 'MyRandomTest_'   # Name of pdfs (MyRandomTest_1.pdf, MyRandomTest_2.pdf, ... )
#' n.test <- 3                # Number of tests to build
#' n.question <- 4            # Number of questions in each test
#'
#' # Get latex example from package
#' f.in <- system.file("extdata", "MyRandomTest.tex", package = "RndTexExams")
#'
#' # Break latex file into a R list
#' list.out <- rte.analize.tex.file(f.in,
#'                                  latex.dir.out = latex.dir.out,
#'                                  pdf.dir.out = pdf.dir.out)
#'
#' # Build pdfs
#' list.build.rdn.exam <- rte.build.rdn.test(list.in = list.out,
#'                                           f.out = f.out,
#'                                           n.test = n.test,
#'                                           n.question = n.question,
#'                                           latex.dir.out = latex.dir.out,
#'                                           pdf.dir.out = pdf.dir.out,
#'                                           do.randomize.questions=TRUE,
#'                                           do.randomize.answers=TRUE,
#'                                           do.clean.up = TRUE)
#'
#' # Grade it!
#' #' # create some (almost) random names
#' my.names <- c('John', 'Max','Marcelo')
#'
#' # version of the test for each student
#' ver.test <- seq(1:length(my.names))
#'
#' # number of simulated questions (same as before)
#' n.questions <- n.question
#'
#' # Get the correct answer sheet from previous code
#' correct.answer.sheet <- list.build.rdn.exam$answer.matrix
#'
#' # create simulated answers from students (cheat a little bit!)
#' q.to.cheat <- 2  # get at least 2 questions right!
#' my.answers <- cbind(correct.answer.sheet[ver.test,1:q.to.cheat],
#'                     matrix(sample(letters[1:5],
#'                                   replace = TRUE,
#'                                   size = length(my.names)*(n.questions-q.to.cheat)),
#'                            ncol = n.questions-q.to.cheat ))
#'
#' # grade exams with rte.grade.exams
#' list.grade <- rte.grade.exams(exam.names = my.names,
#'                               exam.version = ver.test,
#'                               exam.answer.matrix = my.answers,
#'                               list.build.rdn.exam = list.build.rdn.exam)
#'
#' print(list.grade$df.final.score)
#' @export
rte.grade.exams <- function(exam.names,
                            exam.version,
                            exam.answer.matrix,
                            list.build.rdn.exam,
                            question.points = NULL){

  # error checking

  if (length(exam.version)!=length(exam.names)){
    stop('Length of inputs exam.version DO NOT match with the length of exam.names')
  }

  if (nrow(exam.answer.matrix)!=length(exam.names)){
    stop('The number of rows in exam.answer.matrix DO NOT match with the length of exam.names')
  }

  correct.answer.sheet <- list.build.rdn.exam$answer.matrix
  rnd.idx.questions <- list.build.rdn.exam$df.answer.long$rnd.idx.questions
  idx.questions <- list.build.rdn.exam$df.answer.long$n.test
  # check types

  if (class(exam.names)!='character'){
    stop('Input exam.names should be a character class')
  }

  if (!any(class(exam.version)!=c('numeric','integer'))){
    stop('Input exam.version should be a numeric class')
  }

  if (class(exam.answer.matrix)!='matrix'){
    stop('Input exam.answer.matrix should be a matrix')
  }

  if (class(correct.answer.sheet)!='matrix'){
    stop('Input correct.answer.sheet should be a matrix')
  }

  n.question <- ncol(correct.answer.sheet)

  if (is.null(question.points)){
    question.points <- rep(1/n.question ,n.question)
  }

  if (length(question.points)!=n.question){
    stop('Input question.points should have number of elements equal to the number of questions in exam.answer.matrix')
  }


  # check sizes

  size.names <- length(exam.names)
  size.version <- length(exam.version)


  if (size.names!=size.version){
    stop('The number of elements in input exam.names does not match the number of elements in exam.version')
  }

  if (size.names!=nrow(exam.answer.matrix)){
    stop('The number of rows in exam.answer.matrix should match the number of elements in size.names')
  }

  # grade it!


  df.grade<- data.frame()

  for (i.std in seq(1,size.names)){

    name.now <- exam.names[i.std]
    ver.now <- exam.version[i.std]

    answers.now <- exam.answer.matrix[i.std, ]

    # find correct answers by version
    correct.answer.now <- correct.answer.sheet[ver.now,]

    logical.correct <- (answers.now==correct.answer.now)

    names(logical.correct) <- paste0('Q.',seq(1:length(logical.correct)))
    names(logical.correct) <- NULL

    # build temp df for output

    temp.df <- data.frame(exam.names = name.now,
                          exam.ver = ver.now,
                          n.question = seq(1:n.question),
                          rnd.idx.questions = rnd.idx.questions[idx.questions==ver.now],
                          question.score = question.points,
                          grade.logical = logical.correct)

    df.grade <- rbind(df.grade, temp.df)

  }

  df.final.score <- with(df.grade,aggregate(question.score*grade.logical,by = list(exam.names), FUN=sum))

  colnames(df.final.score) <- c('exam.names','final.score')

  df.correction.wide <- data.table::dcast(data = data.table::data.table(df.grade),
                                          formula = exam.names ~ rnd.idx.questions,
                                          fun.aggregate = function(x) return(x), value.var = 'grade.logical', fill = NA )


  grade.l.out <- list(df.grade = df.grade,
                      df.final.score = df.final.score,
                      df.correction.wide = df.correction.wide)

  return(grade.l.out)

}
