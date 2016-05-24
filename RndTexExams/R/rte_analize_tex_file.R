#' Analize a LaTeX file and convert it into a list
#'
#' This function will take as input a LaTeX file and break its components into a
#' single R List.
#'
#' @param f.in The latex file with the exam
#' @param latex.dir.out The name of the folder where the files from the latex
#'   compilation should go (will create if not found)
#' @param pdf.dir.out The name of the folder where the pdf from the latex
#'   compilation should go (will create if not found)
#' @return A list that represents the tex file with preamble, questions, answers
#'   and more. This list is later used by function rte.build.rdn.text \describe{
#'   \item{df.questions}{A data.frame with all questions} \item{df.answers}{A
#'   data.frame with all answers} \item{my.begin.mchoice.line}{text with
#'   beggining of mchoice enviroment} \item{my.preamble}{preamble of tex file,
#'   including everything before the beggining of the multiple choice enviroment
#'   } \item{my.last.part}{All of the tex code after the end of the multiple
#'   choice enviroment } }
#' @examples
#' latex.dir.out <- 'latexOut' # Name of folder where latex files are going
#'                             #(will create if it does not exists)
#'
#' pdf.dir.out <- 'PdfOut'     # Name of folder where resulting pdf files are going
#'
#' # Get latex example from package
#' f.in <- system.file("extdata", "MyRandomTest.tex", package = "RndTexExams")
#'
#' # Break latex file into a R list
#' list.out <- rte.analize.tex.file(f.in,
#'                                 latex.dir.out = latex.dir.out,
#'                                 pdf.dir.out = pdf.dir.out)
#'
#' print(list.out)
#' @export
rte.analize.tex.file <- function(f.in,
                                 latex.dir.out = 'latexOut',
                                 pdf.dir.out = 'PdfOut'){

  cat('\nrte: Changing LaTeX file into dataframe...')

  #require(stringr)
  #require(stringi)

  # create folder for latex output

  if (!file.exists(latex.dir.out)) dir.create(latex.dir.out)

  # create folder for pdf output

  if (!file.exists(pdf.dir.out)) dir.create(pdf.dir.out)

  # clean folders

  my.temp.files <- dir(latex.dir.out, pattern = '*.*', full.names = T)
  if (length(my.temp.files)!=0) file.remove(my.temp.files)

  my.temp.files <- dir(pdf.dir.out, pattern = '*.*', full.names = T)
  if (length(my.temp.files)!=0) file.remove(my.temp.files)

  # error catching (check if .tex is compilable)

  #rte.compile.latex(f.in = f.in, pdf.dir.out = latex.dir.out)

  # clean up after compiling tex file

  #my.temp.files <- dir(latex.dir.out, pattern = '*.*', full.names = T)
  #if (length(my.temp.files)!=0) file.remove(my.temp.files)

  # read tex file

  my.text <- stringi::stri_read_lines(f.in)
  Encoding(my.text) <- 'UTF-8'

  # find idx of text of multiple choice questions

  line.beg.mchoice <- which(stringr::str_detect(string = my.text, stringr::fixed('begin{multiplechoice}')))
  line.end.mchoice   <- which(stringr::str_detect(string = my.text, stringr::fixed('end{multiplechoice}')))

  mchoice.text <- my.text[line.beg.mchoice:line.end.mchoice]

  # find beggining and ending of all multiple choice questions

  idx.beg.mchoices <- which(stringr::str_detect(string = mchoice.text, stringr::fixed('\\begin{question}')))
  idx.end.mchoices <- which(stringr::str_detect(string = mchoice.text, stringr::fixed('\\end{question}')))

  # find text part before and after multiple choice questions

  my.preamble <- paste0(my.text[1:(line.beg.mchoice[1]-1)], collapse = ' \n')
  my.last.part <- paste0(my.text[(line.end.mchoice[length(line.end.mchoice)]+1):length(my.text)], collapse = '\n')

  # Find text with options for \begin{multiplechoice}

  idx.mchoice.comm <- which(stringr::str_detect(string = my.text, stringr::fixed('\\begin{multiplechoice}')))
  my.begin.mchoice.line <- my.text[idx.mchoice.comm]

  # build dataframe with questions to be randomized

  df.questions <- data.frame(n.question = seq(1,length(idx.beg.mchoices)),
                             idx.begin = idx.beg.mchoices,
                             idx.end = idx.end.mchoices)

  # function to get text from m.questions

  locate.questions <- function(x,mchoice.text){
    q.temp <- mchoice.text[x['idx.begin']:(x['idx.end'])]
    q.temp <- paste0(q.temp,collapse = '\n')
    return(q.temp)
  }

  df.questions$q.text <- apply(df.questions,MARGIN = 1 , FUN = locate.questions, mchoice.text=mchoice.text)

  # function to locate all choices of all m.questions

  locate.choices <- function(x){
    q.temp <- x['q.text']
    idx <- stringr::str_locate(q.temp, '\\\\choice')[1]
    q.temp <- stringr::str_sub(q.temp,1,idx-2)
    return(q.temp)
  }

  df.questions$main.text <- apply(df.questions,MARGIN = 1 , FUN = locate.choices)


  df.answers <- data.frame() # Dataframe with all answers to all questions
  for (i.q in df.questions$n.question){

    q.temp <- df.questions$q.text[i.q]
    q.temp <- stringr::str_split(q.temp, pattern = '\\n')[[1]] # break lines

    idx.choice <- which(stringr::str_detect(string = q.temp, pattern = stringr::fixed('\\choice')))



    my.main.text <- as.character(paste0('\n',q.temp[1:(idx.choice[1]-1)], collapse = '\n'))

    my.answers <- as.character(q.temp[idx.choice])
    my.correct.answers <- stringr::str_detect(string = my.answers, '!')

    # Build data.frame for output

    df.answers <- rbind(df.answers, data.frame(n.question = rep(i.q, length(my.answers)),
                                               main.text = rep(my.main.text,length(my.answers)),
                                               text.answer = my.answers,
                                               correct.answer = my.correct.answers, stringsAsFactors=F ))

  }

  # return a list with questions, answers, preamble and last part of tex file

  out <- list(df.questions=df.questions, # df with all questions
              df.answers = df.answers,   # df with all answers
              my.begin.mchoice.line = my.begin.mchoice.line, # text with begging of mchoice enviroment
              my.preamble = my.preamble, # preamble of tex file
              my.last.part = my.last.part) # last part of tex file

  cat(' Done')
  return(out)

}
