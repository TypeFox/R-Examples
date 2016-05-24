#' Performs statistical tests for cheating using \pkg{CopyDetect}
#'
#' Uses function \code{\link[CopyDetect]{CopyDetect1}} from package \pkg{CopyDetect} to test for visual cheating in exams
#' based on correct answers from the students. The statistical tests performed by \code{\link[CopyDetect]{CopyDetect1}} are:
#' \itemize{
#' \item Omega index (Wollack, 1996)
#' \item Generalized Binomial Test ([GBT], van der Linden & Sotaridona (2006)
#' \item K index (Holland, 1996)
#' \item K1 and K2 indices (Sotaridona & Meijer, 2002)
#' \item S1 and S2 indices (Sotaridona & Meijer, 2003) }
#' The function \code{rte.test.cheating} will have as input a dataframe with the names and corrections of students
#' and output a summary of the cheating tests as a list, including suspicious pairs.
#'
#' More details regarding the tests can be found in:
#'
#' Zopluoglu, C. (2013). CopyDetect An R Package for Computing Statistical Indices to Detect Answer Copying on Multiple-Choice Examinations.
#' Applied psychological measurement, 37(1), 93-95.
#'
#' The article can be found \href{http://apm.sagepub.com/content/37/1/93.extract}{here}
#'
#' @param df.grade  A dataframe where first column is the name of students (item exam.names)
#' and the rest of the columns are the correct (TRUE) and incorrect (FALSE) answers. Each column other
#' than exam.names should be a question
#' @param p.level Critial level of statistical testing
#' @param print.suspects Print suspects on screen ? (TRUE or FALSE) (Default=TRUE)
#' @param do.cheat.plot Print plot of cheating tests? (TRUE of FALSE) (Default=TRUE)
#' @param suspicion.threshold Proportion of failed cheating tests that justify suspition, between 0 and 1
#'
#' @return A list with the following items: \describe{ \item{df.pvalue}{A dataframe
#'   with the statistical results for all pairs of students from the upper triangle (1 test for each pair) }
#'   \item{df.suspects}{A dataframe with the suspicious pair of students } }
#'
#' @export
#' @import CopyDetect
#'
#' @examples
#' # number of simulated questions in exam
#' n.sim.questions <- 15
#'
#' base.names <- c('John', 'Marcelo','Ricardo', 'Tarcizio')
#' last.names <- c('Smith', 'Johnson','P.')
#'
#' name.grid <- expand.grid(base.names,last.names)
#'
#' my.names <- paste(name.grid[,1], name.grid[,2])
#' # official names from the university system (will assume it is equal to my.names)
#' # In a practical situation, this list of official names will come from the university system
#' exam.names <- my.names
#'
#' set.seed(10)
#'
#' correction.mat <- matrix(sample(c(TRUE,FALSE),
#'                                 size = length(exam.names)*n.sim.questions,
#'                                 replace = TRUE),nrow = length(exam.names))
#'
#' idx.cheater.1 <- 5 # std 5 and 6 have simillar correct answers
#' idx.cheater.2 <- 6
#' proportion.to.cheat <- 0.6  # proportion of same correct answers
#' q.to.cheat <- floor(proportion.to.cheat*n.sim.questions)
#' correction.mat[idx.cheater.1, ] <-  c(rep(TRUE,q.to.cheat),
#'                                       rep(FALSE,n.sim.questions-q.to.cheat))
#'
#' correction.mat[idx.cheater.2, ] <- correction.mat[idx.cheater.1, ]
#'
#'
#' df.grade <- cbind(data.frame(exam.names),correction.mat)
#'
#'
#' test.cheating.out <- rte.test.cheating(df.grade, do.cheat.plot = FALSE )
rte.test.cheating <- function(df.grade,
                              p.level = 0.05,
                              print.suspects = TRUE,
                              do.cheat.plot = TRUE,
                              suspicion.threshold = 0.5){

  if (is.null(df.grade$exam.names)){
      stop('Input df.grade does not have named element exam.names. Please use the output list from rte.grade.exams')
  }

  if (ncol(df.grade)<2){
    stop('Input df.grade should have more than 2 columns. See help for details.')
  }

  if (any(p.level<0,p.level>1)){
    stop('The values of p.level should be between 0 and 1')
  }

  if (any(!is.logical(print.suspects),!is.logical(do.cheat.plot))){
    stop('Inputs print.suspects and do.cheat.plot should be logical values (TRUE or FALSE)')
  }

  exam.names <- df.grade$exam.names

  data.in <- as.data.frame(df.grade)[,-1]

  n.std <- nrow(data.in)

  my.grid <- expand.grid(seq(n.std),seq(n.std))
  names(my.grid) <- c('Std1','Std2')

  #my.grid <- my.grid[with(my.grid,Std1!=Std2), ]
  #my.grid <- my.grid[with(my.grid,Std1>Std2), ]

  my.grid$exam.names.Std1 <- exam.names[my.grid$Std1]
  my.grid$exam.names.Std2 <- exam.names[my.grid$Std2]

  # try to estimate the parameters of CopyDetect

  est.ipar <- NA
  try(est.ipar <- irtoys::est(data.in, model = "2PL", engine = "ltm")$est,
      silent = TRUE)

  if (is.na(est.ipar[1])){
    stop(paste('The function irtoys::est was not able to estimate the model underlying CopyDetect calculations.',
               'This is usually due to a low number of students/rows in the input dataframe.',
               'Given this problem, the function rte.test.cheating is not able to',
               'move forward with the statistical tests provided in CopyDetect.') )

  }

  # run CopyDetect on all pairs

  df.pvalue <- data.frame()
  for (i.pair in seq(nrow(my.grid))){

    pair.now <- c(my.grid$Std1[i.pair],my.grid$Std2[i.pair])

    # control for test (dont do it for same person and dont do it for i>j ())

    run.test <- any(pair.now[1]==pair.now[2], pair.now[1]>pair.now[2])

    cat(paste0('\n',
               '[',i.pair,'|',nrow(my.grid),'] - ',
               my.grid$exam.names.Std1[i.pair],' | ',my.grid$exam.names.Std2[i.pair]))


    if (run.test){
      temp.df <- data.frame(n.std.1 = pair.now[1],
                            n.std.2 = pair.now[2],
                            Student1 = my.grid$exam.names.Std1[i.pair],
                            Student2 = my.grid$exam.names.Std2[i.pair],
                            W.index.pvalue = NA,
                            GBT.index.pvalue = NA,
                            K.index.pvalue = NA,
                            K1.pvalue = NA,
                            K2.pvalue = NA,
                            S1.pvalue = NA,
                            S2.pvalue = NA)

      df.pvalue <- rbind(df.pvalue, temp.df)

      cat(' - No need for this run (repeated case). Skipping..')

      next()

    } else {


      out.copydetect <- CopyDetect::CopyDetect1(data.in,
                                                item.par = est.ipar,
                                                pair = pair.now)

      if (is.null(out.copydetect$K.variants)){
        out.copydetect$K.index$k.index <- NA
        out.copydetect$K.variants$K1.index <- NA
        out.copydetect$K.variants$K2.index <- NA
        out.copydetect$K.variants$S1.index <- NA
        out.copydetect$K.variants$S2.index <- NA

      }
    }

    temp.df <- data.frame(n.std.1 = pair.now[1],
                          n.std.2 = pair.now[2],
                          Student1 = my.grid$exam.names.Std1[i.pair],
                          Student2 = my.grid$exam.names.Std2[i.pair],
                          W.index.pvalue =out.copydetect$W.index$p.value,
                          GBT.index.pvalue = out.copydetect$GBT.index$p.value,
                          K.index.pvalue = out.copydetect$K.index$k.index,
                          K1.pvalue = out.copydetect$K.variants$K1.index,
                          K2.pvalue = out.copydetect$K.variants$K2.index,
                          S1.pvalue = out.copydetect$K.variants$S1.index,
                          S2.pvalue = out.copydetect$K.variants$S2.index)

    rownames(temp.df) <- NULL
    df.pvalue <- rbind(df.pvalue, temp.df)

    my.tests <- as.matrix(temp.df[ , c(3:8)])
    n.tests <- sum(!is.na(my.tests))
    my.failed.prop <- sum(my.tests<p.level,na.rm = TRUE)/n.tests

    suspicion.flag <- ''
    if (my.failed.prop > suspicion.threshold){
      suspicion.flag <- '  --> Cheater?'
    }

    cat(paste0(' | Proportion of failed tests cheating ',
               format(my.failed.prop*100,digits = 2),
               '% ',
               suspicion.flag))

  }


  tests.mat <- as.matrix(df.pvalue[,c(5:11)])

  n.tests.not.na <- apply(tests.mat,MARGIN = 1,FUN = function(x) return(sum(!is.na(x))))
  n.test.fail <- apply(tests.mat,MARGIN = 1,FUN = function(x) return(sum(x<p.level,na.rm = TRUE)))

  df.pvalue$prop.test.fail <- n.test.fail/n.tests.not.na

  idx <- n.tests.not.na==0
  df.pvalue$prop.test.fail[idx] <- NA

  df.pvalue <- df.pvalue
  idx <- which(df.pvalue$prop.test.fail>suspicion.threshold)
  df.suspects <- df.pvalue[idx, ]
  n.suspects <- nrow(df.suspects)

  if (do.cheat.plot){
    n.std.1 <- n.std.2 <- prop.test.fail <- NULL # fix for CHECK Note ("no visible binding for global variable 'n.std.1")

    p <- ggplot2::ggplot(data = df.pvalue,ggplot2::aes(x = n.std.1, y = n.std.2, size = prop.test.fail))
    p <- p + ggplot2::geom_point()
    p <- p + ggplot2::scale_fill_gradient(low = "white",  high = "black")
    p <- p + ggplot2::labs(x = 'idx Std 1', y = 'idx Std 2')
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))

    if (n.suspects>0){
      p <- p + ggplot2::geom_text(data = df.suspects, ggplot2::aes(x = n.std.1, y = n.std.2, label='X', color='white'),show.legend=F)
    }

    print(p)

  }

  if (print.suspects){

    if (n.suspects>0){
      cat('\nList of Suspects: \n')
      cat(paste0('\n',
                 df.suspects$Student1,
                 ' <-> ',
                 df.suspects$Student2,
                 ' Proportion of failed tests:',
                 format(df.suspects$prop.test.fail*100,digits = 2),
                 '%'))

    } else {
      cat('\nList of Suspects: \n')
      cat('\nNo suspects to print! Nice..')
    }

  }


  return(list(df.pvalue=df.pvalue,
              df.suspects = df.suspects))

}
