fetchSudokuUK <- function(day){
  if(missing(day)){
    th <- url('http://www.sudoku.org.uk/DailySudoku.asp')
  } else {
    th <- url(paste('http://www.sudoku.org.uk/DailySudoku.asp?day=',day,sep=''))
  }

  
  tmp <- readLines(th)

  close(th)

  tmp2 <- grep('InnerTDone',tmp,value=TRUE,fixed=TRUE)

  if(length(tmp2) < 81){
    stop('Unable to download full puzzle, did you specify a correct date?\n')
  }

  tmp3 <- regexpr('.</td>$', tmp2)

  vals <- substr(tmp2,tmp3,tmp3)
  vals <- as.numeric( sub('[^1-9]','0',vals) )

  matrix(vals,9,9,byrow=TRUE)
}
