utils::globalVariables(c("moves", "switches", ".SD"))
#' Function to calculate the number of times a student has changed schools.
#' @description This function calculates the number of times a student has changed 
#' schools, including accounting for gaps in enrollment data. It returns a 
#' \code{\link{data.table}} with the student ID and the number of student moves.
#' @param df a data.frame containing minimally a student identifier, school identifier, enrollment date, and exit date.
#' @param enrollby a date that determines the earliest a student can enroll for 
#' the first time without being credited with having moved at least once.
#' @param exitby  a date that determines the latest a student can exit for the final 
#' time without being credited with having moved at least once.
#' @param gap a number, of days,  that represents the largest gap between an exit date 
#' and the next enrollment date that can occur without indicating the student 
#' moved to a third school not contained within the data set. The default value is 
#' \code{14}.
#' @param sid  a character that indicates the name of the student id attribute 
#' in \code{df}. The default value is \code{sid}.
#' @param schid a character that indicates the name of the school id attribute 
#' in \code{df}. The default value is \code{schid}.
#' @param enroll_date a character that indicates the name of the enrollment date 
#' attribute in \code{df}. The default value is \code{enroll_date}.
#' @param exit_date a character that indicates the name of the student id 
#' attribute in \code{df}. The default value is \code{exit_date}.
#' @details \code{enrollby} and \code{exitby} are specified automatically if not 
#' defined. They are assigned to the default dates of -09-15 and -06-01 of the min 
#' and max year respectively.
#' @author Jason P. Becker
#' @import data.table
#' @return a data.frame
#' @export
#'
#' @examples
#' \dontrun{
#' df <- data.frame(sid = c(rep(1,3), rep(2,4), 3, rep(4,2)), 
#'                  schid = c(1, 2, 2, 2, 3, 1, 1, 1, 3, 1),
#'                  enroll_date = as.Date(c('2004-08-26',
#'                                          '2004-10-01',
#'                                          '2005-05-01',
#'                                          '2004-09-01',
#'                                          '2004-11-03',
#'                                          '2005-01-11',
#'                                          '2005-04-02',
#'                                          '2004-09-26',
#'                                          '2004-09-01',
#'                                          '2005-02-02'), 
#'                                        format='%Y-%m-%d'),
#'                  exit_date = as.Date(c('2004-08-26',
#'                                        '2005-04-10',
#'                                        '2005-06-15',
#'                                        '2004-11-02',
#'                                        '2005-01-10',
#'                                        '2005-03-01',
#'                                        '2005-06-15',
#'                                        '2005-05-30',
#'                                        NA,
#'                                        '2005-06-15'), 
#'                                      format='%Y-%m-%d'))
#' moves <- moves_calc(df)
#' moves
#' moves <- moves_calc(df, enrollby='2004-10-15', gap=22)
#' moves
#' moves <- moves_calc(df, enrollby='2004-10-15', exitby='2005-05-29')
#' moves
#' }
#' 
moves_calc <- function(df, 
                       enrollby,
                       exitby,
                       gap=14,
                       sid='sid', 
                       schid='schid',
                       enroll_date='enroll_date',
                       exit_date='exit_date'){
    # df is a data.frame that minimally contains a student ID (default 'sid'),
  # a school ID (default 'schno'), and two dates an enrollment date and an 
  # exit date for each sid-schid combination.
  # Type checking inputs.
  if (!inherits(df[[enroll_date]], "Date") | !inherits(df[[exit_date]], "Date"))
      stop("Both enroll_date and exit_date must be Date objects")
  # Check if enrollby and exitby arguments are supplied. When they aren't,
  # assign them to default dates of -09-15 and -06-01 of the min and max year.
  # If they are assigned but are not date objects then 
  if(missing(enrollby)){
    enrollby <- as.Date(paste(year(min(df[[enroll_date]], na.rm=TRUE)),
                              '-09-15', sep=''), format='%Y-%m-%d')
  }else{
    if(is.na(as.Date(enrollby, format="%Y-%m-%d"))){
      enrollby <- as.Date(paste(year(min(df[[enroll_date]], na.rm=TRUE)),
                                '-09-15', sep=''), format='%Y-%m-%d')
      warning(paste("enrollby must be a string with format %Y-%m-%d,",
                    "defaulting to", 
                    enrollby, sep=' '))
    }else{
      enrollby <- as.Date(enrollby, format="%Y-%m-%d")
    }
  }
  if(missing(exitby)){
    exitby <- as.Date(paste(year(max(df[[exit_date]], na.rm=TRUE)),
                            '-06-01', sep=''), format='%Y-%m-%d')
  }else{
    if(is.na(as.Date(exitby, format="%Y-%m-%d"))){
      exitby <- as.Date(paste(year(max(df[[exit_date]], na.rm=TRUE)),
                                '-06-01', sep=''), format='%Y-%m-%d')
      warning(paste("exitby must be a string with format %Y-%m-%d,",
                    "defaulting to", 
                    exitby, sep=' '))
    }else{
      exitby <- as.Date(exitby, format="%Y-%m-%d")
    }
  }
  if(!is.numeric(gap)){
    gap <- 14
    warning("gap was not a number, defaulting to 14 days")
  }
  # Generate results data table
  output <- data.frame(id = as.character(unique(df[[sid]])),
                       moves = vector(mode = 'numeric', 
                                      length = length(unique(df[[sid]]))))
  # Students with missing data receive missing moves
  incomplete <- df[!complete.cases(df[, c(enroll_date, exit_date)]), ]
  if(nrow(incomplete)>0){
    output[which(output[['id']] %in% incomplete[[sid]]),][['moves']] <- NA
  }
  output <- data.table(output, key='id')
  
  df <- df[complete.cases(df[, c(enroll_date, exit_date)]), ]
  dt <- data.table(df, key=sid)
  dt[[sid]] <- as.factor(as.character(dt[[sid]]))
  setnames(dt, names(dt)[which(names(dt) %in% enroll_date)], "enroll_date")
  setnames(dt, names(dt)[which(names(dt) %in% exit_date)], "exit_date")
  dt[['moves']] <- 0

  first <- dt[, list(enroll_date=min(enroll_date)), by=sid]
  last <- dt[, list(exit_date=max(exit_date)), by=sid]
  output[id %in% first[enroll_date>enrollby][[sid]], moves:=moves+1L]
  output[id %in% last[exit_date<exitby][[sid]], moves:=moves+1L]
  
  # Select all students who have more than one row. Create a recursive function
  # that checks that rows > 2, selects the min exit_date and difftimes the min # enroll_date thats > the exit_date value. If >gap, add 1 to counter. Remove 
  # observation with min exit_date, then call self with what's left. Break when
  # rows < 2 and return sid and moves counter.
  school_switch <- function(dt, x=0){
    if(dim(dt)[1]<2){
      return(x)
    }else{
      exit <- min(dt[, exit_date])
      exit_school <- dt[exit_date==exit][[schid]]
      rows <- dt[, enroll_date]>exit
      dt <- dt[rows,]
      enroll <- min(dt[, enroll_date])
      enroll_school <- dt[enroll_date==enroll][[schid]]
      if(difftime(min(dt[, enroll_date], na.rm=TRUE), exit)<gap &
         exit_school==enroll_school){
        y = x
      }else if(difftime(min(dt[, enroll_date], na.rm=TRUE), exit)<gap){
        y = x + 1L
      }else{
        y = x + 2L
      }
      school_switch(dt, y)
    }
  }
  dt[, moves:= school_switch(.SD), by=sid]
  dt <- dt[,list(switches=unique(moves)), by=sid]
  # Need to combine dt with output
  output[dt, moves:=moves+switches]
  # Set names properly and return data.frame as supplied.
  setnames(output, 'id', sid)
  return(as.data.frame(output))
}