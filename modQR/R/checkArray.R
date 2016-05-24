checkArray <- function(Subject, IsInt, IsUnique, ValueRange, NRowRange, NColRange, ErrStatus){
#checkArray <- function(Subject, IsInt, IsUnique, ValueRange, NRowRange, NColRange, ErrStatus), output: list(Subject, Status)
#checking an input (numeric array) Subject and making it real or integer if necessary
# Status = 0 if the check was successful, and Status = ErrStatus otherwise
#
#IsInt      ... 1/0           : Subject should/should not be integer-valued
#IsUnique   ... 1/0           : Subject should/should not contain unique elements
#ValueRange ... [p q] or []   : the range for values of Subject
#NRowRange  ... [p q] or []   : the range for the number of rows of Subject
#NColRange  ... [p q] or []   : the range for the number of columns of Subject
#ErrStatus  ... p             : the value of Status when Subject does not meet the conditions
#
#Note: Any NULL ([]) range means no restriction (an arbitrary value)
Status <- 0;

if (!is.numeric(Subject)){                                               Status <- ErrStatus; return(list(Subject, Status))}
NRow <- dim(as.matrix(Subject))[1]
NCol <- dim(as.matrix(Subject))[2]
if (!is.null(NRowRange)&&((NRow < NRowRange[1])||(NRow > NRowRange[2]))){Status <- ErrStatus; return(list(Subject, Status))}
if (!is.null(NColRange)&&((NCol < NColRange[1])||(NCol > NColRange[2]))){Status <- ErrStatus; return(list(Subject, Status))}
if (IsInt) {Subject <- round(Re(Subject))} else {Subject <- Re(Subject)}
if (IsUnique && (any(duplicated(Subject)))){                             Status <- ErrStatus; return(list(Subject, Status))}
if (!is.null(ValueRange)){
    if (any(Subject<ValueRange[1])||any(Subject>ValueRange[2])){         Status <- ErrStatus; return(list(Subject, Status))}
}
return(list(Subject, Status));
}