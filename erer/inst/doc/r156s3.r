# A. New function for adjusting and calculating exam scores
score <- function(y, sort.by = c("Letter", "Exam.A", "Exam.B", "Total"),
  ...) {

  # A1. Inputs
  if (!inherits(y, "data.frame") || ncol(y) != 3) {
    stop("y should be a data fram with three columns.\n")
  }
  sort.by <- match.arg(sort.by)
     
  # A2. Transformation
  y2 <- y
  y2$Total <- y2[, 2] * 0.4 + y2[, 3] * 0.6
  y2$Letter <- ifelse(test = y2$Total >= 90, yes = "A", 
    no = ifelse(test = y2$Total >= 80, yes = "B", 
    no = ifelse(test = y2$Total >= 70, yes = "C", 
    no = ifelse(test = y2$Total >= 60, yes = "D", no = "F"))))  
  y.new <- y2[order(y2[, sort.by], ...), ]
  rownames(y.new) <- 1:nrow(y.new)
      
  # A3. Output
  result <- list(y = y, y.new = y.new)
  class(result) <- c("score", "scoreNum")
  return(result)
}

# B. Define print method
print.score <- function (x, ...) {
  need <- x$y.new
  cat("\n========================================\n")
  cat("    Course Exam Scores\n")
  cat(  "========================================\n")
  print(need)
  invisible(x)
}

summary.scoreNum <- function(object, ...) {
  out <- list(count = table(object$y.new$Letter), 
    average = round(x = mean(object$y.new$Total), ...))
  return(out)
}

# C. Test
args(score); formals(score); body(score)
methods(class = "score"); methods(class = "scoreNum")
methods(generic.function = print)
print.default; getAnywhere(print.acf); stats:::print.acf
print.acf  # Error: object 'print.acf' not found

raw <- data.frame(Student = I(c("Smith, Bob", "Green, Joe", "Hue, Lisa")),
  Exam.A = c(85, 55, 90), Exam.B = c(82, 49, 97))
raw
agg <- score(y = raw, sort.by = "Letter", decreasing = FALSE)
names(agg); agg
summary(agg, digits = 1) 

# D. S3 is convenient but not rigorous / not safe
aa <- 10; class(aa); aa                        # print.default for vector
class(aa) <- c("score", "scoreNum"); class(aa) # class changed
aa                                             # print.score not working
aa <- c(aa, 20); class(aa); aa                 # class removed by c()
aa                                             # print.default works again