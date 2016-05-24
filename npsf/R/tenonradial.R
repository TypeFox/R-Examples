tenonradial <- function(formula, data, subset,
                        rts = c("C", "NI", "V"), 
                        base = c("output", "input"), 
                        ref = NULL, data.ref = NULL, subset.ref = NULL,
                        print.level = 1)
{
 if( !is.null(ref) & is.null(data.ref) ){
  warning("If you use variable names in 'ref', 'data.ref' is required", call. = FALSE)
 }
 
 
 if( is.null(ref) & !is.null(data.ref) ){
  stop("If you use 'data.ref', 'ref' is required", call. = FALSE)
 }
 
 winw <- getOption("width")
 if (winw > 100+5){
  winw <- 100
 }
 else if (winw > 90+5) {
  winw <- 90
 }
 else if (winw > 80+5) {
  winw <- 80
 }
 else if (winw > 70+5) {
  winw <- 70
 }
 else if (winw > 60+5) {
  winw <- 60
 }
 else if (winw > 50+5) {
  winw <- 50
 }
 else {
  winw <- 0
 }
 
 t1 <- .prepareYX(formula = formula, data = data, subset = subset, rts = rts,
                  base = base, ref = ref,	data.ref = data.ref, subset.ref = subset.ref,
                  print.level = print.level, type = "RM", winw = winw)
 
 te <- .teNonrad(t(t1$y), t(t1$x), ncol(t1$y), ncol(t1$x), nrow(t1$y),
              t(t1$y.ref), t(t1$x.ref), nrow(t1$y.ref), t1$myrts, t1$mybase,
              1, print.level = print.level)
 if(print.level >= 3){
  cat("\n")
 }
 if(print.level >= 2){
  cat(paste("",rep("_", (winw-10)/1),"", sep = ""), "\n\n", sep = "")
  cat("Summary of efficiencies:\n\n", sep = "")
  .su(te, print = TRUE)
 }
 tymch <- list(K = nrow(t1$y), M = ncol(t1$y), N = ncol(t1$x), rts = t1$rts.string, base = t1$base.string, te = te, esample = t1$esample, esample.ref = t1$esample_ref)
 class(tymch) <- "npsf"
 return(tymch)
}



#