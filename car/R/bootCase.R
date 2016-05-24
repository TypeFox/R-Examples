# March 9, 2012 modified by SW as suggested by Derek Ogle to return an object
# of class c("bootCase", "matrix").  
# May 2012 added methods for 'bootCase'
# 2012-12-10 replaced .GlobalEnv by car:::.carEnv to suppress warnings
# 2013-01-28 Changed argument f to f.
# 2013-07-08 Changed .carEnv to car:::.carEnv
# 2015-01-27 .carEnv now in global environment. John

nextBoot <- function(object, sample){UseMethod("nextBoot")}
nextBoot.default <- function(object, sample){
   update(object, subset=sample)
   }
nextBoot.lm <- function(object, sample) nextBoot.default(object, sample)
nextBoot.nls <- function(object, sample){
# modify to assure resampling only rows in the original subset 9/1/2005
   update(object,subset=sample,start=coef(object),
    data=data.frame(update(object,model=TRUE)$model))}

bootCase <- function(object, f.=coef, B=999){UseMethod("bootCase")}
bootCase.lm <- function(object, f.=coef, B=999) {
    bootCase.default(object, f., B, names(resid(object)))
#    bootCase.default(update(object, 
#             data=na.omit(model.frame(object))), f, B)
    }
bootCase.glm <- function(object, f.=coef, B=999) {
    bootCase.lm(object, f., B)
    }
bootCase.nls <- function(object, f.=coef, B=999) {
    bootCase.default(object, f., B, seq(length(resid(object))))
    }
bootCase.default <- function (object, f.=coef, B = 999, rows)
{       
    n <- length(resid(object))
    opt<-options(show.error.messages = FALSE)
    on.exit(options(opt))
    pointEstimate <- f.(object)
    coefBoot <- matrix(0, nrow=B, ncol=length(f.(object)))
    colnames(coefBoot) <- names(pointEstimate)  # adds names if they exist
    class(coefBoot) <- c("bootCase", "matrix")
    count.error <- 0
    i <- 0
    while (i < B) {
		assign(".boot.sample", sample(rows, replace=TRUE), envir=.carEnv)
        obj.boot <- try(update(object, subset=get(".boot.sample", envir=.carEnv)))
        if (is.null(class(obj.boot))) {
            count.error <- 0
            i <- i + 1
            coefBoot[i, ] <- f.(obj.boot)
        }
        else {
            if (class(obj.boot)[1] != "try-error") {
                count.error <- 0
                i <- i + 1
                coefBoot[i, ] <- f.(obj.boot)
            }
            else {
                count.error <- count.error + 1
            }
        }
        if (count.error >= 25) {
            options(show.error.messages = TRUE)
            stop("25 consecutive bootstraps did not converge.  Bailing out.")}
    }
	remove(".boot.sample", envir=.carEnv)
	attr(coefBoot, "pointEstimate") <- pointEstimate
    return(coefBoot)
}

