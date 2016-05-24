msr <- function (ms, pSelect = c("adj.r.squared", "bic" ), modelSelect = FALSE) 
{

      pSelect <- match.arg(pSelect)
 
      buildlm <- function(ms){
        nc <- ncol(ms$x)
        nr <- nrow(ms$x)
        lms <- vector("list", length(ms$mins))
        df <- data.frame(y=ms$y, ms$x) 
        if(length(ms$crystalsSize)  == 1) {
          df <- data.frame(y = ms$y, ms$x) 
          l <- lm(y ~ ., data=df)
        }
        else{
          df <- data.frame(y = ms$y, ms$x, cID = as.factor(ms$crystals)) 
          l <- lm(y ~ . + .:cID, data=df)
        }
        if(modelSelect){
          l <- step(l, trace = 0)
        }
        s <- summary(l)
        ar2 <- s$adj.r.squared;
        bic <- BIC(l)
        #Needs stratified sampling to work
        #cv <-  cv.lm(formula(l), df) 
        obj <- list(bic = bic, ar2=ar2, lm = l)
      }


      buildslm <- function(ms){
        obj <- structure(list(ms = ms), class = "msr")
        mm <- model.matrix(obj, ms$x)
        df <- data.frame(y = ms$y, mm)
        slm <- lm(y ~ 0+., df)
        if(modelSelect){
          slm <- step(slm, trace = 0)
        }
        cv <- cv.lm(formula(slm), df)
        obj <- list(lm = slm,  cv = cv$meanSE)
      }

      if(is.null(ms$nLevels)){ 
        l <- buildlm(ms)
        if(class(ms) == "msc.svm"){
          if( ! is.null(ms$svm) ){
            df <- data.frame(cry = as.factor(ms$crystals), obj$x)
            ms$svm <- svm(cry ~ ., df, probability=TRUE, cost=ms$cost, scale = FALSE)
          }
        }
        obj <- structure(list(ms = ms, lm = l, slm = buildslm(ms)), class = "msr")
      }
      else{      
        lms <- c()
        minP <- Inf
        pLevel <- 1
        for(i in 1:ms$nLevels){
          lms[[i]] <- buildlm(ms$mscl[[i]])
          tmp <- 0
          if(pSelect == "cv"){
            tmp <- -lms[[i]]$cv
          }
          else if(pSelect == "bic"){
            tmp <- -lms[[i]]$bic
          }
          else{
            tmp <- lms[[i]]$ar2
          }
          if(minP > tmp){
            minP <- tmp
            pLevel <- i
          }
        }
        if(class(ms) == "msc.svm"){
          if( ! is.null(ms$level[[pLevel]]$svm) ){
            df <- data.frame(cry = as.factor(ms$mscl[[pLevel]]$crystals), obj$mscl[[pLevel]]$x)
            ms$level[[pLevel]]$svm <- svm(cry ~ ., df, probability=TRUE, cost=ms$cost, scale = FALSE)
          }
        }

        slm <- buildslm(ms$level[[pLevel]])
        obj <- structure(list(ms = ms, lms = lms, slm= slm, predictLevel = pLevel), class = "msr")
      }
      obj
}

