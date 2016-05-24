"nuskewCOP" <-
function(cop=NULL, para=NULL, as.sample=FALSE, brute=FALSE, delta=0.002, ...) {
   joeskewCOP(cop=cop, para=para, type="nu",
              brute=brute, delta=delta, as.sample=as.sample, ...)
}
"nustarCOP" <-
function(cop=NULL, para=NULL, as.sample=FALSE, brute=FALSE, delta=0.002, ...) {
   joeskewCOP(cop=cop, para=para, type="nustar",
              brute=brute, delta=delta, as.sample=as.sample, ...)
}
"joeskewCOP" <-
function(cop=NULL, para=NULL, type=c("nu", "nustar"), as.sample=FALSE,
                              brute=FALSE, delta=0.002,...) {

   type = match.arg(type)

   if(as.sample) {
      if(is.null(para)) {
         warning("Sample Joe's Nu-Skew desired but para is NULL, ",
                 "returning NULL")
         return(NULL)
      }
      if(length(names(para)) != 2) {
        warning("para argument must be data.frame having only two columns, ",
                "returning NULL")
        return(NULL)
      }

      n <- length(para[,1]); nn <- n^2
      R <- rank(para[,1]); S <- rank(para[,2])
      samNU <- NA
      if(type == "nu") {
         if(as.sample == -1) message("Sample Nu-Skew after Joe (2015)",
                                     "---CPU intensive!")
         samNU <- sum(sapply(1:n, function(i) {
                     sum(sapply(1:n, function(j) {
                        (j - i)*sum(as.numeric(R <= i & S <= j))
                     } ))/nn
                  } ))
      } else if(type == "nustar") {
         if(as.sample == -1) message("Sample Nu-Skew-Star after Joe (2015)",
                                  "---CPU intensive!")
         samNU <- sum(sapply(1:n, function(i) {
                  sum(sapply(1:n, function(j) {
                     (j + i - n)*sum(as.numeric(R <= i & S <= j))
                  } ))/nn
               } ))
      } else {
         stop("Never should be here in logic")
      }
      samNU <- (6/(nn - 1)) * samNU
      if(type == "nustar") samNU <- samNU - 1/2
      return(samNU)
   }

   if(brute) {
      us <- vs <- seq(.Machine$double.eps, 1-.Machine$double.eps, delta)
      skew <- NA
      if(type == "nu") {
         skew <- sum(sapply(us, function(u) {
                    sum(sapply(vs, function(v) {
                       (v-u)*COP(u,v,cop=cop,para=para, ...)
                    }))
                 }))
      } else if(type == "nustar") {
         skew <- sum(sapply(us, function(u) {
                    sum(sapply(vs, function(v) {
                       (v+u-1)*COP(u,v,cop=cop,para=para, ...)
                    }))
                 }))
      } else {
         stop("Never should be here in logic")
      }
      skew <- 6*skew*delta^2
      if(type == "nustar") skew <- skew - 1/2
      return(skew)
   }

   myint <- NULL
   if(type == "nu") {
      try(myint <- integrate(function(u) {
               sapply(u,function(u) { integrate(function(v) {
                          (v-u)*COP(u,v,cop=cop, para=para,...)
               }, 0, 1)$value })}, 0, 1) )
   } else if(type == "nustar") {
      try(myint <- integrate(function(u) {
               sapply(u,function(u) { integrate(function(v) {
                          (v+u-1)*COP(u,v,cop=cop, para=para,...)
               }, 0, 1)$value })}, 0, 1) )
   } else {
      stop("Never should be here in logic")
   }
   skew <- ifelse(is.null(myint), NA, 6*myint$value)
   if(type == "nustar") skew <- skew - 1/2
   return(skew)
}
