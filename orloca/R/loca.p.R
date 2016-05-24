#
# Class loca.p definition
#
setClass("loca.p",
	representation(x="numeric", y="numeric", w="numeric", label="character")
	)

#
# loca.p Validity method
#
setValidity("loca.p", 
   function(object)
      {
      if(length(object@x)==length(object@y))
         {
	 if (length(object@x)==length(object@w) || length(object@w)==0)
	    {
            if (!any(is.na(object@x)) && !any(is.na(object@y) && !any(is.na(object@w)))) TRUE
            else paste(gettext("NA's values are not allowed", domain = "R-orloca"), sep="")
	    }
	 else paste(gettext("The length of w (", domain = "R-orloca"), length(object@w), gettext(") should be the same as the length of x, and y (", domain = "R-orloca"), length(object@x) ,gettext(") or 0", domain = "R-orloca"))
         }
      else paste(gettext("The length of x and y are different", domain = "R-orloca"), length(object@x), ", ", length(object@y))
      }
   )

#
# loca.p initialize method
#
setMethod("initialize", "loca.p",  
   function(.Object, x, y, w = numeric(0), label="")
      {
      .Object@x <- x
      .Object@y <- y
      if (length(w) == 0) .Object@w <- rep(1,length(x))
      else .Object@w <- w
      .Object@label <- label
      validObject(.Object)
      .Object
      }
)

loca.p <- function(x, y, w = numeric(0), label="") new("loca.p", x, y, w, label)


#
# loca.p summary method
#
setMethod("summary", "loca.p",
   function(object, ...)
          {
            c("label"=object@label, "n"=length(object@x), "xmin"=min(object@x), "xwmean"=weighted.mean(object@x,object@w), "xmax"=max(object@x), "ymin"=min(object@y), "ywmean"=weighted.mean(object@y,object@w), "ymax"=max(object@y))
            }
          )
          

#
# loca.p print method
#
setMethod("print", "loca.p",
   function(x, ...)
      {
      # To ensure that orloca is included in pot file
      gettext("orloca", domain="R-orloca")
      print(summary(x), ...)
      invisible(x)
      }
)
