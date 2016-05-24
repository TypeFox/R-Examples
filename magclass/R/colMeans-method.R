setMethod("colMeans",
          signature(x = "magpie"),
          function (x, na.rm = FALSE, dims = 1, ...) 
          {
            x_array<-as.array(x)
            x_glo<-colMeans(x_array,na.rm=na.rm,...)
            x<-new("magpie",array(x_glo,dim=c(1,dim(x_glo)),dimnames=c("GLO",dimnames(x_glo))))
            return(x)
          }
          )

