index <-
function(
dsgn, data=data, formula=formula
)

{
# definitions
  Terms <- dsgn$Terms
  label <- attr(Terms, "term.labels")
  variables <- as.character(attr(Terms,"variables"))[-c(1)]
  m <- dsgn$m
  int <- dsgn$int
  X <- model.matrix(Terms, m)
  ass <- attr(X,"assign")

# index1 - number of coefficients per covariate
  index1 <- summary(as.factor(ass))
    if (int==0) index1 <- index1[-1]
  names(index1) <- NULL

# more definitions
  # index2 - indicator for special v
  # index3 - indicator for special p
  # index4 - indicator for grouped 
  # index5 - indicator for grouped.fused
  # index6 - indicator for sp/pspline
  # index7 - indicator for SCAD
  # index8 - indicator for elastic
  # index9 - indicator for vspline
    index2 <- index3 <- index4 <- index5 <- index6 <- index7 <- index8 <- index9 <- rep(0, length(index1))
    j <- 1 # label
    r <- 1
  
  # indicator for hierarchical?!
    index.hier <- rep(1,length(index1))
  
  # index.appro - indicates the type of approximation
    index.appro <- rep(NA, length(index1))
    
  # index2b - indicator for b_j in special v()
    index2b <- rep(0,length(index1))  

# index2
  if(!is.null(attr(Terms,"specials")$v)){
      # intercept
      if (index1[1]!=1) {
        u <- model.frame(~eval(parse(text=gsub(")", "", gsub(" ", "",
             strsplit(names(m)[2], "\\,")[[1]][2])))), data)
        erms <- attr(u, "terms")
        if (attr(erms, "dataClasses")[[1]]=="factor") {index2[1] <- c(-1)}
        if (attr(erms, "dataClasses")[[1]]=="ordered"){index2[1] <- c(+1)}
        j <- 2
        r <- 0
      }
  
      # others
      if (length(label)>=j) {
      for (i in j:length(label)) {
        variable <- which(variables==label[i])
        if (length(variable)>0) {if (variable %in% attr(Terms,"specials")$v){
          u <- model.frame(~eval(parse(text=sub(")", "", gsub(" ", "",
               strsplit(label[i], "\\,")[[1]][2])))), data)
          erms <- attr(u, "terms")
          if (attr(erms, "dataClasses")[[1]]=="factor") {index2[i+r] <- -1}
          if (attr(erms, "dataClasses")[[1]]=="ordered"){index2[i+r] <- +1}
          
          index2b[i+r] <- 1
          if(grepl("=F", gsub(pattern=" ", replacement="", x=label[i]))) {index2b[i+r] <- 0}
          
          }}
        } # for
      } # uberhaupt others
  } # uberhaupt


# index3 (p!!) + index.hier
  if(!is.null(attr(Terms,"specials")$p)){
      if (length(label)>=j) {
          for (i in j:length(label)) {
            variable <- which(variables==label[i])
            if (length(variable)>0) {if (variable %in% attr(Terms,"specials")$p){
              u <- model.frame(~eval(parse(text=sub(")", "", gsub(" ", "",
                   strsplit(strsplit(label[i], "p\\(")[[1]][2], ",")[[1]][1])))), data)
              erms <- attr(u, "terms")
              if (attr(erms, "dataClasses")[[1]]=="factor") {index3[i+r] <- -1}
              if (attr(erms, "dataClasses")[[1]]=="ordered"){index3[i+r] <- +1}
              if (attr(erms, "dataClasses")[[1]]=="numeric"){index3[i+r] <- +1}
              }}
    
#           if (length(variable)==0 && grepl(":", label[i]) && !grepl("v\\(", label[i]) &&  # interaktion ohne special*special
#                !grepl("SCAD\\(", label[i]) && !grepl("grouped\\(", label[i]) && !grepl("grouped.fused\\(", label[i]) &&
#                !grepl("sp\\(", label[i]) && !grepl("elastic\\(", label[i])
#               )
#             {
#              eins <- strsplit(label[i], ":")[[1]][1]
#              zwei <- strsplit(label[i], ":")[[1]][2]
#              for (m in c(eins, zwei)){
#              welche <- which(variables==m)
#                if (length(welche)>0) {if (welche %in% attr(Terms,"specials")$p){
#                  u <- model.frame(~eval(parse(text=sub(")", "", gsub(" ", "",
#                       strsplit(strsplit(label[i], "p\\(")[[1]][2], ",")[[1]][1])))), data)
#                  erms <- attr(u, "terms")
#                  if (attr(erms, "dataClasses")[[1]]=="factor") {index3[i+r] <- -1}
#                  if (attr(erms, "dataClasses")[[1]]=="ordered"){index3[i+r] <- +1}
#                  if (attr(erms, "dataClasses")[[1]]=="numeric"){index3[i+r] <- +1}
#    
#                  if (variables[welche] %in% label && attr(erms, "dataClasses")[[1]]!="numeric") # nur int mit p(cat) in index.hier schreiben!!
#                     {index.hier[i+r] <- which(label==variables[welche])+r}
#                  }}
#             }}
# interaktionen ganz verboten, weil p()*numeric entspricht v() und weil p()*cat Unsinn fuer penalty.
    
            } # for
      } # uberhaupt others
  } # uberhaupt

# index4 - indicator for grouped
  if(!is.null(attr(Terms,"specials")$grouped)){
      if (length(label)>=j) {
          for (i in j:length(label)) {
            variable <- which(variables==label[i])
            if (length(variable)>0) {if (variable %in% attr(Terms,"specials")$grouped){
              u <- model.frame(~eval(parse(text=sub(")", "", gsub(" ", "",
                   strsplit(label[i], "grouped\\(")[[1]][2])))), data)
              erms <- attr(u, "terms")
              if (attr(erms, "dataClasses")[[1]]=="factor") {index4[i+r] <- -1} # spielt keine Rolle!!
              if (attr(erms, "dataClasses")[[1]]=="ordered"){index4[i+r] <- +1}
              if (attr(erms, "dataClasses")[[1]]=="numeric"){index4[i+r] <- +1}
              }}
        
            } # for
      } # uberhaupt others
  } # uberhaupt

# index5 - indicator for grouped.fused
  if(!is.null(attr(Terms,"specials")$grouped.fused)){
      if (length(label)>=j) {
          for (i in j:length(label)) {
            variable <- which(variables==label[i])
            if (length(variable)>0) {if (variable %in% attr(Terms,"specials")$grouped.fused){
#              u <- model.frame(~eval(parse(text=sub(")", "", gsub(" ", "",
#                   strsplit(label[i], "grouped.fused\\(")[[1]][2])))), data)                   
              u <- model.frame(~eval(parse(text=sub(")", "", gsub(" ", "",
                   strsplit(label[i], "\\,")[[1]][2])))), data)                                  
              erms <- attr(u, "terms")
              if (attr(erms, "dataClasses")[[1]]=="factor") {index5[i+r] <- -1} 
              if (attr(erms, "dataClasses")[[1]]=="ordered"){index5[i+r] <- +1}
              }}
    
            } # for
      } # uberhaupt others
  } # uberhaupt

# index6 - indicator for spline
  if(!is.null(attr(Terms,"specials")$sp)){
      if (length(label)>=j) {
          for (i in j:length(label)) {
            variable <- which(variables==label[i])
            if (length(variable)>0) {if (variable %in% attr(Terms,"specials")$sp){
              index6[i+r] <- 1
              }}
    
            } # for
      } # uberhaupt others
  } # uberhaupt

# index7 - SCAD - wie p!!
  if(!is.null(attr(Terms,"specials")$SCAD)){
      if (length(label)>=j) {
          for (i in j:length(label)) {
            variable <- which(variables==label[i])
            if (length(variable)>0) {if (variable %in% attr(Terms,"specials")$SCAD){
              u <- model.frame(~eval(parse(text=sub(")", "", gsub(" ", "",
                   strsplit(strsplit(label[i], "SCAD\\(")[[1]][2], ",")[[1]][1])))), data)
              erms <- attr(u, "terms")
              if (attr(erms, "dataClasses")[[1]]=="factor") {index7[i+r] <- -1}
              if (attr(erms, "dataClasses")[[1]]=="ordered"){index7[i+r] <- +1}
              if (attr(erms, "dataClasses")[[1]]=="numeric"){index7[i+r] <- +1}
              }}
    
            } # for
      } # uberhaupt others
  } # uberhaupt

# index8 - elastic - wie p!!
  if(!is.null(attr(Terms,"specials")$elastic)){
      if (length(label)>=j) {
          for (i in j:length(label)) {
            variable <- which(variables==label[i])
            if (length(variable)>0) {if (variable %in% attr(Terms,"specials")$elastic){
              u <- model.frame(~eval(parse(text=sub(")", "", gsub(" ", "",
                   strsplit(strsplit(label[i], "elastic\\(")[[1]][2], ",")[[1]][1])))), data)
              erms <- attr(u, "terms")
              if (attr(erms, "dataClasses")[[1]]=="factor") {index8[i+r] <- -1}
              if (attr(erms, "dataClasses")[[1]]=="ordered"){index8[i+r] <- +1}
              if (attr(erms, "dataClasses")[[1]]=="numeric"){index8[i+r] <- +1}
              }}
    
            } # for
      } # uberhaupt others
  } # uberhaupt
  
# index 9 - vspline
  if(!is.null(attr(Terms,"specials")$vspline)){
      if (length(label)>=j) {
      for (i in j:length(label)) {
        variable <- which(variables==label[i])
        if (length(variable)>0) {if (variable %in% attr(Terms,"specials")$vspline){
          u <- model.frame(~eval(parse(text=sub(")", "", gsub(" ", "",
               strsplit(label[i], "\\,")[[1]][2])))), data)
          erms <- attr(u, "terms")
          index2b[i+r] <- ncol(model.matrix(erms, u)) # usual spline, number levels u
          if (attr(erms, "dataClasses")[[1]]=="factor") {index9[i+r] <- -1}
          if (attr(erms, "dataClasses")[[1]]=="ordered"){index9[i+r] <- +1}
          
          }}
        } # for
      } # uberhaupt others
  } # uberhaupt  

# index.appro
for (i in 1:(length(label))) {
  sp <- strsplit(label[[i]], "(", fixed=TRUE)[[1]][1]
  if(sp %in% c("p", "v")) {
     index.appro[i+int] <- if (grepl("\"", label[[i]])) strsplit(label[[i]], "\"")[[1]][2] else "L1"
     } else {
              if(sp == "sp") {
                 index.appro[i+int] <- if (grepl("\"", label[[i]]))  strsplit(label[[i]], "\"")[[1]][2] else "L2"
                 } else  {index.appro[i+int] <- sp}
     }
}


# return
output <- rbind(
   index1,
   index2,
   index2b,
   index3,
   index4,
   index5,
   index6,
   index7,
   index8,
   index9
)
rownames(output) <- c("index1", "index2", "index2b", "index3", "index4", 
                      "index5", "index6", "index7", "index8", "index9")
colnames(output) <- index.appro
return(output)

}
