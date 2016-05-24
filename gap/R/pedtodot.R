pedtodot <- function(pedfile,makeped=FALSE,sink=TRUE,page="B5",
            url="http://www.mrc-epid.cam.ac.uk",height=0.5,width=0.75,rotate=0,dir="none")
{
  if (makeped) ped <- pedfile[,-c(5,6,7,9)]
  else ped <- pedfile
  pedigree <- ped[,1]
  member <- ped[,2]
  father <- ped[,3]
  mother <- ped[,4]
  sex <- ped[,5]
  aff <- ped[,6]
  page.int <- charmatch(page,c("A4","A5","B5","Legal","Letter","Executive"))
  pagesize <- c("8.2677165,11.692913",
                "5.83,8.27",
                "7.17,10.12",
                "8.5,14",
                "8.5,11",
                "7.25,10.5") 
  ashape <- matrix(c(
             "m","box,regular=1",
             "1","box,regular=1",
             "f","circle",
             "2","circle"),ncol=2,byrow=T)
  ashade <- matrix(c(
             "y","style=filled,color=grey",
             "2","style=filled,color=grey",
             "n","style=\"setlinewidth(2)\"",
             "1","style=\"setlinewidth(2)\"",
             "x","green",
             "0","green"),ncol=2,byrow=T)
  ssize <- dim(ped)[1]
  shape <- shade <- rep('1',ssize)
  for (s in 1:ssize) {
      for (t in 1:4) if (sex[s]==ashape[t,1]) shape[s] <- ashape[t,2]
      for (t in 1:6) if (aff[s]==ashade[t,1]) shade[s] <- ashade[t,2]
  }
  uid <- unique(pedigree)
  for (j in 1:length(uid))
  {
    if(sink) cat(paste("[",uid[j],"]",sep=""))
    if(sink) sink(paste(uid[j],".dot",sep=""))
    cat(paste("digraph ped_",uid[j],sep=""),"{\n")
    if (page!="") {
       if (is.na(page.int)) cat(paste("page=\"", page, "\"",sep="")," ;\n")
       else if (page.int>0) cat(paste("page=\"", pagesize[page.int], "\"",sep="")," ;\n")
    }
    cat("ratio=\"auto\" ;\n")
    cat("mincross = 2.0 ;\n")
    cat("label=\"pedigree",uid[j],"\" ;\n")
    cat(paste("rotate=",rotate,sep="")," ;\n")
    if(url!="") cat(paste("URL=\"",url,"\"",sep="")," ;\n")
    selected <- pedigree==uid[j]
    id.j <- member[selected]
    dad.j <- father[selected]
    mom.j <- mother[selected]
    sex.j <- sex[selected]
    aff.j <- aff[selected]
  # Pedigree diagrams via kinship:
  # ped.j <- pedigree(id=id.j, dadid=dad.j, momid=mom.j, sex=sex.j, affected=aff.j)
  # plot(ped.j)
    shape.j <- shape[selected]
    shade.j <- shade[selected]
    n <- length(id.j)
    for (s in 1:n) cat(paste("\"", id.j[s], "\" [shape=", sep=""), shape.j[s], 
                      ",height=", height, ",width=",width, shade.j[s], "] ;\n")
    fid <- match(dad.j,id.j)
    mid <- match(mom.j,id.j)
    fid <- fid[!is.na(fid)]
    mid <- mid[!is.na(mid)]
    marriage <- matrix(rep(0,3*n*(n+1)/2),ncol=3)
    child <- array(rep('0',n*n*(n+1)/2+2),dim=c(n*(n+1)/2,n+2))
    k <- 1
    for (s in 1:n) {
       s1 <- fid[k]
       s2 <- mid[k]
       l <- min(s1,s2)
       u <- max(s1,s2)
       if (dad.j[s]!="x" && dad.j[s]!="0") {
          loc <- u*(u-1)/2 + l
          marriage[loc,1] <- s1
          marriage[loc,2] <- s2
          marriage[loc,3] <- marriage[loc,3] + 1
        # child[loc,1] <- s1
        # child[loc,2] <- s2
          child[loc,marriage[loc,3]+2] <- id.j[s]
          k <- k + 1
       }
    }
    marriage <- as.data.frame(marriage)
    child <- as.data.frame(child)
    married <- marriage[marriage[,3]>0,]
    n <- dim(married)[1]
    for (m in 1:n) {
        s1 <- married[m,1]
        s2 <- married[m,2]
        l <- min(s1,s2)
        u <- max(s1,s2)
        loc <- u*(u-1)/2 + l
        s1 <- id.j[s1]
        s2 <- id.j[s2]
        mating <- paste("\"", s1, "x", s2, "\"",sep="")
        cat(mating, "[shape=diamond,style=filled,label=\"\",height=.1,width=.1] ;\n")
        cat(paste("\"", s1, "\"",sep="")," -> ", mating, paste(" [dir=",dir, ",weight=1]",sep="")," ;\n")
        cat(paste("\"", s2, "\"",sep="")," -> ", mating, paste(" [dir=",dir, ",weight=1]",sep="")," ;\n")
        for (k in 1:married[m,3]) {
            cat(mating, " -> ",paste("\"", child[loc,k+2], "\"",sep=""), paste("[dir=",dir, ",weight=2]",sep="")," ;\n")
        }
    }
    cat("}\n")
    if(sink) sink()
  }
  cat("\n")
}

# History
# 02/01/2005 start experiment
# 03/01/2005 keep pedtodot in .Rd file with further work
# 04/01/2005 success with the use of 2-d array in no need of sort, put to gap
# 10/03/2008 use of setlinewidth(2) for unaffected
#
# The R/S program for GAW14 was originally written for Xiaoyan Liu's data from Stata
# the data with the following variables as required in genassoc by David Clayton
# pedigree, member, father, mother, sex, affected
# pedigrees 27, 106 have loops which can be broken via individuals 4, 5
# i.e., pedigrees 10051, 10065, individuals 10000529, 10001161
#
