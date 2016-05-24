sendFA <- function(FAseq, database, hl_size=20, program="blastn", filter="L", expect=10, email, method="POST"){
 if(method=="PUT"){
   blastRun <- scan(paste("http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?QUERY=",FAseq,"&DATABASE=",database,"&HITLIST_SIZE=",hl_size,"&FILTER=",filter,"&EXPECT=",expect,"&PROGRAM=",program,"&CLIENT=web&SERVICE=plain&NCBI_GI=on&PAGE=Nucleotides&CMD=Put&EMAIL=",email,"&TOOL=hoardeR" ,sep=""),what="raw")
   RIDNo <- which((blastRun=="RID")==TRUE)
   RIDNo <- blastRun[RIDNo+2]
 } else if(method=="POST"){  
   tries <- 0
   newError <- TRUE
   while(newError & tries < 6){
     newError <- FALSE
     tryCatch(post <- POST(url="http://www.ncbi.nlm.nih.gov/blast/Blast.cgi",
                           body=(list(QUERY=FAseq,
                                      DATABASE=database,
                                      HITLIST_SIZE=hl_size,
                                      FILTER=filter,
                                      EXPECT=expect,
                                      PROGRAM=program,
                                      CLIENT="web",
                                      SERVICE="plain",
                                      NCBI_GI="on",
                                      PAGE="Nucleotides",
                                      CMD="Put",
                                      EMAIL=email,
                                      TOOL="R"))
                   )
               , error = function(e){
                         cat("An error occured, try",tries,"\n")  
                         Sys.sleep(10)
                         newError <<- TRUE
                         tries <<- tries + 1
                         }
             )
   }
   blastRun <- read.csv(text=rawToChar(post$content),header=T,sep="\t", stringsAsFactor=FALSE)
   RIDNo <- blastRun[grepl("RID =",blastRun[,1]),1]
   RIDNo <- strsplit(RIDNo," = ")[[1]][2]
 }
 RIDNo
}