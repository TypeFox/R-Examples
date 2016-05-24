raw.alignment <-function(s,Syllabic=5,Place=40,Stop=50,Voice=10,Nasal=10,Retroflex=10,Lateral=10,Aspirated=5,Long=1,High=5,Back=5,Round=5,sk=10)
  {
    
    ### Input Arguments ###  
    features<-as.integer(c(Syllabic,Place,Stop,Voice,Nasal,Retroflex,Lateral,Aspirated,Long,High,Back,Round))
    ### C++ interface ####
    z<-.Call("exchange",s,features,as.integer(sk))
    m<-list(s,z[[3]],z[[1]],z[[2]])
    names(m)<-c("word pair","similarity score","alignment1","alignment2")
    return(m)    
  }
