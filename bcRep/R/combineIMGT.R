## Julia Bischof
## 15-10-2015

combineIMGT<-function(folders=NULL,name=NULL){
  if(length(folders)<2){
    stop("--> At least 2 folders are required")
  }
  if(length(name)!=1){
    stop("--> Project name for combined files is required")
  }
  
  print(paste("Reading ",folders[1]," (",1,"/",length(folders),")",sep=""))
  if(!file.exists(as.character(name))){
    dir.create(file.path(as.character(name)))
  }
  
  if(length(grep("unix",.Platform$OS.type))>0){
    system(paste("cat ",folders[1],"/",list.files(folders[1],pattern = "1_Summary")," > ",name,"/1_Summary.txt",sep=""))
    system(paste("cat ",folders[1],"/",list.files(folders[1],pattern = "2_IMGT-gapped-nt-sequences")," > ",name,"/2_IMGT-gapped-nt-sequences.txt",sep=""))
    system(paste("cat ",folders[1],"/",list.files(folders[1],pattern = "3_Nt-sequences")," > ",name,"/3_Nt-sequences.txt",sep=""))
    system(paste("cat ",folders[1],"/",list.files(folders[1],pattern = "4_IMGT-gapped-AA-sequences")," > ",name,"/4_IMGT-gapped-AA-sequences.txt",sep=""))
    system(paste("cat ",folders[1],"/",list.files(folders[1],pattern = "5_AA-sequences")," > ",name,"/5_AA-sequences.txt",sep=""))
    system(paste("cat ",folders[1],"/",list.files(folders[1],pattern = "6_Junction")," > ",name,"/6_Junction.txt",sep=""))
    system(paste("cat ",folders[1],"/",list.files(folders[1],pattern = "7_V-REGION-mutation-and-AA-change-table")," > ",name,"/7_V-REGION-mutation-and-AA-change-table.txt",sep=""))
    system(paste("cat ",folders[1],"/",list.files(folders[1],pattern = "8_V-REGION-nt-mutation-statistics")," > ",name,"/8_V-REGION-nt-mutation-statistics.txt",sep=""))
    system(paste("cat ",folders[1],"/",list.files(folders[1],pattern = "9_V-REGION-AA-change-statistics")," > ",name,"/9_V-REGION-AA-change-statistics.txt",sep=""))
    system(paste("cat ",folders[1],"/",list.files(folders[1],pattern = "10_V-REGION-mutation-hotspots")," > ",name,"/10_V-REGION-mutation-hotspots.txt",sep=""))
    
    for(i in 2:length(folders)){
      print(paste("Reading ",folders[i]," (",i,"/",length(folders),")",sep=""))    
      system(paste("cat ",folders[i],"/",list.files(folders[i],pattern = "1_Summary")," | tail -n+2 >> ",name,"/1_Summary.txt",sep=""))
      system(paste("cat ",folders[i],"/",list.files(folders[i],pattern = "2_IMGT-gapped-nt-sequences")," | tail -n+2 >> ",name,"/2_IMGT-gapped-nt-sequences.txt",sep=""))
      system(paste("cat ",folders[i],"/",list.files(folders[i],pattern = "3_Nt-sequences")," | tail -n+2 >> ",name,"/3_Nt-sequences.txt",sep=""))
      system(paste("cat ",folders[i],"/",list.files(folders[i],pattern = "4_IMGT-gapped-AA-sequences")," | tail -n+2 >> ",name,"/4_IMGT-gapped-AA-sequences.txt",sep=""))
      system(paste("cat ",folders[i],"/",list.files(folders[i],pattern = "5_AA-sequences")," | tail -n+2 >> ",name,"/5_AA-sequences.txt",sep=""))
      system(paste("cat ",folders[i],"/",list.files(folders[i],pattern = "6_Junction")," | tail -n+2 >> ",name,"/6_Junction.txt",sep=""))
      system(paste("cat ",folders[i],"/",list.files(folders[i],pattern = "7_V-REGION-mutation-and-AA-change-table")," | tail -n+2 >> ",name,"/7_V-REGION-mutation-and-AA-change-table.txt",sep=""))
      system(paste("cat ",folders[i],"/",list.files(folders[i],pattern = "8_V-REGION-nt-mutation-statistics")," | tail -n+2 >> ",name,"/8_V-REGION-nt-mutation-statistics.txt",sep=""))
      system(paste("cat ",folders[i],"/",list.files(folders[i],pattern = "9_V-REGION-AA-change-statistics")," | tail -n+2 >> ",name,"/9_V-REGION-AA-change-statistics.txt",sep=""))
      system(paste("cat ",folders[i],"/",list.files(folders[i],pattern = "10_V-REGION-mutation-hotspots")," | tail -n+2 >> ",name,"/10_V-REGION-mutation-hotspots.txt",sep=""))
    }
  }else if(length(grep("windows",.Platform$OS.type))>0){
    shell(paste("type ",folders[1],"\\",list.files(folders[1],pattern = "1_Summary")," > ",name,"\\1_Summary.txt",sep=""))
    shell(paste("type ",folders[1],"\\",list.files(folders[1],pattern = "2_IMGT-gapped-nt-sequences")," > ",name,"\\2_IMGT-gapped-nt-sequences.txt",sep=""))
    shell(paste("type ",folders[1],"\\",list.files(folders[1],pattern = "3_Nt-sequences")," > ",name,"\\3_Nt-sequences.txt",sep=""))
    shell(paste("type ",folders[1],"\\",list.files(folders[1],pattern = "4_IMGT-gapped-AA-sequences")," > ",name,"\\4_IMGT-gapped-AA-sequences.txt",sep=""))
    shell(paste("type ",folders[1],"\\",list.files(folders[1],pattern = "5_AA-sequences")," > ",name,"\\5_AA-sequences.txt",sep=""))
    shell(paste("type ",folders[1],"\\",list.files(folders[1],pattern = "6_Junction")," > ",name,"\\6_Junction.txt",sep=""))
    shell(paste("type ",folders[1],"\\",list.files(folders[1],pattern = "7_V-REGION-mutation-and-AA-change-table")," > ",name,"\\7_V-REGION-mutation-and-AA-change-table.txt",sep=""))
    shell(paste("type ",folders[1],"\\",list.files(folders[1],pattern = "8_V-REGION-nt-mutation-statistics")," > ",name,"\\8_V-REGION-nt-mutation-statistics.txt",sep=""))
    shell(paste("type ",folders[1],"\\",list.files(folders[1],pattern = "9_V-REGION-AA-change-statistics")," > ",name,"\\9_V-REGION-AA-change-statistics.txt",sep=""))
    shell(paste("type ",folders[1],"\\",list.files(folders[1],pattern = "10_V-REGION-mutation-hotspots")," > ",name,"\\10_V-REGION-mutation-hotspots.txt",sep=""))
    
    for(i in 2:length(folders)){
      print(paste("Reading ",folders[i]," (",i,"/",length(folders),")",sep=""))    
      shell(paste("more +1 ",folders[i],"\\",list.files(folders[i],pattern = "1_Summary")," >> ",name,"\\1_Summary.txt",sep=""))
      shell(paste("more +1 ",folders[i],"\\",list.files(folders[i],pattern = "2_IMGT-gapped-nt-sequences")," >> ",name,"\\2_IMGT-gapped-nt-sequences.txt",sep=""))
      shell(paste("more +1 ",folders[i],"\\",list.files(folders[i],pattern = "3_Nt-sequences")," >> ",name,"\\3_Nt-sequences.txt",sep=""))
      shell(paste("more +1 ",folders[i],"\\",list.files(folders[i],pattern = "4_IMGT-gapped-AA-sequences")," >> ",name,"\\4_IMGT-gapped-AA-sequences.txt",sep=""))
      shell(paste("more +1 ",folders[i],"\\",list.files(folders[i],pattern = "5_AA-sequences")," >> ",name,"\\5_AA-sequences.txt",sep=""))
      shell(paste("more +1 ",folders[i],"\\",list.files(folders[i],pattern = "6_Junction")," >> ",name,"\\6_Junction.txt",sep=""))
      shell(paste("more +1 ",folders[i],"\\",list.files(folders[i],pattern = "7_V-REGION-mutation-and-AA-change-table")," >> ",name,"\\7_V-REGION-mutation-and-AA-change-table.txt",sep=""))
      shell(paste("more +1 ",folders[i],"\\",list.files(folders[i],pattern = "8_V-REGION-nt-mutation-statistics")," >> ",name,"\\8_V-REGION-nt-mutation-statistics.txt",sep=""))
      shell(paste("more +1 ",folders[i],"\\",list.files(folders[i],pattern = "9_V-REGION-AA-change-statistics")," >> ",name,"\\9_V-REGION-AA-change-statistics.txt",sep=""))
      shell(paste("more +1 ",folders[i],"\\",list.files(folders[i],pattern = "10_V-REGION-mutation-hotspots")," >> ",name,"\\10_V-REGION-mutation-hotspots.txt",sep=""))
    }
  }
}

