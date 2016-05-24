flemming.table<- function(lang=lang){
  tab.leg<-matrix(nrow=29,ncol=2)
  tab.leg<-as.data.frame(tab.leg)
  code<-c("S","A-I","A-II","","B-I","B-II","B-III","B-IV","","C-I","C-II",
          "C-III","C-IV","C-V","C-VI","","D-I","D-II","D-III","D-IV","D-V",
          "D-VI","","E-I","E-II","E-III","E-IV","E-V","E-VI")
  
  if (lang=="en-US" | lang=="en-GR"| lang=="eng"| lang=="e"){
    colnames(tab.leg)<-c("Code","Textural class")
    text.class<-c("Sand","Slightly silty sand","Slightly clayey sand",
                  "","Very silty sand","Silty sand","Clayey sand","Very clayey sand",
                  "","Extremely silty sandy mud","Very silty sandy mud",
                  "Silty sandy mud","Clayey sandy mud","Very clayey sandy mud",
                  "Extremely clayey sandy mud","","Extremely silty slightly sandy mud",
                  "Very silty slightly sandy mud","Silty slightly sandy mud",
                  "Clayey slightly sandy mud","Very clayey slightly sandy mud",
                  "Extremely clayey slightly sandy mud","","Silt","Slightly clayey silt",
                  "Clayey silt","Silty clay","Slightly silty clay","Clay")
  }
  if (lang=="pt-BR" | lang=="pt-PT"| lang=="port"| lang=="p"){
    colnames(tab.leg)<-c("C\u00F3digo","Classe Textural")
    text.class<-c("Areia","Areia ligeiramente s\u00EDltica","Areia ligeiramente argilosa",
                  "","Areia muito s\u00EDltica","Areia s\u00EDltica","Areia argilosa","Areia muito argilosa",
                  "","Lama arenosa extremamente s\u00EDltica","Lama arenosa muito s\u00EDltica",
                  "Lama arenosa s\u00EDltica","Lama arenosa argilosa","Lama arenosa muito argilosa",
                  "Lama arenosa extremamente argilosa","","Lama extremamente s\u00EDltica ligeiramente arenosa",
                  "Lama muito s\u00EDltica ligeiramente arenosa","Lama s\u00EDltica ligeiramente arenosa",
                  "Lama argilosa ligeiramente arenosa","Lama muito argilosa ligeiramente arenosa",
                  "Lama extremamente argilosa ligeiramente arenosa",""," Silte","Silte ligeiramente argiloso",
                  "Silte argiloso","Argila s\u00EDltica","Argila ligeiramente s\u00EDltica","Argila")
  }
  tab.leg[,1]<-code
  tab.leg[,2]<-text.class
  return(tab.leg)
  }



