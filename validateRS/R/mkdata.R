p.0<-c(0.0002, 0.0007, 0.0022, 0.0086, 0.0428, 0.2685)
names(p.0)<-paste("rc.", 1:length(p.0), sep="")


# SP.sizes<-rbind( 
#   c(374,  1330,  1637,  1047,  1471,  154), #Global
#   c(100,   563,  1084,   836,  1277,	131), #Nonfinancials
#   c(274,   767,   553,   211,   194,	 23), #Financials
#   c(148,   387,   188,    48,    27,	  6), #Insurance
#   c(126,   380,   365,	 163,	  167,   17)  #Financial Institution
# )

sizes<-rbind( 
  c(374,  1330,  1637,  1047,  1471,  154), #Global
  c(100,   563,  1084,   836,  1277,	131), #Nonfinancials
  c(148,   387,   188,    48,    27,	  6), #Insurance
  c(100,   100,   100,   100,   100,  100), #small
  c(10000,10000,10000,10000,10000,10000)   #large
)

rownames(sizes)<-c("Global (BENCHMARK)", 
                      "Nonfinancials (NFC)", 
                      "Insurance (INS)", 
                      "small", 
                      "large")
colnames(sizes)<-paste("rc.", 1:ncol(sizes), sep="")


# eba.data<-t(read.csv(file="S:\\AM\\ICAS\\Modelling\\ECB Paper\\packages\\eba_based.csv"))
# colnames(eba.data)<-paste("CQS", 1:ncol(eba.data), sep="")
# 
# other.sizes<-rbind(rep(100, times=6), 
#                    rep(10000, times=6))
# rownames(other.sizes)<-paste("scen.", 1:nrow(other.sizes), sep="")
# colnames(other.sizes)<-paste("CQS", 1:ncol(other.sizes), sep="")


ratingData<-list(p.0=p.0, 
                 sizes=sizes)

save(ratingData, file="ratingData.rdata")
