library(d3Network)

Nodes <- c(
  "Head", "Sternum", "Shoulder.r", "Shoulder.l", # 0 1 2 3
  "Elbow.r", "Elbow.l", "Wrist.r", "Wrist.l", # 4 5 6 7
  "rib.r", "hip.c", "hip.r", "hip.l", # 8 9 10 11
  "knee.r", "knee.l", "ankle.r", "ankle.l", # 12 13 14 15
  "rib.l") # 16

Nodes.df <- data.frame(name=Nodes, group=rep(1,length(Nodes)))

Con <- rbind(
  c(0, 1, 300),
  c(1, 2, 50),
  c(1, 3, 50),
  c(2, 4, 20),
  c(3, 5, 20),
  c(4, 6, 1),
  c(5, 7, 1),
  c(8, 9, 50),
  c(9, 10, 20),
  c(9, 11, 20),
  c(10, 12, 10),
  c(11, 13, 10),
  c(12, 14, 5),
  c(13, 15, 5),
  c(2, 8, 50),
  c(3, 16, 50),
  c(1, 8, 50),
  c(1, 16, 50),
  c(8, 16, 50),
  c(16, 9, 50))

colnames(Con) <- c("source", "target", "value")
Con <- as.data.frame(Con)

d3ForceNetwork(Links = Con, Nodes = Nodes.df, Source = "source",
               Target = "target", Value = "value", NodeID = "name",
               Group = "group", opacity = 0.4, file="human.html")
