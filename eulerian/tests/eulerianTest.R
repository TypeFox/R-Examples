require(graph)
require(eulerian)

g <- new("graphNEL", nodes=as.character(1:10), edgemode="directed")
g <- addEdge(graph=g, from="1", to="10")
g <- addEdge(graph=g, from="2", to="1")
g <- addEdge(graph=g, from="2", to="6")
g <- addEdge(graph=g, from="3", to="2")
g <- addEdge(graph=g, from="4", to="2")
g <- addEdge(graph=g, from="5", to="4")
g <- addEdge(graph=g, from="6", to="5")
g <- addEdge(graph=g, from="6", to="8")
g <- addEdge(graph=g, from="7", to="9")
g <- addEdge(graph=g, from="8", to="7")
g <- addEdge(graph=g, from="9", to="6")
g <- addEdge(graph=g, from="10", to="3")


testNum <- 1
cat("Test-", testNum, ": ", sep="")
has <- hasEulerianPath(g)
msg <- ifelse(has==TRUE, "passed", "failed!!!")
cat(msg, "\n")


cat("Test-", testNum <- testNum + 1, ": ", sep="")
has <- hasEulerianCycle(g)
msg <- ifelse(has==TRUE, "passed", "failed!!!")
cat(msg, "\n")


cat("Test-", testNum <- testNum + 1, ": ", sep="")
epath <- eulerian(g)
msg <- ifelse(length(epath)==numEdges(g)+1, "passed", "failed!!!")
cat(msg, "\n")

cat("##############\n")

###########################################
g <- addEdge(graph=g, from="5", to="6")
##########################################


cat("Test-", testNum <- testNum + 1, ": ", sep="")
has <- hasEulerianPath(g)
msg <- ifelse(has==TRUE, "passed", "failed!!!")
cat(msg, "\n")

cat("Test-", testNum <- testNum + 1, ": ", sep="")
has <- hasEulerianPath(g, "5")
msg <- ifelse(has==TRUE, "passed", "failed!!!")
cat(msg, "\n")

cat("Test-", testNum <- testNum + 1, ": ", sep="")
has <- hasEulerianPath(g, "6")
msg <- ifelse(has==FALSE, "passed", "failed!!!")
cat(msg, "\n")

cat("Test-", testNum <- testNum + 1, ": ", sep="")
has <- hasEulerianCycle(g)
msg <- ifelse(has==FALSE, "passed", "failed!!!")
cat(msg, "\n")

cat("Test-", testNum <- testNum + 1, ": ", sep="")
epath <- eulerian(g, "5")
msg <- ifelse(length(epath)==numEdges(g)+1, "passed", "failed!!!")
cat(msg, "\n")

cat("Test-", testNum <- testNum + 1, ": ", sep="")
epath <- tryCatch(eulerian(g, "7"), error = function(e) NA);
msg <- ifelse(is.na(epath), "passed", "failed!!!")
cat(msg, "\n")

cat("##############\n")

########################################
g <- new("graphNEL", nodes=LETTERS[6:1], edgemode="undirected")
g <- addEdge(graph=g, from=c("A","B","B","C","D"), to=c("B","C","D","E","E"))
########################################


cat("Test-", testNum <- testNum + 1, ": ", sep="")
has <- hasEulerianPath(g)
msg <- ifelse(has==TRUE, "passed", "failed!!!")
cat(msg, "\n")

cat("Test-", testNum <- testNum + 1, ": ", sep="")
has <- hasEulerianPath(g, "A")
msg <- ifelse(has==TRUE, "passed", "failed!!!")
cat(msg, "\n")

cat("Test-", testNum <- testNum + 1, ": ", sep="")
has <- hasEulerianPath(g, "B")
msg <- ifelse(has==TRUE, "passed", "failed!!!")
cat(msg, "\n")

cat("Test-", testNum <- testNum + 1, ": ", sep="")
has <- hasEulerianPath(g, "C")
msg <- ifelse(has==FALSE, "passed", "failed!!!")
cat(msg, "\n")

cat("Test-", testNum <- testNum + 1, ": ", sep="")
has <- hasEulerianPath(g, "F")
msg <- ifelse(has==FALSE, "passed", "failed!!!")
cat(msg, "\n")

cat("Test-", testNum <- testNum + 1, ": ", sep="")
has <- hasEulerianCycle(g)
msg <- ifelse(has==FALSE, "passed", "failed!!!")
cat(msg, "\n")

cat("Test-", testNum <- testNum + 1, ": ", sep="")
epath <- eulerian(g, "B")
msg <- ifelse(length(epath)==numEdges(g)+1 && epath[1]=="B", "passed", "failed!!!")
cat(msg, "\n")

cat("Test-", testNum <- testNum + 1, ": ", sep="")
epath <- tryCatch(eulerian(g, "C"), error = function(e) NA);
msg <- ifelse(is.na(epath), "passed", "failed!!!")
cat(msg, "\n")

cat("###############\n")
########################################
g <- new("graphNEL", nodes=LETTERS[1:4], edgemode="undirected")
g <- addEdge(graph=g, from=c("A","B","C","D"), to=c("B","A","D","C"))
########################################

cat("Test-", testNum <- testNum + 1, ": ", sep="")
has <- hasEulerianPath(g)
msg <- ifelse(has==FALSE, "passed", "failed!!!")
cat(msg, "\n")

cat("Test-", testNum <- testNum + 1, ": ", sep="")
has <- hasEulerianCycle(g)
msg <- ifelse(has==FALSE, "passed", "failed!!!")
cat(msg, "\n")

cat("Test-", testNum <- testNum + 1, ": ", sep="")
epath <- tryCatch(eulerian(g, "C"), error = function(e) NA);
msg <- ifelse(is.na(epath), "passed", "failed!!!")
cat(msg, "\n")


cat("###############\n")
########################################
g <- new("graphNEL", nodes=LETTERS[1:5], edgemode="undirected")
g <- addEdge(graph=g, from=c("A","B","B","C", "C", "D", "C","C"), to=c("B","C","D","E", "C", "E","C","C"))
########################################
cat("Test-", testNum <- testNum + 1, ": ", sep="")
has <- hasEulerianPath(g)
msg <- ifelse(has==TRUE, "passed", "failed!!!")
cat(msg, "\n")

cat("Test-", testNum <- testNum + 1, ": ", sep="")
has <- hasEulerianCycle(g)
msg <- ifelse(has==FALSE, "passed", "failed!!!")
cat(msg, "\n")

cat("Test-", testNum <- testNum + 1, ": ", sep="")
epath <- tryCatch(eulerian(g, "B"), error = function(e) NA);
msg <- ifelse(length(epath)==numEdges(g)+1 && epath[1]=="B", "passed", "failed!!!")
cat(msg, "\n")