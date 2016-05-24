### R code from vignette source 'rEMM.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: rEMM.Rnw:127-130
###################################################
options(width = 70, prompt="R> ", digits=4)
### for sampling
set.seed(1234)


###################################################
### code chunk number 2: rEMM.Rnw:848-851
###################################################
library("rEMM")
data(EMMTraffic)
EMMTraffic


###################################################
### code chunk number 3: rEMM.Rnw:863-867
###################################################
emm <- EMM(threshold=0.2, measure="eJaccard")
build(emm, EMMTraffic)
size(emm)
ntransitions(emm)


###################################################
### code chunk number 4: rEMM.Rnw:877-878
###################################################
cluster_counts(emm)


###################################################
### code chunk number 5: rEMM.Rnw:883-884
###################################################
cluster_centers(emm)


###################################################
### code chunk number 6: Traffic_graph
###################################################
plot(emm, method="graph")


###################################################
### code chunk number 7: rEMM.Rnw:915-916
###################################################
transition_matrix(emm)


###################################################
### code chunk number 8: rEMM.Rnw:920-922
###################################################
transition_matrix(emm, type="counts")
#transition_matrix(emm, type="log_odds")


###################################################
### code chunk number 9: rEMM.Rnw:934-935
###################################################
transition(emm, "2", "1", type="probability")


###################################################
### code chunk number 10: rEMM.Rnw:942-943
###################################################
predict(emm, n=2, current="2")


###################################################
### code chunk number 11: rEMM.Rnw:948-949
###################################################
predict(emm, n=2, current="2", probabilities=TRUE)


###################################################
### code chunk number 12: Traffic_r3
###################################################
emm_3removed <- remove_clusters(emm, "3")
plot(emm_3removed, method="graph")


###################################################
### code chunk number 13: Traffic_rt52
###################################################
emm_52removed <- remove_transitions(emm, "5", "2")
plot(emm_52removed, method="graph")


###################################################
### code chunk number 14: Traffic_m25
###################################################
emm_25merged <- merge_clusters(emm, c("2","5"))
plot(emm_25merged, method="graph")


###################################################
### code chunk number 15: Traffic_l
###################################################
emm_fading <- EMM(threshold=0.2, measure="eJaccard", lambda = 1)
build(emm_fading, EMMTraffic)
plot(emm_fading, method="graph")


###################################################
### code chunk number 16: Traffic_lp
###################################################
emm_fading_pruned <- prune(emm_fading, count_threshold=0.1,
    clusters=TRUE, transitions=TRUE)
plot(emm_fading_pruned, method="graph")


###################################################
### code chunk number 17: rEMM.Rnw:1106-1107
###################################################
data("EMMsim")


###################################################
### code chunk number 18: sim_data
###################################################
plot(EMMsim_train, col="gray", pch=EMMsim_sequence_train)
lines(EMMsim_test, col ="gray")
points(EMMsim_test, col="red", pch=5)
text(EMMsim_test, labels=1:nrow(EMMsim_test), pos=3)


###################################################
### code chunk number 19: sim_graph
###################################################
emm <- EMM(threshold=0.1, measure="euclidean")
build(emm, EMMsim_train)
plot(emm)


###################################################
### code chunk number 20: sim_graphviz
###################################################
plot(emm, method="graph")


###################################################
### code chunk number 21: sim_MDS
###################################################
plot(emm, method="MDS")


###################################################
### code chunk number 22: sim_MDS2
###################################################
plot(emm, method = "MDS", data=EMMsim_train)


###################################################
### code chunk number 23: simil
###################################################
x <- seq(0,5, length.out=50)
plot(x, rEMM:::.simil_weight(x,1), type="l", xlab="d(x,s)/t", ylab="simil(x,s)")


###################################################
### code chunk number 24: rEMM.Rnw:1342-1347
###################################################
score(emm, EMMsim_test, method="log_loss")
score(emm, EMMsim_test, method="likelihood")
score(emm, EMMsim_test, method="product")
score(emm, EMMsim_test, method="sum")
score(emm, EMMsim_test, method="supported_transitions")


###################################################
### code chunk number 25: rEMM.Rnw:1362-1363
###################################################
transition_table(emm, EMMsim_test)


###################################################
### code chunk number 26: rEMM.Rnw:1382-1384
###################################################
score(emm, EMMsim_test, method="product", match_cluster="nn")
score(emm, EMMsim_test, method="product", match_cluster="weighted")


###################################################
### code chunk number 27: rEMM.Rnw:1397-1398
###################################################
score(emm, EMMsim_test, method="supported_transitions", match_cluster=1.1)


###################################################
### code chunk number 28: rEMM.Rnw:1413-1416
###################################################
methods <- c("product", "sum", "log_loss", "likelihood")
sapply(methods, FUN = function(m) 
    score(emm, EMMsim_test, method=m, match="weighted"))


###################################################
### code chunk number 29: sim_hc
###################################################
## find best predicting model (clustering)
k <- 2:10
emmc <- recluster_hclust(emm, k=k, method ="average") 
plot(attr(emmc, "cluster_info")$dendrogram)


###################################################
### code chunk number 30: rEMM.Rnw:1453-1456
###################################################
sc <- sapply(emmc, score, EMMsim_test, "log_likelihood")
names(sc) <- k
sc


###################################################
### code chunk number 31: sim_optc_graph
###################################################
plot(emmc[[which.max(sc)]], method="MDS")


###################################################
### code chunk number 32: sim_optc_MDS
###################################################
plot(emmc[[which.max(sc)]], method="MDS", data=EMMsim_train)


###################################################
### code chunk number 33: rEMM.Rnw:1508-1510
###################################################
data(Derwent)
summary(Derwent)


###################################################
### code chunk number 34: Derwent1
###################################################
plot(Derwent[,1], type="l", ylab="Gauged flow", 
main=colnames(Derwent)[1])


###################################################
### code chunk number 35: Derwent_cluster_counts
###################################################
Derwent_scaled <- scale(Derwent)
emm <- EMM(measure="euclidean", threshold=3)
build(emm, Derwent_scaled)
cluster_counts(emm)
cluster_centers(emm)
plot(emm, method = "cluster_counts", log="y")


###################################################
### code chunk number 36: Derwent_EMM1
###################################################
plot(emm, method="MDS")


###################################################
### code chunk number 37: Derwent_EMM2
###################################################
rare_threshold <- sum(cluster_counts(emm))*0.005
rare_threshold
plot(prune(emm, rare_threshold), method="MDS")


###################################################
### code chunk number 38: Derwent2
###################################################
catchment <- 1 
plot(Derwent[,catchment], type="l", ylab="Gauged flows", 
main=colnames(Derwent)[catchment])
state_sequence <- find_clusters(emm, Derwent_scaled)

mark_states <- function(states, state_sequence, ys, col=0, label=NULL, ...) {
    x <- which(state_sequence %in% states)
    points(x, ys[x], col=col, ...)
    if(!is.null(label)) text(x, ys[x], label, pos=4, col=col)
}

mark_states("11", state_sequence, Derwent[,catchment], col="blue", label="11")
mark_states("12", state_sequence, Derwent[,catchment], col="red", label="12")


###################################################
### code chunk number 39: Derwent3
###################################################
catchment <- 6 
plot(Derwent[,catchment], type="l", ylab="Gauged flow", 
main=colnames(Derwent)[catchment])

mark_states("11", state_sequence, Derwent[,catchment], col="blue", label="11")
mark_states("12", state_sequence, Derwent[,catchment], col="red", label="12")


###################################################
### code chunk number 40: rEMM.Rnw:1746-1750
###################################################
data("16S")

emm <- EMM(threshold=0.1, "Kullback")
build(emm, Mollicutes16S+1)


###################################################
### code chunk number 41: Mollicutes_graph
###################################################
plot(emm, method = "graph")
## start state for sequences have an initial state probability >0
it <- initial_transition(emm)
it[it>0]


