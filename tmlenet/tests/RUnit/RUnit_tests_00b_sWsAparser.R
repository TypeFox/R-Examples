`%+%` <- function(a, b) paste0(a, b)

#***************************************************************************************
# Testing new parser DefineSummariesClass:
#***************************************************************************************
data(df_netKmax6) # load observed data
head(df_netKmax6)
data(NetInd_mat_Kmax6)  # load the network ID matrix
netind_cl <- simcausal:::NetIndClass$new(nobs = nrow(df_netKmax6), Kmax = 6)
netind_cl$NetInd <- NetInd_mat_Kmax6
head(netind_cl$nF)

# --------------------------------------------------------------------------------------
# test 1 (simple vars, expressions unnamed):
# --------------------------------------------------------------------------------------
# as expressions:
def_sW <- def.sW(W1, W2, W3)

mat1 <- def_sW$eval.nodeforms(data.df = df_netKmax6, netind_cl = netind_cl)
def_sW$sVar.names.map
head(mat1)
checkTrue(all.equal(mat1[,"W1"], df_netKmax6[,"W1"]))
checkTrue(all.equal(mat1[,"W2"], df_netKmax6[,"W2"]))
checkTrue(all.equal(mat1[,"W3"], df_netKmax6[,"W3"]))
checkTrue(all.equal(mat1[,"nF"], netind_cl$nF))

# old parser (deprecated):
# mat1b <- def_sW$get.mat.sVar(data.df = df_netKmax6, netind_cl = netind_cl)
# head(mat1b)
# checkTrue(all.equal(mat1, mat1b)) # passed

def_sW <- def.sW(W1[[0]],W2,W3)

mat2 <- def_sW$eval.nodeforms(data.df = df_netKmax6, netind_cl = netind_cl)
def_sW$sVar.names.map
head(mat2)
checkTrue(all.equal(mat1, mat2))

# old parser (deprecated):
# mat2b <- def_sW$get.mat.sVar(data.df = df_netKmax6, netind_cl = netind_cl)
# head(mat2b)
# checkTrue(all.equal(mat2, mat2b)) # passed

# as character strings:
def_sW <- def.sW("W1", "W2", "W3")

mat2 <- def_sW$eval.nodeforms(data.df = df_netKmax6, netind_cl = netind_cl)
head(mat2)
def_sW$sVar.names.map
checkTrue(all.equal(mat1, mat2))

# old parser (deprecated):
# mat2b <- def_sW$get.mat.sVar(data.df = df_netKmax6, netind_cl = netind_cl)
# checkTrue(all.equal(mat2, mat2b)) # passed

def_sW <- def.sW("W1[[0]]", "W2", "W3")

mat2 <- def_sW$eval.nodeforms(data.df = df_netKmax6, netind_cl = netind_cl)
def_sW$sVar.names.map
head(mat2)
checkTrue(all.equal(mat1, mat2))

# old parser (deprecated):
# mat2b <- def_sW$get.mat.sVar(data.df = df_netKmax6, netind_cl = netind_cl)
# def_sW$sVar.names.map
# head(mat2b)
# (all.equal(mat2, mat2b)) # passed


# --------------------------------------------------------------------------------------
# test 2 (matrix result for named expression):
# --------------------------------------------------------------------------------------
def_sW <- def.sW(netW2 = W2[[0:Kmax]], W3 = W3[[0]])

mat1 <- def_sW$eval.nodeforms(data.df = df_netKmax6, netind_cl = netind_cl)
(map1 <- def_sW$sVar.names.map)
head(mat1)

checkTrue(all.equal(names(map1)[1], "netW2"))
checkTrue(all.equal(mat1[,"W2"], df_netKmax6[,"W2"]))

# old parser (deprecated):
# mat2 <- def_sW$get.mat.sVar(data.df = df_netKmax6, netind_cl = netind_cl, addnFnode = "nF")
# (map2 <- def_sW$sVar.names.map)
# head(mat2)
# colnames(mat2)

# --------------------------------------------------------------------------------------
# test 3 (matrix result with unnamed expression):
# --------------------------------------------------------------------------------------
def_sW <- def.sW(W2[[0:Kmax]], W3 = W3[[0]])

mat1 <- def_sW$eval.nodeforms(data.df = df_netKmax6, netind_cl = netind_cl)
(map1 <- def_sW$sVar.names.map)
head(mat1)

checkTrue(all.equal(colnames(head(mat1)), as.vector(unlist(map1))))
checkTrue(all.equal(names(map1)[1], "W2"))
checkTrue(all.equal(mat1[,"W2"], df_netKmax6[,"W2"]))


# --------------------------------------------------------------------------------------
# test 4 (vector result, complex expression with one parent):
# --------------------------------------------------------------------------------------
# named expression:
def_sW <- def.sW(sum.netW3 = sum(W3[[1:Kmax]]), replaceNAw0 = TRUE)

mat1a <- def_sW$eval.nodeforms(data.df = df_netKmax6, netind_cl = netind_cl)
(map1a <- def_sW$sVar.names.map)
head(mat1a)
checkTrue(all.equal(names(map1a)[1], "sum.netW3"))

# the same unnamed expression:
def_sW <- def.sW(sum(W3[[1:Kmax]]), replaceNAw0 = TRUE)

mat1b <- def_sW$eval.nodeforms(data.df = df_netKmax6, netind_cl = netind_cl)
(map1b <- def_sW$sVar.names.map)
head(mat1b)
checkTrue(all.equal(names(map1b)[1], "W3"))
checkTrue(all.equal(mat1a[,"sum.netW3"], mat1b[,"W3"]))

# --------------------------------------------------------------------------------------
# test 5 (vector result, complex expression with more than one parent):
# --------------------------------------------------------------------------------------
# named expression:
def_sW <- def.sW(sum.netW2W3 = sum(W3[[1:Kmax]]*W2[[1:Kmax]]), replaceNAw0 = TRUE)

mat1a <- def_sW$eval.nodeforms(data.df = df_netKmax6, netind_cl = netind_cl)
(map1a <- def_sW$sVar.names.map)
head(mat1a)
checkTrue(all.equal(names(map1a)[1], "sum.netW2W3"))
checkTrue(all.equal(colnames(mat1a)[1], names(map1a)[1]))

# the same unnamed expression (should throw an exception):
def_sW <- def.sW(sum(W3[[1:Kmax]]*W2[[1:Kmax]]), replaceNAw0 = TRUE)

checkException(mat1b <- def_sW$eval.nodeforms(data.df = df_netKmax6, netind_cl = netind_cl))

# --------------------------------------------------------------------------------------
# test 6 (matrix result, complex expression with more than one parent):
# --------------------------------------------------------------------------------------
# when more than one parent is present, new naming convention for resulting column names:
# sVar.name%+%c(1:ncol)
def_sW <- def.sW(sum.netW2W3 = W3[[1:Kmax]]*W2[[1:Kmax]])

mat1a <- def_sW$eval.nodeforms(data.df = df_netKmax6, netind_cl = netind_cl)
(map1a <- def_sW$sVar.names.map)
head(mat1a)
checkTrue(all.equal(map1a$sum.netW2W3, "sum.netW2W3."%+%c(1:6)))

# the same as unnamed expression (should throw an exception):
def_sW <- def.sW(W3[[1:Kmax]]*W2[[1:Kmax]])
checkException(mat1b <- def_sW$eval.nodeforms(data.df = df_netKmax6, netind_cl = netind_cl))

# --------------------------------------------------------------------------------------
# test 7 (removing duplicate summary measure names):
# --------------------------------------------------------------------------------------
def_sW <- def.sW(netW1 = W1[[1:Kmax]], W1, W2, W3) + def.sW(netW1 = W2[[1:Kmax]])
mat1 <- def_sW$eval.nodeforms(data.df = df_netKmax6, netind_cl = netind_cl)
def_sW$sVar.names.map
head(mat1)

# check that only one summary measure names netW1 exists:
checkTrue(sum(names(def_sW$sVar.names.map)%in%"netW1")==1L)
# check that columns from the first summary measure (W1_netF1, ..., W1_netF6) were removed:
checkTrue(!"W1_netF1"%in%colnames(mat1))
checkTrue(!"W1_netF2"%in%colnames(mat1))
checkTrue(!"W1_netF3"%in%colnames(mat1))
checkTrue(!"W1_netF4"%in%colnames(mat1))
checkTrue(!"W1_netF5"%in%colnames(mat1))
checkTrue(!"W1_netF6"%in%colnames(mat1))

# --------------------------------------------------------------------------------------
# test 8 (removing duplicate column names from evaluation matrix):
# --------------------------------------------------------------------------------------
def_sW <- def.sW(netW1 = W1[[1:Kmax]], W1, W2, W3) + def.sW(netW2 = W1[[1:Kmax]])
mat1 <- def_sW$eval.nodeforms(data.df = df_netKmax6, netind_cl = netind_cl)
def_sW$sVar.names.map
head(mat1)

# check that only one summary measure names netW1 exists:
checkTrue(sum(colnames(mat1) %in% "W1_netF1")==1L)
checkTrue(sum(colnames(mat1) %in% "W1_netF2")==1L)
checkTrue(sum(colnames(mat1) %in% "W1_netF3")==1L)
checkTrue(sum(colnames(mat1) %in% "W1_netF4")==1L)
checkTrue(sum(colnames(mat1) %in% "W1_netF5")==1L)
checkTrue(sum(colnames(mat1) %in% "W1_netF6")==1L)



