############################################################################
#     MLwiN User Manual
#
# 18  Modelling Cross-classified Data . . . . . . . . . . . . . . . . . .271
#
#     Rasbash, J., Steele, F., Browne, W. J. and Goldstein, H. (2012).
#     A User's Guide to MLwiN, v2.26. Centre for Multilevel Modelling,
#     University of Bristol.
############################################################################
#     R script to replicate all analyses using R2MLwiN
#
#     Zhang, Z., Charlton, C., Parker, R, Leckie, G., and Browne, W.J.
#     Centre for Multilevel Modelling, 2012
#     http://www.bristol.ac.uk/cmm/software/R2MLwiN/
############################################################################

library(R2MLwiN)
# MLwiN folder
mlwin <- getOption("MLwiN_path")
while (!file.access(mlwin, mode = 1) == 0) {
  cat("Please specify the root MLwiN folder or the full path to the MLwiN executable:\n")
  mlwin <- scan(what = character(0), sep = "\n")
  mlwin <- gsub("\\", "/", mlwin, fixed = TRUE)
}
options(MLwiN_path = mlwin)


# 18.1 An introduction to cross-classification . . . . . . . . . . . . . 271

# 18.2 How cross-classified models are implemented in MLwiN . . . . . . .273

# 18.3 Some computational considerations . . . . . . . . . . . . . . . . 273



# 18.4 Modelling a two-way classification: An example . . . . . . . . . .275


data(xc, package = "R2MLwiN")
summary(xc)

sid_dummy <- model.matrix(~factor(xc$sid) - 1)
sid_dummy_names <- paste0("s", 1:19)
colnames(sid_dummy) <- sid_dummy_names
xc <- cbind(xc, sid_dummy)

covmatrix <- NULL
for (i in 1:19) {
  for (j in 1:i) {
    if (i != j) {
      covmatrix <- cbind(covmatrix, c(3, sid_dummy_names[j], sid_dummy_names[i]))
    }
  }
}

random.ui <- matrix(0, 21, 18)
random.ui[1, ] <- 1
random.ui[2:19, ] <- diag(18) * -1
random.ci <- rep(0, 18)

(mymode1 <- runMLwiN(attain ~ 1 + (s1 + s2 + s3 + s4 + s5 +
s6 + s7 + s8 + s9 + s10 + s11 + s12 + s13 + s14 +
  s15 + s16 + s17 + s18 + s19 | cons) + (1 | pid) + (1 | pupil), estoptions = list(clre = covmatrix, constraints = list(random.ui = random.ui,
  random.ci = random.ci)), data = xc))


# 18.5 Other aspects of the SETX command . . . . . . . . . . . . . . . . 277

sid_inter <- sid_dummy * xc$vrq
sid_inter_names <- paste0("s",1:19, "Xvrq")
colnames(sid_inter) <- sid_inter_names
xc <- cbind(xc, sid_inter)

sid_names <- c(sid_dummy_names, sid_inter_names)

covmatrix <- NULL
for (i in 1:38) {
  for (j in 1:i) {
    if (i != j) {
      covmatrix <- cbind(covmatrix, c(3, sid_names[j], sid_names[i]))
    }
  }
}

random.ui <- matrix(0, 40, 36)
random.ui[1, 1:18] <- 1
random.ui[2:19, 1:18] <- diag(18) * -1
random.ui[20, 19:36] <- 1
random.ui[21:38, 19:36] <- diag(18) * -1
random.ci <- rep(0, 36)

(mymodel2 <- runMLwiN(attain ~ 1 + (s1 + s2 + s3 + s4 + s5 +
s6 + s7 + s8 + s9 + s10 + s11 + s12 + s13 + s14 +
  s15 + s16 + s17 + s18 + s19 + s1Xvrq + s2Xvrq + s3Xvrq + s4Xvrq + s5Xvrq + s6Xvrq + s7Xvrq + s8Xvrq + s9Xvrq +
  s10Xvrq + s11Xvrq + s12Xvrq + s13Xvrq + s14Xvrq + s15Xvrq + s16Xvrq + s17Xvrq + s18Xvrq + s19Xvrq | cons) + (1 | pid) +
  (1 | pupil), estoptions = list(clre = covmatrix, constraints = list(random.ui = random.ui, random.ci = random.ci)),
  data = xc))

# Note: The final models in this section of the manual are for demonstration only.
# The models presented in the manual do not converge with the current data.
# We have therefore not given the R2MLwiN commands for these models.

# 18.6 Reducing storage overhead by grouping . . . . . . . . . . . . . . 279

findClust <- function(data, var1, var2) {
  grplist <- by(xc, xc$sid, function(x) x)
  moreclust <- TRUE
  while (moreclust) {
    found <- FALSE
    for (i in (length(grplist) - 1):1) {
      if (length(grplist) == 1) {
        moreclust <- FALSE
        break
      }
      for (j in length(grplist):(i + 1)) {
        if (length(intersect(grplist[[i]][[var2]], grplist[[j]][[var2]])) != 0) {
          grplist[[i]] <- rbind(grplist[[i]], grplist[[j]])
          grplist[[j]] <- NULL
          found <- TRUE
        }
      }
    }
    if (found) {
      break
    } else {
      moreclust <- FALSE
    }
  }
  
  cat(paste0("Number of clusters: ", length(grplist), "\n"))
  
  ids <- NULL
  for (i in 1:length(grplist)) {
    ids <- rbind(ids, cbind(unique(grplist[[i]][[var1]]), i))
  }
  colnames(ids) <- c(var1, "id")
  merge(data[, c(var1, var2)], ids, sort = FALSE)$id
}

if (!require(doBy)) install.packages("doBy")
library(doBy)

data(xc, package = "R2MLwiN")

xc$region <- findClust(xc, "sid", "pid")
xc$region <- NULL

numchildren <- summaryBy(cons ~ sid + pid, FUN = length, data = xc)
colnames(numchildren) <- c("sid", "pid", "numchildren")
xc <- merge(xc, numchildren, sort = FALSE)

xc <- xc[order(xc$sid, xc$pid), ]

xc <- xc[xc$numchildren > 2, ]

xc$region <- findClust(xc, "sid", "pid")

xc <- xc[order(xc$region, xc$sid), ]

xc$rsid <- rep(1, nrow(xc))

rgrp <- xc$region[1]
sgrp <- xc$sid[1]
id <- 1
for (i in 2:nrow(xc)) {
  if (xc$region[i] != rgrp) {
    rgrp <- xc$region[i]
    sgrp <- xc$sid[i]
    id <- 1
  }
  if (xc$sid[i] != sgrp) {
    sgrp <- xc$sid[i]
    id <- id + 1
  }
  xc$rsid[i] <- id
}

rs_dummy <- model.matrix(~factor(xc$rsid) - 1)
rs_dummy_names <- paste0("rs", 1:8)
colnames(rs_dummy) <- rs_dummy_names
xc <- cbind(xc, rs_dummy)

covmatrix <- NULL
for (i in 1:8) {
  for (j in 1:i) {
    if (i != j) {
      covmatrix <- cbind(covmatrix, c(3, rs_dummy_names[j], rs_dummy_names[i]))
    }
  }
}

random.ui <- matrix(0, 10, 7)
random.ui[1, ] <- 1
random.ui[2:8, ] <- diag(7) * -1
random.ci <- rep(0, 7)

xc <- xc[order(xc$region, xc$pid, xc$pupil), ]

(mymode1 <- runMLwiN(attain ~ 1 + (rs1 + rs2 + rs3 + rs4 + rs5 + rs6 + rs7 + rs8 | region) + (1 | pid) + (1 | pupil), 
  estoptions = list(clre = covmatrix, constraints = list(random.ui = random.ui, random.ci = random.ci)), data = xc))

# 18.7 Modelling a multi-way cross-classification . . . . . . . . . . . .280

# 18.8 MLwiN commands for cross-classifications . . . . . . . . . . . . .281

# Chapter learning outcomes . . . . . . . . . . . . . . . . . . . . . . .282

############################################################################
