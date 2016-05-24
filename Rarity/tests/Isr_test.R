library(Rarity)
data(spid.occ)

# Preparation of rarity weights
rarity.weights <- rWeights(spid.occ, extended = TRUE)

# Generation of an assemblage matrix
assemblages.matrix <- cbind(assemblage.1 = sample(c(0, 1), 50, replace = TRUE),
                            assemblage.2 = sample(c(0, 1), 50, replace = TRUE),
                            assemblage.3 = sample(c(0, 1), 50, replace = TRUE),
                            assemblage.4 = sample(c(0, 1), 50, replace = TRUE),
                            assemblage.5 = sample(c(0, 1), 50, replace = TRUE))
# Random attribution of names to the sampled species
rownames(assemblages.matrix) <- sample(rownames(spid.occ), 50, replace = FALSE)


# Test 1: A vector - W vector
curW <- rarity.weights$W
names(curW) <- rownames(rarity.weights)
curA <- assemblages.matrix[, 1]

result <- Isr(W = curW, assemblages = curA)
result
tot.res <- result[1]
if(any(is.na(result)) | any(result[grep("Isr", names(result))] < 0))
{
  stop("Error in the test")
}
# Testing NAs
curW[2] <- NA
result <- Isr(W = curW, assemblages = curA)
result
if(any(is.na(result)) | any(result[grep("Isr", names(result))] < 0))
{
  stop("Error in the test")
}


# Test 2: A vector - W matrix
curW <- as.matrix(rarity.weights)
curA <- assemblages.matrix[, 1]

result <- Isr(W = curW, assemblages = curA)
result
if(any(is.na(result)) | any(result[grep("Isr", names(result))] < 0))
{
  stop("Error in the test")
}
tot.res <- c(tot.res,
             result["Isr_W"])

# Testing NAs
curW[1, 2] <- NA
curW[10, 2] <- NA
curW[5, 3] <- NA
result <- Isr(W = curW, assemblages = curA)
result
if(any(is.na(result)) | any(result[grep("Isr", names(result))] < 0))
{
  stop("Error in the test")
}


# Test 3: A matrix - W vector
curW <- rarity.weights$W
names(curW) <- rownames(rarity.weights)
curA <- assemblages.matrix

result <- Isr(W = curW, assemblages = curA)
result
tot.res <- c(tot.res,
             result[1, 1])
if(any(is.na(result)) | any(result[grep("Isr", names(result))] < 0))
{
  stop("Error in the test")
}

# Testing NAs
curW[2] <- NA
result <- Isr(W = curW, assemblages = curA)
result
if(any(is.na(result)) | any(result[grep("Isr", names(result))] < 0))
{
  stop("Error in the test")
}

# Test 4: A matrix - W matrix
curW <- rarity.weights
curA <- assemblages.matrix

result <- Isr(W = curW, assemblages = curA)
result
tot.res <- c(tot.res,
             result[1, 3])
if(any(is.na(result)) | any(result[grep("Isr", names(result))] < 0))
{
  stop("Error in the test")
}

# Testing NAs
curW[1, 2] <- NA
curW[10, 2] <- NA
curW[5, 3] <- NA
result <- Isr(W = curW, assemblages = curA)
result
if(any(is.na(result)) | any(result[grep("Isr", names(result))] < 0))
{
  stop("Error in the test")
}


tot.res




#### TEST WITH ABUNDANCES ####

# Generation of an assemblage matrix
assemblages.matrix <- cbind(assemblage.1 = sample(0:25, 50, replace = TRUE),
                            assemblage.2 = sample(0:25, 50, replace = TRUE),
                            assemblage.3 = sample(0:25, 50, replace = TRUE),
                            assemblage.4 = sample(0:25, 50, replace = TRUE),
                            assemblage.5 = sample(0:25, 50, replace = TRUE))
# Random attribution of names to the sampled species
rownames(assemblages.matrix) <- sample(rownames(spid.occ), 50, replace = FALSE)


# Test 1: A vector - W vector
curW <- rarity.weights$W
names(curW) <- rownames(rarity.weights)
curA <- assemblages.matrix[, 1]

result <- Isr(W = curW, assemblages = curA, abundance = TRUE)
result
tot.res <- result[1]
if(any(is.na(result)) | any(result[grep("Isr", names(result))] < 0))
{
  stop("Error in the test")
}
# Testing NAs
curW[2] <- NA
result <- Isr(W = curW, assemblages = curA, abundance = TRUE)
result
if(any(is.na(result)) | any(result[grep("Isr", names(result))] < 0))
{
  stop("Error in the test")
}


# Test 2: A vector - W matrix
curW <- as.matrix(rarity.weights)
curA <- assemblages.matrix[, 1]

result <- Isr(W = curW, assemblages = curA, abundance = TRUE)
result
if(any(is.na(result)) | any(result[grep("Isr", names(result))] < 0))
{
  stop("Error in the test")
}
tot.res <- c(tot.res,
             result["Isr_W"])

# Testing NAs
curW[1, 2] <- NA
curW[10, 2] <- NA
curW[5, 3] <- NA
result <- Isr(W = curW, assemblages = curA)
result
if(any(is.na(result)) | any(result[grep("Isr", names(result))] < 0))
{
  stop("Error in the test")
}


# Test 3: A matrix - W vector
curW <- rarity.weights$W
names(curW) <- rownames(rarity.weights)
curA <- assemblages.matrix

result <- Isr(W = curW, assemblages = curA, abundance = TRUE)
result
tot.res <- c(tot.res,
             result[1, 1])
if(any(is.na(result)) | any(result[grep("Isr", names(result))] < 0))
{
  stop("Error in the test")
}

# Testing NAs
curW[2] <- NA
result <- Isr(W = curW, assemblages = curA, abundance = TRUE)
result
if(any(is.na(result)) | any(result[grep("Isr", names(result))] < 0))
{
  stop("Error in the test")
}

# Test 4: A matrix - W matrix
curW <- rarity.weights
curA <- assemblages.matrix

result <- Isr(W = curW, assemblages = curA, abundance = TRUE)
result
tot.res <- c(tot.res,
             result[1, 3])
if(any(is.na(result)) | any(result[grep("Isr", names(result))] < 0))
{
  stop("Error in the test")
}

# Testing NAs
curW[1, 2] <- NA
curW[10, 2] <- NA
curW[5, 3] <- NA
result <- Isr(W = curW, assemblages = curA, abundance = TRUE)
result
if(any(is.na(result)) | any(result[grep("Isr", names(result))] < 0))
{
  stop("Error in the test")
}


tot.res



#### TEST WITH NORMALISATION ####

# Generation of an assemblage matrix
assemblages.matrix <- cbind(assemblage.1 = sample(0:25, 50, replace = TRUE),
                            assemblage.2 = sample(0:25, 50, replace = TRUE),
                            assemblage.3 = sample(0:25, 50, replace = TRUE),
                            assemblage.4 = sample(0:25, 50, replace = TRUE),
                            assemblage.5 = sample(0:25, 50, replace = TRUE))
# Random attribution of names to the sampled species
rownames(assemblages.matrix) <- sample(rownames(spid.occ), 50, replace = FALSE)


# Test 1: A vector - W vector
curW <- rarity.weights$W
names(curW) <- rownames(rarity.weights)
curA <- assemblages.matrix[, 1]

result <- Isr(W = curW, assemblages = curA, normalise = TRUE)
result
tot.res <- result[1]
if(any(is.na(result)) | any(result[grep("Isr", names(result))] < 0))
{
  stop("Error in the test")
}
# Testing NAs
curW[2] <- NA
result <- Isr(W = curW, assemblages = curA, normalise = TRUE)
result
if(any(is.na(result)) | any(result[grep("Isr", names(result))] < 0))
{
  stop("Error in the test")
}


# Test 2: A vector - W matrix
curW <- as.matrix(rarity.weights)
curA <- assemblages.matrix[, 1]

result <- Isr(W = curW, assemblages = curA, normalise = TRUE)
result
if(any(is.na(result)) | any(result[grep("Isr", names(result))] < 0))
{
  stop("Error in the test")
}
tot.res <- c(tot.res,
             result["Isr_W"])

# Testing NAs
curW[1, 2] <- NA
curW[10, 2] <- NA
curW[5, 3] <- NA
result <- Isr(W = curW, assemblages = curA, normalise = TRUE)
result
if(any(is.na(result)) | any(result[grep("Isr", names(result))] < 0))
{
  stop("Error in the test")
}


# Test 3: A matrix - W vector
curW <- rarity.weights$W
names(curW) <- rownames(rarity.weights)
curA <- assemblages.matrix

result <- Isr(W = curW, assemblages = curA, normalise = TRUE)
result
tot.res <- c(tot.res,
             result[1, 1])
if(any(is.na(result)) | any(result[grep("Isr", names(result))] < 0))
{
  stop("Error in the test")
}

# Testing NAs
curW[2] <- NA
result <- Isr(W = curW, assemblages = curA, normalise = TRUE)
result
if(any(is.na(result)) | any(result[grep("Isr", names(result))] < 0))
{
  stop("Error in the test")
}

# Test 4: A matrix - W matrix
curW <- rarity.weights
curA <- assemblages.matrix

result <- Isr(W = curW, assemblages = curA, normalise = TRUE)
result
tot.res <- c(tot.res,
             result[1, 3])
if(any(is.na(result)) | any(result[grep("Isr", names(result))] < 0))
{
  stop("Error in the test")
}

# Testing NAs
curW[1, 2] <- NA
curW[10, 2] <- NA
curW[5, 3] <- NA
result <- Isr(W = curW, assemblages = curA, normalise = TRUE)
result
if(any(is.na(result)) | any(result[grep("Isr", names(result))] < 0))
{
  stop("Error in the test")
}


if(any(diff(tot.res)) > 0)
{
  stop("Different values found for the same indices")
}