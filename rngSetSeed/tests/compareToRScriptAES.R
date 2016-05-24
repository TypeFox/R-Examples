library(rngSetSeed)

# Comparison of the output of C-level functions in rngSetSeed
# with AES implemented in a (slow) R script below prepared using
# the description from
# http://en.wikipedia.org/wiki/Advanced_Encryption_Standard
# Petr Savicky 2012

# Prepare tables for AES

m <- 8
M <- 2^m

vect <- unname(as.matrix(rev(expand.grid(rep(list(0:1), times=m)))))
weight <- 2^((m - 1):0)
xor <- matrix(nrow=M, ncol=M)
for (i in seq.int(length.out=M)) {
	for (j in seq.int(length.out=M)) {
		xor[i, j] <- sum(((vect[i, ] + vect[j, ]) %% 2)*weight)
	}
}

h <- c("0","1","2","3","4","5","6","7","8","9","a","b","c","d","e","f")
hh <- rev(expand.grid(h, h))
hexa <- 0:255
names(hexa) <- paste(hh[, 1], hh[, 2], sep="")

mult.x <- unname(hexa[c(
  "00","02","04","06","08","0a","0c","0e","10","12","14","16","18","1a","1c","1e",
  "20","22","24","26","28","2a","2c","2e","30","32","34","36","38","3a","3c","3e",
  "40","42","44","46","48","4a","4c","4e","50","52","54","56","58","5a","5c","5e",
  "60","62","64","66","68","6a","6c","6e","70","72","74","76","78","7a","7c","7e",
  "80","82","84","86","88","8a","8c","8e","90","92","94","96","98","9a","9c","9e",
  "a0","a2","a4","a6","a8","aa","ac","ae","b0","b2","b4","b6","b8","ba","bc","be",
  "c0","c2","c4","c6","c8","ca","cc","ce","d0","d2","d4","d6","d8","da","dc","de",
  "e0","e2","e4","e6","e8","ea","ec","ee","f0","f2","f4","f6","f8","fa","fc","fe",
  "1b","19","1f","1d","13","11","17","15","0b","09","0f","0d","03","01","07","05",
  "3b","39","3f","3d","33","31","37","35","2b","29","2f","2d","23","21","27","25",
  "5b","59","5f","5d","53","51","57","55","4b","49","4f","4d","43","41","47","45",
  "7b","79","7f","7d","73","71","77","75","6b","69","6f","6d","63","61","67","65",
  "9b","99","9f","9d","93","91","97","95","8b","89","8f","8d","83","81","87","85",
  "bb","b9","bf","bd","b3","b1","b7","b5","ab","a9","af","ad","a3","a1","a7","a5",
  "db","d9","df","dd","d3","d1","d7","d5","cb","c9","cf","cd","c3","c1","c7","c5",
  "fb","f9","ff","fd","f3","f1","f7","f5","eb","e9","ef","ed","e3","e1","e7","e5")])

s.box <- unname(hexa[c(
  "63","7c","77","7b","f2","6b","6f","c5","30","01","67","2b","fe","d7","ab","76",
  "ca","82","c9","7d","fa","59","47","f0","ad","d4","a2","af","9c","a4","72","c0",
  "b7","fd","93","26","36","3f","f7","cc","34","a5","e5","f1","71","d8","31","15",
  "04","c7","23","c3","18","96","05","9a","07","12","80","e2","eb","27","b2","75",
  "09","83","2c","1a","1b","6e","5a","a0","52","3b","d6","b3","29","e3","2f","84",
  "53","d1","00","ed","20","fc","b1","5b","6a","cb","be","39","4a","4c","58","cf",
  "d0","ef","aa","fb","43","4d","33","85","45","f9","02","7f","50","3c","9f","a8",
  "51","a3","40","8f","92","9d","38","f5","bc","b6","da","21","10","ff","f3","d2",
  "cd","0c","13","ec","5f","97","44","17","c4","a7","7e","3d","64","5d","19","73",
  "60","81","4f","dc","22","2a","90","88","46","ee","b8","14","de","5e","0b","db",
  "e0","32","3a","0a","49","06","24","5c","c2","d3","ac","62","91","95","e4","79",
  "e7","c8","37","6d","8d","d5","4e","a9","6c","56","f4","ea","65","7a","ae","08",
  "ba","78","25","2e","1c","a6","b4","c6","e8","dd","74","1f","4b","bd","8b","8a",
  "70","3e","b5","66","48","03","f6","0e","61","35","57","b9","86","c1","1d","9e",
  "e1","f8","98","11","69","d9","8e","94","9b","1e","87","e9","ce","55","28","df",
  "8c","a1","89","0d","bf","e6","42","68","41","99","2d","0f","b0","54","bb","16")])

rcon <- unname(hexa[c(
  "01","02","04","08","10","20","40","80","1b","36","6c","d8","ab","4d","9a","2f",
  "5e","bc","63","c6","97","35","6a","d4","b3","7d","fa","ef","c5","91","39","72",
  "e4","d3","bd","61","c2","9f","25","4a","94","33","66","cc","83","1d","3a","74",
  "e8","cb","8d","01","02","04","08","10","20","40","80","1b","36","6c","d8","ab",
  "4d","9a","2f","5e","bc","63","c6","97","35","6a","d4","b3","7d","fa","ef","c5",
  "91","39","72","e4","d3","bd","61","c2","9f","25","4a","94","33","66","cc","83",
  "1d","3a","74","e8","cb","8d","01","02","04","08","10","20","40","80","1b","36",
  "6c","d8","ab","4d","9a","2f","5e","bc","63","c6","97","35","6a","d4","b3","7d",
  "fa","ef","c5","91","39","72","e4","d3","bd","61","c2","9f","25","4a","94","33",
  "66","cc","83","1d","3a","74","e8","cb","8d","01","02","04","08","10","20","40",
  "80","1b","36","6c","d8","ab","4d","9a","2f","5e","bc","63","c6","97","35","6a",
  "d4","b3","7d","fa","ef","c5","91","39","72","e4","d3","bd","61","c2","9f","25",
  "4a","94","33","66","cc","83","1d","3a","74","e8","cb","8d","01","02","04","08",
  "10","20","40","80","1b","36","6c","d8","ab","4d","9a","2f","5e","bc","63","c6",
  "97","35","6a","d4","b3","7d","fa","ef","c5","91","39","72","e4","d3","bd","61",
  "c2","9f","25","4a","94","33","66","cc","83","1d","3a","74","e8","cb")])

# Functions implementing AES

sched.core <- function(t,i)
{
	t <- t[c(2,3,4,1)]
	t <- s.box[t+1]
	t[1] <- xor[t[1]+1,rcon[i]+1]
	t
}

mix.col <- function(a)
{
	c1 <- mult.x[a+1]
	c2 <- a[c(4,1:3)]
	c3 <- a[c(3,4,1,2)]
	c4 <- a[c(2:4,1)]
	c5 <- mult.x[a[c(2:4,1)]+1]
	c45 <- xor[cbind(c4,c5) + 1]
	c23 <- xor[cbind(c2,c3) + 1]
	c123 <- xor[cbind(c1,c23) + 1]
	xor[cbind(c123,c45) + 1]
}

key.schedule <- function(key.len,key)
{
	if (key.len == 128) {
		n <- 16
		b <- 176
	} else if(key.len == 192) {
		n <- 24
		b <- 208
	} else if(key.len == 256) {
		n <- 32
		b <- 240
	} else {
		stop("chyba key.len")
	}
	stopifnot(length(key)==n)
	curr <- n
	i <- 1
	while (curr < b) {
		t <- key[(curr-3):curr]
		t <- sched.core(t,i)
		i <- i+1
		s <- key[(curr+1-n):(curr+4-n)]
		key[(curr+1):(curr+4)] <- xor[cbind(t,s) + 1]
		curr <- curr+4
		for (k in 1:3) {
			t <- key[(curr-3):curr]
			s <- key[(curr+1-n):(curr+4-n)]
			key[(curr+1):(curr+4)] <- xor[cbind(t,s) + 1]
			curr <- curr+4
		}
		if (key.len == 256) {
			t <- key[(curr-3):curr]
			t <- s.box[t+1]
			s <- key[(curr+1-n):(curr+4-n)]
			key[(curr+1):(curr+4)] <- xor[cbind(t, s) + 1]
			curr <- curr+4
		}
		if (key.len == 128) {
			kk <- c()
		} else if (key.len == 192) {
			kk <- 1:2
		} else if (key.len == 256) {
			kk <- 1:3
		} else stop()
		for (k in kk) {
			t <- key[(curr-3):curr]
			s <- key[(curr+1-n):(curr+4-n)]
			key[(curr+1):(curr+4)] <- xor[cbind(t, s) + 1]
			curr <- curr+4
		}
	}
	key[1:b]
}

encrypt <- function(exp.key,input)
{
	stopifnot(length(input)==16)
	key <- matrix(exp.key,ncol=16,byrow=TRUE)
	rounds <- nrow(key)
	# Initial step
	input <- xor[cbind(input, key[1, ]) + 1]
	# Cycle
	for (k in 2:(rounds-1)) {
		input <- s.box[input+1]
		input <- input[c(1,6,11,16,5,10,15,4,9,14,3,8,13,2,7,12)]
		for (j in 1:4) {
			input[4*(j-1)+(1:4)] <- mix.col(input[4*(j-1)+(1:4)])
		}
		input <- xor[cbind(input, key[k, ]) + 1]
	}
	# Final round
	input <- s.box[input+1]
	input <- input[c(1,6,11,16,5,10,15,4,9,14,3,8,13,2,7,12)]
	xor[cbind(input, key[rounds, ]) + 1]
}

# R script implementation of the algorithm used by C-level functions
# in generateInitialization(vseed) and the comparison between the results of
# C and R implementation.

N <- 624
NN <- ceiling(N/4)

for (n in c(1, 2, 7, 8, 9, 47, 48)) {
    vseed <- seq.int(from=floor(0.35*n), length.out=n)
    cat("length(vseed) =", length(vseed), "\n")
    v <- vseed
    v <- c(v, length(v))
    s <- numeric(8*ceiling(length(v)/8))
    s[seq.int(along.with=v)] <- v
    key <- rep(NA, times=32)
    plainText <- rep(0, times=16)
    outBytes <- rep(0, times=16*NN)
    j <- 0
    for (i in seq.int(along.with=s)) { # length(s) is a multiple of 8
        j <- j + 1
        key[j] <- (s[i] %/% 2^24) %% 2^8
        j <- j + 1
        key[j] <- (s[i] %/% 2^16) %% 2^8
        j <- j + 1
        key[j] <- (s[i] %/% 2^8) %% 2^8
        j <- j + 1
        key[j] <- s[i] %% 2^8
        if (j == 32) {
            expandedKey <- key.schedule(256, key)
            u <- (i - 1) %/% 8
            plainText[1] <- (u %/% 2^24) %% 2^8
            plainText[2] <- (u %/% 2^16) %% 2^8
            plainText[3] <- (u %/% 2^8) %% 2^8
            plainText[4] <-  u %% 2^8
            k <- 0
            while (k < length(outBytes)) { # length(outBytes) is a multiple of 16
                u <- k %/% 16
                plainText[5] <- (u %/% 2^24) %% 2^8
                plainText[6] <- (u %/% 2^16) %% 2^8
                plainText[7] <- (u %/% 2^8) %% 2^8
                plainText[8] <-  u %% 2^8
                cipher <- encrypt(expandedKey, plainText)
                outBytes[k + 1:16] <- xor[cbind(outBytes[k + 1:16], cipher) + 1]
                k <- k + 16
            }
            j <- 0
        }
    }
    stopifnot(j == 0)
    out <- c(rbind(2^c(24, 16, 8, 0)) %*% matrix(outBytes, nrow=4, ncol=4*NN))
    out <- out[seq.int(length.out=N)]
    out[out >= 2^31] <- out[out >= 2^31] - 2^32
    OK <- identical(as.integer(out), generateInitialization(vseed, N))
    print(OK)
    stopifnot(OK)
    print(out[1:5])
}

