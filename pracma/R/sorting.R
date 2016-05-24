##
##  s o r t i n g . R  Sorting Routines
##


.comp <- function(u, v) {
	# strictly increasing order
	if (u < v) TRUE
	else FALSE
}

bubbleSort <- function(a) {
    if (!is.numeric(a))
        stop("Argument 'a' must be a non-empty numeric vector.")

    n <- length(a)
    if (n <= 1) return(a)
    for (i in 1:n) {
        for (j in 2:n) {
            b <- a[j]
            if (.comp(a[j], a[j-1]) ) {
                a[j] <- a[j-1]
                a[j-1] <- b
            }
        }
    }
    return(a)
}

insertionSort <- function(a) {
    if (!is.numeric(a))
        stop("Argument 'a' must be a non-empty numeric vector.")

	n <- length(a)
	if (n <= 1) return(a)
	for (i in 2:n) {
		t <- a[i]
		j <- i
		while (j >=2 && a[j-1] > t) {
			a[j] <- a[j-1]
			j <- j-1
		}
		a[j] <- t
	}
	return(a)
}

selectionSort <- function(a) {
    if (!is.numeric(a))
        stop("Argument 'a' must be a non-empty numeric vector.")

    n <- length(a)
    if (n <= 1) return(a)
    for (i in 1:(n-1)) {
        min <- i
        for (j in (i+1):n) {
            if (a[j] < a[min]) min <- j
        }
        t <- a[min]; a[min] <- a[i]; a[i] <- t
    }
    return(a)
}

shellSort <- function(a, f = 2.3) {
    if (!is.numeric(a))
        stop("Argument 'a' must be a non-empty numeric vector.")
    if (!is.numeric(f) || length(f) != 1 || f <= 1)
        stop("The retracting factor 'f' must be a numeric scalar > 1.0")

	n <- length(a)
	if (n <= 1) return(a)
	h <- n %/% 2
	while ( h >= 1) {
		for (i in (h+1):n) {
			t <- a[i]
			j <- i
			while (a[j-h] > t) {
				a[j] <- a[j-h]
				j <- j - h
				if (j <= h) break
			}
			a[j] <- t
		}
		h <- round(h/f)
	}
	return(a)
}

heapSort <- function(a) {
    if (!is.numeric(a))
        stop("Argument 'a' must be a non-empty numeric vector.")
    warning("Function 'heapSort' not yet implemented.")

    n <- length(a)
    if (n <= 1) return(a)
	return(a)
}

mergeSort <- function(a, m = 10) {
    stopifnot(is.numeric(m), length(m) == 1, floor(m) == ceiling(m))
    if (!is.numeric(a))
        stop("Argument 'a' must be a non-empty numeric vector.")

	n <- length(a)
	#if (n <= 1) return(a)
	if (n <= m) return(insertionSort(a))
	m <- n %/% 2
	left <- mergeSort(a[1:m])
	right <- mergeSort(a[(m+1):n])
	return(mergeOrdered(left, right))
}

mergeOrdered <- function(a, b){
	na <- length(a); nb <- length(b)
	ab <- numeric(na + nb)
	i <- j <- 1
	repeat {
		if (a[i] <= b[j]) {
			ab[i+j-1] <- a[i]
			i <- i + 1
			if (i > na) {
				ab[(i+j-1):(na+nb)] <- b[j:nb]
				break
			}
		} else {
			ab[i+j-1] <- b[j]
			j <- j + 1
			if (j > nb) {
				ab[(i+j-1):(na+nb)] <- a[i:na]
				break
			}
		}
	}
	return(ab)
}

quickSort <- function(a, m = 3) {  # m = 3..30
	n <- length(a)
	if (n <= m) {
		return( insertionSort(a) )
	} else {
		v <- (a[1]+a[2])/2
		return ( c(quickSort(a[a <= v], m=m), quickSort(a[a > v], m=m)))
	}
}

quickSortx <- function(a, m = 25) {  # m=20..40, m=25 is favoured
	n <- length(a)
	if (n <= m) return(insertionSort(a))
	i <- 0; j <- n
	v <- a[n]
	repeat {
		while (i < n && a[i+1]  < v) i <- i+1
		while (j > 1 && a[j-1] >= v) j <- j-1
		if (i >= j-1) break
		t <- a[i+1]; a[i+1] <- a[j-1]; a[j-1] <- t
	}
	if (i == 0) return( c(a[n], quickSortx(a[1:(n-1)])) )
	if (j == n) return( c(quickSortx(a[1:(n-1)]), a[n]) )
	return( c(quickSortx(a[1:i]), a[n], quickSortx(a[(i+1):(n-1)])) )
}

is.sorted <- function(a) !is.unsorted(a)

testSort <- function(n = 1000) {
	if (n >= 1e5) warning("n quite large: This will take some time!")
	x <- runif(n)
	
	cat("Test Begin...\n")
	cat("Do not test bubble sort (too slow).\n\n")

	# elapsed <- system.time(y <- bubbleSort(x))[3]
	# cat("Bubble sort: Elapsed time = ", elapsed, "\n")
	# if (is.sorted(y)) cat("Bubble sort successful.\n\n")
	# else cat("Bubble sort test FAILED.\n\n")
	# flush(stdout())

	elapsed <- system.time(y <- insertionSort(x))[3]
	cat("Insertion sort: Elapsed time = ", elapsed, "secs\n")
	if (is.sorted(y)) cat("Insertion sort successful.\n\n")
	else cat("Insertion sort test FAILED.\n\n")
	flush(stdout())

	elapsed <- system.time(y <- selectionSort(x))[3]
	cat("Selection sort: Elapsed time = ", elapsed, "secs\n")
	if (is.sorted(y)) cat("Selection sort successful.\n\n")
	else cat("Selection sort test FAILED.\n\n")
	flush(stdout())

	elapsed <- system.time(y <- shellSort(x))[3]
	cat("Shellsort: Elapsed time = ", elapsed, "secs\n")
	if (is.sorted(y))  { cat("Shell sort successful.\n\n")
	} else { cat("Shell sort test FAILED.\n\n") }
	flush(stdout())

	elapsed <- system.time(y <- mergeSort(x))[3]
	cat("Mergesort: Elapsed time = ", elapsed, "secs\n")
	if (is.sorted(y))  { cat("Merge sort successful.\n\n")
	} else { cat("Merge sort test FAILED.\n\n") }
	flush(stdout())

	elapsed <- system.time(y <- heapSort(x))[3]
	cat("Heapsort: Elapsed time = ", elapsed, "secs\n")
	if (is.sorted(y))  { cat("heap sort successful.\n\n")
	} else { cat("Heap sort test FAILED.\n\n") }
	flush(stdout())

	elapsed <- system.time(y <- quickSort(x))[3]
	cat("Quicksort: Elapsed time = ", elapsed, "secs\n")
	if (is.sorted(y))  { cat("Quicksort successful.\n\n")
	} else { cat("Quicksort test FAILED.\n\n") }

	cat("Test End\n")
}

#-- HwB (C) 2010 -------------------------------------------------------
