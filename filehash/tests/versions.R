## Test databases

suppressMessages(library(filehash))

testdblist <- dir(pattern = glob2rx("testdb-v*"))

for(testname in testdblist) {
        msg <- sprintf("DATABASE: %s\n", testname)
        cat(paste(rep("=", nchar(msg)), collapse = ""), "\n")
        cat(msg)
        cat(paste(rep("=", nchar(msg)), collapse = ""), "\n")
        db <- dbInit(testname, "DB1")
        keys <- dbList(db)
        print(keys)

        for(k in keys) {
                cat("key:", k, "\n")
                val <- dbFetch(db, k)
                print(val)
                cat("\n")
        }
}
