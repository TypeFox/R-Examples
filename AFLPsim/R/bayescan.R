bayescan <-
function (mat, filename, nbp = 20, pilot = 5000, burn = 50000, exec=NULL) 
{
    sim2bayescan(mat, filename)
    prefix <- c("-nbp", "-pilot", "-burn")
    suffix <- c(nbp, pilot, burn)
    opts <- paste(prefix, suffix, sep = " ", collapse = " ")
       os <- Sys.info()[1]
    if (is.null(exec)) {
        if (os == "Linux") 
            exec= "bayescan_2.1"
        if (os == "Darwin") 
            exec= "bayescan_2.1"
        if (os == "Windows") 
            shortPathName("C:/Program Files/BayeScan2.1/binaries/BayeScan2.1_win32bits_cmd_line.exe")
    }
    
    system(paste(exec, opts, filename))
    out <- paste(filename, "_fst.txt", sep = "")
    if (file.exists(out)==FALSE){
      print("Bayescan_2.1 not found")
      return(NULL)
    } 
    else {
    resultfst <- read.table(out, header = T, row.names = 1)
    result <- resultfst
    result}
}
