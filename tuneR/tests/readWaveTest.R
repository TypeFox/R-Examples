library("tuneR")
d <- list.files("Testfiles", pattern="\\.wav$")
garbage <- lapply(d, 
    function(x) {
        wav <- readWave(file.path("Testfiles", x))
        cat(x, "\n")
        print(wav)
        if(is(wav, "Wave")){
            print(wav@left[1:10])
            if(wav@stereo)
                print(wav@right[1:10])
        } else {
            print(wav@.Data[1:10,])
        }
    }
)
