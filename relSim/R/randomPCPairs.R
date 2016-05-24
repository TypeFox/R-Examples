randomPCPairs = function(Freqs, BlockSize = 1){
    nLoci = length(Freqs$loci)
    Profile = vector(mode = "list", length = BlockSize)

    Parent = randomProfiles(Freqs$freqs, BlockSize)
    Child = randomChildren(Parent, Freqs$freqs, BlockSize)

    for(b in 1:BlockSize){
        i1 = (b - 1) * 2 * nLoci + 1
        i2 =  b * 2 * nLoci

        Profile[[b]] = vector(mode = "list", length = 2)
        names(Profile[[b]]) = c("parent", "child")

        Profile[[b]]$parent = Parent[i1:i2]
        class(Profile[[b]]$parent) = "profile"

        Profile[[b]]$child = Child[i1:i2]
        class(Profile[[b]]$child) = "profile"
    }

    if(BlockSize == 1){
        return(Profile[[1]])
    }else{
        return(Profile)
    }
}
