library(bitops)

i7 <- 0:127
ri7 <- bitFlip(i7)
stopifnot(identical(bitAnd(i7,ri7),rep(0,length(i7))),
          ri7+i7 == 2^32-1,
          ## flipping the bits twice should be the identity (modulo overflow):
          identical(i7, as.integer(bitFlip(ri7))),
          bitAnd(i7, ri7) == 0,

          bitAnd(15,17) ==  1,
          bitOr (15,17) == 31,
          bitXor(15,17) == 30
          )

for(N in 1:200) {
    j7 <- sample(i7)
    ## Commutative Law:
    stopifnot(identical(bitOr (i7, j7), bitOr (j7, i7)),
              identical(bitAnd(i7, j7), bitAnd(j7, i7)),
              identical(bitXor(i7, j7), bitXor(j7, i7)))
    ## Xor "+" And  == Or :
    stopifnot(identical(bitOr(i7, j7),
                        bitOr(bitXor(i7,j7), bitAnd(i7,j7))))
    ## Logic:  !(A & B)  <->  (!A) | (!B)
    stopifnot(identical(bitFlip(bitAnd(i7, j7)),
                        bitOr(bitFlip(i7), bitFlip(j7))))
    ##         !(A | B)  <->  (!A) & (!B)
    stopifnot(identical(bitFlip(bitOr(i7, j7)),
                        bitAnd(bitFlip(i7), bitFlip(j7))))
    ## Associative Law:
    k7 <- sample(j7)
    stopifnot(identical(bitOr(bitOr(i7, j7), k7),
                        bitOr(i7, bitOr(j7, k7))),
              identical(bitAnd(bitAnd(i7, j7), k7),
                        bitAnd(i7, bitAnd(j7, k7))),
              identical(bitXor(bitXor(i7, j7), k7),
                        bitXor(i7, bitXor(j7, k7))))
}

# 
# verify cksum 
#
b<-sapply(1:92, FUN=function(a,b=" !#$%&()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[]^_`abcdefghijklmnopqrstuvwxyz{|}~") {
	paste(substring(b,a),substring(b,1,a-1),collapse="",sep="")
})

d<-c(2442416856, 1669542993, 313999433, 178729808, 3394733963, 2164389741,
        3871734349, 3789449038, 40636212, 1452746146, 541480198, 
        2979936832, 2923782422, 792265197, 3640409291, 1202696403, 
        4011398543, 2699207183, 2985612474, 1439186030, 1508213684, 
        1865388774, 2380454843, 454855490, 1019166481, 924244674, 
        1406204380, 2429078660, 1046223291, 1230078089, 1548993556, 
        280855472, 421066716, 2967223269, 1100914587, 886676022, 
        1657109189, 843923270, 620178494, 1609552402, 1787171819, 
        4006198310, 1023859819, 1411671880, 513493423, 2495633464, 
        1866449535, 4291277827, 3301230818, 381214501, 2497598429, 
        675736398, 3735311659, 2170409126, 3731386467, 1015853879, 
        4060922207, 1023658490, 2980477601, 350747207, 2650042644, 
        600967562, 4254175774, 1970787970, 4065204194, 1521286262, 
        3589949651, 879070207, 1152896007, 2418807455, 2666637124, 
        2577590065, 4208759298, 3274144307, 1957580223, 3095930811, 
        3625810032, 126832280, 1912362968, 515865842, 3876027886, 
        304043927, 785523686, 3840974701, 2587165204, 1710947718, 
        2356035548, 430213333, 3484582166, 885948210, 1348073033, 
        2652440189)

stopifnot( all.equal(cksum(b),d))

# verify bit shifts:

stopifnot( identical(2^(0:31), bitShiftL(1,0:31)) )
stopifnot( identical(2^(31:0),bitShiftR(2^31,0:31)) )

# test boundary value behavior:  +/- Inf, NA, NaN, 2^32:


a<-round(runif(500)*2^33)
b<-which(a<4294967296)

stopifnot(identical(bitAnd(a,a)[-b],as.numeric(rep(NA,length(a))[-b])))
stopifnot(identical(bitAnd(a,a)[b],a[b]))
stopifnot(identical(bitOr(a,a)[-b],as.numeric(rep(NA,length(a))[-b])))
stopifnot(identical(bitOr(a,a)[b],a[b]))
stopifnot(identical(bitXor(a,0)[-b],as.numeric(rep(NA,length(a))[-b])))
stopifnot(identical(bitXor(a,0)[b],a[b]))
stopifnot(identical(bitXor(0,a)[-b],as.numeric(rep(NA,length(a))[-b])))
stopifnot(identical(bitXor(0,a)[b],a[b]))
stopifnot(identical(bitFlip(bitFlip(a))[-b],as.numeric(rep(NA,length(a))[-b])))
stopifnot(identical(bitFlip(bitFlip(a))[b],bitAnd(a,2^32-1)[b]))
stopifnot(identical(bitShiftR(a,runif(10)*32)[-b],as.numeric(rep(NA,length(a))[-b])))
stopifnot(identical(bitShiftL(a,runif(10)*32)[-b],as.numeric(rep(NA,length(a))[-b])))

a[-b]<-1/0
stopifnot(identical(bitAnd(a,a)[-b],as.numeric(rep(NA,length(a))[-b])))
stopifnot(identical(bitAnd(a,a)[b],a[b]))
stopifnot(identical(bitOr(a,a)[-b],as.numeric(rep(NA,length(a))[-b])))
stopifnot(identical(bitOr(a,a)[b],a[b]))
stopifnot(identical(bitXor(a,0)[-b],as.numeric(rep(NA,length(a))[-b])))
stopifnot(identical(bitXor(a,0)[b],a[b]))
stopifnot(identical(bitXor(0,a)[-b],as.numeric(rep(NA,length(a))[-b])))
stopifnot(identical(bitXor(0,a)[b],a[b]))
stopifnot(identical(bitFlip(bitFlip(a))[-b],as.numeric(rep(NA,length(a))[-b])))
stopifnot(identical(bitFlip(bitFlip(a))[b],bitAnd(a,2^32-1)[b]))
stopifnot(identical(bitShiftR(a,runif(10)*32)[-b],as.numeric(rep(NA,length(a))[- b])))
stopifnot(identical(bitShiftL(a,runif(10)*32)[-b],as.numeric(rep(NA,length(a))[- b])))

a[-b]<--1/0
stopifnot(identical(bitAnd(a,a)[-b],as.numeric(rep(NA,length(a))[-b])))
stopifnot(identical(bitAnd(a,a)[b],a[b]))
stopifnot(identical(bitOr(a,a)[-b],as.numeric(rep(NA,length(a))[-b])))
stopifnot(identical(bitOr(a,a)[b],a[b]))
stopifnot(identical(bitXor(a,0)[-b],as.numeric(rep(NA,length(a))[-b])))
stopifnot(identical(bitXor(a,0)[b],a[b]))
stopifnot(identical(bitXor(0,a)[-b],as.numeric(rep(NA,length(a))[-b])))
stopifnot(identical(bitXor(0,a)[b],a[b]))
stopifnot(identical(bitFlip(bitFlip(a))[-b],as.numeric(rep(NA,length(a))[-b])))
stopifnot(identical(bitFlip(bitFlip(a))[b],bitAnd(a,2^32-1)[b]))
stopifnot(identical(bitShiftR(a,runif(10)*32)[-b],as.numeric(rep(NA,length(a))[- b])))
stopifnot(identical(bitShiftL(a,runif(10)*32)[-b],as.numeric(rep(NA,length(a))[- b])))

a[-b]<-suppressWarnings(sqrt(-1))
stopifnot(identical(bitAnd(a,a)[-b],as.numeric(rep(NA,length(a))[-b])))
stopifnot(identical(bitAnd(a,a)[b],a[b]))
stopifnot(identical(bitOr(a,a)[-b],as.numeric(rep(NA,length(a))[-b])))
stopifnot(identical(bitOr(a,a)[b],a[b]))
stopifnot(identical(bitXor(a,0)[-b],as.numeric(rep(NA,length(a))[-b])))
stopifnot(identical(bitXor(a,0)[b],a[b]))
stopifnot(identical(bitXor(0,a)[-b],as.numeric(rep(NA,length(a))[-b])))
stopifnot(identical(bitXor(0,a)[b],a[b]))
stopifnot(identical(bitFlip(bitFlip(a))[-b],as.numeric(rep(NA,length(a))[-b])))
stopifnot(identical(bitFlip(bitFlip(a))[b],bitAnd(a,2^32-1)[b]))
stopifnot(identical(bitShiftR(a,runif(10)*32)[-b],as.numeric(rep(NA,length(a))[- b])))
stopifnot(identical(bitShiftL(a,runif(10)*32)[-b],as.numeric(rep(NA,length(a))[- b])))

a[-b]<-NA
stopifnot(identical(bitAnd(a,a)[-b],as.numeric(rep(NA,length(a))[-b])))
stopifnot(identical(bitAnd(a,a)[b],a[b]))
stopifnot(identical(bitOr(a,a)[-b],as.numeric(rep(NA,length(a))[-b])))
stopifnot(identical(bitOr(a,a)[b],a[b]))
stopifnot(identical(bitXor(a,0)[-b],as.numeric(rep(NA,length(a))[-b])))
stopifnot(identical(bitXor(a,0)[b],a[b]))
stopifnot(identical(bitXor(0,a)[-b],as.numeric(rep(NA,length(a))[-b])))
stopifnot(identical(bitXor(0,a)[b],a[b]))
stopifnot(identical(bitFlip(bitFlip(a))[-b],as.numeric(rep(NA,length(a))[-b])))
stopifnot(identical(bitFlip(bitFlip(a))[b],bitAnd(a,2^32-1)[b]))
stopifnot(identical(bitShiftR(a,runif(10)*32)[-b],as.numeric(rep(NA,length(a))[-b])))
stopifnot(identical(bitShiftL(a,runif(10)*32)[-b],as.numeric(rep(NA,length(a))[-b])))

