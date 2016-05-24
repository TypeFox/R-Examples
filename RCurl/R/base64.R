# There is some support for this in caTools but 
# that didn't install on some of my machines 
# and so we add it here since it is already available
# from libcurl and so is a natural facility in RCurl.
# I don't like duplicating functionality in other packages
# and discourage it. 


base64 =
function(txt, encode = !inherits(txt, "base64"), mode = "character")
{
   asRaw = (as.character(mode) == "raw")

   encode
   
   if(typeof(txt) != "raw")
      txt = as.character(txt)

   if(encode) {
      ans = .Call(R_base64_encode, txt, asRaw)
      class(ans) <- "base64"
      ans
   } else
      .Call(R_base64_decode, txt, asRaw)

}


base64Encode = 
function(txt, mode = "character")
  base64(txt, TRUE, mode)

base64Decode = 
function(txt, mode = "character")
  base64(txt, FALSE, mode)

