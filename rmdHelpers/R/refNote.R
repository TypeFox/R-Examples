refNote <-
function(text = "This is a test note", number = "*"){
  out <- paste('<span class="ref"><span class="refnum">[',
               number,
               ']</span><span class="refbody">',
               text,
               '</span></span>',
               sep = "")
  return(out)
}
