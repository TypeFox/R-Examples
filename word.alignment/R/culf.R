culf <-
function(x, n = 1, first = TRUE, second = FALSE , lower = FALSE)
{
if (first)  x = paste (tolower (substr (x, 1, n)), substring (x, n+1), sep = '')
       
if (second) x [substring (x, 2, 2) %in% LETTERS] = tolower (x [substring (x, 2, 2) %in% LETTERS])
      
if (lower)  x = tolower (x)

return (x)
}