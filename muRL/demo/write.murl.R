data(murljobs)

## Create mailmerge.tex required for LaTeX import
write.murl(murljobs)

## Specify a file containing the letters' body text
## write.murl(murljobs, letter.file = "mybodytext.txt")

## Specify a string containing the letters' body text
write.murl(murljobs, letter.text = "This is the whole body of my letters.")

## Specify salutation, valediction options (overwrites previous mailmerge.tex)
write.murl(murljobs, file.name = "mailmerge.tex", salutation = "Greetings", 
           sal.punct = ",", valediction = "Truly Yours,", include.opening = FALSE)
           
## Specify opening line also (overwrites previous mailmerge.tex)
write.murl(murljobs, file.name = "mailmerge.tex", salutation = "Greetings", 
           sal.punct = ",", valediction = "Truly Yours,", 
           opening = "I am applying for the job as", include.opening = TRUE)