filename.clean <-
function# Clean up filenames, replacing spaces and bad characters with "_" and "-"
### internal function for sisus
(filename
### internal variable
)
{
  ##details<<
  ## interal function for sisus.run()

    #bad.char.list = "!@#$%^&*+=/`'\"{}[]:;<>,.?|";  # represents list of bad characters

  filename = gsub(" ", "_", filename);                                        # replace space by underscore "_"
  #filename = gsub("[\\\|\[\{\^\$\*\+\?!@#%&=/`'\"}:;<>,]", "-", filename);  # replace bad char by hyphen "-", all but "]" and \.
    # in R2.6.1, doesn't like the 8 escapes removed below  1/2/2008 3:30PM
  filename = gsub("[\\!@#%&=/`'\"}:;<>,]", "-", filename);  # replace bad char by hyphen "-", all but "]" and \.

  return(filename);

  ### internal variable
}
