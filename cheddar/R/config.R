# Everything to do with configuration
.UnnamedString <- function()
{
    # A string for unnamed objects
    return ('<unnamed>')
}

.maxQueueOption <- function()
{
    return (getOption('cheddarMaxQueue', 1e7))
}
