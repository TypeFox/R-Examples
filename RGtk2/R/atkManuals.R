atkComponentAddFocusHandler <- atkObjectConnectPropertyChangeHandler <-
function(object, handler)
{
	.notimplemented("does not have user data for the handler")
}
atkAddFocusTracker <-
function(focus.tracker)
{
    .notimplemented("does not have user data for the tracker")
}
atkFocusTrackerInit <-
function(add.function)
{
	.notimplemented("does not have user data for init function. Besides, it's only for ATK implementations")
}
atkAddGlobalEventListener <-
function(listener, event.type)
{
    .notimplemented("does not have user data for the listener")
}
atkStreamableContentGetStream <-
function (object, mime.type) 
{
	.notimplemented("is not yet written. It will support retrieving an R connection from AtkStreamableContent")
}
