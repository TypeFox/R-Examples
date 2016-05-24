library("R.rsp")

rcat("A random integer in [1,100]: <%=sample(1:100, size=1)%>\n")

# Passing arguments
rcat("A random integer in [1,<%=K%>]: <%=sample(1:K, size=1)%>\n", args=list(K=50))

text <- 'The <%=n <- length(letters)%> letters in the English alphabet are:
<% for (i in 1:n) { %>
<%=letters[i]%>/<%=LETTERS[i]%><%=if(i < n) ", "-%>
<% } %>.\n'
rcat(text)


# Informative syntax error messages
text <- '<%={
10 print "Hello world!"
20 goto 10
}%>\n'
tryCatch({ rcat(text) }, error=function(ex) cat(ex$message))
