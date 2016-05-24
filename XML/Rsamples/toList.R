library(XML)
tt = '<x>
         <a>text</a>
         <b foo="1"/>
         <c bar="me">
           <d>a phrase</d>
           <d>a second phrase</d>
          <d>a conclusion</d>
         </c>
       </x>'

xmlToList(tt)

tt = '<x>
         <c bar="me">
           <d>a phrase</d>
           <d>a second phrase</d>
          <d>a conclusion</d>
         </c>
       </x>'

tt = '<x>
         <c>
           <d>a phrase</d>
           <d>a second phrase</d>
          <d>a conclusion</d>
         </c>
       </x>'

