## Helper function for 'setData'

bindCombo <-
function(cbd1, cbd2, cbd3, cbd4, cbd5, cbd6, cbd7, cbd8,
         dst1, dst2, dst3, dst4, dst5, dst6, dst7, dst8,
         cbs1, cbs2, cbs3, cbs4, cbs5, cbs6, cbs7, cbs8,
         str1, str2, str3, str4, str5, str6, str7, str8,
         tbl1, tbl2, tbl3, tbl4, tbl5, tbl6, tbl7, tbl8,
         var1, var2, var3, var4, var5, var6, var7, var8){
  tkbind(cbd1, "<FocusIn>", function() setDist(dst1, str1, tbl1, var1))
  tkbind(cbd2, "<FocusIn>", function() setDist(dst2, str2, tbl2, var2))
  tkbind(cbd3, "<FocusIn>", function() setDist(dst3, str3, tbl3, var3))
  tkbind(cbd4, "<FocusIn>", function() setDist(dst4, str4, tbl4, var4))
  tkbind(cbd5, "<FocusIn>", function() setDist(dst5, str5, tbl5, var5))
  tkbind(cbd6, "<FocusIn>", function() setDist(dst6, str6, tbl6, var6))
  tkbind(cbd7, "<FocusIn>", function() setDist(dst7, str7, tbl7, var7))
  tkbind(cbd8, "<FocusIn>", function() setDist(dst8, str8, tbl8, var8))

  tkbind(cbs1, "<FocusIn>", function() setDist(dst1, str1, tbl1, var1))
  tkbind(cbs2, "<FocusIn>", function() setDist(dst2, str2, tbl2, var2))
  tkbind(cbs3, "<FocusIn>", function() setDist(dst3, str3, tbl3, var3))
  tkbind(cbs4, "<FocusIn>", function() setDist(dst4, str4, tbl4, var4))
  tkbind(cbs5, "<FocusIn>", function() setDist(dst5, str5, tbl5, var5))
  tkbind(cbs6, "<FocusIn>", function() setDist(dst6, str6, tbl6, var6))
  tkbind(cbs7, "<FocusIn>", function() setDist(dst7, str7, tbl7, var7))
  tkbind(cbs8, "<FocusIn>", function() setDist(dst8, str8, tbl8, var8))

  setDist(dst1, str1, tbl1, var1);  setDist(dst5, str5, tbl5, var5)
  setDist(dst2, str2, tbl2, var2);  setDist(dst6, str6, tbl6, var6)
  setDist(dst3, str3, tbl3, var3);  setDist(dst7, str7, tbl7, var7)
  setDist(dst4, str4, tbl4, var4);  setDist(dst8, str8, tbl8, var8)
}

