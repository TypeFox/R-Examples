ll<-function(data,data.type,tmother,tfather,model.type){
  if(data.type=='sample'){
    cnt=strata.cnt(data)
    if ((dim(data)[2]==4)
        &(model.type=='grri'|model.type=='domi'|model.type=='gdi')){
      cat('Data is without parent-of-origin information,
          cannot fit model with imprinting effect.',fill=T)
    } else if ((dim(data)[2]==5)
              &(model.type=='grr'|model.type=='dom'|model.type=='gd')){
      cnt=c(cnt[1:8],cnt[9]+cnt[10],cnt[11:16])
    }
  }else if (data.type=='counts'){
    cnt=data
    if ((length(data)==15)
        &(model.type=='grri'|model.type=='domi'|model.type=='gdi')){
      cat('Data is without parent-of-origin information,
          cannot fit model with imprinting effect.',fill=T)
    }else if ((length(data)==16)
          &(model.type=='grr'|model.type=='dom'|model.type=='gd')){
      cnt=c(cnt[1:8],cnt[9]+cnt[10],cnt[11:16])
    }
  }


  tm=tmother
  tf=tfather

  if(model.type=='grr'){
    r=tm
    stat=ll.grr(cnt,r)
  }else if (model.type=='dom'){
    r=tm
    stat=ll.dom(cnt,r)
  }else if (model.type=='gd'){
    r=tm
    stat=ll.gd(cnt,r)
  }else if(model.type=='grri'){
    stat=ll.grri(cnt,tm,tf)
  }else if (model.type=='domi'){
    stat=ll.domi(cnt,tm,tf)
  }else if (model.type=='gdi'){
    stat=ll.gdi(cnt,tm,tf)
  } else {
    cat('Model not supported.',fill=T)
  }

  out=out.stats(stat,model.type)
}
