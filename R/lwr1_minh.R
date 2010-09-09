lwr1_minh <-
function(form,window=0,bandwidth=0,kern="tcub",method="gcv",print=TRUE,alldata=FALSE,data=NULL) {
  attach(data,warn.conflicts=FALSE)

  nw = length(window)
  nb = length(bandwidth)
  if ((nw==1)&(nb==1)) { cat("Must provide a vector of window or bandwidth values","\n") }
  if ((nw>1)&(nb>1))   { cat("Both window and bandwidth vectors specified; will use window vector","\n") }

  minval = 9999999
  minterm = 0
  icross = ifelse(method=="gcv",FALSE,TRUE)

  if (nw>1) {
    for (iw in window) {
      fit1 <- lwr1(form,window=iw,bandwidth=0,kern=kern,alldata=alldata,predx=0,graph.yhat=FALSE,graph.predx=FALSE)
      hval = ifelse(icross==TRUE,fit1$cv,fit1$gcv)
      if (print==TRUE) {print(c(iw,hval))}
      if (hval<minval) {
        minh = iw
        minval = hval
        outfit <- fit1 
      }
     }
   }

  if (nb>1) {
    for (ib in bandwidth) {
      fit1 <- lwr1(form,window=0,bandwidth=ib,kern=kern,alldata=alldata,predx=0,graph.yhat=FALSE,graph.predx=FALSE)
      hval = ifelse(icross==TRUE,fit1$cv,fit1$gcv)
      if (print==TRUE) {print(c(ib,hval))}
      if (hval<minval) {
        minh = ib
        minval = hval
        outfit <- fit1 
      }
    }
   }

   cat("\n","h = ",minh,"Function Value = ",minval,"\n") 
 
   return(minh)
}

