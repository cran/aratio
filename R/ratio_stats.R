ratio_stats <-
function(av,price,nboot=0,print=TRUE,data=NULL) {
  attach(data,warn.conflicts=FALSE)

  sampobs <- !is.na(av)&!is.na(price)
  av <- av[sampobs]
  price <- price[sampobs]
  ratio <- av/price
  rdata <- data.frame(av,price,ratio)
  n = length(ratio)
  i <- seq(1:n)

  basef <- function(rdata,i) {
    bdata <- rdata[i,]
    rmean = mean(bdata$ratio)
    rmedian = median(bdata$ratio)
    abar = mean(bdata$av)
    pbar = mean(bdata$price)
    wmean = abar/pbar
    prd = rmean/wmean
    x <- abs(bdata$ratio-rmedian)
    cod = 100*mean(x)/rmedian
    x <- (bdata$price/pbar)*x
    wcod = 100*mean(x)/rmedian
    
    avect <- c(rmean,rmedian,wmean,cod,wcod,prd)
    return(avect)
  }
  
  if (nboot==0) {
    bvect <- basef(rdata,i)
    names(bvect) <- c("mean","median","wgt mean","cod","wgt cod","prd")
  }
  if (nboot>0) {
    library(boot)
    rboot <- boot(rdata,basef,R=nboot)
    bvect <- cbind(rboot$t0,sqrt(diag(var(rboot$t))), 0,0,0,0)
    bvect[,3] <- bvect[,1] - 1.96*bvect[,2]
    bvect[,4] <- bvect[,1] + 1.96*bvect[,2]
    for (i in seq(1:6)) {
      bvect[i,5] = quantile(rboot$t[,i],.025)
      bvect[i,6] = quantile(rboot$t[,i],.975)
    }
    rownames(bvect) <- c("mean","median","wgt mean","cod","wgt cod","prd")
    colnames(bvect) <- c("observed","Std Err","Normal-lo","Normal-hi","Percentile-lo","Percentile-hi")
  }

  if(print==TRUE) print(bvect)
  return(bvect)
}

