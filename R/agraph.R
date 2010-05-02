agraph <-
function(av,price,cftitle="95% CI",title="Assessment Ratios",width=.01,
    legloc="none",freq=FALSE,normdens=FALSE,kdens=FALSE,cfint=c(FALSE,FALSE),statute=FALSE,file="none") {
  if (file!="none") png(file)

  ratio <- av/price

  statute <- ifelse(statute==FALSE,0,statute)
  nclass = ceiling(max(ratio)/width)
  h <- hist(ratio,nclass=nclass,plot=FALSE)
  ylabel = ifelse(freq==TRUE,"Frequency","Percent")
  plot(h,freq=TRUE,axes=FALSE,xlab="Assessment Ratio",ylab=ylabel,main=title,col="khaki",border="brown",xlim=c(min(ratio),max(max(ratio),statute)))
  axis(1)
  if (statute!=FALSE) {abline(v=statute,lty="dashed",lwd=2)}
  if (!identical(cfint,c(FALSE,FALSE))) {
    abline(v=cfint[1],lty="solid",col="blue",lwd=2)
    abline(v=cfint[2],lty="solid",col="blue",lwd=2)
  }

  if (freq==TRUE) axis(2)
  if (freq==FALSE) {
    ynorm <- 100*h$density/sum(h$density)
    ytick <- seq(0,max(h$count),length=6)
    ylabel <- round(seq(0,max(ynorm),length=6),digits=1)
    axis(2,at=ytick,labels=ylabel)
  }
  
  if (normdens==TRUE) {
    k <- dnorm((h$breaks-mean(ratio))/sd(ratio))
    sumk <- sum(k)
    normy <- length(ratio)*k/sumk
    lines(h$breaks,normy,col="red",lwd=2)
  }
  if (kdens==TRUE) {
    n = length(h$breaks)
    k <- density(ratio,from=min(h$breaks),to=max(h$breaks), n=length(h$breaks))
    sumk = sum(k$y)
    ky <- length(ratio)*k$y/sumk
    lines(h$breaks,ky,col="green",lwd=2)
  }

  if (legloc!="none") {
    vect1 <- c("Normal Density","Kernel Density",cftitle,"Statutory Rate")
    vect2 <- c("red","green","blue","black")
    vect3 <- c("solid","solid","solid","dashed")
    x <- c(1*normdens,2*kdens,3*(cfint[1]!=FALSE),4*(statute!=FALSE))
    x <- x[x!=0]
    vect1 <- vect1[x]
    vect2 <- vect2[x]
    vect3 <- vect3[x]

    legend(legloc,vect1,col=vect2,lty=vect3,lwd=1,bty="n")
   }

  if (file!="none") dev.off()

}

