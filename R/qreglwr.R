qreglwr <-
function(form,window=.50,bandwidth=0,kern="tcub",taumat=c(.1,.25,.5,.75,.9),
    alldata=FALSE,predx=0,graph.yhat=FALSE,graph.predx=FALSE,data=NULL) {
  attach(data,warn.conflicts=FALSE)

  mat <- model.frame(form)
  y <- mat[,1]
  x <- mat[,2]
  sdx = sd(x)
  n = length(y)

  if (kern=="rect")  { wgt <- function(psi) {1 } }
  if (kern=="tria")  { wgt <- function(psi) {1 - abs(psi) } }
  if (kern=="epan")  { wgt <- function(psi) { 1-psi^2 } }
  if (kern=="bisq")  { wgt <- function(psi) { (1-psi^2)^2 } }
  if (kern=="tcub")  { wgt <- function(psi) { (1 - abs(psi)^3)^3 } }
  if (kern=="trwt")  { wgt <- function(psi) { (1 - psi^2)^3 } }
  if (kern=="gauss") { wgt <- function(psi) { exp(-((2.5*psi)^2)/2) } }

  if (bandwidth>0) {window = 0}

  if (alldata==FALSE) {
    if (window>0)    {fit <- locfit(~lp(x,nn=window,deg=1)) }
    if (bandwidth>0) {fit <- locfit(~lp(x,h=2*bandwidth,deg=1)) }
    target <- lfeval(fit)$xev
    nt = length(target)
  }
  if (alldata==TRUE) {
    target <- x
    nt = n 
  }

  namevect <- names(mat)
  yname = namevect[1]
  xname = namevect[2]
  ntau = length(taumat)
  ytarget <- array(0,dim=c(nt,ntau))
  yhat    <- array(0,dim=c(n,ntau))
  
  for (i in seq(1:nt)) {
    xtarget <- x - target[i]
    dist <- abs(xtarget)/sdx
    if (window>0) {h = quantile(dist,window)}
    if (bandwidth>0) {h = bandwidth}
    samp <- dist<=h
    if (kern=="gauss") {samp <- dist<=max(dist)}
    wx <- wgt(dist[samp]/h)
    for (j in seq(1:ntau)) {
      fit <- rq(y[samp]~xtarget[samp],weights=wx,tau=taumat[j],ci=FALSE)
      ytarget[i,j] <- fit$coef[1] 
    }
  }

  if (alldata==FALSE) {
    for (j in seq(1,ntau)) {
      yhat[,j] <- aspline(target,ytarget[,j],x)$y
    }
  }
  if (alldata==TRUE) {yhat <- ytarget}

  predx_yhat <- NULL
  if (!identical(predx,0)) {
    predx_yhat <- yhat
    for (j in seq(1,ntau)) {
      predx_yhat[,j] <- aspline(x,yhat[,j],predx)$y
    }
  }
    
  o <- order(x)
  if (graph.yhat==TRUE) {
    ymin = min(yhat)
    ymax = max(yhat)
    plot(x[o],yhat[o,1],type="l",xlab=xname,ylab=yname,ylim=c(ymin,ymax))
    if (ntau>1) {
      for (j in seq(2,ntau)) {
        lines(x[o],yhat[o,j])
      }
    }
  }
  if (graph.predx==TRUE) {
    ymin = min(predx_yhat)
    ymax = max(predx_yhat)
    plot(predx,predx_yhat[,1],type="l",xlab=xname,ylab=yname,ylim=c(ymin,ymax))
    if (ntau>1) {
      for (j in seq(2,ntau)) {
        lines(predx,predx_yhat[,j])
      }
    }
  }

  out <- list(yhat,predx_yhat,taumat,window,bandwidth)
  names(out) <- c("yhat","predx_yhat","taumat","window","bandwidth")
  return(out)
}

