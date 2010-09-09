lwr1 <-
function(form,window=.25,bandwidth=0,kern="tcub",search.method="gcv",search.print=TRUE,alldata=FALSE,predx=0,graph.yhat=FALSE,
  graph.predx=FALSE,data=NULL) {
  attach(data,warn.conflicts=FALSE)

  nw = length(window)
  nb = length(bandwidth)
  if ((nw>1)&(nb>1))   { cat("Both window and bandwidth vectors specified; will use window vector","\n") }
  if (nw>1) {
    h = lwr1_minh(form,window=window,bandwidth=0,kern=kern,method=search.method,print=search.print,alldata=alldata)
    window = h
  }
  if (nb>1&nw==1) {
    h = lwr1_minh(form,window=0,bandwidth=bandwidth,kern=kern,method=search.method,print=search.print,alldata=alldata)
    bandwidth = h
  }
  
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
  
  ytarget     <- array(0,dim=nt)
  ytarget.se  <- array(0,dim=nt)
  dtarget    <- array(0,dim=nt)
  dtarget.se <- array(0,dim=nt)
  df1target   <- array(0,dim=nt)
  df2target   <- array(0,dim=nt)

  for (i in seq(1:nt)) {
    dist <- abs(x-target[i])/sdx
    if (window>0) {h = quantile(dist,window) }
    if (bandwidth>0) {h = bandwidth}
    samp <- dist<=h
    if (kern=="gauss") {samp <- dist<=max(dist)}
    
    xmat1 <- cbind(1,x[samp]-target[i])
    k <- wgt(dist[samp]/h)
    xmat2 <- k*xmat1
    xx <- solve(crossprod(xmat1,xmat2))
    xmat1 <- xx%*%t(xmat2)
    bmat <- xmat1%*%y[samp]
    ytarget[i] = bmat[1,1]
    dtarget[i] = bmat[2,1]
    df1target[i] = xx[1,1]
    df2target[i] = sum(xmat1[1,]^2)
    vmat <- tcrossprod(xmat1)
    ytarget.se[i] = sqrt(vmat[1,1])
    dtarget.se[i] = sqrt(vmat[2,2])
  }

  if (alldata==FALSE) {
    hat <- aspline(target,ytarget,x)
    yhat <- hat$y
    hat <- aspline(target,dtarget,x)
    dhat <- hat$y
    hat <- aspline(target,df1target,x)
    infl <- hat$y
    df1 = sum(infl)
    hat <- aspline(target,df2target,x)
    df2 = sum(hat$y)
    hat <- aspline(target,ytarget.se,x)
    yhat.se <- hat$y
    hat <- aspline(target,dtarget.se,x)
    dhat.se <- hat$y
  }

  if (alldata==TRUE) {
    yhat <- ytarget
    dhat <- dtarget
    yhat.se  <- ytarget.se
    dhat.se <- dtarget.se
    infl <- df1target
    df1 = sum(infl)
    df2 = sum(df2target)
  }

  rss = sum((y-yhat)^2)
  sig2 = rss/(n-2*df1 + df2)
  cv = mean(((y-yhat)/(1-infl))^2)
  gcv = n*rss/((n-df1)^2)
  ytarget.se  <- sqrt(sig2)*ytarget.se
  dtarget.se <- sqrt(sig2)*dtarget.se
  yhat.se  <- sqrt(sig2)*yhat.se
  dhat.se <- sqrt(sig2)*dhat.se

# predicted values 
  predx_yhat <- NULL
  predx_dhat <- NULL
  predx_yhat.se <- NULL
  predx_dhat.se <- NULL
  if (!identical(predx,0)) {
    hat <- aspline(x, yhat, predx)
    predx_yhat <- hat$y
    hat <- aspline(x, dhat, predx)
    predx_dhat <- hat$y
    hat <- aspline(x, yhat.se, predx)
    predx_yhat.se <- hat$y
    hat <- aspline(x, dhat.se, predx)
    predx_dhat.se <- hat$y
  }
  
  if (graph.yhat==TRUE) {
    o <- order(x)
    plot(x[o],yhat[o],type="l",xlab="x",ylab="Predicted y")
  }
  if (graph.predx==TRUE) {plot(predx,predx_yhat,type="l",xlab="Target x",ylab="Predicted y") }

 out <- list(target, ytarget, dtarget, ytarget.se, dtarget.se, yhat, dhat, yhat.se, dhat.se, df1, df2, sig2, cv, gcv, infl,
   predx_yhat, predx_dhat, predx_yhat.se, predx_dhat.se,window,bandwidth)
 names(out) <- c("target","ytarget","dtarget","ytarget.se","dtarget.se","yhat","dhat","yhat.se","dhat.se","df1","df2","sig2","cv","gcv","infl",
   "predx_yhat", "predx_dhat", "predx_yhat.se", "predx_dhat.se","window","bandwidth")

 return(out) 
}

