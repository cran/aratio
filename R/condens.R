condens <-
function(x,y,window=.7,kern="tcub",mingrid.x=min(x),maxgrid.x=max(x),mingrid.y=min(y),maxgrid.y=max(y),ngrid=50,
  contour=TRUE,level=FALSE,wire=FALSE,dens=FALSE,namex="Sales Price",namey="Assessment Ratio",targetx.dens=NULL,
  quantile.dens=c(.10,.25,.50,.75,.90),data=NULL) {
  attach(data,warn.conflicts=FALSE)

  n = length(x)
  sdx = 1.06*min(sd(x), (quantile(x,.75)-quantile(x,.25))/1.349)
  sdy = 1.06*min(sd(y), (quantile(y,.75)-quantile(y,.25))/1.349)

  fitx  <- locfit(~lp(x/sdx,nn=window),kern=kern)
  fitxy <- locfit(~lp(x/sdx,y/sdy,nn=window),kern=kern)
  fxhat <- fitted(fitx)/sdx
  fxyhat <- fitted(fitxy)/(sdx*sdy)
  fxy <- fxyhat/fxhat

  grid.x <- seq(mingrid.x/sdx,maxgrid.x/sdx,length=ngrid) 
  grid.y <- seq(mingrid.y/sdy,maxgrid.y/sdy,length=ngrid)
  xy <- expand.grid(grid.x,grid.y)
  grid.fxy <- predict(fitxy,as.matrix(xy))/predict(fitx,xy[,1])
  grid.fxy <- grid.fxy/sdy
  grid.x <- xy[,1]*sdx
  grid.y <- xy[,2]*sdy
  gridmat <- cbind(grid.x,grid.y,grid.fxy)

  if (contour==TRUE) {print(contourplot(grid.fxy~grid.x*grid.y,xlab=namex,ylab=namey))}
  if (level==TRUE)     {print(levelplot(grid.fxy~grid.x*grid.y,xlab=namex,ylab=namey))}
  if (wire==TRUE)      {print(wireframe(grid.fxy~grid.x*grid.y,xlab=namex,ylab=namey,zlab="Density"))}


  if (identical(targetx.dens,NULL)) {targetx.dens <- quantile(x,quantile.dens) }
  if (length(targetx.dens)>5) {
    cat("Target x for 2d density graph > 5;  will only use first 5","\n")
    targetx.dens <- targetx.dens[1:5]
  }
  nq <- length(targetx.dens)
  densmat <- array(0,dim=c(n,nq))
  for (j in seq(1:nq)) {
    densmat[,j] <- predict(fitxy,cbind(targetx.dens[j]/sdx,y/sdy))/(predict(fitx,targetx.dens[j]/sdx)*sdy)
  }
  if (dens==TRUE) {
    colmat <- c("black","blue","red","green","orange")
    o <- order(y)
    plot(y[o],densmat[o,1],xlab=namey,ylab="Conditional Density",type="l",ylim=c(min(densmat),max(densmat)))
    for (j in seq(2,nq)) {
      lines(y[o],densmat[o,j],col=colmat[j])
    }
    legend("topright",as.character(targetx.dens),col=colmat[1:nq],lwd=1) 
  }

  out <- list(fxy, gridmat, densmat)
  names(out) <- c("fxy","gridmat","densmat")
  return(out)
}

