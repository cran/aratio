\name{condens}
\alias{condens}
\title{
Conditional density estimation
}
\description{condens estimates conditional density functions of the form f(y| x) = f(x,y)/f(x).  
Kernel density estimators are used to estimate f(x,y) and f(x).  
The conditional density function can be plotted as a three-dimensional surface or as contour map.  
Alternatively, the conditional density of y can be graphed for as many as five target values of x.  
Though by default x is assumed to represent sales price and y is the assessment ratio, 
condens can be used to estimate conditional density functions for any variables x and y.
}

\usage{
condens(x, y, window = 0.7, kern = "tcub", mingrid.x = min(x), 
  maxgrid.x = max(x), mingrid.y = min(y), maxgrid.y = max(y), 
  ngrid = 50, contour = TRUE, level = FALSE, wire = FALSE, 
  dens = FALSE, namex = "Sales Price", namey = "Assessment Ratio", 
  targetx.dens = NULL, quantile.dens = c(0.1, 0.25, 0.5, 0.75, 0.9), 
  data=NULL)
}

\arguments{
  \item{x}{Price variable.}
  \item{y}{Ratio variable.}
  \item{window}{Window size for kernel density functions.  Default is .7.}
  \item{kern}{Kernel weighting functions.  Default is the tri-cube. Options include  "rect", "tria", "epan", "bisq", "tcub", "trwt", and "gauss".}
  \item{mingrid.x, maxgrid.x, mingrid.y, maxgrid.y,
ngrid}{The mingrid and maxgrid values are the boundaries for the ngrid x ngrid lattice used in the graphs produced by condens.  
By default, mingrid.x = min(x), maxgrid.x = max(x), mingrid.y = min(y), maxgrid.y = max(y), and ngrid=50.  
}
  \item{contour}{If contour=T, produces a two-dimensional contour plot of the conditional density estimates.  
Evaluated at the ngrid x ngrid lattice.  Default is contour=T.}
  \item{level}{If level=T, produces a two-dimensional level plot of the conditional density estimates.  
Evaluated at the ngrid x ngrid lattice.  Default is level=F.}
  \item{wire}{If wire=T, produces a three-dimensional plot of the conditional density estimates.  Evaluated at the ngrid x ngrid lattice.  Default is wire=T.}
  \item{dens}{If dens=T, produces a plot showing how f(y|x) varies over y given target values of x.  
Target values of x are provided using the targetx.dens or quantile.dens options.  Default is dens=F.}
  \item{namex}{Label for x axis in contour, level, and wire plots.  Default is "Sales Price".}
  \item{namey}{Label for y axis in contour, level, wire, and dens plots.  Default is "Assessment Ratio".}
  \item{targetx.dens}{Target values for x in density plots, e.g, targetx.dens = c(200,400,600).  Maximum number of entries is 5.  
If targetx.dens has more than 5 entries, only the first 5 will be used.  Default is targetx.dens =NULL, 
meaning that the target values for x are determined by the quantile.dens option.}
  \item{quantile.dens}{Quantiles for the target values for x in density plots, e.g, quantile.dens = c(.25,.50,.75).  Maximum number of entries is 5.  
If quantile.dens has  more than 5 entries, only the first 5 will be used.  Default is q}
  \item{data }{A data frame containing the data. Default:  use data in the current working directory.}
}


\value{

\item{fxy}{The values of f(y|x), one for each data point.}

\item{gridmat}{An (ngrid*ngrid)x3 matrix used to produce the contour, level, and wire maps.  
The first column contains the lattice values for x, the second column contains the lattice values for y, and the third column is the estimated 
value of f(y|x) at the target values for x and y.}
\item{densmat}{The estimated values of f(y|x) for the two-dimensional dens graphs.  
If the number of observations in the call to condens is n and the number of entries in quantile.dens is nq, then densmat is an n x nq matrix.}
}

\references{
Li, Oi and Jeffrey Scott Racine.  \emph{Nonparametric Econometrics:  Theory and Practice.}  Princeton, NJ:  Princeton University Press, 2007.

Loader, Clive.  \emph{Local Regression and Likelihood.}  New York:  Springer, 1999.

Pagan, Adrian and Aman Ullah.  \emph{Nonparametric Econometrics.} New York:  Cambridge University Press, 1999.
}

\author{
Daniel McMillen
}

\examples{
data(dupage99)
dupage99$ratio <- dupage99$av/dupage99$price
dropobs1 <- nptrim_obs(dupage99$ratio,k=3)
dropobs2 <- (dupage99$price<quantile(dupage99$price,.01))|
  	(dupage99$price>quantile(dupage99$price,.99))
dupage99 <- dupage99[dropobs1==FALSE&dropobs2==FALSE,]
o <- order(dupage99$price)
dupage99 <- dupage99[o,]
attach(dupage99)

# show all plots
fit <- condens(price, ratio, window=.5, contour=TRUE, level=TRUE, wire=TRUE, 
  dens=TRUE, targetx.dens=c(100,200,300,400), ngrid=40)

# skip plots in base command; add them afterward
fit <- condens(price, ratio, window=.5, contour=FALSE, level=FALSE, wire=FALSE, 
  dens=FALSE, targetx.dens=c(100,200,300,400), ngrid=40)

contourplot(fit$gridmat[,3]~fit$gridmat[,1]*fit$gridmat[,2], xlab="Price", 
ylab="Assessment Ratio")
levelplot(fit$gridmat[,3]~fit$gridmat[,1]*fit$gridmat[,2], xlab="Price", 
ylab="Assessment Ratio")
wireframe(fit$gridmat[,3]~fit$gridmat[,1]*fit$gridmat[,2], xlab="Price", 
ylab="Assessment Ratio", zlab="Density")
}


\details{
Uses the locfit package to estimate f(x,y) and f(x) and then forms f(y|x) = f(x,y)/f(x).  
The procedure is very fast because locfit evaluates the f(x,y) and f(x) and carefully chosen 
target values and then interpolates to either the full data set or to a set of arbitrary alternative points.  
Condens uses the following commands to produce the conditional density estimates:

n = length(x) \cr
sdx = 1.06*min(sd(x), (quantile(x,.75)-quantile(x,.25))/1.349) \cr
sdy = 1.06*min(sd(y), (quantile(y,.75)-quantile(y,.25))/1.349) \cr
fitx  <- locfit(~lp(x/sdx,nn=window),kern=kern) \cr
fitxy <- locfit(~lp(x/sdx,y/sdy,nn=window),kern=kern) \cr
fxhat <- fitted(fitx)/sdx \cr
fxyhat <- fitted(fitxy)/(sdx*sdy) \cr
fxy <- fxyhat/fxhat \cr

The locfit command "predict" is then used to interpolate to other points.  For example, the following commands produce estimates of the conditional 
density function at the grid of values used to produce the contour, level, and wire plots:

  grid.x <- seq(mingrid.x/sdx,maxgrid.x/sdx,length=ngrid)  \cr
  grid.y <- seq(mingrid.y/sdy,maxgrid.y/sdy,length=ngrid) \cr
  xy <- expand.grid(grid.x,grid.y) \cr
  grid.fxy <- predict(fitxy,as.matrix(xy))/predict(fitx,xy[,1]) \cr
  grid.fxy <- grid.fxy/sdy \cr
  grid.x <- xy[,1]*sdx \cr
  grid.y <- xy[,2]*sdy \cr
  gridmat <- cbind(grid.x,grid.y,grid.fxy) \cr


Available kernel weighting functions include the following:


\tabular{lll}{
Kernel \tab Call abbreviation \tab Kernel function K(z) \cr
Rectangular \tab ``rect'' \tab \eqn{\frac{1}{2} I(|z| <1)}{1/2 * I(|z|<1)} \cr
Triangular \tab ``tria'' \tab \eqn{(1-|z|)I(|z|<1)}{(1-|z|) * I(|z|<1)}\cr
Epanechnikov \tab ``epan'' \tab \eqn{\frac{3}{4} (1-z^2) * I(|z| <1)}{3/4 * (1-z^2)*I(|z| < 1)} \cr
Bi-Square \tab ``bisq'' \tab \eqn{\frac{15}{16} (1-z^2)^2 * I(|z| <1)}{15/16 * (1-z^2)^2 * I(|z| < 1)} \cr
Tri-Cube \tab ``tcub'' \tab \eqn{\frac{70}{81} (1-|z|^3)^3 * I(|z| <1)}{70/81 * (1-|z|^3)^3 * I(|z| < 1)} \cr
Tri-Weight \tab ``trwt'' \tab \eqn{\frac{35}{32} (1-z^2)^3 * I(|z| <1)}{35/32 * (1-z^2)^3 * I(|z| < 1)} \cr
Gaussian \tab ``gauss'' \tab \eqn{(2\pi)^.5 e^{z^2/2}}{2pi^.5 exp(z^2/2)} \cr }
}

