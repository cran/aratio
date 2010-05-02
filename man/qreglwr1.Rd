\name{qreglwr1}
\alias{qreglwr1}
\title{Nonparametric quantiles based on weighted quantile regressions }
\description{
Estimates nonparametric quantile functions by applying a weight function to the standard quantreg estimator. }
\usage{
qreglwr1(form, window = .50, bandwidth = 0, kern = "tcub", 
  taumat = c(0.1, 0.25, 0.5, 0.75, 0.9), alldata = FALSE, 
  predx = 0, graph.yhat = FALSE, graph.predx = FALSE, data=NULL)

}

\arguments{
  \item{form}{ A formula of the type y~x. }
  \item{window}{Window size.  Default:  0.50}
  \item{bandwidth}{Bandwidth.  Default:  not specified.  If bandwidth is specified, the \emph{window} option is ignored.
  A common rule of thumb is bandwidth = 1.06*min(sd(x), (quantile(y,.75)-quantile(x,.25))/1.349)*(n^(-.2)). }
  \item{kern}{Kernel weighting functions.  Default is the tri-cube. Options include "rect", "tria", "epan", "bisq", "tcub", "trwt", and "gauss".}
  \item{taumat}{Vector of quantiles.   Default:  taumat=c(.10, .25, .50, .75, .90). }
  \item{alldata}{ If alldata=T, each observation is used as a target value for x.  
When alldata=F, the function is estimated at a set of points chosen by the locfit program using an adaptive decision tree approach, 
and the \emph{akima} program is used to interpolate to the full set of observations.  Specifying alldata=T can lead to long estimation times.}
  \item{predx}{
A vector of arbitrary values to form predictions for the quantile functions. 
Helpful when comparing predictions across similar models with different numbers of observations. 
Default is predx=0, meaning predictions are only formed for the full data set rather than an alternative set of arbitrary points.
}
  \item{graph.yhat}{
If graph.yhat=T, graph of results is produced using all observations. Default:  graph.yhat=F.
}
  \item{graph.predx}{
If graph.predx=T, graph of results is produced using predx data points. Default:  graph.predx=F.
}
  \item{data }{A data frame containing the data. Default:  use data in the current working directory.}
}

\value{

\item{yhat}{Matrix of quantile predictions. Dimension is n x length(taumat)}
\item{predx_yhat}{Matrix of quantile predictions at predx values.  Dimension is length(predx) x length(taumat).}
\item{taumat}{Vector of target quantile values.}
\item{window}{Window size.}
\item{bandwidth}{Bandwidth.}
}

\references{
Koenker, Roger.  \emph{Quantile Regression.}  New York:  Cambridge University Press, 2005.
Chapter 7 and Appendix A.9.
}
\author{
Daniel McMillen
}
\seealso{
 \link{condens1}
}

\examples{
set.seed(23849103)
n = 1000
x <- runif(n,0,2*pi)
x <- sort(x)
ybase <- x - .1*(x^2) + sin(x) - cos(x) -.5*sin(2*x) + .5*cos(2*x)
sig = sd(ybase)/2
y <- ybase + rnorm(n,0,sig)
fit <- qreglwr1(y~x, graph.yhat=TRUE)

# Plot results after estimation
data(dupage99)
dupage99$ratio <- dupage99$av/dupage99$price
dropobs1 <- nptrim_obs(dupage99$ratio,k=3)
dropobs2 <- (dupage99$price<quantile(dupage99$price,.01))|
  (dupage99$price>quantile(dupage99$price,.99))
dupage99 <- dupage99[dropobs1==FALSE&dropobs2==FALSE,]
o <- order(dupage99$price)
dupage99 <- dupage99[o,]
attach(dupage99)
price <- price/1000

fit <- qreglwr1(ratio~price,window=.3)
ymin = min(fit$yhat)
ymax = max(fit$yhat)
plot(price, fit$yhat[,1],type="l",xlab="Sales Price (1000s)",ylab="Assessment Ratio",
  ylim=c(ymin,ymax),main="Nonparametric Quantile Regression")
for (j in seq(2,5)) {
  lines(price,fit$yhat[,j])
}
fit$window
fit$bandwidth
}


\details{
Serves as an interface to the quantreg package.  
Uses a kernel weight function in quantreg's "weight" option to estimate quantile regressions at a series of target values of x.  
The target values are found using locfits adaptive decision tree approach.   
The predictions are then interpolated to the full set of x values using the akima package.  
If alldata=T, the procedure is applied to every value of x rather than a set of target points.  
If a vector is listed for predx, akima is used to interpolate the estimates to the points listed in predx.

Fixed bandwidths may prove too small if there are regions where x is sparse.  A nearest neighbor approach is generally preferable (e.g, window=.50).  


Available kernel weighting functions include the following:  

\tabular{lll}{
Kernel \tab Call abbreviation  \tab Kernel function K(z) \cr
Rectangular \tab ``rect'' \tab \eqn{\frac{1}{2} I(|z|<1)}{1/2 * I(|z|<1)} \cr
Triangular \tab ``tria'' \tab \eqn{(1-|z|)I(|z|<1)}{(1-|z|) * I(|z|<1)}\cr
Epanechnikov \tab ``epan'' \tab \eqn{\frac{3}{4} (1-z^2) * I(|z| <1)}{3/4 * (1-z^2)*I(|z| < 1)} \cr
Bi-Square \tab ``bisq'' \tab \eqn{\frac{15}{16} (1-z^2)^2 * I(|z| <1)}{15/16 * (1-z^2)^2 * I(|z| < 1)} \cr
Tri-Cube \tab ``tcub'' \tab \eqn{\frac{70}{81} (1-|z|^3)^3 * I(|z| <1)}{70/81 * (1-|z|^3)^3 * I(|z| < 1)} \cr
Tri-Weight \tab ``trwt'' \tab \eqn{\frac{35}{32} (1-z^2)^3 * I(|z| <1)}{35/32 * (1-z^2)^3 * I(|z| < 1)} \cr
Gaussian \tab ``gauss'' \tab \eqn{(2\pi)^.5 e^{z^2/2}}{2pi^.5 exp(z^2/2)} \cr }
}




