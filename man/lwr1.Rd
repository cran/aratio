\name{lwr1}
\alias{lwr1}
\title{Locally Weighted Regression}
\description{
Estimates a model of the form y = f(x) using locally weighted regression.  x must be a single variable.  
Returns estimated values, derivatives, and standard errors for both the estimated values of f(x) and df(x)/dx.}
\usage{
lwr1(form, window = .25, bandwidth = 0, kern = "tcub", 
  search.method = "gcv", search.print = TRUE, alldata = FALSE, 
  predx = 0, graph.yhat = FALSE, graph.predx = FALSE, data=NULL)

}
\arguments{
  \item{form}{Formula of the type y~x.}
  \item{window}{
Window size.   Default:  h = .25.  Estimated using the indicated window size if specified as a scalar.  
If specified as a vector (e.g., window = seq(from=.10, to=.9, by=.10)), lwr1 chooses the value that minimizes a criterion function.  
The criterion is specified using \emph{search.method.}
}
  \item{bandwidth}{

Bandwidth.   Default: not used.  If bandwidth is specified, the \emph{window} option is ignored.
A common rule of thumb is bandwidth = \eqn{1.06*sdx*(n^{-.2})} where sdx = 
min(sd(x) , (quantile(x,.75) - quantile(x,.25))/1.349 ).   
The model is estimated using the indicated bandwidth if specified as a scalar (e.g., bandwidth=.25).  
If specified as a vector (e.g., bandwidth = seq(from=.10, to=.9, by=.10)), lwr1 chooses the value that minimizes a criterion function.  
The criterion is specified using \emph{search.method.}
}
  \item{kern}{
Kernel weighting functions.  Default is the tri-cube. Options include  "rect", "tria", "epan", "bisq", "tcub", "trwt", and "gauss".
}
  \item{search.method}{
The selection criterion when searching for the optimal value for the bandwidth or window.  Must be in quotes.  
The default is the generalized cross-validation criterion, or "gcv".  The alternative is the cross-validation criterion, or "cv".  
See the "values" section for definitions.
}
  \item{search.print}{
If search.print=T, a summary of the search results is printed when a vector of possible values is provided for window or bandwidth.  Default is search.print=T.
}
  \item{alldata}{
If alldata=T, each observation is used as a target value for x.  
When alldata=F, the function is estimated at a set of points chosen by the locfit program using an adaptive decision tree approach, 
and the akima program is used to interpolate to the full set of observations.  Specifying alldata=T can lead to long estimation times.
}
  \item{predx}{
A vector of arbitrary values to form predictions for y.  Helpful when comparing predictions across similar models with different numbers of observations.  
Default is predx=0, meaning predictions are only formed for the full data set rather than for a set of arbitrary points.
}
  \item{graph.yhat}{If graph.yhat=T, produces a plot of the predicted values at all data points.  Default is graph.yhat=F.
}
  \item{graph.predx}{
If graph.predx=T and a vector is provided for \emph{predx}, produces a plot of the predicted values the \emph{predx} points.  Default is graph.predx=F.
}
  \item{data }{A data frame containing the data. Default:  use data in the current working directory.}
}

\value{
\item{target}{The target points for the original estimation of the function.}
\item{ytarget}{The predicted values of y at the original target points.}
\item{dtarget}{The estimated derivatives dy/dx at the target points.}
\item{ytarget.se}{Standard errors for the predicted values of y at the target points.}
\item{dtarget.se}{Standard errors for the derivatives dy/dx at the target points.}
\item{yhat}{The predicted values of y for the full data set.}
\item{dhat}{The estimated derivatives dy/dx for the full data set.}
\item{yhat.se}{Standard errors for the predicted values of y, full data set.}
\item{dhat.se}{Standard errors for estimated derivatives dy/dx for the full data set.}
\item{df1}{tr(L), a measure of the degrees of freedom used in estimation.}
\item{df2}{tr(L'L), an alternative measure of the degrees of freedom used in estimation.}
\item{sig2}{Estimated residual variance, sig2 = rss/(n-2*df1+df2).}
\item{cv}{Cross-validation measure.  cv = mean(((y-yhat)/(1-infl))^2), where yhat is vector of predicted values for y and infl is the vector of diagonal terms for L.}
\item{gcv}{gcv = n*(n*sig2)/((n-nreg)^2), where sig2 is the estimated residual variance and nreg = 2*df1 - 2f2.}
\item{infl}{A vector containing the diagonal elements of L.}
\item{predx_yhat}{The predicted values of y at the predx values.}
\item{predx_dhat}{The estimated derivatives dy/dx at the predx values.}
\item{predx_yhat.se}{Estimated standard errors of the predicted values for y at the predx values.}
\item{predx_dhat.se}{Estimated standard errors for dy/dx at the predx values of x.}
\item{window}{Window size used in estimation.}
\item{bandwidth}{Bandwidth used in estimation.}
}



\references{

Cleveland, William S. and Susan J. Devlin, ``Locally Weighted Regression:  An Approach to Regression Analysis by Local Fitting,'' \emph{Journal of the American Statistical Association} 83 (1988), 596-610.

Loader, Clive.  \emph{Local Regression and Likelihood.}  New York:  Springer, 1999.

McMillen, Daniel P., ``Issues in Spatial Data Analysis,'' \emph{Journal of Regional Science} 50 (2010), 119-141. 

McMillen, Daniel P., ``Employment Densities, Spatial Autocorrelation, and Subcenters in Large Metropolitan Areas,'' \emph{Journal of Regional Science} 44 (2004), 225-243.

McMillen, Daniel P. and John F. McDonald, ``A Nonparametric Analysis of Employment Density in a Polycentric City,'' \emph{Journal of Regional Science} 37 (1997), 591-612.

McMillen, Daniel P. and Christian Redfearn, ``Estimation and Hypothesis Testing for Nonparametric Hedonic House Price Functions,'' 
\emph{Journal of Regional Science} 50 (2010), 712-733.

Pagan, Adrian and Aman Ullah.  \emph{Nonparametric Econometrics.}  New York:  Cambridge University Press, 1999.

Silverman, B.W.  \emph{Density Estimation for Statistics and Data Analysis.}  London:  Chapman and Hall, 1986. 
}


\author{
Daniel McMillen
}
\seealso{
 \link{cubespline1}

 \link{fourier1}
}



\details{
The estimated value of y at a target value x0 is
\deqn{\hat{\theta}(x_0)=(\sum_{i=1}^n K(\phi_i)Z_i Z_i^\top)^{-1} \sum_{i=1}^n K(\phi_i) Z_i ^\top y_i}{theta hat(x_0)=(sum_{i=1 to n} K(phi_i) Z_i Z_i')^-1 * sum_{i=1 to n} K(phi_i) Z_i' y_i,} 
where \eqn{Z_i = (1,  x_i)} and theta is a vector holding the corresponding coefficients \eqn{\theta=(\alpha, \beta)}{theta=(alpha, beta)}. 

\eqn{K(\phi_i)}{K(phi_i)} is a kernel weighting function, with \eqn{\phi_i=(x_i-x_0)/h}{phi_i=(x_i-x_0)/h}. 
h is specified by the bandwidth or window options. 
The intercept, \eqn{\alpha}{alpha} provides an estimate of y at \eqn{x_0} and \eqn{\beta}{beta} provides an estimate of the slope, dy/dx at \eqn{x_0}. 
When alldata=T, each data point in turn is used as a target point. 
Since each estimate is a linear function of all n values for y, the full set of estimates takes the form \eqn{\hat{Y}=LY}{Yhat=LY}, where L is an nxn matrix. 
Loader (1999) suggests two measures of the number of degrees of freedom used in estimation: df1=tr(L) and df2=tr(L'L), both of which are stored by \emph{lwr1}. 
The diagonal elements of tr(L) are stored in the array \emph{infl}.
 Again, following Loader (1999), the degrees of freedom correction used to estimate the error variance \eqn{\sigma^2}{sigma2} is df=2*df1-df2. 
Let \eqn{e_i} represent the residual for observation i, or \eqn{e_i=y_i - \hat{y_i}}{e_i=y_i - y_i hat}, 

the estimated variance is \eqn{\hat{\sigma}^2 = \frac{1}{n-df} \sum_{i=1}^n e_i^2}{sigmahat^2 = 1/(n-df) * sum_{i=1 to n} e_i^2}. 
The covariance matrix estimate for 
\eqn{\hat{\theta}(x_0)}{thetahat(x_0)} is \deqn{\hat{\sigma}^2(\sum_{i=1}^n Z_i K(\phi_i) Z_i^\top)^{-1}(\sum_{i=1}^n Z_i (K(\phi_i))^2 Z_i^\top )(\sum_{i=1}^n Z_i K(\phi_i) Z_i^\top)^{-1}.}{sigmahat^2(\sum_{i=1 to n} Z_i K(phi_i) Z_i')^{-1}(\sum_{i=1 to n} Z_i (K(phi_i))^2 Z_i')(\sum_{i=1 to n} Z_i K(phi_i) Z_i')^{-1}.}

Estimation can be very slow when alldata=T. When alldata=F, the package locfit is used to find a good set of target points at which to evaluate the function. 
See Loader (1999, section 12.2) for a description of the algorithm used to determine the target points. 
The akima package is then used to interpolate the coefficient estimates, the standard errors,and the values used to form df1 and df2. 

h can be specified to be either a fixed bandwidth or a window size set to a percentage of the sample size. 
Optionally, a vector of values can be specified for h with lwr1 picking the one that minimizes a criterion function. 
Fixed bandwidths may prove too small if there are regions where x is sparse.  A nearest neighbor approach is generally preferable (e.g, window=.50).  

Available kernel weighting functions include the following: 

\tabular{lll}{
Kernel \tab Call abbreviation \tab Kernel function K(z) \cr
Rectangular \tab ``rect'' \tab \eqn{\frac{1}{2} I(|z|<1)}{1/2 * I(|z|<1)} \cr
Triangular \tab ``tria'' \tab \eqn{(1-|z|)I(|z|<1)}{(1-|z|) * I(|z|<1)}\cr
Epanechnikov \tab ``epan'' \tab \eqn{\frac{3}{4} (1-z^2) * I(|z| <1)}{3/4 * (1-z^2)*I(|z| < 1)} \cr
Bi-Square \tab ``bisq'' \tab \eqn{\frac{15}{16} (1-z^2)^2 * I(|z| <1)}{15/16 * (1-z^2)^2 * I(|z| < 1)} \cr
Tri-Cube \tab ``tcub'' \tab \eqn{\frac{70}{81} (1-|z|^3)^3 * I(|z| <1)}{70/81 * (1-|z|^3)^3 * I(|z| < 1)} \cr
Tri-Weight \tab ``trwt'' \tab \eqn{\frac{35}{32} (1-z^2)^3 * I(|z| <1)}{35/32 * (1-z^2)^3 * I(|z| < 1)} \cr
Gaussian \tab ``gauss'' \tab \eqn{(2\pi)^.5 e^{z^2/2}}{2pi^.5 exp(z^2/2)} \cr }
}

\examples{
	set.seed(23849103)
	n = 1000
	x <- runif(n,0,2*pi)
	x <- sort(x)
	ybase <- x - .1*(x^2) + sin(x) - cos(x) -.5*sin(2*x) + .5*cos(2*x)
	sig = sd(ybase)/2
	y <- ybase + rnorm(n,0,sig)
        fit <- lwr1(y~x, window=seq(from=.05, to=.15, by = .01))
        plot(x,y,xlab="x",ylab="y")
        lines(x,fit$yhat,col="red")

}
