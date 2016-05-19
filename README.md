
fmlogit: Estimate Fractional Multinomial Logit Models
=====================================================
Description
-----------
Used to estimate fractional multinomial logit models using quasi-maximum likelihood estimations following Papke and Wooldridge(1996).

Usage
-----
fmlogit(y, X, beta0 = NULL, MLEmethod = "CG", maxit = 50000,
  abstol = 1e-04, ...)
Arguments
---------
y	
the dependent variable (N*J). Can be a matrix or a named data frame. The first column of the matrix is automatically treated as the baseline.

X	
independent variable (N*K). Can be a matrix or a named data frame. If there is no intercept term in the X, then an intercept term is automatically added.

beta0	
Initial value for beta used in optimization. Uses a 1*K(J-1) vector. Default to a vector of zeros.

MLEmethod	
Method of optimization. Goes into optim(method=MLEmethod)). Default to "CG", the conjugate gradients method.

maxit	
Maximum number of iteration.

abstol	
Tolerence.

...	
additional parameters that goes into optim()

Details
-------
The fractional multinomial model is the expansion of the multinomial logit to fractional responses. Unlike standard multinomial logit models, which only considers 0-1 respones, fractional multinomial model considers the case where the response variable is fractions that sums up to one. Examples of these type of data are, percentages of budget spent in education, defense, public health; fractions of a population that have middle school, high school, college, or post college education, etc.

This function follows Papke and Wooldridge(1996)'s paper, in which they proposed a quasi-maximum likelihood estimator for fractional response data. The likelihood function used here is a standard multinomial likelihood function, see http://maartenbuis.nl/software/likelihoodFmlogit.pdf for the likelihood used here. Robust standard errors are provided following Papke and Wooldridge(1996), in which they proposed an asymptotically consistent estimator of variance.

Maximization is done by calling optim. It seems that difference choice of maximization method may yield drastically different parameter estimates. Using the conjugate gradients method yields the same results as using Stata's fmlogit package by Maarten Buis.

Value
-----
The function returns a list of elements below:

$estimates A list of matrices containing parameter estimates, standard errors, and hypothesis testing results.

$likelihood The likelihood value

$convergence Convergence diagnostics code

$count Provides dataset information

References
----------
Papke, L. E. and Wooldridge, J. M. (1996), Econometric methods for fractional response variables with an application to 401(k) plan participation rates. J. Appl. Econ., 11: 619-632.

Examples
--------
```Rcode
require(foreign)
data = read.dta("http://fmwww.bc.edu/repec/bocode/c/citybudget.dta")
X = data[,2:5]
y = data[,6:11]
results1 = fmlogit(y,X)
```
