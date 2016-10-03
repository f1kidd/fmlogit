---
title: "The fmlogit Package: A Document"
author: "Xinde James Ji" 
date: "Thursday, May 19, 2016"
output: pdf_document
---

This document provides documentations for the fmlogit package in R. Updates will be published at [my github site](https://github.com/f1kidd/fmlogit). Any suggestions or concerns are welcomed. 

# What is fractional multinomial logit?
Fractional multinomial responses arises naturally in various occasions. For example, a municipality allocates its budgets to multiple departments, and we are interested in the proportion of the budgets that each department receives. Or, there are multiple candidates in a presendential election, and we are interested in the percentage of support for each candidate in each state. 

# So, why do we even need fmlogit in R? 
Don't we already have an fmlogit module in Stata? Yes, and you are very welcome to [check that out](http://maartenbuis.nl/software/fmlogit.html). 

However, this package offers several advantages over Stata's fmlogit module, namely,
## 1. Integration with the R Platform
Implementating the model in R offers the opportunity to integrate the whole empirical process. With the help of numerous R packages, everything can be accomplished in R from data processing, estimation, post-estimation, to final manuscript writing. This is a huge advantage over stata. 

## 2. Post-estimation Commands
The marginal effect estimation in this package is much faster than Stata's fmlogit package. In this package user can specify which variable(s), and what effect to be calculated. This results in a huge gain in running time for the post-estimation commands. 

Also, this package allows the calculation of marginal and discrete effect at the individual level, which is crucial to detect possible heterogeneities among observations. 

## 3. Estimation Flexibility
This package allows factor variable inputs, and automatically transform it into dummy variables. This is not (explicitly) allowed in Stata. 

## 4. Extensions
This package also allows the user to calculate and infer from a vector of "willingness to pay" measures. The user can calculate the mean and the standard error of the effect on overall WTP increase, which is implicitly derived from the choice dynamics. 

# How the estimator works?
The estimator used here is an extension of Papke and Wooldridge(1996)'s paper, in which they proposed a quasi-maximum likelihood(QMLE) estimator for fractional response variables. As their approach applies to binary response variables, here we expand it to a multinomial response variables with fractional structure. 

The basic step of the estimator is the following: 
## Step 1 Construct the multinomial logit likelihood, which takes the form: 
$$ ln(L_i) = \sum_{j=1}^J y_{ij} ln(G(x_i\beta_j))$$
where G(.) is the multinomial logit function with form: 
$$G(z_j) = \frac{exp(z_j)}{\sum_{k=1}^J exp(z_k)}$$
for which $i \in {1,...,N}$ denotes N individual observation, $j \in {1,...,J}$ denotes J available choices, and $k \in {1,...,K}$ denotes K explanatory variables. 
## Step 2 Maximize the sum of the log likelihood function
Generally speaking, R is not the most efficient scientific computing machine that exists, and that is the tradeoff we have to face. Here the program offers several maximization methods provided in the *maxLik* package. The recommended algorithm is either conjugate gradients (CG), or Berndt-Hall-Hall-Hausman (BHHH). For a large dataset it may take a while (running for one hour is entirely possible, so don't terminate the program just yet).
## step 3 Calculate robust standard error
Here the program follows Papke & Wooldridge (1996)'s paper, and construct the robust standard error estimator for the parameters. The program also offers a simple z-test for parameters based on the standard error. 

# How post-estimation commands works?

Interpreting marginal and discrete effects for limited dependent variables can be tricky, and this is especially true for multinomial logit models. The coefficients obtained in the regression model represents the logit-transformed odds ratio for that specific choice against the baseline choice. This is not intuitive at all in terms of what are the actual effects on that specific choice. The bottom line is, the coefficients and standard errors obtained in the original models are not the basis for evaluating hypotheses. 

# Marginal Effects

Instead, researchers need to compute what's called the "marginal effect", as we usually do in linear models. The marginal effect of the multinomial models actually has a very distinctive form: 
 
$$ME_{jk} = \frac{\partial p_j}{x_k} = p_j(\beta_{kj} - \bar{\beta}_i)$$
 
where $p_j$ is an 1*N vector of predicted probabilities for choice j, and $\bar{\beta}_i = \sum_{m=1}^J \beta_{km} p_m$ is the probability weighted average of $\beta{km}$. This shows that the marginal effects among different individuals are actually different given different predicted probabilities of choice j. 
 
Typically, two types of summary measures are used to illustrate the global average marginal effects. The first one is called marginal effects at the mean (MEM) in the code. In algebric form, this is represented as
 
$$MEM_{jk} = \bar{p}_j(\beta_{kj} - \bar{\beta}_i)$$
 
where $\bar{p}_j$ is the predicted value of choice j at the mean of all X covariates. Centering observations around the mean simplifies the calculation, however it ignores the potential heterogeneity in marginal effects, especially at the extreme values. 
 
Another measure is called average marginal effects (AME). This can be written as:
 
$$AME_{jk} = \frac{1}{N}\sum_{i=1}^N p_j(\beta_{kj} - \bar{\beta}_i)$$
 
However, according to Greene(2003), there is no agreement as to which one is prefered. A more inclusive approach will be to plot the marginal effect of interest across all individuals. This is not provided in the function, but can certainly be implemented in a straightforward way in R. 
 
## Discrete Effects 
 
Discrete effect is a little bit different from marginal effects. Instead of calculating the slope of the 
coefficients, discrete effect considers the impact of a discrete change in one covariates on the predicted outcome variables. This is especially useful for dummy variables, where calculating marginal effects does not make much sense. 

The discrete effect has a straight-forward form. Consider a discrete change of a dummy variable k from 0 to 1. This is just
$$D_j = Pr(y=j|\mathbf{x}_{x_k=1}) - Pr(y=j|\mathbf{x}_{x_k=0})$$
The change in predicted value by setting $x_k=1$ and $x_k=0$. 

Similar to the marginal effect case, we can calculate discrete effects at the mean (DEM) by predicting the outcome when all other covariates at the mean, or average discrete effects, which averages the predicted difference across all observations. 

## Standard Errors

Here we adopt the Krinsky-Robb method to compute standard errors for marginal and discrete effects. As oppose to the delta method commonly used in other programs such as Stata, Krinsky-Robb is a simulation-based method. The idea of Krinsky-Robb is that, to calculate the variance of a function $Var(f(\mathbf{\beta}))$, we do the following step:

For i in 1:N, where N is a very large number, 

1) Sample from the known distribution of $\mathbf{\beta}$

2) For each of the sample, calculate $f(\mathbf{\beta})$

3) Take the empirical variance of $f(\mathbf{\beta})$. 

And after sufficiently large sample size, the empirical variance converges to the theoretical variance. 

Note that the Krinsky-Robb sampling is done by setting all other covariates to the mean, i.e., the MEM and the DEM approach. 

# Practical Concerns
## Optimization Method
This function calls *maxLik()* in package *maxLik* to maximize the quasi-likelihood function. The *maxLik* function is a wrapper which provides several different maximization methods, including most *optim()* methods in the base package, as well as other useful methods such as BHHH(Berndt-Hall-Hall-Hausman). The choice of optimization method can create vastly different parameter estimates. Here it is recommended that either conjugate gradients(CG), or Berndt-Hall-Hall-Hausman(BHHH) to ensure convergence. In limited testing scenarios, BHHH typically has the best performace in terms of convergence for large datasets, while CG is faster in computation speed for smaller, easy to converge datasets.

## Robust Standard Error
It is worth noting that the robust standard error created in this function is consistently lower than that created in Stata's fmlogit package, typically by about 20\%. However, the robust SE here is a consistent estimator following Pakpe and Wooldridge(1996)'s $\hat{A_j}^{-1}\hat{B_j}\hat{A_j}^{-1}$ estimator, so it is recommended that the number should be used with causion. 

#References
Papke, L. E. and Wooldridge, J. M. (1996), Econometric methods for fractional response variables with an application to 401(k) plan participation rates. J. Appl. Econ., 11: 619-632.

Wulff, Jesper N. "Interpreting Results From the Multinomial Logit Model Demonstrated by Foreign Market Entry." Organizational Research Methods (2014): 1094428114560024.






