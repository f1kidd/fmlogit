---
title: "The fmlogit Package: A Light Document"
author: "Xinde James Ji" 
date: "Oct 10, 2016"
output: pdf_document
---

This document provides an overview of the fmlogit package in R. Updates will be published at [my github site](https://github.com/f1kidd/fmlogit). Any suggestions or concerns are welcome. 

# What is the fractional multinomial logit model?
Fractional multinomial logit models estimate fractional responses by modelling the dependent variables as fractions using multinomial logits. It is the preferred model when the true data generation process is indeed fractions of multiple choices. Fractional responses arise naturally in various settings. For example, a municipality allocates its budget across multiple departments, and we are interested in the proportion of the budget that each department receives. Or, there are multiple candidates in a presendential election, and we are interested in explaining the percentage of support for each candidate in each state. 

The model is distinct in that: 1) each of the responses lies between 0 and 1, and 2) the share of all responses adds up to one. The fmlogit model uses these two distinct factors, and models them explicitly. 

# How to install fmlogit
Type the following code into your R console:
```R
require(devtools)
install_github("f1kidd/fmlogit")
library(fmlogit)
```

# Why do we need fmlogit in R? 
Don't we already have an fmlogit module in Stata? Yes, and you are very welcome to [check that out](http://maartenbuis.nl/software/fmlogit.html) if you can afford a Stata license. 

However, this package offers several advantages over Stata's fmlogit module, namely:
### 1. Integration with the R Platform
Implementating the model in R offers the opportunity to integrate the whole empirical process within a free, open-source platform. With the help of numerous R packages, everything can be accomplished in a single environment including data processing, estimation, post-estimation, and final manuscript writing. This is a huge advantage over stata. 

### 2. Post-estimation improvements
The marginal effect estimation in this package is much faster than Stata's fmlogit package. In this package user can specify which variable(s), and what type of partial effect to be calculated. This results in a huge gain in running time for the post-estimation commands. 

Also, this package allows hypothesis testing for marginal and discrete effects while Stata does not. The standard error is calculated via Krinsky-Robb method, which allows empirical hypothesis testing without knowing the underlying distribution of the effects. 

### 3. Estimation flexibility
This package allows factor variable inputs, and automatically transform it into dummy variables. This is not (explicitly) allowed in Stata. 

### 4. Extensions
This package also allows the user to easily calculate and infer the "average aggregate partial effect" given a user-specified weight scheme. This is done through linearly aggregating the attribute of each choice (e.g., expected profit/utility of each choice) with the calculated APE. 

# How does the estimator work?
The estimator used here is an extension of that used in [Papke and Wooldridge (1996)](http://onlinelibrary.wiley.com.ezproxy.lib.utexas.edu/doi/10.1002/(SICI)1099-1255(199611)11:6%3C619::AID-JAE418%3E3.0.CO;2-1/abstract). There, they proposed a quasi-maximum likelihood(QMLE) estimator for fractional response variables. As their approach applies to binary response variables, here we expand it to a multinomial response variables with fractional structure. 

The steps involved in calculating the estimator are as follows: 
## Step 1. Construct the multinomial logit likelihood
This step is straightforward. A simple multinomial logit transformation will do the job. For detailed derivations and formula, please see the technical document [here](https://github.com/f1kidd/fmlogit/blob/master/Documentation/fmlogit_docs.pdf) where I explain the econometric steps in detail.  
## Step 2. Maximize the sum of the log likelihood function
Generally, R is not the most efficient scientific computing machine that exists, and that is the tradeoff we have to face. Here, the program offers several maximization methods provided in the *maxLik* package. The recommended algorithm is either conjugate gradients (CG), or Berndt-Hall-Hall-Hausman (BHHH). For a large dataset it may take a while (running for one hour is entirely possible, so don't terminate the program prematurely).
## step 3. Calculate robust standard error
Here the program follows Papke & Wooldridge (1996), and construcst the robust standard error estimator for the parameters. The program also offers a simple z-test for parameters based on the standard error. 

# How do the post-estimation commands work?
Calculating partial effects for limited dependent variables can be tricky, and this is especially true for multinomial logit models. The coefficients obtained in the regression model represent the logit-transformed odds ratio for that specific choice against the baseline choice. This is not intuitive at all in terms of actual effects on that specific choice. The bottom line is, the coefficients and standard errors obtained in the original models are not the basis for evaluating hypotheses. 

##  Marginal and discrete effects
Instead, researchers need to compute what are called the "partial effects", as we usually do in linear models. However, the partial effect in logit-type models is tricky because the effects are heterogeneous across different observations. In other word, each unique observation have a different set of partial effects.
 
We provide two types of partial effects: marginal and discrete. The marginal effect represents how a unit change in one variable k changes the value in choice j, i.e., $\frac{\partial x_k}{\partial y_j}$. The discrete effect represents how a discrete change in variable k, usually from the minimum to the maximum, changes the value in choice j, i.e., $\hat{y}_{j,x_k=1}-\hat{y}_{j,x_k=0}$. 

Typically, two types of aggregation measures are used to illustrate the global APE: one is the partial effects at the mean (PEM), which is the partial effect of variable k when every other variables are set at their mean ; and the other is partial effect of the average (PEA), which is the average of partial effect for all observations. We allow both of the two options to be specified. 

A more inclusive approach will be to plot the marginal effect of interest across all individuals. This is not provided in the function, but can certainly be implemented in future developments. Another possibility will be to calculate the so-called locally averaged treatment effect (LATE), where the effect of interest will be centered around a certain range of values.  

## Standard Errors for APEs
Here we adopt the simulation-based Krinsky-Robb method to compute standard errors for marginal and discrete effects as opposed to the empirical delta method used in Stata. These two methods should be asymptotically equivalent. However, using Krinsky-Robb allow us to perform hypothesis testing on the effects, while Delta method cannot accomplish that.

Hypothesis testing is done using the standard normal z-test by treating the APE estimates as normally distributed. The approach is very simple: say we test $H_0: D_j=0$. We just need to compare 0 with our N draws, and see if it falls out of the 95% mass. This is a major advantage we provide here compared to Stata's fmlogit module. 

# Practical Concerns
One of the concerns for the package is the computation speed of the estimation process. The maximization process can take somewhere from 20 seconds to 1 hour, depending on how large the dataset is. This is certainly a limitation. This is the inherent drawback of R's computation speed, and I can do nothing about that. 

However, the loss in estimation will certainly be compensated in the post-estimation process. Stata's dfmlogit command is very slow (takes somewhere between 5-60 minutes), while here the effects calculation takes seconds to complete. 

# References
Papke, L. E. and Wooldridge, J. M. (1996), Econometric methods for fractional response variables with an application to 401(k) plan participation rates. J. Appl. Econ., 11: 619-632.

Wulff, Jesper N. "Interpreting Results From the Multinomial Logit Model Demonstrated by Foreign Market Entry." Organizational Research Methods (2014): 1094428114560024.

Mullahy, J., 2015. Multivariate fractional regression estimation of econometric share models. Journal of Econometric Methods 4(1), 71-100.






