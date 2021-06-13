# teg_regr_tests
Function in R for tests of linear constraints on regression coefficients.

The function (a practice project after reading Bingham & Fry, but a handy helper function potentially) tests linear hypotheses on regression coefficients beta_i of the form Constraints * betas = constants.

The function requires predictor and outcome matrices X and y (an intercept column must be added to X explicitly if required), plus a Hypothesis: a list of a Constraints matrix and a constants vector that defines the constraints to test. Constraints can, e.g., test the effect of removing a predictor or equalizing two predictors etc. The function runs the F-test for the hypothesis (significance meaning the hypothesis is rejected). If the constraints are equivalent to testing a nested model, the AIC is provided for the full model and constrained model.

Data from a data.frame will need to be converted to matrices, e.g., via data.matrix().

The script includes example code for various constraints. E.g.,

H$Constraints = matrix(c(1, 0, 0, 0), nrow=1)
H$Constraints = rbind(H$Constraints, matrix(c(0, 0, 0, 1), nrow=1))
H$constants = matrix(c(0, 0), ncol=1)
O <- teg_regr_tests(X, y, H)

tests against the hypothesis that predictors X1 and X4 have a weight of 0; or, equivalently, it tests whether adding the set of predictors X1 and X4 to a model without them results in a significant increase in explained variance. The function (unless told to suppress output) will provide output like this, with the final two lines providing the AIC comparison (for nested models only) and F-test:

$Constraints
     [,1] [,2] [,3] [,4]
[1,]    1    0    0    0
[2,]    0    0    0    1

$constants
     [,1]
[1,]    0
[2,]    0

# # # Full model:

Call:
lm(formula = y ~ X - 1)

Residuals:
    Min      1Q  Median      3Q     Max 
-8.4956 -2.0026 -0.0373  2.3558  6.0999 

Coefficients:
   Estimate Std. Error t value Pr(>|t|)    
X1  0.49579    0.34448   1.439  0.15425    
X2  0.05733    0.40657   0.141  0.88824    
X3  0.46371    0.32892   1.410  0.16274    
X4  1.25658    0.40567   3.098  0.00275 ** 
X5  9.68838    0.37376  25.921  < 2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 3.289 on 75 degrees of freedom
Multiple R-squared:  0.9075,	Adjusted R-squared:  0.9014 
F-statistic: 147.2 on 5 and 75 DF,  p-value: < 2.2e-16

          Df Sum Sq Mean Sq F value Pr(>F)    
X          5   7961  1592.3   147.2 <2e-16 ***
Residuals 75    811    10.8                   
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
AIC =  424 

# # # Reduced model:

Call:
lm(formula = y ~ X_reduced - 1)

Residuals:
    Min      1Q  Median      3Q     Max 
-8.1720 -2.3229 -0.0206  2.7177  6.4043 

Coefficients:
           Estimate Std. Error t value Pr(>|t|)    
X_reduced1   0.2102     0.4304   0.488    0.627    
X_reduced2   0.4833     0.3514   1.375    0.173    
X_reduced3   9.9010     0.3936  25.157   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 3.515 on 77 degrees of freedom
Multiple R-squared:  0.8915,	Adjusted R-squared:  0.8873 
F-statistic:   211 on 3 and 77 DF,  p-value: < 2.2e-16

          Df Sum Sq Mean Sq F value Pr(>F)    
X_reduced  3   7821  2607.1     211 <2e-16 ***
Residuals 77    951    12.4                   
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Reduced model AIC =  433 , versus full-fit AIC =  424 

F-test of linear constraint: F(2, 75) = 6.47, p = 0.00255
