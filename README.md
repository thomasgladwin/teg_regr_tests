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

tests against the hypothesis that predictors X1 and X4 have a weight of 0; or, equivalently, it tests whether adding the set of predictors X1 and X4 to the model without them results in a significant increase in explained variance. The function (unless told to suppress output) will provide the usual regression output from lm(), with the final two lines providing the AIC comparison (for nested models only) and F-test for the hypothesis:

Reduced model AIC =  433 , versus full-fit AIC =  424 

F-test of linear constraint: F(2, 75) = 6.47, p = 0.00255
