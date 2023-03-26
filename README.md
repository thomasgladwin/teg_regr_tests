# teg_regr_tests
Function in R and Python for tests of linear constraints on regression coefficients.

The function (a practice project after reading Bingham & Fry, but a handy helper function potentially) tests linear hypotheses on regression coefficients beta_i of the form Constraints * betas = constants.

The function requires predictor and outcome matrices X and y (an intercept column must be added to X explicitly if required), plus a Hypothesis: a list of a Constraints matrix and a constants vector that defines the constraints to test. Constraints can, e.g., test the effect of removing a predictor or equalizing two predictors etc. The function runs the F-test for the hypothesis (significance meaning the hypothesis is rejected) and provides the AIC's for the full model and constrained model.

# For R:

To run a typical test aimed at the predictors, where the intercept isn't of interest, center the variables and run without an explicit intercept (otherwise you'd have to specify you're testing against the constraint that all parameters are 0 except the intercept).

Data from a data.frame will need to be converted to matrices, e.g., via data.matrix().

The script includes example code for various constraints. E.g.,

H <- c()

H$Constraints = matrix(c(1, 0, 0, 0), nrow=1)

H$Constraints = rbind(H$Constraints, matrix(c(0, 0, 0, 1), nrow=1))

H$constants = matrix(c(0, 0), ncol=1)

O <- teg_regr_tests(X, y, H)

tests against the hypothesis that predictors X1 and X4 have a weight of 0; or, equivalently, it tells you whether adding the set of predictors {X1 and X4} to the model without them results in a significant increase in explained variance. The function (unless told to suppress output) will provide the usual regression output from lm() and the tests of the linear hypothesis:

F-test of linear constraint: F(2, 395) = 0.639, p = 0.529  (p < .05 rejects constrained model.)

Constrained model AIC =  -633 

Full model AIC =  -630 

Difference =  -2.71  (negative supports constrained model).

So, in this case, the reduced model is better than the full model in terms of AIC and the F-test agrees, as removing the predictors does not result in a significant increase in unexplained variance.

# For Python
The usage is illustrated, with simulated data, in test_teg_regression.py.

By default, do not add an explicit intercept. The intercept will be appended as an additional predictor.

The Constraints argument is a dictionary. The example below shows the Constraints setup to set a specific predictor-coefficient to 0.

pred_to_test = 1

Constraints = {}

Constraints['coefficients'] = np.array([0 for a in range(X.shape[1])]).reshape(1, X.shape[1])

Constraints['coefficients'][0][pred_to_test] = 1

Constraints['constants'] = np.array([0])

O = teg_regression.run_regression(X, y, Constraints)


[![DOI](https://zenodo.org/badge/376601604.svg)](https://zenodo.org/badge/latestdoi/376601604)

