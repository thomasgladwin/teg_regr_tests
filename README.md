# teg_regr_tests
Function in R for tests of linear constraints on regression coefficients.

The function (a practice project after reading Bingham & Fry, but a handy helper function potentially) tests linear hypotheses on regression coefficients beta_i of the form Constraints * betas = constants.

The function requires matrices X and y as usual (an intercept column must be added to X explicitly if required), plus a Hypothesis: a list of a Constraints matrix and a constants vector that defines the constraints to test. Constraints can, e.g., test the effect of removing a predictor or equalizing two predictors etc. The function runs the F-test for the hypothesis (significance meaning the hypothesis is rejected). If the constraints are equivalent to testing a nested model, the AIC is provided for the full model and constrained model.

The file includes example code for various constraints.
