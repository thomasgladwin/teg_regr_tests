teg_regr_tests <- function(X, y, H, verbose0 = 1, names0 = c()) {
  # Function to run a regression and test linear hypotheses.
  # Method from Bingham & Fry.
  #
  # X and y must be matrices in the correct orientation.
  # X = n x k, y = n x 1.
  # Names is a list of labels for predictors in X.
  # The intercept must be explicitly included as a column in X.
  # (This is for consistency, as a linear constraint could involve the intercept.)
  #
  # The hypothesis representation H contains fields:
  # Constraints (a matrix, uppercase)
  # and
  # constants (a vector, lowercase)
  # that define linear constraints on coefficients, such that:
  # Constraints * fitted_coefficients = constants
  # So Constraints = h x k, constants = h x 1.
  # Each row of B and c represents a different hypothesis.
  # If H is c(), only the ordinary regression is performed.

  # Prep
  Output <- c()
  N <- dim(X)[1]
  nPred <- dim(X)[2]
  
  # Show label mapping
  if (verbose0 == 1) {
    if (!is.null(names0)) {
      for (n in 1:nPred) {
        cat(c('X', n, ' = ', names0[n], '\n'))
      }
    }
    print(H)
  }
  nDigits <- 3
  
  # Fit the full model
  fit_full <- lm(y ~ X - 1)
  Output$fit_full <- fit_full
  beta_full <- matrix(fit_full$coefficients, ncol=1)
  
  # Print results
  if (verbose0 == 1) {
    cat("# # # Full model:\n")
    print(summary(fit_full))
    print(summary.aov(fit_full))
    cat("AIC = ", signif(AIC(fit_full), digits=nDigits), "\n\n")
  }

  # If not testing a linear constraint, just return the 
  # regression fit.
  if (is.null(H)) {
    return(Output)
  }

  # Linear constraint test
  hyp_B = H$Constraints
  hyp_B = cbind(hyp_B, rep(0,dim(hyp_B)[1]))
  hyp_c = H$constants
  C = t(X) %*% X
  beta_constr <- beta_full - solve(C) %*% t(hyp_B) %*% solve(hyp_B %*% solve(C) %*% t(hyp_B)) %*% (hyp_B %*% beta_full - hyp_c)
  SSH <- t(beta_full - beta_constr) %*% C %*% (beta_full - beta_constr)
  nConstr = length(hyp_c)
  SSE <- sum(fit_full$residuals ** 2)
  F <- (SSH/nConstr) / (SSE/(N - nPred))
  p <- 1 - pf(F, nConstr, (N - nPred))

  # If linear constraint consists of exclusions:
  H_Constraint_vals <- unique(as.vector(H$Constraints))
  H_constants_vals <- unique(H$constants)
  if (all(H_Constraint_vals %in% c(0, 1)) && all(H_constants_vals == 0)) {
    X_reduced = X[, -which(colSums(H$Constraints) == 1)]
    fit_reduced <- lm(y ~ X_reduced - 1)
    Output$fit_reduced <- fit_reduced
    if (verbose0 == 1) {
      cat("# # # Reduced model:\n")
      
      if (!is.null(names0)) {
        retained <- 1:nPred
        retained <- retained[which(colSums(H$Constraints) == 0)]
        m = 1
        for (n in retained) {
          cat(c('X', m, ' = ', names0[n], '\n'))
          m = m + 1
        }
      }

      print(summary(fit_reduced))
      print(summary.aov(fit_reduced))
      cat("Reduced model AIC = ", signif(AIC(fit_reduced), digits=nDigits), ", versus full-fit AIC = ", signif(AIC(fit_full), digits=nDigits), "\n\n")
    }    
  }
  
  if (verbose0 == 1) {
    if (!is.null(names0)) {
      for (n in 1:nPred) {
        cat(c('b_constrained', n, ' = ', signif(beta_constr[n], digits = nDigits), '\n'))
      }
    }
    cat(paste('\nF-test of linear constraint: F', '(', nConstr, ', ', (N - nPred), ') = ', signif(F, digits = nDigits), ', p = ', signif(p, digits = nDigits), sep=""))
  }

  Output$beta_constrained <- beta_constr
  Output$F <- F
  Output$df1 <- nConstr
  Output$df2 <- (N - nPred)
  Output$p <- p
  
}

# Example code

# Create test data
N <- 280
nPred <- 4 # Intercept must be added to these as an explicit column in X
names0 <- paste(rep('v', nPred - 1), 1:(nPred - 1), sep="")
names0 <- c(names0, 'Intercept')
X <- matrix(rnorm(N * nPred), ncol = nPred)
X <- cbind(X, rep(1,N)) # Add intercept explicitly
print(dim(X))
beta_true = matrix(c(3 * (1:nPred), 100), ncol=1) # True coefficients increase by three, and set Intercept to 100
beta_true[2] = 0 # Set Beta_2 to zero for testing
beta_true[3] = beta_true[4] # Set Beta_3 == Beta_4 for testing
print(beta_true)
e <- 0.1 * rnorm(N)
y <- X %*% beta_true + e

# Run basic regression, where H = c()
O <- teg_regr_tests(X, y, c(), 1, names0)

# Linear constraint tests
# Test whether the data contradict the specified hypothesis involving coefficients,
# e.g., that Beta_1 == 0 or Beta_3 == Beta_4.

## Examples where one or more weights are set to 0, removing the predictor(s).
# Example: set X1 to 0
H$Constraints = matrix(c(1, 0, 0, 0), nrow=1)
H$constants = matrix(c(0), ncol=1)
O <- teg_regr_tests(X, y, H, 1, names0)
# Example: set X2 to 0
H$Constraints = matrix(c(0, 1, 0, 0), nrow=1)
H$constants = matrix(c(0), ncol=1)
O <- teg_regr_tests(X, y, H, 1, names0)
# Example: Set X2 and X4 to 0
H$Constraints = matrix(c(1, 0, 0, 0), nrow=1)
H$Constraints = rbind(H$Constraints, matrix(c(0, 0, 0, 1), nrow=1))
H$constants = matrix(c(0, 0), ncol=1)
O <- teg_regr_tests(X, y, H, 1, names0)

## Examples of equating coefficients
# Example: set X1 == X4
H$Constraints = matrix(c(1, 0, 0, -1), nrow=1)
H$constants = matrix(c(0), ncol=1)
O <- teg_regr_tests(X, y, H, 1, names0)
# Example: set X3 == X4
H$Constraints = matrix(c(0, 0, -1, 1), nrow=1)
H$constants = matrix(c(0), ncol=1)
O <- teg_regr_tests(X, y, H, 1, names0)

# Example of multi-constraints: set X1 == X2 and X3 to 5
H$Constraints = matrix(c(1, -1, 0, 0), nrow=1)
H$Constraints = rbind(H$Constraints, matrix(c(0, 0, 1, 0), nrow=1))
H$constants = matrix(c(0, 5), ncol=1)
O <- teg_regr_tests(X, y, H, 1, names0)
