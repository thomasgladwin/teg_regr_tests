teg_regr_betas <- function(X, y) {
  C <- t(X) %*% X
  betas <- solve(C) %*% t(X) %*% y
  return(betas)
}

teg_logL <- function(y_err, p) {
  # https://www.le.ac.uk/users/dsgp1/COURSES/MATHSTAT/13mlreg.pdf
  # https://www.quantstart.com/articles/Maximum-Likelihood-Estimation-for-Linear-Regression/
  N <- length(y_err)
  s2 <- sum(y_err ** 2)/N
  s <- sqrt(s2)
  logL <- -N * log(s) - (N/2)*log(2*pi) - (1/(2*s2)) * sum(y_err ** 2)
}

teg_AIC <- function(logL, p) {
  return(-2*logL + 2 * (p + 1))
}

teg_regr <- function(X, y, betas = c(), hyp_B = c()) {
  O = c()
  if (is.null(betas)) {
    O$betas <- teg_regr_betas(X, y)
    p <- dim(X)[2]
  } else {
    O$betas = betas
    p <- dim(X)[2] - dim(hyp_B)[1]
  }
  O$y_pred <- X %*% O$betas
  O$y_err <- y - O$y_pred
  O$logL <- teg_logL(O$y_err, p)
  O$AIC <- teg_AIC(O$logL, p)
  return(O)
}

teg_regr_tests <- function(X, y, H, verbose0 = 1, names0 = c()) {
  # Function to test linear hypotheses on regression coefficients.
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
        cat(c('X', n, ' = ', names0[n], '\n'), sep="")
      }
    }
    cat('Constraints:\n')
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
  hyp_c = H$constants
  C = t(X) %*% X
  beta_constr <- beta_full - solve(C) %*% t(hyp_B) %*% solve(hyp_B %*% solve(C) %*% t(hyp_B)) %*% (hyp_B %*% beta_full - hyp_c)
  SSH <- t(beta_full - beta_constr) %*% C %*% (beta_full - beta_constr)
  nConstr = length(hyp_c)
  SSE <- sum(fit_full$residuals ** 2)
  F <- (SSH/nConstr) / (SSE/(N - nPred))
  p <- 1 - pf(F, nConstr, (N - nPred))

  if (verbose0 == 1) {
    if (!is.null(names0)) {
      for (n in 1:nPred) {
        cat(c('\nb_constrained', n, ' = ', signif(beta_constr[n], digits = nDigits)))
      }
      cat(paste("\n\n"))
    }
    cat(paste('F-test of linear constraint: F', '(', nConstr, ', ', (N - nPred), ') = ', signif(F, digits = nDigits), ', p = ', signif(p, digits = nDigits), sep=""), " (p < .05 rejects constrained model.)")
  }

  if (verbose0 == 1) {
    O_constr <- teg_regr(X, y, beta_constr, hyp_B)
    O_full <- teg_regr(X, y)
    cat("\nConstrained model AIC = ", signif(O_constr$AIC, digits = 3), "")
    cat("\nFull model AIC = ", signif(O_full$AIC, digits = 3), "")
    cat("\nDifference = ", signif(O_constr$AIC - O_full$AIC, digits=3), " (negative supports constrained model).")
  }
  
  Output$beta_constrained <- beta_constr
  Output$F <- F
  Output$df1 <- nConstr
  Output$df2 <- (N - nPred)
  Output$p <- p
}

# Example code

# Create test data
N <- 400
nPred <- 4 # Intercept must be added to these as an explicit column in X
names0 <- paste(rep('v', nPred), 1:nPred, sep="")
names0 <- c(names0, 'Intercept')
X <- matrix(rnorm(N * nPred), ncol = nPred)
X <- cbind(X, rep(1,N)) # Add intercept explicitly
print(dim(X))
beta_true = 1 + 1 * matrix(c(3 * (1:nPred), 0), ncol=1) # True coefficients increase by three, and set Intercept to 100
beta_true[2] = 0 # Set Beta_2 to zero for testing
beta_true[3] = 0 # Set Beta_2 to zero for testing
beta_true[1] = beta_true[4] # Set Beta_3 == Beta_4 for testing
print(beta_true)
e <- 0.1 * rnorm(N)
y <- X %*% beta_true + e

# Set to 0
H <- c()
H$Constraints = matrix(c(0, 1, 0, 0, 0), nrow=1)
H$Constraints = rbind(H$Constraints, matrix(c(0, 0, 1, 0, 0), nrow=1))
H$constants = matrix(c(0, 0), ncol=1)
O <- teg_regr_tests(X, y, H, 1, names0)

# To run basic regression, where H = c()
O <- teg_regr_tests(X, y, c(), 1, names0)

# Set X1 == X2 and X3 to 5
H <- c()
H$Constraints = matrix(c(1, -1, 0, 0, 0), nrow=1)
H$Constraints = rbind(H$Constraints, matrix(c(0, 0, 1, 0, 0), nrow=1))
H$constants = matrix(c(0, 5), ncol=1)
O <- teg_regr_tests(X, y, H, 1, names0)

# Set X2 == X3 == X4
H <- c()
H$Constraints = matrix(c(-1, 0, 1, 0, 0), nrow=1)
H$Constraints = rbind(H$Constraints, matrix(c(0, 0, -1, 1, 0), nrow=1))
H$constants = matrix(c(0, 0), ncol=1)
O <- teg_regr_tests(X, y, H, 1, names0)

# Check consistency with R functions
O <- teg_regr(X, y); print(O$logL); print(O$AIC)
fit0 <- lm(y ~ X - 1); print(logLik(fit0)); print(AIC(fit0))
