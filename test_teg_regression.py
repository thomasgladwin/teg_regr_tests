import numpy as np
import teg_regression

nObs = 300
nPred = 5
fix_coeffs = {0: 1, 3: 2}
#fix_coeffs = {}  # Set to empty dict to simulate the null model.
X, y = teg_regression.sim_data(nObs, nPred, fix_coeffs=fix_coeffs, fix_intercept=20)
O = teg_regression.run_regression(X, y)

# Constraints example: Set two given predictors to 0
pred_to_test = [1, 2]
Constraints = {}
Constraints['coefficients'] = np.array([[0 for a in range(X.shape[1])] for newrow in range(2)]).reshape(2, X.shape[1])
Constraints['coefficients'][0][pred_to_test[0]] = 1
Constraints['coefficients'][1][pred_to_test[1]] = 1
Constraints['constants'] = np.array([0, 0])
O = teg_regression.run_regression(X, y, Constraints)
