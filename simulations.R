# Sample size
n = 1000
# Length of function support
j = 1000
# Number of basis functions
nbasis = 25
# Number of PLS components
npls = 5

set.seed (12345)
# Generate the data
sim_dat = data_generation(n = n, j = j)
X = sim_dat$X
Y = sim_dat$Y
beta_true = sim_dat$Beta

# Model fits
simpls_model = reg_fun(Y = Y, X = X, nbasis = nbasis, npls = npls, rangeval = c(0,1), method = "simpls")
nipals_model = reg_fun(Y = Y, X = X, nbasis = nbasis, npls = npls, rangeval = c(0,1), method = "nipals")
bidiag1_model = reg_fun(Y = Y, X = X, nbasis = nbasis, npls = npls, rangeval = c(0,1), method = "bidiag1")
bidiag2_model = reg_fun(Y = Y, X = X, nbasis = nbasis, npls = npls, rangeval = c(0,1), method = "bidiag2")

# MISE
mean((beta_true - simpls_model$bhat)^2) # 0.708649
mean((beta_true - nipals_model$bhat)^2) # 0.7086455
mean((beta_true - bidiag1_model$bhat)^2) # 0.7086455
mean((beta_true - bidiag2_model$bhat)^2) # 0.7086455

# MSE
mean((Y - simpls_model$fits)^2) # 1.064305
mean((Y - nipals_model$fits)^2) # 1.007794
mean((Y - bidiag1_model$fits)^2) # 1.007794
mean((Y - bidiag2_model$fits)^2) # 1.007794
