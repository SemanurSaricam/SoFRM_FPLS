#___________________________
#__________Packages_________
#___________________________
library('fda')
library('expm')
library('pracma')
library('matrixStats')
library('fda.usc')
#___________________________
#___________________________

#________________________________________________________________________________________________________
#______________________________________________FPLS_projection_________________________________________
#________________________________________________________________________________________________________

get_pls_mat = function(data, nbasis, rangeval){
  n = dim(data)[1]
  p = dim(data)[2]
  dimnames(data)=list(as.character(1:n), as.character(1:p))
  grid_points = seq(rangeval[1], rangeval[2], length.out = p)
  bs_basis = create.bspline.basis(rangeval, nbasis = nbasis)
  evalbase = eval.basis(grid_points, bs_basis)
  innp_mat = inprod(bs_basis, bs_basis)
  sqrt_innp_mat = sqrtm(innp_mat)$B
  fdobj = fdPar(bs_basis, int2Lfd(2), lambda=0)
  pcaobj = smooth.basisPar(grid_points, t(data), bs_basis, Lfdobj=NULL, lambda=0)$fd
  data_mat = t(pcaobj$coefs) %*% sqrt_innp_mat
 
  return(list(data_mat = data_mat, sqrt_innp_mat = sqrt_innp_mat,
              eval_base = evalbase))
}
#________________________________________________________________________________________________________
#________________________________________________________________________________________________________


#___________________________________________________________________________________________________________________
#________________________________________________Regressin_functions________________________________________________
#___________________________________________________________________________________________________________________
reg_fun = function(Y, X, nbasis, npls, rangeval, method = c("simpls", "nipals", "bidiag1", "bidiag2")){
  Y = as.matrix(Y)
  n = dim(Y)[1]
 
  dmat = get_pls_mat(data = X, nbasis = nbasis, rangeval = rangeval)
  matX = dmat$data_mat
 
 
  if(method == "nipals"){
    model_pls = nipalsN(X = matX, Y = Y, a = npls)
  }else if(method == "simpls"){
    model_pls = simpls(X = matX, Y = Y, a = npls)
  }else if(method == "bidiag1"){
    model_pls = bidiag1(X = matX, Y = Y, a = npls)
  }else if(method == "bidiag2"){
    model_pls = bidiag2(X = matX, Y = Y, a = npls)
  }
 
  coef_pls = model_pls$beta
 
  fits = matX %*% coef_pls
 
  bhat = dmat$eval_base %*% (solve(dmat$sqrt_innp_mat) %*% coef_pls)
 
  return(list(fits = fits, bhat = bhat))
}
#___________________________________________________________________________________________________________________
#___________________________________________________________________________________________________________________



#_____________________________________________________________________________________________
#___________________________________NIPALS_with_normalization_________________________________
#_____________________________________________________________________________________________
nipalsN = function(X,Y,a){
  W = NULL
  T = NULL
  P = NULL
  q = NULL
  for(i in 1:a){
    w = t(X) %*% Y
    w = w / norm(w, type = "2")
    W = cbind(W, w)
    t = X %*% w
    t = t / norm(t, type = "2")
    T = cbind(T, t)
    P = cbind(P, t(X) %*% t)
    X = X - t %*% t(P[,i])
    q = c(q, t(Y) %*% t)
    Y = Y - q[i] * t
  }
  beta = rowCumsums(t(t(mrdivide(W, triu(t(P) %*% W))) * q))[,a]
  return(list(beta=beta))
}
#_____________________________________________________________________________________________
#_____________________________________________________________________________________________


#_____________________________________________________________________________________________
#_________________________________Bidiag2_with_reorthogonalization____________________________
#_____________________________________________________________________________________________
bidiag2 = function(X, Y, a){
  B = matrix(0, nrow = a, ncol = 2)
  w = t(X) %*% Y
  w = w / norm(w, type = "2")
  W = w
  t = X %*% w
  rho = norm(t, type = "2")
  t = t / rho
  T = t
  B[1,1] = rho
  d = w / rho
  beta = d %*% (t(t) %*% Y)
  for(i in 2:a){
    w = t(X) %*% t - rho * w
    w = w - W %*% (t(W) %*% w)
    theta = norm(w, type = "2")
    w = w / theta
    W = cbind(W, w)
    t = X %*% w - theta * t
    t = t - T %*% (t(T) %*% t)
    rho = norm(t, type = "2")
    t = t / rho
    T = cbind(T, t)
    B[(i-1),2] = theta
    B[i,1] = rho
    d = (w - theta * d) / rho
    beta = cbind(beta, beta[,(i-1)] + d %*% (t(t) %*% Y))
  }
  return(list(beta=beta[,a]))
}
#_____________________________________________________________________________________________
#_____________________________________________________________________________________________



#_____________________________________________________________________________________________
#_________________________________Bidiag1_with_reorthogonalization____________________________
#_____________________________________________________________________________________________
bidiag1 = function(X, Y, a){
  beta = NULL
  B = matrix(0, nrow = a, ncol = 2)
  gamma = norm(Y, type = "2")
  t = Y / gamma
  T = t
  w = t(X) %*% t
  alpha = norm(w, type = "2")
  w = w / alpha
  W = w
  d = w
  reg = 0
  rhobar = alpha
  phibar = gamma
 
  for(i in 1:a){
    t = X %*% w - alpha * t
    t = t - T %*% (t(T) %*% t)
    gamma = norm(t, type = "2")
    t = t /gamma
    T = cbind(T, t)
    w = t(X) %*% t - gamma * w
    w = w - W %*% (t(W) %*% w)
    alpha = norm(w, type = "2")
    w = w / alpha
    W = cbind(W, w)
    rho = norm(c(rhobar,gamma), type = "2")
    cos = rhobar / rho
    sin = gamma / rho
    theta = sin * alpha
    rhobar = -cos * alpha
    phi = cos * phibar
    phibar = sin * phibar
    G = matrix(c(cos, sin,sin,-cos), nrow = 2, ncol = 2, byrow = TRUE)
    T[,(i:(i+1))] = T[,(i:(i+1))] %*% G
    B[(i-1), 2] = theta
    B[i, 1] = rho
    reg = reg + d * (phi / rho)
    beta = cbind(beta, reg)
    d = w - d * (theta / rho)
  }
  return(list(beta=beta[,a]))
}
#_____________________________________________________________________________________________
#_____________________________________________________________________________________________


#_____________________________________________________________________________________________
#_____________________________________________SIMPLS___________________________________________
#_____________________________________________________________________________________________

simpls <- function(X, Y, a) {
  Y = as.matrix(Y)
  n = nrow(X)
  k = ncol(X)
  m = ncol(Y)
 
  Ps = matrix(0, k, a)
  Cs = matrix(0, m, a)
  Rs = matrix(0, k, a)
  Ts = matrix(0, n, a)
 
  mx = apply(X, 2, mean)
  sdx = apply(X, 2, sd)
  X = sapply(1:k, function(i) (X[,i]-mx[i]))
  my = apply(Y, 2, mean)
  sdy = apply(Y, 2, sd)
  Y = sapply(1:m, function(i) (Y[,i]-my[i]))
  S = t(X)%*%Y
 
  Snew = S
  for (i in 1:a){    
    rs = svd(Snew)$u[,1,drop=FALSE]
    rs = rs/norm(rs,type="F")
    ts = X %*% rs
    tsn = ts/norm(ts,type="F")
    ps = t(X) %*% tsn
    cs = t(Y) %*% tsn
    Rs[,i] = rs
    Ts[,i] = ts
    Ps[,i] = ps
    Cs[,i] = cs
    Snew = Snew-Ps[,1:i] %*% solve(t(Ps[,1:i]) %*% Ps[,1:i]) %*% t(Ps[,1:i]) %*% Snew
  }
  beta = Rs %*% solve(t(Ps) %*% Rs) %*% t(Cs)
 
  return(list(beta = beta))
}
#_____________________________________________________________________________________________
#_____________________________________________________________________________________________
