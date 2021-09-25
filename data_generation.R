data_generation = function(n, j){
  s = seq(0, 1, length.out = j)
 
  ksi = list()
  for(ik in 1:5){
    ksi[[ik]] = rnorm(n, 4, sd = (3*ik^(-1/2)))
  }
 
  phi = list()
  for(ik in 1:5){
    phi[[ik]] = 4*sin(ik * pi * s^2 ) + cos(4*ik * pi * s^2)
  }
 
  fX = Reduce("+", lapply(1:5, function(k){ksi[[k]] %*% t(phi[[k]])}))
 
  Beta_fun = 2*cos(pi * s)
 
  fX = fdata(fX, argvals = s)
  Beta_fun = fdata(Beta_fun, argvals = s)
 
  err = rnorm(n, mean=0, sd=1)
 
  Y = inprod.fdata(fX, Beta_fun)
 
  Y = Y + err
 
  return(list(Y = Y, X = fX$data, Beta = c(Beta_fun$data)))
 
}
