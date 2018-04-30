IS_NMF_SMOOTH <- function(V, iters, W, H, lambda){
  f = dim(V)[1]
  n = dim(V)[2]
  k = dim(W)[2]
  cost = numeric(length = iters)
  cost_pen = numeric(length = iters)
  V_appx = W %*% H
  cost[1] = sum(V/V_appx - log(V/V_appx)) - f * n
  r = H[, 1:(n-1)]/H[, 2:n]
  cost_pen[1] = cost[1] + lambda * (sum(r - log(r)) - k * (n-1))
  for (i in 2:iters) {
    Gp = t(W) %*% (1/V_appx)
    Gn = t(W) %*% (V * (V_appx^(-2)))
    Ht = H
    
    p2 = matrix(Gp[,1],ncol = 1) + lambda/matrix(H[,2],ncol = 1)
    p1 = -lambda
    p0 = - matrix(Gn[,1] * (Ht[,1]^2), ncol = 1)
    H[,1] = (sqrt(p1^2-4*p2*p0) - p1)/(2*p2)
    
    for (j in 2:(n-1)) {
      H[,j] = sqrt((Gn[,j] * (Ht[,j]^2) + lambda * H[,(j-1)])/(Gp[,j] + lambda/H[,j+1]))
    }
    
    p2 = matrix(Gp[,n], ncol = 1)
    p1 = lambda
    p0 = - matrix(Gn[,n] * (Ht[,n]^2) + lambda * H[,(n-1)], ncol = 1)
    H[,n] = (sqrt(p1^2-4*p2*p0) - p1)/(2*p2)
    
    W = W * ((V * V_appx^(-2)) %*% t(H))/((1/V_appx) %*% t(H))
    V_appx = W %*% H
    
    total = apply(W, 2, sum)
    W = W * matrix(rep(1/total, f), nrow = f, byrow = TRUE)
    H = H * matrix(rep(total, n), ncol = n, byrow = FALSE)
    
    cost[i] = sum(V/V_appx - log(V/V_appx)) - f * n
    
    r = H[, 1:(n-1)]/H[, 2:n]
    cost_pen[i] = cost[i] + lambda * (sum(r - log(r)) - k * (n-1))
  }
  list(Matrix_W = W, Matrix_H = H, cost = cost, cost_pen = cost_pen)
}

IS_NMF_EM <- function(V, iters, W, H){
  f = dim(V)[1]
  n = dim(V)[2]
  cost = numeric(length = iters)
  V_appx = W %*% H
  cost[1] <- sum(V/V_appx - log(V/V_appx)) - f * n
  for (i in 2:iters) {
    for (j in 1:k) {
      PowC_j <- matrix(W[,j], ncol = 1) %*% matrix(H[j,], nrow = 1)
      PowR_j <- V_appx - PowC_j
      G_j <- PowC_j/V_appx
      V_j <- G_j * (G_j * V + PowR_j)
      H[j,] <- t(1/matrix(W[,j], ncol = 1)) %*% V_j / f
      W[,j] <- V_j %*% t(1/matrix(H[j,], nrow = 1)) / n
      multiplier <- sqrt(sum(W[,j]^2))
      W[,j] <- W[,j] / multiplier
      H[j,] <- H[j,] * multiplier
      V_appx <- PowR_j + matrix(W[,j], ncol = 1) %*% matrix(H[j,], nrow = 1)
    }
    cost[i] <- sum(V/V_appx - log(V/V_appx)) - f * n
  }
  list(Matrix_W = W, Matrix_H = H, cost = cost)
}

IS_NMF_MU <- function(V, iters, W, H){
  f = dim(V)[1]
  n = dim(V)[2]
  cost = numeric(length = iters)
  V_appx = W %*% H
  cost[1] <- sum(V/V_appx - log(V/V_appx)) - f * n
  for (i in 2:iters) {
    W = W * ((V * (V_appx^(-2))) %*% t(H))/((1/V_appx) %*% t(H))
    V_appx = W %*% H
    
    H = H * (t(W) %*% (V * (V_appx^(-2))))/(t(W) %*% (1/V_appx))
    V_appx = W %*% H
    
    sqroot = sqrt(apply(W^2, 2, sum))
    W = W * matrix(rep(1/sqroot, f), nrow = f, byrow = TRUE)
    H = H * matrix(rep(sqroot, n), ncol = n, byrow = FALSE)
    
    cost[i] <- sum(V/V_appx - log(V/V_appx)) - f * n
  }
  list(Matrix_W = W, Matrix_H = H, cost = cost)
}