##### Continuous X sample generator
#sample_model1 <- function(N, delta1) {
#  ## delta1 = 0 for no interaction for X1
#  sample_data_linear_binary(N, 4, 1, 3, delta1, -2, 4)
#}

assign_to_data <- function(pi, data) {
  S = data[["S"]]
  profiles = data[["profiles"]]
  N = length(S)
  
  Y1 = data[["Y1"]]
  Y0 = data[["Y0"]]
  
  A_profolio = random.schemes(N, pi, S, profiles)
  
  SRS_A = A_profolio$SRS_A
  SBR_A = A_profolio$SBR_A
  SBCD_A = A_profolio$SBCD_A
  PS_A = A_profolio$PS_A
  
  SRS_Y = Y1*SRS_A + Y0*(1-SRS_A)
  SBR_Y = Y1*SBR_A + Y0*(1-SBR_A)
  SBCD_Y = Y1*SBCD_A + Y0*(1-SBCD_A)
  PS_Y = Y1*PS_A + Y0*(1-PS_A)
  
  data$SRS_A = SRS_A
  data$SBR_A = SBR_A
  data$SBCD_A = SBCD_A
  data$PS_A = PS_A
 
  data$SRS_Y = SRS_Y
  data$SBR_Y = SBR_Y
  data$SBCD_Y = SBCD_Y
  data$PS_Y = PS_Y

  return(data)
}

sample_model1 <- function(N, pi, delta1, transformed = FALSE) {
  ## delta1 = 0 for no interaction for X1
  data = sample_data_linear(N, 4, 1, 3, delta1, -2, 4)
  data$S = stratify(data[["X_d"]], data[["W_d"]])
  data$profiles = cbind(data[["X_d"]], data[["W_d"]])
  data <- assign_to_data(pi, data)
  
  S = data[["S"]]
  X = data[["X"]]
  Z = data[["W"]]
  
  SRS_A = data$SRS_A
  SBR_A = data$SBR_A
  SBCD_A = data$SBCD_A
  PS_A = data$PS_A
  
  SRS_Y = data$SRS_Y
  SBR_Y = data$SBR_Y
  SBCD_Y = data$SBCD_Y
  PS_Y = data$PS_Y
  
  SRS_fitted_model = model1.Ya.fit(SRS_Y, SRS_A, S, X, Z)
  SBR_fitted_model = model1.Ya.fit(SBR_Y, SBR_A, S, X, Z)
  SBCD_fitted_model = model1.Ya.fit(SBCD_Y, SBCD_A, S, X, Z)
  PS_fitted_model = model1.Ya.fit(PS_Y, PS_A, S, X, Z)
  
  data$SRS_Y1_fit = SRS_fitted_model$Y1_fit
  data$SRS_Y0_fit = SRS_fitted_model$Y0_fit
  data$SBR_Y1_fit = SBR_fitted_model$Y1_fit
  data$SBR_Y0_fit = SBR_fitted_model$Y0_fit
  data$SBCD_Y1_fit = SBCD_fitted_model$Y1_fit
  data$SBCD_Y0_fit = SBCD_fitted_model$Y0_fit
  data$PS_Y1_fit = PS_fitted_model$Y1_fit
  data$PS_Y0_fit = PS_fitted_model$Y0_fit
  
  if (transformed) {
    X = 2 - 3 + delta1 *X
    data$X = X + 1 / sqrt(N) * rnorm(N) / 5 * ifelse(delta1 == 0, 1, 0)
  }
  
  return(data)
}

sample_model2 <- function(N, pi, delta1, transformed = FALSE) {
  ## delta1 = 0 for no interaction for X
  data = sample_data_exp(N, 5, 4, 0.5, delta1, 2,6)
  S = stratify(data[["X_d"]], data[["W_d"]])
  data$S = S
  data$profiles = cbind(data[["X_d"]], data[["W_d"]])
  data <- assign_to_data(pi, data)
  
  S = data[["S"]]
  X = data[["X"]]
  Z = data[["W"]]
  
  SRS_A = data$SRS_A
  SBR_A = data$SBR_A
  SBCD_A = data$SBCD_A
  PS_A = data$PS_A
  
  SRS_Y = data$SRS_Y
  SBR_Y = data$SBR_Y
  SBCD_Y = data$SBCD_Y
  PS_Y = data$PS_Y
  
  SRS_fitted_model = model2.Ya.fit(SRS_Y, SRS_A, S, X, Z)
  SBR_fitted_model = model2.Ya.fit(SBR_Y, SBR_A, S, X, Z)
  SBCD_fitted_model = model2.Ya.fit(SBCD_Y, SBCD_A, S, X, Z)
  PS_fitted_model = model2.Ya.fit(PS_Y, PS_A, S, X, Z)
  
  data$SRS_Y1_fit = SRS_fitted_model$Y1_fit
  data$SRS_Y0_fit = SRS_fitted_model$Y0_fit
  data$SBR_Y1_fit = SBR_fitted_model$Y1_fit
  data$SBR_Y0_fit = SBR_fitted_model$Y0_fit
  data$SBCD_Y1_fit = SBCD_fitted_model$Y1_fit
  data$SBCD_Y0_fit = SBCD_fitted_model$Y0_fit
  data$PS_Y1_fit = PS_fitted_model$Y1_fit
  data$PS_Y0_fit = PS_fitted_model$Y0_fit
  
  if (transformed) {
    X = 5 - 4 + exp((0.5 + delta1)*X) - exp(0.5*X)
    data$X = X + 1 / sqrt(N) * rnorm(N) / 5 * ifelse(delta1 == 0, 1, 0)
  }
  
  return(data)
}
## N, mu1, mu0, beta1, alpha1,alpha2, gamma1, gamma2, delta1

sample_model3 <- function(N, pi, delta1, transformed = FALSE) {
  ## delta1 = 0 for no interaction for X
  data = sample_data_log(N, 2, -1, 1, -1, 2, 4,-2, delta1)
  #S = stratify(data[["W_d"]])
  S = stratify(data[["W_d"]], data[["X_d"]])
  data$S = S
  data$profiles = cbind(data[["X_d"]], data[["W_d"]])
  data <- assign_to_data(pi, data)
  
  S = data[["S"]]
  X = data[["X"]]
  Z = data[["W"]]
  
  SRS_A = data$SRS_A
  SBR_A = data$SBR_A
  SBCD_A = data$SBCD_A
  PS_A = data$PS_A
  
  SRS_Y = data$SRS_Y
  SBR_Y = data$SBR_Y
  SBCD_Y = data$SBCD_Y
  PS_Y = data$PS_Y
  
  SRS_fitted_model = model3.Ya.fit(SRS_Y, SRS_A, S, X, Z)
  SBR_fitted_model = model3.Ya.fit(SBR_Y, SBR_A, S, X, Z)
  SBCD_fitted_model = model3.Ya.fit(SBCD_Y, SBCD_A, S, X, Z)
  PS_fitted_model = model3.Ya.fit(PS_Y, PS_A, S, X, Z)
  
  data$SRS_Y1_fit = SRS_fitted_model$Y1_fit
  data$SRS_Y0_fit = SRS_fitted_model$Y0_fit
  data$SBR_Y1_fit = SBR_fitted_model$Y1_fit
  data$SBR_Y0_fit = SBR_fitted_model$Y0_fit
  data$SBCD_Y1_fit = SBCD_fitted_model$Y1_fit
  data$SBCD_Y0_fit = SBCD_fitted_model$Y0_fit
  data$PS_Y1_fit = PS_fitted_model$Y1_fit
  data$PS_Y0_fit = PS_fitted_model$Y0_fit
  
  if (transformed) {
    X =  log((1 + delta1)*X + 4) - log(1*X + 4)
    data$X = X + 1 / sqrt(N) * rnorm(N) / 5 * ifelse(delta1 == 0, 1, 0)
  }
  
  
  return(data)
}

sample_model4 <- function(N, pi, delta1, transformed = FALSE) {
  ## delta1 = 0 for no interaction for X
  data = sample_data_log_high(N, 2, -1, 1, -1, 2, 4,-2, delta1)
  S = stratify(data[["W_d"]], data[["X_d"]])
  data$S = S
  data$profiles = cbind(data[["X_d"]], data[["W_d"]])
  data <- assign_to_data(pi, data)
  
  S = data[["S"]]
  X = data[["X"]]
  Z = cbind(data[["W"]], data[["U"]], data[["high_covariates"]])
  
  SRS_A = data$SRS_A
  SBR_A = data$SBR_A
  SBCD_A = data$SBCD_A
  PS_A = data$PS_A
  
  SRS_Y = data$SRS_Y
  SBR_Y = data$SBR_Y
  SBCD_Y = data$SBCD_Y
  PS_Y = data$PS_Y
  
  SRS_fitted_model = model4.Ya.fit(SRS_Y, SRS_A, S, X, Z)
  SBR_fitted_model = model4.Ya.fit(SBR_Y, SBR_A, S, X, Z)
  SBCD_fitted_model = model4.Ya.fit(SBCD_Y, SBCD_A, S, X, Z)
  PS_fitted_model = model4.Ya.fit(PS_Y, PS_A, S, X, Z)
  
  data$SRS_Y1_fit = SRS_fitted_model$Y1_fit
  data$SRS_Y0_fit = SRS_fitted_model$Y0_fit
  data$SBR_Y1_fit = SBR_fitted_model$Y1_fit
  data$SBR_Y0_fit = SBR_fitted_model$Y0_fit
  data$SBCD_Y1_fit = SBCD_fitted_model$Y1_fit
  data$SBCD_Y0_fit = SBCD_fitted_model$Y0_fit
  data$PS_Y1_fit = PS_fitted_model$Y1_fit
  data$PS_Y0_fit = PS_fitted_model$Y0_fit
  
  if (transformed) {
    X =  log((1 + delta1)*X + 4) - log(1*X + 4)
    data$X = X + 1 / sqrt(N) * rnorm(N) / 5 * ifelse(delta1 == 0, 1, 0)
  }
  
  
  return(data)
}

sample_model_small <- function(N, pi, delta1) {
  ## delta1 = 0 for no interaction for X
  data = sample_data_linear_small(N, 4, 1, 3, 0, delta1, -2, 4)
  #data = sample_data_log(N, 2, -1, 1, -1, 2, 4,-2, delta1)
  S = stratify(data[["W_d"]])
  data$S = S
  data$profiles = cbind(data[["W_d"]])
  data <- assign_to_data(pi, data)
  
  S = data[["S"]]
  X = data[["X"]]
  Z = data[["W"]]
  
  SRS_A = data$SRS_A
  SBR_A = data$SBR_A
  SBCD_A = data$SBCD_A
  PS_A = data$PS_A
  
  SRS_Y = data$SRS_Y
  SBR_Y = data$SBR_Y
  SBCD_Y = data$SBCD_Y
  PS_Y = data$PS_Y
  
  SRS_fitted_model = model.small.Ya.fit(SRS_Y, SRS_A, S, X, Z)
  SBR_fitted_model = model.small.Ya.fit(SBR_Y, SBR_A, S, X, Z)
  SBCD_fitted_model = model.small.Ya.fit(SBCD_Y, SBCD_A, S, X, Z)
  PS_fitted_model = model.small.Ya.fit(PS_Y, PS_A, S, X, Z)
  
  data$SRS_Y1_fit = SRS_fitted_model$Y1_fit
  data$SRS_Y0_fit = SRS_fitted_model$Y0_fit
  data$SBR_Y1_fit = SBR_fitted_model$Y1_fit
  data$SBR_Y0_fit = SBR_fitted_model$Y0_fit
  data$SBCD_Y1_fit = SBCD_fitted_model$Y1_fit
  data$SBCD_Y0_fit = SBCD_fitted_model$Y0_fit
  data$PS_Y1_fit = PS_fitted_model$Y1_fit
  data$PS_Y0_fit = PS_fitted_model$Y0_fit
  
  return(data)
}



#### cross-fitting
model1.Ya.fit <- function(Y, A, S, X, Z) {
  model1and2.cross.fit_strata(Y, A, S, X, Z, M=4)
}

model2.Ya.fit <- function(Y, A, S, X, Z) {
  model1and2.cross.fit_strata(Y, A, S, X, Z, M=4)
}

model3.Ya.fit <- function(Y, A, S, X, Z) {
  model3.cross.fit_strata(Y, A, S, X, Z, M=4)
}

model4.Ya.fit <- function(Y, A, S, X, Z) {
  model4.cross.fit_strata(Y, A, S, X, Z, M=4)
}

model.small.Ya.fit <- function(Y, A, S, X, Z) {
  model.small.cross.fit_strata(Y, A, S, X, Z, M=4)
}

model.syn.Ya.fit <- function(Y, A, S, X, Z) {
  model.syn.cross.fit_strata(Y, A, S, X, Z, M=4)
}

model.cross.fit_strata_old <- function(Y, A, S, X, Z, M=2) {
  n = length(Y)
  nStrata = max(S)
  fold = create.fold(n, M)
  
  total_data = data.frame(Y = Y, A = A, S = S, X = X, Z = Z)
  
  Y1_fit = numeric(n)
  Y0_fit = numeric(n)
  Y_fit = numeric(n)
  res = numeric(n)
  
  for (m in 1:M) {
    for (i in 1:nStrata) {
      index_in_fold = which(fold == m & S == i)
      index_out_fold = which(fold != m & S == i)
      train_data = total_data[index_out_fold,]
      train_data_1 = train_data %>% subset(A==1)
      train_data_0 = train_data %>% subset(A==0)
      
      bws1 = c(0.1256611,  0.1750254)
      bws0 = c(0.3156482, 0.09994133)
      
      kernel_model_1 = npregbw(Y ~ X + Z, bandwidth.compute = F, 
                               bws = bws1, data = train_data_1, regtype = "ll")
      kernel_model_0 = npregbw(Y ~ X + Z, bandwidth.compute = F, 
                               bws = bws0, data = train_data_0, regtype = "ll")
      Y1_fit[index_in_fold] = npreg(kernel_model_1, exdat = cbind(X, Z)[index_in_fold,])$mean
      Y0_fit[index_in_fold] = npreg(kernel_model_0, exdat = cbind(X, Z)[index_in_fold,])$mean
      Y_fit[index_in_fold] = A[index_in_fold] * Y1_fit[index_in_fold] + (1-A)[index_in_fold] * Y0_fit[index_in_fold]
      res[index_in_fold] = Y[index_in_fold] - Y_fit[index_in_fold]
    }
  }
  return(data.frame(Y_fit = Y_fit,res = res, Y1_fit = Y1_fit, Y0_fit = Y0_fit ))
}

model1and2.cross.fit_strata <- function(Y, A, S, X, Z, M=5) {
  model1.Y1.fit <- function(train_data_1) {
    bws1 = c(0.13,  0.18)
    #bws1 = c(0.26,  0.36)
    kernel_model_1 = npregbw(Y ~ X + Z, bandwidth.compute = F, 
                             bws = bws1, data = train_data_1, regtype = "ll")
    return(kernel_model_1)
  }
  
  model1.Y0.fit <- function(train_data_0) {
    bws0 = c(0.32, 0.10)
    kernel_model_0 = npregbw(Y ~ X + Z, bandwidth.compute = F, 
                             bws = bws0, data = train_data_0, regtype = "ll")
    return(kernel_model_0)
  }
  
  model1.Ya.predict <- function(model_a, covariates_in_fold) {
    Ya_fit = npreg(model_a, exdat = covariates_in_fold)$mean
    return(Ya_fit)
  }
  
  return(model.cross.fit_strata(Y, A, S, X, Z, M, 
                                model1.Y1.fit, model1.Y0.fit, model1.Ya.predict))
}

model3.cross.fit_strata <- function(Y, A, S, X, Z, M=5) {
  model.Y1.fit <- function(train_data_1) {
    bws1 = c(0.8,  0.4)
    #bws1 = c(1.6,  0.8)
    kernel_model_1 = npregbw(Y ~ X + Z, bandwidth.compute = F, 
                             bws = bws1, data = train_data_1, regtype = "ll")
    return(kernel_model_1)
  }
  
  model.Y0.fit <- function(train_data_0) {
    #bws0 = c(0.3, 0.4)
    bws0 = c(0.6, 0.3)
    kernel_model_0 = npregbw(Y ~ X + Z, bandwidth.compute = F, 
                             bws = bws0, data = train_data_0, regtype = "ll")
    return(kernel_model_0)
  }
  
  model.Ya.predict <- function(model_a, covariates_in_fold) {
    Ya_fit = npreg(model_a, exdat = covariates_in_fold)$mean
    return(Ya_fit)
  }
  
  return(model.cross.fit_strata(Y, A, S, X, Z, M, 
                                model.Y1.fit, model.Y0.fit, model.Ya.predict))
}

model4.cross.fit_strata <- function(Y, A, S, X, Z, M=5) {
  model.Y1.fit <- function(train_data_1) {
    rf_model_1 = ranger(Y ~ ., num.trees = 500, mtry = 5, oob.error = F,
                        data = train_data_1[,-c(2,3)], replace = F, importance = "none")
    #rf_model_1 = randomForest(Y ~ ., num.trees = 50, mtry = 5, data = train_data_1[,-c(2,3)])
    return(rf_model_1)
  }
  
  model.Y0.fit <- function(train_data_0) {
    rf_model_0 = ranger(Y ~ ., num.trees = 500, mtry = 5, oob.error = F,
                        data = train_data_0[,-c(2,3)], replace = F, importance = "none")
    #rf_model_0 = randomForest(Y ~ ., num.trees = 50, mtry = 3, data = train_data_0[,-c(2,3)])
    return(rf_model_0)
  }
  
  model.Ya.predict <- function(model_a, covariates_in_fold) {
    Ya_fit = predict(model_a, covariates_in_fold)$prediction
    #Ya_fit = predict(model_a, covariates_in_fold)
    return(Ya_fit)
  }
  
  return(model.cross.fit_strata(Y, A, S, X, Z, M, 
                                model.Y1.fit, model.Y0.fit, model.Ya.predict))
}

model.small.cross.fit_strata <- function(Y, A, S, X, Z, M=5) {
  model.Y1.fit <- function(train_data_1) {
    #ctrl <- trainControl(method = "cv", number = 2)
    ctrl <- trainControl(method = "none")
    #  Manually set k (e.g., k = 7)
    k = floor(sqrt(nrow(train_data_1)))
    knn_model_1 <- train(
      Y ~ X + Z,                # regression target
      data = train_data_1,
      method = "knn",
      trControl = ctrl,
      preProcess = c("center", "scale"),  # important for kNN
      tuneGrid = data.frame(k = k)        # manually specify k
    )
    return(knn_model_1)
  }
  # ---- Fit Y=0 model ----
  model.Y0.fit <- function(train_data_0) {
    #ctrl <- trainControl(method = "cv", number = 2)
    ctrl <- trainControl(method = "none")
    #  Manually set k (e.g., k = 7)
    k = floor(sqrt(nrow(train_data_0)))
    knn_model_0 <- train(
      Y ~ X + Z,                # regression target
      data = train_data_0,
      method = "knn",
      trControl = ctrl,
      preProcess = c("center", "scale"),  # important for kNN
      tuneGrid = data.frame(k = k)        # manually specify k
    )
    return(knn_model_0)
  }
  
  model.Ya.predict <- function(model_a, covariates_in_fold) {
    Ya_fit = predict(model_a, covariates_in_fold)
    #Ya_fit = predict(model_a, covariates_in_fold)
    return(Ya_fit)
}
  return(model.cross.fit_strata(Y, A, S, X, Z, M, 
                                model.Y1.fit, model.Y0.fit, model.Ya.predict))
}


model.syn.cross.fit_strata <- function(Y, A, S, X, Z, M=5) {
  model.Y1.fit <- function(train_data_1) {
    rf_model_1 = ranger(Y ~ ., num.trees = 500, mtry = 5, oob.error = F,
                        data = train_data_1[,-c(2,3)], replace = F, importance = "none")
    #rf_model_1 = randomForest(Y ~ ., num.trees = 50, mtry = 5, data = train_data_1[,-c(2,3)])
    return(rf_model_1)
  }
  
  model.Y0.fit <- function(train_data_0) {
    rf_model_0 = ranger(Y ~ ., num.trees = 500, mtry = 5, oob.error = F,
                        data = train_data_0[,-c(2,3)], replace = F, importance = "none")
    #rf_model_0 = randomForest(Y ~ ., num.trees = 50, mtry = 3, data = train_data_0[,-c(2,3)])
    return(rf_model_0)
  }
  
  model.Ya.predict <- function(model_a, covariates_in_fold) {
    Ya_fit = predict(model_a, covariates_in_fold)$prediction
    #Ya_fit = predict(model_a, covariates_in_fold)
    return(Ya_fit)
  }
  
  return(model.cross.fit_strata(Y, A, S, X, Z, M, 
                                model.Y1.fit, model.Y0.fit, model.Ya.predict))
}

model.cross.fit_strata <- function(Y, A, S, X, Z, M=5, model.Y1.fit, model.Y0.fit, model.Ya.predict) {
  n = length(Y)
  nStrata = max(S)
  fold = create.fold(n, M)
  
  total_data = data.frame(Y = Y, A = A, S = S, X = X, Z = Z)
  
  Y1_fit = numeric(n)
  Y0_fit = numeric(n)
  Y_fit = numeric(n)
  res = numeric(n)
  
  for (m in 1:M) {
    for (i in 1:nStrata) {
      index_in_fold = which(fold == m & S == i)
      index_out_fold = which(fold != m & S == i)
      train_data = total_data[index_out_fold,]
      train_data_1 = train_data %>% subset(A==1)
      train_data_0 = train_data %>% subset(A==0)
      
      model_1 = model.Y1.fit(train_data_1)
      model_0 = model.Y0.fit(train_data_0)
      Y1_fit[index_in_fold] = model.Ya.predict(model_1, data.frame(X = X, Z = Z)[index_in_fold,])
      Y0_fit[index_in_fold] = model.Ya.predict(model_0, data.frame(X = X, Z = Z)[index_in_fold,])
      
      Y_fit[index_in_fold] = A[index_in_fold] * Y1_fit[index_in_fold] + (1-A)[index_in_fold] * Y0_fit[index_in_fold]
      res[index_in_fold] = Y[index_in_fold] - Y_fit[index_in_fold]
    }
  }
  return(data.frame(Y_fit = Y_fit,res = res, Y1_fit = Y1_fit, Y0_fit = Y0_fit ))
}

##### Continuous X model
sample_data_linear <- function(N, mu1, mu0, beta1,   
                               delta1, alpha1, gamma1) {
  ## N:      simulated sample size
  ## mu1:    base potential outcome for the treatment group
  ## mu0:    base potential outcome for the control group
  ## beta1:  coefficient for X1
  ## gamma1: coefficient for Z1
  sigma1 = 1
  sigma0 = 0.5
  
  X1 = runif(N, 0, 2)
  X_d = ifelse(X1 > 1, 1, 0)
  W = runif(N, -1.5, 1.5)
  W_d = ifelse(W > 0, 1, 0)
  
  base = beta1 * X1  + alpha1 * W
  error1 = sigma1 * rnorm(N)
  error0 = sigma0 * rnorm(N)
  
  Y1 = mu1 + base + delta1 * X1  + gamma1 * W * X1 + error1
  Y0 = mu0 + base + error0
  return(list(Y1=Y1, Y0=Y0, X=X1, X_d = X_d,W_d=W_d, W = W))
}

sample_data_linear_small <- function(N, mu1, mu0, beta1, beta2, 
                                     delta1, alpha1, gamma1) {
  ## N:      simulated sample size
  ## mu1:    base potential outcome for the treatment group
  ## mu0:    base potential outcome for the control group
  ## beta1:  coefficient for X1
  ## gamma1: coefficient for Z1
  sigma1 = 1
  sigma0 = 0.5
  
  X1 = runif(N, 0, 2)
  X_d = ifelse(X1 > 1, 1, 0)
  W = runif(N, -1.5, 1.5)
  W_d = ifelse(W > 0, 1, 0)
  
  base = beta1 * X1 + beta2 * X1*X1 + alpha1 * W
  error1 = sigma1 * rnorm(N)
  error0 = sigma0 * rnorm(N)
  
  Y1 = mu1 + base + delta1 * X1  + gamma1 * W * X1 + error1
  Y0 = mu0 + base + error0
  return(list(Y1=Y1, Y0=Y0, X=X1, X_d = X_d,W_d=W_d, W = W))
}

sample_data_linear_binary <- function(N, mu1, mu0, beta1, delta1, alpha1, gamma1) {  
  data = sample_data_linear(N, mu1, mu0, beta1, delta1, alpha1, gamma1)
  data$X = data$X_d
  return(data)
}
#N=1000;pi=0.5;mu1=6;mu0=1;beta1=3;alpha1 = 2;gamma1 = 4;delta1 = 0
#df = sample_data_linear(N, mu1, mu0, beta1, delta1,alpha1, gamma1)

sample_data_exp <- function(N, mu1, mu0, beta1,   
                            delta1, alpha1, gamma1) {
  ## N:      simulated sample size
  ## mu1:    base potential outcome for the treatment group
  ## mu0:    base potential outcome for the control group
  
  X_star = runif(N, -1, 1)
  #X_star = rnorm(N, 0, 0.5)
  X_d = ifelse(X_star > 0, 1, 0)
  W = runif(N, -1.5, 1.5)
  W_d = ifelse(W > 0, 1, 0)
  
  error1 = exp(0.5*X_star)*rnorm(N)
  error0 = 0.5*exp(0.5*X_star)*rnorm(N)
  error = rnorm(N)
  
  #Y1 = mu1 + exp((beta1 + delta1)*X_star) + alpha1*W+ gamma1 * W *X_star + error1
  #Y0 = mu0 +  exp(beta1*X_star) + alpha1*W + error0
  Y1 = mu1 + exp((beta1 + delta1)*X_star) + alpha1*W+ gamma1 * W *X_star + error
  Y0 = mu0 +  exp(beta1*X_star) + alpha1*W + error
  return(list(Y1=Y1, Y0=Y0, X=X_star, X_d = X_d, W_d=W_d, W = W))
}
#N=1000;pi=0.5;mu1=5;mu0=4;beta1=0.5;alpha1 = 0.5;delta1 = 0
#df = sample_data_exp(N, mu1, mu0, beta1, delta1,alpha1)
sample_data_exp_binary <- function(N, mu1, mu0, beta1,   
                            delta1, alpha1, gamma1) {
  data = sample_data_exp(N, mu1, mu0, beta1, delta1, alpha1, gamma1)
  data$X = data$X_d
  return(data)
}

sample_data_log <- function(N, mu1, mu0, beta1, alpha1,
                             alpha2, gamma1, gamma2, delta1) {
  ## N:      simulated sample size
  ## mu1:    base potential outcome for the treatment group
  ## mu0:    base potential outcome for the control group
  
  X = runif(N, -1, 1)
  X_d = ifelse(X > 0, 1, 0)
  W = runif(N, -1.5, 1.5)
  W_d = ifelse(W > 0, 1, 0)
  U = rbeta(N, 2, 4) - 1/3
  
  sigma1 = 0.5
  sigma0 = 1
  error1 = sigma1*rnorm(N)
  error0 = sigma0*rnorm(N)
  error = rnorm(N)
  
  Y1 = mu1 + log((beta1 + delta1)*X + 4) + alpha1*W + alpha2*U^2 + gamma1*W*X + gamma2*U*X + error1
  Y0 = mu0 +  log(beta1*X + 4) + alpha1*sin(W) + error0
  
  ## high dimension setting
  return(list(Y1=Y1, Y0=Y0, X=X, X_d=X_d, W_d=W_d, W = W, U = U))
}

sample_data_log_high <- function(N, mu1, mu0, beta1, alpha1,
                            alpha2, gamma1, gamma2, delta1) {
  ## N:      simulated sample size
  ## mu1:    base potential outcome for the treatment group
  ## mu0:    base potential outcome for the control group
  data = sample_data_log(N, mu1, mu0, beta1, alpha1,
                         alpha2, gamma1, gamma2, delta1)
  
  ## high dimension setting
  p <- 3
  # zeta <- c(1, -0.8, 0.6)
    # Define the covariance matrix
    cov_matrix <- matrix(0.2, nrow = p, ncol = p)
    diag(cov_matrix) <- 1
    # Generate iid standard normal variables
    standard_normals <- matrix(rnorm(N * p), nrow = N, ncol = p)
    # Cholesky decomposition of the covariance matrix
    chol_decomp <- chol(cov_matrix)
    # Transform standard normal samples to target multivariate normal distribution
    #high_covariates <- standard_normals %*% chol_decomp 
    high_covariates <- cbind(data[["X"]]*data[["W"]], data[["X"]]*data[["U"]], standard_normals %*% chol_decomp )
    #high_covariates[,1] <- high_covariates[,1] * data[["W"]]
    #high_covariates[,2] <- high_covariates[,2] * data[["U"]]
    
  return(list(Y1 = data[["Y1"]], Y0 = data[["Y0"]], X = data[["X"]], X_d = data[["X_d"]], W_d = data[["W_d"]],
              W  = data[["W"]], U = data[["U"]], high_covariates = high_covariates))
}
