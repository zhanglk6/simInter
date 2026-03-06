#library(bigsplines)
#library(xtable)
load(file = "cont_subset.RData")
load(file = "spline0.RData")
load(file = "spline1.RData")
#source("randomization_methods.R")
#source("sample_model.R")
#source("test_formula.R")
#library(dplyr)
#library(MASS)
#library(doSNOW)
#library(parallel)
#library(xtable)
source("randomization_methods.R")
source("sample_model.R")
source("test_formula.R")
library(ranger)
library(bigsplines)

syn<-function(cont_subset,ssint0,ssint1,pi, n_boot = 600, Iternum=100){
  seeds = 439
  set.seed(seeds)
  n = nrow(cont_subset)
  syndata_origin = cont_subset[,c("A2","GENDER","Mstatus2","TreatPD","AGE",
                                  "HAMD17","HAMD24","HAMD_COGNID","number.of.depressive.episode","HAMA_SOMATI")]
  
  boot_ind = sample(1:n,n_boot,replace = TRUE)
  syndata = syndata_origin[boot_ind,]
  n = n_boot
  #syndata = syndata_origin
  
  #HAMD17_S = rep(1, n)
  #HAMD17_S[which(syndata$HAMD17<18)]=0
  #HAMD17_S[which(syndata$HAMD17>21)]=2
  
  AGE_D = rep(1,n)
  AGE_D[which(syndata$AGE<40)]=0
  AGE_D[which(syndata$AGE>49)]=2
  
  #S = stratify(syndata$GENDER,HAMD17_S)
  S = stratify(syndata$GENDER)
  
  #2 ssint
  str_ind0 = which(syndata$GENDER == 0)
  str_ind1 = which(syndata$GENDER == 1)
  
  SRS_Y = numeric(n)
  SBR_Y = numeric(n)
  SBCD_Y = numeric(n)
  WEI_Y = numeric(n)
  
  X1 = syndata$HAMA_SOMATI
  X2 = syndata$AGE
 
  Z = syndata[,-c(1,10)]
  
  numCores <- 8
  cl <- makeCluster(numCores)
  registerDoSNOW(cl)
  
  result <- foreach(i = 1:Iternum, .combine = 'rbind', .packages = c("sandwich","np","MASS", "dplyr", "caratINT")) %dopar% {
    set.seed(seeds + i)
    #for (i in 1:Iternum) {
    source("randomization_methods.R")
    source("sample_model.R")
    source("test_formula.R")
    library(ranger)
    library(bigsplines)
    SRS_A = SRS(n, pi)
    SBR_A = SBR(S, pi, bsize = 6)
    SBCD_A = SBCD(S, pi)
    #WEI_A = WEI(S, pi)
    
    syndata$A2 = SRS_A
    SRS_Y[str_ind0] = predict(ssint0,syndata[str_ind0,])
    SRS_Y[str_ind1] = predict(ssint1,syndata[str_ind1,])
    
    syndata$A2 = SBR_A
    SBR_Y[str_ind0] = predict(ssint0,syndata[str_ind0,])
    SBR_Y[str_ind1] = predict(ssint1,syndata[str_ind1,])
    
    syndata$A2 = SBCD_A
    SBCD_Y[str_ind0] = predict(ssint0,syndata[str_ind0,])
    SBCD_Y[str_ind1] = predict(ssint1,syndata[str_ind1,])
    
    SRS_fitted_model = model.syn.Ya.fit(SRS_Y, SRS_A, S, X1, Z)
    
    
    Y1_fit = SRS_fitted_model$Y1_fit
    Y0_fit = SRS_fitted_model$Y0_fit
    
    c(
    c(ols.test.cont(SRS_Y, SRS_A, X1), 
      usual.test.cont(SRS_Y, SRS_A, X1), 
      mod.test.cont(SRS_Y, SRS_A, S, X1, pi, pi*(1-pi)),
      eff.test.cont(SRS_Y, SRS_A, S, X1, pi, Y1_fit, Y0_fit),
      ols.test.cont(SBR_Y, SBR_A, X1), 
      usual.test.cont(SBR_Y, SBR_A, X1), 
      mod.test.cont(SBR_Y, SBR_A, S, X1, pi, 0),
      eff.test.cont(SBR_Y, SBR_A, S, X1, pi, Y1_fit, Y0_fit),
      ols.test.cont(SBCD_Y, SBCD_A, X1), 
      usual.test.cont(SBCD_Y, SBCD_A, X1), 
      mod.test.cont(SBCD_Y, SBCD_A, S, X1, pi, 0),
      eff.test.cont(SBCD_Y, SBCD_A, S, X1, pi, Y1_fit, Y0_fit)),
    c(ols.test.cont(SRS_Y, SRS_A, X2), 
      usual.test.cont(SRS_Y, SRS_A, X2), 
      mod.test.cont(SRS_Y, SRS_A, S, X2, pi, pi*(1-pi)),
      eff.test.cont(SRS_Y, SRS_A, S, X2, pi, Y1_fit, Y0_fit),
     ols.test.cont(SBR_Y, SBR_A, X2), 
        usual.test.cont(SBR_Y, SBR_A, X2), 
      mod.test.cont(SBR_Y, SBR_A, S, X2, pi, 0),
      eff.test.cont(SBR_Y, SBR_A, S, X2, pi, Y1_fit, Y0_fit),
     ols.test.cont(SBCD_Y, SBCD_A, X2), 
        usual.test.cont(SBCD_Y, SBCD_A, X2), 
      mod.test.cont(SBCD_Y, SBCD_A, S, X2, pi, 0),
      eff.test.cont(SBCD_Y, SBCD_A, S, X2, pi, Y1_fit, Y0_fit))
    )
  }
  stopCluster(cl)
  return(apply(result, 2, mean))
  #return(result)
}

syn.effect.plot <- function(cont_subset,ssint0,ssint1, n_boot = 600) {
  seeds = 439
  set.seed(seeds)
  n = nrow(cont_subset)
  syndata_origin = cont_subset[,c("A2","GENDER","Mstatus2","TreatPD","AGE",
                                  "HAMD17","HAMD24","HAMD_COGNID","number.of.depressive.episode","HAMA_SOMATI")]
  
  #boot_ind = sample(1:n,n_boot,replace = TRUE)
  #syndata = syndata_origin[boot_ind,]
  #n = n_boot
  
  syndata = syndata_origin
  
  #2 ssint
  str_ind0 = which(syndata$GENDER == 0)
  str_ind1 = which(syndata$GENDER == 1)
  X1 = syndata$HAMA_SOMATI
  X2 = syndata$AGE
  
  syn_Y1 = numeric(n)
  syn_Y0 = numeric(n)
  
  syndata$A2 = 1
  syn_Y1[str_ind0] = predict(ssint0,syndata[str_ind0,])
  syn_Y1[str_ind1] = predict(ssint1,syndata[str_ind1,])
  syndata$A2 = 0
  syn_Y0[str_ind0] = predict(ssint0,syndata[str_ind0,])
  syn_Y0[str_ind1] = predict(ssint1,syndata[str_ind1,])
  
  par(mfrow = c(1,2))
  
  plot(X1, syn_Y1 - syn_Y0, xlab = "HAMA_SOMATI", ylab = "Treatment Effect", cex.lab=1)
  spline.fit = smooth.spline(X1, syn_Y1 - syn_Y0, cv=T)
  lines(spline.fit, col = "red", lwd = 2)
  mtext("A", side = 3, line = -2, adj = 0.05, font = 1) 
  
  plot(X2, syn_Y1 - syn_Y0, xlab = "AGE", ylab = "Treatment Effect",cex.lab=1)
  spline.fit = smooth.spline(X2, syn_Y1 - syn_Y0, cv=T)
  lines(spline.fit,col = "red", lwd = 2)
  mtext("B", side = 3, line = -2, adj = 0.05, font = 1) 
  
}



