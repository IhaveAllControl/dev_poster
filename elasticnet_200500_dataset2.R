library(elasticnet)
library(foreach)
library(doParallel)
library(glmnet)
library(MASS)
cl <- makeCluster(11)
registerDoParallel(cl)


simul_data2 = function(n, p){
  A = matrix(0,p,p)
  diag(A) = 1
  for(i in 1:(p-1)) for(j in (i+1):p) if(abs(i-j)<=2) {A[i,j]=A[j,i]=0.4}
  X = mvrnorm(n = n, mu = rep(0,p) , Sigma = A)
  #B = c(rep(1, 0.05*p), rep(-1, 0.05*p), rep(2,0.05*p),rep(0, 0.85*p))
  c1_be <- rep(0,p)
  c1_be[1:20] = 2
  c1_be[21:25] = -1
  c1_be[26:50] = 3
  
  B = c1_be
  e = rnorm(n, 0, 1)
  y = X %*% B + e
  
  test_X = mvrnorm(n = n, mu = rep(0,p) , Sigma = A)
  test_e = rnorm(n, 0, 1)
  test_y = test_X %*% B + test_e
  return(list(X = scale(X,T,F), y=scale(y,T,F), testx = scale(test_X,T,F), testy = scale(test_y,T,F)))
}


elastic3type =function(n, p, alpha0,lambda_grid,fold_number){
  set.seed(2021)
  type1traincv = vector()
  type1starerr = vector()
  
  type2traincv = vector()
  type2starerr = vector()
  
  type3traincv = vector()
  type3starerr = vector()
  alpha_star = seq(2, 0.5, -0.01)
  
  for (i in 1:50) {
    data1 = simul_data2(n, p)
    
    num.fol = floor(dim(data1$X)[1]/fold_number)
    res.vector = vector()
    for (fold in 1:fold_number) {
      
      if (fold!=fold_number) {
        i.except.train = (1 + (fold-1) * num.fol):(fold * num.fol)
      }
      else{
        i.except.train = (1 + (fold-1) * num.fol):dim(data1$X)[1]
      }
      glmnet_cv_train = glmnet(data1$X[-i.except.train,], data1$y[-i.except.train], alpha = alpha0, 
                               lambda = lambda_grid/(dim(data1$X)[1]-length(i.except.train))
                               ,  intercept = F, family = 'gaussian')
      glmnet_cv_pred = predict(glmnet_cv_train, newx = matrix(data1$X[i.except.train,], length(i.except.train), p))
      apply_glmnet = sapply(1:length(lambda_grid), function(sx) mean( (data1$y[i.except.train] - glmnet_cv_pred[,sx])^2 ))
      
      res.vector = cbind(res.vector, apply_glmnet)
      print(paste0(fold,"/", fold_number,"----fold in type1"))
      
    }
   
    cvout = apply(res.vector, 1, mean)
    n_train = dim(data1$X)[1]-length(i.except.train) # ¼öÁ¤ ÇÊ¿ä
    
    ela1_2 = glmnet(data1$X, data1$y, alpha = alpha0, lambda = (lambda_grid[which.min(cvout)]/n_train) * alpha_star, intercept = F
                    , family = 'gaussian')
    yhat1_test2 = predict(ela1_2, newx = data1$testx)
    star_test_error = sapply(1:length(alpha_star), 
                             function(xx2) mean( (data1$testy - yhat1_test2[,xx2])^2 ))
    type1traincv = cbind(type1traincv, cvout)
    type1starerr = cbind(type1starerr, star_test_error)
    
    ######################################################################################
    
    num.fol = floor(dim(data1$X)[1]/fold_number)
    res.vector = vector()
    for (fold in 1:fold_number) {
      
      if (fold!=fold_number) {
        i.except.train <- (1 + (fold-1) * num.fol):(fold * num.fol)
      }
      else{
        i.except.train <- (1 + (fold-1) * num.fol):dim(data1$X)[1]
      }
      ela2 = glmnet(data1$X[-i.except.train,], data1$y[-i.except.train], alpha = alpha0,
                    lambda = lambda_grid/(dim(data1$X)[1]-length(i.except.train)),
                     intercept = F, family = 'gaussian')
      # alpha=0 --------> foreach too slow.
      if (alpha0==0) {
        pred_alpha0 = predict(ela2, newx = matrix(data1$X[i.except.train,], length(i.except.train), p))
        cv2 = sapply(1:length(lambda_grid), function(aa) mean( (data1$y[i.except.train] - pred_alpha0[,aa])^2 ))
      }
      else{
        fr = foreach(lv = 1:length(lambda_grid), .combine = rbind, .packages = 'glmnet') %dopar% {
          
          m = as.vector(ela2$beta[,lv]!=0)
          if (sum(m)>=2) {
            ela2_2 = glmnet(data1$X[-i.except.train, m], data1$y[-i.except.train],
                            alpha = 0, lambda = seq(lambda_grid[length(lambda_grid)], lambda_grid[lv],length.out = 100)/(dim(data1$X)[1]-length(i.except.train)),
                             intercept = F, family = 'gaussian')
            pred = predict(ela2_2, newx = matrix(data1$X[i.except.train, m], length(i.except.train), sum(m)))[, 100]
            
            mean( (data1$y[i.except.train] - pred)^2 )
            
            
          }else if (sum(m)==0){
            mean(data1$y[i.except.train]^2)
          }
          else{
            nX = as.matrix(data1$X[-i.except.train, m])
            beta_hat = 1/(t(nX) %*% nX + lambda_grid[lv]/(dim(data1$X)[1]-length(i.except.train))) * t(nX) %*% data1$y[-i.except.train]
            
            mean( (data1$y[i.except.train] - as.vector(data1$X[i.except.train, m]) %*%beta_hat)^2 ) 
          }
          
        }
        res.vector = cbind(res.vector,fr[,1])
        print(paste0(fold,"/",fold_number,"--------------------type2_1"))
      }
      
    }
    
    if (alpha0!=0) {
      cv2 = apply(res.vector, 1, mean)
      print("alpha isn't zero")  
    }
    lambdamin = lambda_grid[which.min(cv2)]
    
    n_train = n_train = dim(data1$X)[1]-length(i.except.train)
    ela2_test_err = glmnet(data1$X, data1$y, alpha = alpha0,
                           lambda = alpha_star*lambdamin/n_train,
                           intercept = F, family = 'gaussian')
    fr2 = foreach(lv = 1:length(alpha_star), .combine = rbind, .packages = 'glmnet') %dopar% {
      
      
      m = as.vector(ela2_test_err$beta[,lv]!=0)
      if (sum(m)>=2) {
        ela2_2_test_err = glmnet(data1$X[, m], data1$y,
                                 alpha = 0, lambda = alpha_star*lambdamin/n_train,
                                intercept = F, family = 'gaussian')
        pred = predict(ela2_2_test_err, newx = data1$testx[, m])[,lv]
        
        mean( (data1$testy - pred)^2 )
        
        
      }else if (sum(m)==0){
        mean(data1$testy^2)
      }
      else{
        nX = as.matrix(data1$X[, m])
        beta_hat = 1/(t(nX) %*% nX + lambdamin*alpha_star[lv]/n_train) * t(nX) %*% data1$y
        
        mean( (data1$testy - as.vector(data1$testx[, m]) %*%beta_hat)^2 ) 
      }
      
    }
    
    type2traincv = cbind(type2traincv, cv2)
    type2starerr = cbind(type2starerr, fr2[,1])
    
    
    
    print("end type2")
    ########################################################################################
    
    
    num.fol = floor(dim(data1$X)[1]/fold_number)
    res.vector = vector()
    if (alpha0==0) {
      next
      
    }
    
    for (fold in 1:fold_number) {
      
      if (fold!=fold_number) {
        i.except.train <- (1 + (fold-1) * num.fol):(fold * num.fol)
      }
      else{
        i.except.train <- (1 + (fold-1) * num.fol):dim(data1$X)[1]
      }
      ela3 = glmnet(data1$X[-i.except.train,], data1$y[-i.except.train], alpha = alpha0,
                    lambda = lambda_grid/(dim(data1$X)[1]-length(i.except.train)),
                     intercept = F, family = 'gaussian')
      
      
      fr = foreach(lv = 1:length(lambda_grid), .combine = rbind, .packages = "glmnet") %dopar% {
        m = as.vector(ela3$beta[,lv]!=0)
        
        if ((sum(m) <= dim(data1$X)[1] - length(i.except.train)) & sum(m)>0) {
          lX = as.data.frame(data1$X[-i.except.train, m])
          lX["y"] = data1$y[-i.except.train]
          lmout = lm(y~.+0, data = lX)
          pred = matrix(data1$X[i.except.train, m], length(i.except.train), sum(m)) %*% coef(lmout)
          mean( (data1$y[i.except.train] - pred)^2 )
          
        }
        else if(sum(m)==0){
          mean( (data1$y[i.except.train])^2 )
        }
        else{
          NA
        }
      }
      
      res.vector = cbind(res.vector, fr[,1])
      print(paste0(fold,"/",fold_number,"--------------------type3"))
    }
    cv3 = apply(res.vector, 1, mean,na.rm=T)
    lambdamin = lambda_grid[which.min(cv3)]
    
    ela3_star = glmnet(data1$X, data1$y, alpha = alpha0,
                       lambda = lambdamin * alpha_star/(dim(data1$X)[1] - length(i.except.train)),
                       intercept = F, family = 'gaussian')
    fr2 = foreach(lv = 1:length(alpha_star), .combine = rbind, .packages = 'glmnet') %dopar% {
      m = as.vector(ela3_star$beta[,lv]!=0)
      if ( (sum(m) <= dim(data1$X)[1]) & sum(m)>0  ) {
        lX = as.data.frame(data1$X[,m])
        lX["y"] = data1$y
        lmout = lm(y~.+0, data = lX)
        pred = matrix(data1$testx[,m], dim(data1$testx)[1], sum(m)) %*% coef(lmout)
        mean( (data1$testy - pred)^2 )
      }
      else if(sum(m)==0) {
        mean(data1$testy^2)
      }
      else{
        NA
      }
    }
    
    
    print(paste0(i,"-----iteration"))
    type3traincv = cbind(type3traincv, cv3)
    type3starerr = cbind(type3starerr, fr2[,1])
    

    
  }
  
  return(list(type1traincv, type1starerr,
              type2traincv, type2starerr,
              type3traincv, type3starerr))
  
}

e2_05ch = elastic3type(200, 500, 0.5, lambda_grid = c(seq(8000,6020,-20), seq(6000,0,-2)), fold_number = 2)
e5_01 = elastic3type(200, 500, 0.1, lambda_grid = c(seq(8000,6020,-20), seq(6000,0,-2)), fold_number = 5)
e10_01 = elastic3type(200, 500, 0.1, lambda_grid = c(seq(8000,6020,-20), seq(6000,0,-2)), fold_number = 10)
e20_01 = elastic3type(200, 500, 0.1, lambda_grid = c(seq(8000,6020,-20), seq(6000,0,-2)), fold_number = 20)
e25_01 = elastic3type(200, 500, 0.1, lambda_grid = c(seq(8000,6020,-20), seq(6000,0,-2)), fold_number = 25)
e40_01 = elastic3type(200, 500, 0.1, lambda_grid = c(seq(8000,6020,-20), seq(6000,0,-2)), fold_number = 40)
e50_01 = elastic3type(200, 500, 0.1, lambda_grid = c(seq(8000,6020,-20), seq(6000,0,-2)), fold_number = 50)
e100_01 = elastic3type(200, 500, 0.1, lambda_grid = c(seq(8000,6020,-20), seq(6000,0,-2)), fold_number = 100)
e200_01 = elastic3type(200, 500, 0.1, lambda_grid = c(seq(8000,6020,-20), seq(6000,0,-2)), fold_number = 200)


save(list=ls(), file = "C:/Users/pc/OneDrive/ë°”íƒ• ?™”ë©?/cv/May/elastic200500_01.RData")

