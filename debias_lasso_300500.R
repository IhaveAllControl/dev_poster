library(CDLasso)
library(foreach)
library(doParallel)
cl <- makeCluster(11)
registerDoParallel(cl)

simul_data <- function(n, p,sd=1){
  X = matrix(rnorm(n*p, 0, 1), n, p)
  #B = c(rep(1, 0.05*p), rep(-1, 0.05*p), rep(2,0.05*p),rep(0, 0.85*p))
  c1_be <- rep(0,p)
  c1_be[1:20] = 2
  c1_be[21:25] = -1
  c1_be[26:50] = 3
  
  B = c1_be
  e = rnorm(n, 0, sd)
  y = X %*% B + e
  
  test_X = matrix(rnorm(n*p, 0, 1), n, p)
  test_e = rnorm(n, 0, sd)
  test_y = test_X %*% B + test_e
  return(list(X = scale(X,T,F), y=scale(y,T,F), testx = scale(test_X,T,F), testy = scale(test_y,T,F)))
  
}



setwd("D:/Down/0115_feedback")



cvcompare_no_bias_v2 <- function(n, p, X_data, y_data, X.test, y.test, k, using_data, st=5,len=2,lv=NULL) {
  tol=1e-6;
  cd_iter = 1e6
  niter = 100
  init=rep(0,p)
  if(is.null(lv)) lv <- max(abs(t(dat$X)%*%dat$y))
  
  num.fol <- floor(dim(X_data)[1]/k)
  #lv <- max(abs(t(dat$X)%*%dat$y))
  lam.vec <- seq(from=st, to=lv, by=len)
  nl = length(lam.vec)
  
  res.matrix <- vector()
  for (fold in 1:k) {
    
    if (fold!=k) {
      i.except.train <- (1 + (fold-1) * num.fol):(fold * num.fol)
    }
    else{
      i.except.train <- (1 + (fold-1) * num.fol):dim(X_data)[1]
    }
    
    
    if(k==1)
    {
      print("k=1")
      for(j in 1:nl)
      {
        lam = lam.vec[j]
        if(j==1){
          init = rep(0,p)
        } else {init = temp}
        
        l <-  cd_lasso(X_data,y_data,lam,tol,cd_iter,init)
        temp = l$est
        nz_est = sum(temp!=0)
        if (nz_est == 0) {
          res.matrix[j,k] = res.matrix[j,k]+sum((y.test)^2)
        } else if(nz_est <= n){
          lX <- as.data.frame(X_data[, temp!=0])
          lX[,"y"] <- y_data
          llm <- lm(y~.+0, data = lX)
          res.matrix[j,k] = res.matrix[j,k]+sum((y.test - as.matrix(X.test[,temp!=0]) %*% coef(llm))^2 )
        }
        else {
          res.matrix[j,k] = NA
        }
      }
    } 
    if(using_data == "test") 
    {
      print(paste0("test","_",fold,"/",k))
      for(j in 1:nl)
      {
        lam = lam.vec[j]
        if(j==1){
          init = rep(0,p)
        } else {init = temp}
        
        l <-  cd_lasso(X_data[-i.except.train,],y_data[-i.except.train],lam,tol,cd_iter,init)
        temp = l$est
        nz_est = sum(temp!=0)
        if (nz_est == 0) {
          res.matrix[j,k] = res.matrix[j,k]+sum((y.test)^2)
        } else if(nz_est <= n){
          lX <- as.data.frame(X_data[, temp!=0])
          lX[,"y"] <- y_data
          llm <- lm(y~.+0, data = lX)
          res.matrix[j,k] = res.matrix[j,k]+sum((y.test - as.matrix(X.test[,temp!=0]) %*% coef(llm))^2 )
        }
        else {
          res.matrix[j,k] = NA
        }
      }      
    }
    if(using_data == 'train')
    {
      
      print(paste0("train","_",fold,"/",k))
      fr = foreach(lam = lam.vec, .combine = rbind )%dopar% {
        cd_lasso <- function(X,y,lambda,tol,niter,init,loc="D:/Down/0115_feedback")
        { 
          n <- as.integer(nrow(X))
          p <- as.integer(ncol(X))
          initer <- as.integer(niter)
          iriter <- as.integer(0)
          
          dX <- as.double(X)
          dy <- as.double(y)
          dlam <- as.double(lambda)
          dbe <- as.double(rep(0,p))
          dtol <- as.double(tol)
          ddiff <- as.double(0)
          dres <- as.double(rep(0,n))
          dinit <- as.double(init)
          
          dyn.load(paste0(loc,"/cd_lasso",.Platform$dynlib.ext))
          out <- .C("cd_lasso",n=n,p=p,y=dy,X=dX,beta=dbe,lam=dlam,niter=initer,tol=dtol,riter=iriter,diff=ddiff,resid=dres,init=dinit)
          dyn.unload(paste0(loc,"/cd_lasso",.Platform$dynlib.ext))
          
          return(list(est=out$beta,resid=out$resid,max_diff=out$diff,riter=out$riter))
          
        }
        
        l <-  cd_lasso(X_data[-i.except.train,],y_data[-i.except.train],lam,1e-6,1e6,rep(0,p))
        temp = l$est
        nz_est = sum(temp!=0)
        if (nz_est == 0) {
          sum((y_data[i.except.train])^2)
        } else if(nz_est <= (n - length(i.except.train))){
          lX <- as.data.frame(X_data[-i.except.train, temp!=0])
          lX[,"y"] <- y_data[-i.except.train]
          llm <- lm(y~.+0, data = lX)
          sum((y_data[i.except.train] -matrix(X_data[i.except.train,temp!=0], length(i.except.train), nz_est) %*% coef(llm))^2)
        }
        else {
          NA
        }
      }
    }
    res.matrix <- cbind(res.matrix, fr[,1]) 
  }
  rownames(res.matrix) <- lam.vec
  #res.matrix <- na.omit(res.matrix)
  res.matrix <- apply(res.matrix, 1, mean, na.rm=T)
  best.fit.lam <- as.double(names(res.matrix)[which.min(res.matrix)])
  return(list(lam_s=best.fit.lam, res_mat=res.matrix,lam_vec=lam.vec))
  
}

make_iter = function(n, p, k, using_data,
                     st,len,lv, iteration) {
  cand_alpha = seq(0.5,2,by=0.01)
  result = 0
  result2 = vector()
  set.seed(2021)
  for (i in 1:iteration) {
    dat = simul_data(n,p,sd=0.5)
    X = dat$X
    y = dat$y
    tX = dat$testx
    ty = dat$testy
    # out = cvcomparelm(n, X_data = X, y_data = y, X.test = tX, y.test = ty,
    #                   k, using_data = "train")
    # 
    out = cvcompare_no_bias_v2(n, p, X_data=X, y_data=y, X.test=tX, y.test=ty,
                               k, using_data='train',st=st,len=len,lv=lv)
    result = cbind(result, out[[2]])
    fr = foreach(LXA = cand_alpha * out$lam_s, .combine = rbind) %dopar% {
      cd_lasso <- function(X,y,lambda,tol,niter,init,loc="D:/Down/0115_feedback")
      { 
        n <- as.integer(nrow(X))
        p <- as.integer(ncol(X))
        initer <- as.integer(niter)
        iriter <- as.integer(0)
        
        dX <- as.double(X)
        dy <- as.double(y)
        dlam <- as.double(lambda)
        dbe <- as.double(rep(0,p))
        dtol <- as.double(tol)
        ddiff <- as.double(0)
        dres <- as.double(rep(0,n))
        dinit <- as.double(init)
        
        dyn.load(paste0(loc,"/cd_lasso",.Platform$dynlib.ext))
        out <- .C("cd_lasso",n=n,p=p,y=dy,X=dX,beta=dbe,lam=dlam,niter=initer,tol=dtol,riter=iriter,diff=ddiff,resid=dres,init=dinit)
        dyn.unload(paste0(loc,"/cd_lasso",.Platform$dynlib.ext))
        
        return(list(est=out$beta,resid=out$resid,max_diff=out$diff,riter=out$riter))
        
      }
      l <-  cd_lasso(X,y,LXA,1e-6,1e6,init=rep(0,p))
      temp = l$est
      lX <- as.data.frame(X[, temp!=0])
      lX[,"y"] <- y
      llm <- lm(y~.+0, data = lX)
      sum((ty - as.matrix(tX[,temp!=0]) %*% coef(llm))^2) 
      
    }
    result2 = cbind(result2, fr[,1])
    
    print(paste0(i,"--iteration"))  
  }
  return(list(result[,-1], result2))
}


m2 = make_iter(300,500,2,"train", st = 5, len =0.4, lv =500, iteration = 50)
m5 = make_iter(300,500,5,"train", st = 5, len =0.4, lv =500, iteration = 50)
m10 = make_iter(300,500,10,"train", st = 5, len =0.4, lv =500, iteration = 50)
m15 = make_iter(300,500,15,"train", st = 5, len =0.4, lv =500, iteration = 50)
m20 = make_iter(300,500,20,"train", st = 5, len =0.4, lv =500, iteration = 50)
m25 = make_iter(300,500,25,"train", st = 5, len =0.4, lv =500, iteration = 50)
m30 = make_iter(300,500,30,"train", st = 5, len =0.4, lv =500, iteration = 50)
m50 = make_iter(300,500,50,"train", st = 5, len =0.4, lv =500, iteration = 50)
m60 = make_iter(300,500, 60,"train", st = 5, len =0.4, lv =500, iteration = 50)
m100 = make_iter(300,500,100,"train", st = 5, len =0.4, lv =500, iteration = 50)
m150 = make_iter(300,500,150,"train", st = 5, len =0.4, lv =500, iteration = 50)
m300 = make_iter(300,500,300,"train", st = 5, len =0.4, lv =500, iteration = 50)



# save(list = ls()[c(3:8,10:12, 14,15,16)],file = "C:/Users/pc/OneDrive/¹ÙÅÁ È­¸é/cv/April/debias_lasso_obj300500.RData")
t2 = paste0("2fold:",round(as.numeric(names(which.min(apply(m2[[1]],1,mean)))),2))
t5 = paste0("5fold:",round(as.numeric(names(which.min(apply(m5[[1]],1,mean)))),2))
t10 = paste0("10fold:",round(as.numeric(names(which.min(apply(m10[[1]],1,mean)))),2))
t20 = paste0("20fold:",round(as.numeric(names(which.min(apply(m20[[1]],1,mean)))),2))
t25 = paste0("25fold:",round(as.numeric(names(which.min(apply(m25[[1]],1,mean)))),2))
t50 = paste0("50fold:",round(as.numeric(names(which.min(apply(m50[[1]],1,mean)))),2))
t100 = paste0("100fold:",round(as.numeric(names(which.min(apply(m100[[1]],1,mean)))),2))

plot(seq(5, 500, 0.4),apply(m2[[1]],1,mean), type = 'l', yaxt = 'n', ylab = "", lwd =2, lty =2, 
     xlab = expression(lambda))
par(new = T)
plot(seq(5, 500, 0.4),apply(m5[[1]],1,mean), type = 'l', yaxt = 'n', ylab = "", lwd =2, lty =2, 
     xlab = expression(lambda), xaxt = 'n', col ="red")
par(new = T)

plot(seq(5, 500, 0.4),apply(m10[[1]],1,mean), type = 'l', yaxt = 'n', ylab = "", lwd =2, lty =2, 
     xlab = expression(lambda), xaxt = 'n',col = "blue")
par(new = T)

plot(seq(5, 500, 0.4),apply(m20[[1]],1,mean), type = 'l', yaxt = 'n', ylab = "", lwd =2, lty =2, 
     xlab = expression(lambda), xaxt = 'n', col = "orange")
par(new = T)

plot(seq(5, 500, 0.4),apply(m25[[1]],1,mean), type = 'l', yaxt = 'n', ylab = "", lwd =2, lty =2,
     xlab = expression(lambda), xaxt = 'n', col = "purple")
par(new = T)
plot(seq(5, 500, 0.4),apply(m50[[1]],1,mean), type = 'l', yaxt = 'n', ylab = "", lwd =2, lty =2, 
     xlab = expression(lambda), xaxt = 'n', col = "hotpink")
par(new = T)
plot(seq(5, 500, 0.4),apply(m100[[1]],1,mean), type = 'l', yaxt = 'n', ylab = "", lwd =2, lty =2,
     xlab = expression(lambda), xaxt = 'n', col = "gray40")
legend("bottomright",legend = paste0(c(2,5,10,20,25, 50, 100),"fold"),
       fill =c("black", "red" , "blue", "orange", "purple", "hotpink",
               "gray40"), bty = 'n')



minalpha2= seq(0.5, 2, 0.01)[which.min(apply(m2[[2]],1,mean))]
minalpha5= seq(0.5, 2, 0.01)[which.min(apply(m5[[2]],1,mean))]
minalpha10= seq(0.5, 2, 0.01)[which.min(apply(m10[[2]],1,mean))]
minalpha20= seq(0.5, 2, 0.01)[which.min(apply(m20[[2]],1,mean))]
minalpha25= seq(0.5, 2, 0.01)[which.min(apply(m25[[2]],1,mean))]
minalpha50= seq(0.5, 2, 0.01)[which.min(apply(m50[[2]],1,mean))]
minalpha100= seq(0.5, 2, 0.01)[which.min(apply(m100[[2]],1,mean))]


plot(seq(0.5, 2, 0.01),apply(m2[[2]],1,mean), type ='l', ylab = "test_error(using lm)",
     xlab = expression(alpha), lwd =2, lty = 2, yaxt = 'n')
par(new = T)
plot(seq(0.5, 2, 0.01),apply(m5[[2]],1,mean), type ='l', 
     xaxt = 'n', xlab = "", ylab = "", lwd =2, lty = 2, yaxt ='n', col ='red')
par(new =T)
plot(seq(0.5, 2, 0.01),apply(m10[[2]],1,mean), type ='l', xaxt = 'n', 
     xlab = "", ylab = "", lwd =2, lty = 2, yaxt ='n', col = "blue")
par(new =T)
plot(seq(0.5, 2, 0.01),apply(m20[[2]],1,mean), type ='l', xaxt = 'n', 
     xlab = "", ylab = "", lwd =2, lty = 2, yaxt ='n', col = "orange")
par(new =T)
plot(seq(0.5, 2, 0.01),apply(m25[[2]],1,mean), type ='l', xaxt = 'n', 
     xlab = "", ylab = "", lwd =2, lty = 2, yaxt ='n', col = "purple")
par(new =T)
plot(seq(0.5, 2, 0.01),apply(m50[[2]],1,mean), type ='l', xaxt = 'n', 
     xlab = "", ylab = "", lwd =2, lty = 2, yaxt ='n', col = "hotpink")
par(new =T)
plot(seq(0.5, 2, 0.01),apply(m100[[2]],1,mean), type ='l', xaxt = 'n', 
     xlab = "", ylab = "", lwd =2, lty = 2, yaxt ='n', col = "gray40")
legend(0.5,26900,legend = paste0(c(2, 5, 10, 20, 25, 50,100),"fold"),
       fill =c("black", "red" , "blue", "orange", "purple", "hotpink",
               "gray40", "cadetblue3"), bty = 'n')

plot(1:7, c(minalpha2, minalpha5, minalpha10, minalpha20, minalpha25, minalpha50, 
            minalpha100),
     pch = 1, cex =2, type = 'b', ylab = expression(alpha), xlab = "", 
     xaxt = 'n',
     main = "min alpha of each fold ")
axis(1, at=1:7, labels = paste0(c(2,5,10,20,25,50,100),"fold"))
text(1:7, c(minalpha2, minalpha5, minalpha10, minalpha20, minalpha25, minalpha50, 
            minalpha100)+0.005,
     labels = c(minalpha2, minalpha5, minalpha10, minalpha20, minalpha25, minalpha50, 
                minalpha100))
legend(0.5,0.7, legend = c("data_num_ratio","2fold : 2", "5fold : 1.25", "10fold : 1.11",
                           "20fold : 1.052", "25fold : 1.041",
                           "50fold : 1.02", "100fold : 1.01"), bty = 'n')


mbl2 = round(mean(apply(m2[[2]],2,function(x) seq(0.5, 2, 0.01)[which.min(x)])),2)
mbl5 = round(mean(apply(m5[[2]],2,function(x) seq(0.5, 2, 0.01)[which.min(x)])),2)
mbl10 = round(mean(apply(m10[[2]],2,function(x) seq(0.5, 2, 0.01)[which.min(x)])),2)
mbl20 = round(mean(apply(m20[[2]],2,function(x) seq(0.5, 2, 0.01)[which.min(x)])),2)
mbl25 = round(mean(apply(m25[[2]],2,function(x) seq(0.5, 2, 0.01)[which.min(x)])),2)
mbl50 = round(mean(apply(m50[[2]],2,function(x) seq(0.5, 2, 0.01)[which.min(x)])),2)
mbl100 = round(mean(apply(m100[[2]],2,function(x) seq(0.5, 2, 0.01)[which.min(x)])),2)

plot(1:7, c(mbl2, mbl5, mbl10, mbl20, mbl25, mbl50, mbl100),
     pch = 16, cex =2, type = 'b', ylab = expression(alpha), xlab = "", 
     xaxt = 'n',
     main = "min alpha of each fold for columns apply data")
axis(1, at=1:7, labels = paste0(c(2,5,10,20,25,50,100),"fold"))
legend(2,1, legend = c("data_num_ratio","2fold : 2", "5fold : 1.25", "10fold : 1.11",
                       "20fold : 1.053", "25fold : 1.041",
                       "50fold : 1.02", "100fold : 1.01"), bty = 'n')
text(1:7,c(mbl2, mbl5, mbl10, mbl20, mbl25, mbl50, mbl100) +0.009,
     labels = c(mbl2, mbl5, mbl10, mbl20, mbl25, mbl50, mbl100),col = 'blue')



box_dat = as.data.frame(cbind(apply(m2[[2]],2,function(x) seq(0.5, 2, 0.01)[which.min(x)]),
                              apply(m5[[2]],2,function(x) seq(0.5, 2, 0.01)[which.min(x)]),
                              apply(m10[[2]],2,function(x) seq(0.5, 2, 0.01)[which.min(x)]),
                              apply(m20[[2]],2,function(x) seq(0.5, 2, 0.01)[which.min(x)]),
                              apply(m25[[2]],2,function(x) seq(0.5, 2, 0.01)[which.min(x)]),
                              apply(m50[[2]],2,function(x) seq(0.5, 2, 0.01)[which.min(x)]),
                              apply(m100[[2]],2,function(x) seq(0.5, 2, 0.01)[which.min(x)])))

boxplot(box_dat, xaxt = 'n', main = "min alphas for 50 iterations")
text(1:7,c(mbl2, mbl5, mbl10, mbl20, mbl25, mbl50, mbl100),
     labels = paste0(c(mbl2, mbl5, mbl10, mbl20, mbl25, mbl50, mbl100), "(mean)"),col = 'blue')
par(xpd = TRUE)
text(1:7,0.4,
     labels = c("2fold : 2", "5fold : 1.25", "10fold : 1.11",
                "20fold : 1.053", "25fold : 1.041",
                "50fold : 1.02", "100fold : 1.01")
     ,col = 'blue')

par(xpd = FALSE)

