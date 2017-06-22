# Bayesian models to canlculate LVB parameters
# June 8 2017
# BLM

usePackage("R2jags")
usePackage("snow")
usePackage("dclone")

# ============================= Models for no CWT to Scale Age Calculation ============================== #####
# CWT == FALSE

# =================== Ludwig von Bertalanffy Length growth model ========================

Run_LVB <- function(Data){
  attach(Data)
  
  # Model Specs
  n.chains=2
  n.burnin=5000
  n.iter=10000
  n.thin=4
  
  
  LVBmodel <- function(){
    for (i in 1:N) {
      Lmm[i] ~ dnorm(eta[i], prec)
      eta[i] <- Linf * (1-exp(-k * (Age[i] - t0)))
      Rep[i] ~ dnorm(eta[i], prec)
    }
    
    # Prior for the error precision (=1/Sigma^2)
    prec ~ dgamma(1.0E-4, 1.0E-4)
    sigma2 <- 1 /prec
    
    Linf ~ dunif(0, 10000)
    k     ~ dunif(0, 100)     
    t0    ~ dunif(-1000, 1000)
  }
  
  
  # Run vbStarts for get reasonable initial values
  svTypical <- as.vector(unlist(vbStarts( Lmm~Age )))
  
  print(svTypical)
  
  Inits <- function(Inits=svTypical){
    list(Linf=Inits[1], k=Inits[2], t0=Inits[3], prec=1)
    list(Linf=Inits[1]-100, k=Inits[2]-0.1, t0=-Inits[3], prec=1)
  }
  
  data <- c("N", "Age", "Lmm")
  
  parameters = c("Linf","k","t0","sigma2")  
  
  
  jags.LVB = jags.parallel(data,Inits,parameters,LVBmodel,n.chains=n.chains,n.iter=n.iter,n.burnin=n.burnin,n.thin=n.thin,DIC=T)    
  detach(Data)
  return(jags.LVB)
  
} # End Jags code


# =================== LOG ERROR ***** Ludwig von Bertalanffy Length growth model ========================

Run_LOG_LVB <- function(Data,Inits){
  Data$logLmm <- log(Data$Lmm)
  attach(Data)
  
  # Model Specs
  n.chains=2
  n.burnin=5000
  n.iter=10000
  n.thin=4
  N=length(Data$Lmm)
  
  LVBmodel <- function(){
    for (i in 1:N) {
      logLmm[i] ~ dnorm(eta[i], prec)
      eta[i] <- log(Linf * (1-exp(-k * (Age[i] - t0))))
      # Rep[i] ~ dnorm(eta[i], prec)
    }
    
    # Prior for the error precision (=1/Sigma^2)
    prec ~ dgamma(1.0E-4, 1.0E-4)
    sigma2 <- 1 /prec
    
    Linf ~ dunif(1, 10000)
    k     ~ dunif(0.001, 100)     
    t0    ~ dunif(-1000, 2)
  }
  
  
  # Run vbStarts for get reasonable initial values
  svTypical <- as.vector(unlist(vbStarts( Lmm~Age )))
  
  print(svTypical)
  
  Inits <- function(Inits=svTypical){
    list(Linf=Inits[1], k=Inits[2], t0=Inits[3], prec=1)
    list(Linf=Inits[1]-100, k=Inits[2]-0.1, t0=-Inits[3], prec=1)
  }
  #   
  data <- c("N", "Age", "logLmm")
  
  parameters = c("Linf","k","t0","sigma2")  
  
  
  jags.LVB = jags.parallel(data,Inits,parameters,LVBmodel,n.chains=n.chains,n.iter=n.iter,n.burnin=n.burnin,n.thin=n.thin,DIC=T)    
  detach(Data)
  return(jags.LVB)
  
} # End Jags code


# ========================== Normal LVB Model with stock index =======================

Run_Stock_LVB <- function(Data){
  nPop = length(unique(Data$Pop))
  attach(Data)
  N <- N[1]
  
  # Model Specs
  n.chains=2
  n.burnin=5000
  n.iter=10000
  n.thin=4
  
  
  LVBmodel <- function(){
    for (i in 1:N) {
      Lmm[i] ~ dnorm(eta[i], prec[POP_ID[i]])
      eta[i] <- Linf[POP_ID[i]] * (1-exp(-k[POP_ID[i]] * (Age[i] - t0[POP_ID[i]])))
      Rep[i] ~ dnorm(eta[i], prec[POP_ID[i]])
    }
    
    for(i in 1:nPop)
    {
      Linf[i] ~ dunif(0, 10000)
      k[i]    ~ dunif(0, 100) 
      t0[i]   ~ dunif(-1000, 1000)
      prec[i] ~ dgamma(1.0E-4, 1.0E-4)
      sigma2[i]<-1/prec[i]
    } 
  }
  
  # Run vbStarts for get reasonable initial values
  svTypical <<- as.vector(unlist(vbStarts( Lmm~Age )))
  
  print(svTypical)
  
  Inits <- function(svTypical.=svTypical){
    list(
      Linf=rep(unlist(svTypical[[1]]), times=nPop), 
      k=rep(unlist(svTypical[[2]]),times=nPop), 
      t0=rep(unlist(svTypical[[3]]),times=nPop), 
      prec=rep(1,times=nPop)
    )
    list(
      Linf=rep(unlist(svTypical[[1]]), times=nPop), 
      k=rep(unlist(svTypical[[2]]),times=nPop), 
      t0=rep(unlist(svTypical[[3]]),times=nPop), 
      prec=rep(1,times=nPop)
    )
  }
  
  data <- c("N", "Age", "Lmm", "nPop", "POP_ID")
  
  parameters = c("Linf","k","t0","sigma2")  
  
  jags.LVB = jags(data,Inits,parameters,LVBmodel,n.chains=n.chains,n.iter=n.iter,n.burnin=n.burnin,n.thin=n.thin,DIC=T) 
  #jags.LVB = do.call(jags.parallel, 
  #                   list(data=data,inits=Inits,parameters.to.save=parameters,model.file=LVBmodel,
  #                        n.chains=n.chains,n.iter=n.iter,n.burnin=n.burnin,n.thin=n.thin,DIC=T)
  #)
  detach(Data)
  return(jags.LVB)
  
} # End Jags code


# ========================== Normal LVB Model with stock index =======================

Run_Stock_LOG_LVB <- function(Data){
  nPop = length(unique(Data$Pop))
  Data$logLmm <- log(Data$Lmm)
  attach(Data)
  N <- length(logLmm)
  
  # Model Specs
  n.chains=2
  # n.burnin=5000
  n.adapt=250
  n.update=250 
  n.iter=2000
  n.thin=4
  
  
  LVBmodel <- function(){
    for (i in 1:N) {
      logLmm[i] ~ dnorm(eta[i], prec[POP_ID[i]])
      eta[i] <- log(Linf[POP_ID[i]] * (1-exp(-k[POP_ID[i]] * (Age[i] - t0[POP_ID[i]]))))
      #Rep[i] ~ dnorm(eta[i], prec[POP_ID[i]])
    }
    
    for(i in 1:nPop)
    {
      Linf[i] ~ dunif(1, 10000)
      k[i]    ~ dunif(0.001, 100) 
      t0[i]   ~ dunif(-1000, 1.9)
      prec[i] ~ dgamma(1.0E-4, 1.0E-4)
      sigma2[i]<-1/prec[i]    
    } 
  }

  # Run vbStarts for get reasonable initial values
  sv <<- as.vector(unlist(vbStarts( Lmm~Age )))
  
  print(svTypical)
  
  Inits <- function(svTypical=sv){
    list(
      Linf=rep(unlist(svTypical[[1]]), times=nPop), 
      k=rep(unlist(svTypical[[2]]),times=nPop), 
      t0=rep(unlist(svTypical[[3]]),times=nPop), 
      prec=rep(1,times=nPop)
    )
    list(
      Linf=rep(unlist(svTypical[[1]]), times=nPop), 
      k=rep(unlist(svTypical[[2]]),times=nPop), 
      t0=rep(unlist(svTypical[[3]]),times=nPop), 
      prec=rep(1,times=nPop)
    )
  }
  
  data <- list(N=N, Age=Age, logLmm=logLmm, nPop=nPop, POP_ID=POP_ID)
  
  parameters = c("Linf","k","t0","sigma2")  

  # using dclone
  cl <- makeCluster(n.chains, "SOCK")
  jags.LVB <- jags.parfit(cl, data = data, params = parameters, model = LVBmodel,
                    inits = Inits, n.adapt =n.adapt, n.update = n.update,
                    n.iter = n.iter, thin = n.thin, n.chains = n.chains)
  stopCluster(cl)
  
 # jags.LVB = jags(data,Inits,parameters,LVBmodel,n.chains=n.chains,n.iter=n.iter,n.burnin=n.burnin,n.thin=n.thin,DIC=TRUE) 
  #jags.LVB = do.call(jags.parallel, 
  #                   list(data=data,inits=Inits,parameters.to.save=parameters,model.file=LVBmodel,
  #                        n.chains=n.chains,n.iter=n.iter,n.burnin=n.burnin,n.thin=n.thin,DIC=T)
  #)
  detach(Data)
  return(jags.LVB)
  
} # End Jags code



# ======================= End models for no CWT to Scale Age Calculation





# ==================== Models WITH CWT to Scale Age Calculation ======================================== ####
# CWT == TRUE
# Calculating the thetas

# =================== Ludwig von Bertalanffy Length growth model ========================

Run_LVB_Age_Fit <- function(Data){
  attach(Data)
  
  # Model Specs
  n.chains=2
  n.adapt=50
  n.update=50   # n.adapt+n.update=burnin
  n.iter=100
  n.thin=4
  
  N=length(Data$Lmm)
  nAge <- nK <- max( length(1:max(unique(True),na.rm=T)), length(1:max(unique(Age),na.rm=T)) )
  alpha <- rep(1, length=nK)
  
  LVBmodelAge <- function(){
    for (i in 1:N) {
      Lmm[i] ~ dnorm(eta[i], prec)
      eta[i] <- Linf * (1-exp(-k * (True[i] - t0)))
      True[i] ~ dcat(theta[1:nK,(Age[i])]) 
      #Rep[i] ~ dnorm(eta[i], prec)
    }
    
    for(k in 1:nAge){  
      theta[1:nK,k]~ ddirch(alpha[])          #prior for age model
    } 
    # Prior for the error precision (=1/Sigma^2)
    prec ~ dgamma(1.0E-4, 1.0E-4)
    sigma2 <- 1 /prec
    
    Linf ~ dunif(1, 10000)
    k ~ dunif(0.1, 100) 
    t0 ~ dunif(-1000, 0)
  }
  
  
  # Run vbStarts for get reasonable initial values
  svTypical <- as.vector(unlist(vbStarts( Lmm~Age )))
  
  print(svTypical)
  
  Inits <- function(Inits=svTypical){
    list(Linf=Inits[1], k=Inits[2], t0=Inits[3], prec=1, theta=structure(.Data=rep(1/nK, length=nK*nK),.Dim=c(nK,nK)))
    list(Linf=Inits[1]-100, k=Inits[2]-0.1, t0=-Inits[3], prec=1, theta=structure(.Data=rep(1/nK, length=nK*nK),.Dim=c(nK,nK)))
  }
  
  # data <- c("N", "Age", "Lmm", "True", "nAge", "nK","alpha")
  
  data <- list(N=N, Age=Age, Lmm=Lmm, True=True, nAge=nAge, nK=nK,alpha=alpha)
  
  parameters = c("Linf","k","t0","theta","sigma2")  #True
  
  # using dclone
  cl <- makeCluster(n.chains, "SOCK")
  m2 <- jags.parfit(cl, data = data, params = parameters, model = LVBmodelAge,
                    inits = Inits, n.adapt =n.adapt, n.update = n.update,
                    n.iter = n.iter, thin = n.thin, n.chains = n.chains)
  stopCluster(cl)
  
  #Run n.chains for n.iter, with one chain on each of n.chains cores
  #Need a wrapper for the fitting process. 
  #1) The snow function 'clusterApply' requires functions with just one argument. 
  #2) JAGS requires that both jags.model and coda.samples be run under a single call to clusterApply.
  ##Note that it's important to specify a different random number seed for each worker.
  # coda.samples.wrapper <- function(j)
  # { 
  #   temp.model = jags.model("model.txt", 
  #                           inits=Inits, 
  #                           data=data, n.chains=n.chains, quiet=FALSE)
  #   coda.samples(temp.model, parameters, n.iter=n.iter, thin=n.thin) 
  # }
  # 
  
  # inits2 <- jags.fit(data, parameters, model=LVBmodelAge, inits=Inits, n.chains, 
  #                    n.adapt = 0, n.update = 0, n.iter = 0)$state(internal = TRUE)   # set the updates to 0 for this
  # clusterEvalQ(cl, library(dclone))
  # clusterEvalQ(cl, setwd(getwd()))
  # 
  # filename <- write.jags.model(LVBmodelAge, overwrite=TRUE)
  # 
  # cldata <- list(data=data, params=parameters, model=filename, inits=inits2)
  # clusterExport(cl, "cldata")
  # 
  # jagsparallel <- function(i, ...) {
  #   jags.fit(data = cldata$data, params = cldata$params, model = cldata$model,
  #            inits = cldata$inits[[i]], n.chains = 1, updated.model = FALSE, ...)
  #   }
  # res <- parLapply(cl, 1:n.chains, jagsparallel, n.adapt = 10, n.update = 10, n.iter = 30, thin = thin)
  # 
  # res <- as.mcmc.list(lapply(res, as.mcmc))
  # 
  # list(inits2[[1]])
  
  ##Make sure the rjags library is loaded in each worker
  # clusterEvalQ(cl, library(rjags))
  # ##Send data to workers, then fit models. One disadvantage of this parallelization is that you lose the ability to watch the progress bar.
  # clusterExport(cl, list(data,"n.iter","n.thin"))
  # par.samples = clusterApply(cl, 1:n.chains, coda.samples.wrapper)
  # ##Reorganize 'par.samples' so that it is recognizeable as an 'mcmc.list' object
  # for(i in 1:length(par.samples)) { par.samples[[i]] <- par.samples[[i]][[1]] }
  # class(par.samples) <- "mcmc.list"
  # stopCluster(cl)
  # snow.end.time = proc.time()
  # snow.dtime = snow.end.time - snow.start.time
  
  
  
  #jags.LVB = jags(data,Inits,parameters,LVBmodelAge,n.chains=n.chains,n.iter=n.iter,n.burnin=n.burnin,n.thin=n.thin,DIC=T)    
  detach(Data)
  return(m2)
  
} # End Jags code


# =================== LOG ERROR ***** Ludwig von Bertalanffy Length growth model ========================

Run_LOG_LVB_Age_Fit <- function(Data,Inits){
  Data$logLmm <- log(Data$Lmm)
  attach(Data)
  
  N=length(Data$Lmm)
  nAge <- nK <- max( length(1:max(unique(True),na.rm=T)), length(1:max(unique(Age),na.rm=T)) )
  alpha <- rep(1, length=nK)
  
  # Model Specs
  n.chains=2
  # n.burnin=5000
  n.adapt=50
  n.update=50   # n.adapt+n.update=burnin
  n.iter=100
  n.thin=4
  
  N=length(Data$Lmm)
  nAge <- nK <- max( length(1:max(unique(True),na.rm=T)), length(1:max(unique(Age),na.rm=T)) )
  alpha <- rep(1, length=nK)
  
  logLVBmodelAge <- function(){
    for (i in 1:N) {
      logLmm[i] ~ dnorm(eta[i], prec)
      eta[i] <- log(Linf * (1-exp(-k * (True[i] - t0))))
      True[i] ~ dcat(theta[1:nK,(Age[i])]) 
      # Rep[i] ~ dnorm(eta[i], prec)
    }
    for(k in 1:nAge){  
      theta[1:nK,k]~ ddirch(alpha[])          #prior for age model
    } 
    # Prior for the error precision (=1/Sigma^2)
    prec ~ dgamma(1.0E-4, 1.0E-4)
    # sigma2 <- 1/prec
    
    Linf ~ dunif(1, 10000)
    k     ~ dunif(0.001, 100)     
    t0    ~ dunif(-1000, 2)
  }
  
  # Run vbStarts for get reasonable initial values
  svTypical <- as.vector(unlist(vbStarts( Lmm~Age )))
  
  print(svTypical)
  
  Inits <- function(Inits=svTypical){
    list(Linf=Inits[1], k=Inits[2], t0=Inits[3], prec=.001, theta=structure(.Data=rep(1/nK, length=nK*nK),.Dim=c(nK,nK)))
    list(Linf=Inits[1]-100, k=Inits[2]-0.1, t0=-Inits[3], prec=.0001, theta=structure(.Data=rep(1/nK, length=nK*nK),.Dim=c(nK,nK)))
  }
  
  # data <- c("N", "Age", "logLmm", "True", "nAge", "nK","alpha")
  
  data <- list(N=N, Age=Age, logLmm=logLmm, True=True, nAge=nAge, nK=nK,alpha=alpha)
  
  parameters = c("Linf","k","t0","theta","sigma2")  #True
  
  # using dclone
  cl <- makeCluster(n.chains, "SOCK")
  m2 <- jags.parfit(cl, data = data, params = parameters, model = logLVBmodelAge,
                    inits = Inits, n.adapt =n.adapt, n.update = n.update,
                    n.iter = n.iter, thin = n.thin, n.chains = n.chains)
  stopCluster(cl)
  
  #jags.LVB = jags(data,Inits,parameters,LVBmodel,n.chains=n.chains,n.iter=n.iter,n.burnin=n.burnin,n.thin=n.thin,DIC=T)    
  detach(Data)
  return(m2)
  
} # End Jags code


# ========================== Normal LVB Model with stock index =======================

Run_Stock_LVB_Age_Fit <- function(Data){
  nPop = length(unique(Data$Pop))
  N=length(Data$Lmm)
  
  nAge <- nK <-   max(length(1:max(unique(Data$True),na.rm=T)), length(1:max(unique(Data$Age),na.rm=T)) )
  alpha <- rep(1, length=nK)
  # for(Pop in 1:nPop){
  #   nAge[,,Pop] <- max(length(1:max(unique(Data$True[Data$POP_ID==Pop]),na.rm=T)), length(1:max(unique(Data$Age[Data$POP_ID==Pop]),na.rm=T)) )
  # }
  
  attach(Data)
  
  #alpha <- array(apply(nAge,3, rep,x=1), dim=c(1,1,1))  
  
  # Model Specs
  n.chains=2
  # n.burnin=5000
  n.adapt=50
  n.update=50 
  n.iter=10000
  n.thin=4
  
  
  LVBmodelPOP <- function(){
    for (i in 1:N) {
      Lmm[i] ~ dnorm(eta[i], prec[POP_ID[i]])
      eta[i] <- Linf[POP_ID[i]] * (1-exp(-k[POP_ID[i]] * (True[i] - t0[POP_ID[i]])))
      #Rep[i] ~ dnorm(eta[i], prec[POP_ID[i]])
      True[i] ~ dcat(theta[1:nK,(Age[i]),POP_ID[i]])
    }
    
    for(i in 1:nPop)
    {
      for(k in 1:nAge){  
        theta[1:nK,k,i]~ ddirch(alpha[])          #prior for age model
      } 
      Linf[i] ~ dunif(0, 10000)
      k[i]    ~ dunif(0, 100) 
      t0[i]   ~ dunif(-1000, 1000)
      prec[i] ~ dgamma(1.0E-4, 1.0E-4)
      sigma2[i]<-1/prec[i]
    } 
  }
  
  # Run vbStarts for get reasonable initial values
  svTypical <<- as.vector(unlist(vbStarts( Lmm~Age )))
  
  print(svTypical)
  
  Inits <- function(svTypical.=svTypical){
    list(
      Linf=rep(unlist(svTypical[[1]]), times=nPop), 
      k=rep(unlist(svTypical[[2]]),times=nPop), 
      t0=rep(unlist(svTypical[[3]]),times=nPop), 
      prec=rep(1,times=nPop),
      theta=structure(.Data=rep(1/nK,length=nK*nK*nPop),.Dim=c(nK,nK,nPop))
    )
    list(
      Linf=rep(unlist(svTypical[[1]]), times=nPop), 
      k=rep(unlist(svTypical[[2]]),times=nPop), 
      t0=rep(unlist(svTypical[[3]]),times=nPop), 
      prec=rep(1,times=nPop),
      theta=structure(.Data=rep(1/nK,length=nK*nK*nPop),.Dim=c(nK,nK,nPop))
      #theta=array(apply(nAge,3, function(nK){structure(.Data=rep(1/nK,each=nK*nK),.Dim=c(nK,nK))}))
    )
  }
  
  # data <- c("N", "Age", "Lmm", "True", "nAge", "nK","alpha", "nPop", "POP_ID")
  
  data <- list(N=N, Age=Age, Lmm=Lmm, True=True, nAge=nAge, nK=nK, alpha=alpha, nPop=nPop, POP_ID=POP_ID)
  
  parameters = c("Linf","k","t0","theta","sigma2")  #True
  
  # using dclone
  cl <- makeCluster(n.chains, "SOCK")
  m2 <- jags.parfit(cl, data = data, params = parameters, model = LVBmodelPOP,
                    inits = Inits, n.adapt =n.adapt, n.update = n.update,
                    n.iter = n.iter, thin = n.thin, n.chains = n.chains)
  stopCluster(cl)
  
  #jags.LVB = jags(data,Inits,parameters,LVBmodel,n.chains=n.chains,n.iter=n.iter,n.burnin=n.burnin,n.thin=n.thin,DIC=T) 
  #jags.LVB = do.call(jags.parallel, 
  #                   list(data=data,inits=Inits,parameters.to.save=parameters,model.file=LVBmodel,
  #                        n.chains=n.chains,n.iter=n.iter,n.burnin=n.burnin,n.thin=n.thin,DIC=T)
  #)
  detach(Data)
  #return(jags.LVB)
  return(m2)
} # End Jags code


# ========================== FITTED LOGNormal LVB Model with stock index =======================

Run_Stock_LOG_LVB_Age_Fit <- function(Data){
  nPop = length(unique(Data$POP_ID))
  Data$logLmm <- log(Data$Lmm)
  N=length(Data$logLmm)
  
  nAge <- nK <-   max(length(1:max(unique(Data$True),na.rm=T)), length(1:max(unique(Data$Age),na.rm=T)) )
  alpha <- rep(1, length=nK)
  attach(Data)
  
  
  # Model Specs
  n.chains=2
  # n.burnin=5000
  n.adapt=250
  n.update=250
  n.iter=2000
  n.thin=4
  
  
  LVBLOGPOPmodel <- function(){
    for (i in 1:N) {
      logLmm[i] ~ dnorm(eta[i], prec[POP_ID[i]])
      eta[i] <- log(Linf[POP_ID[i]] * (1-exp(-k[POP_ID[i]] * (Age[i] - t0[POP_ID[i]]))))
      Rep[i] ~ dnorm(eta[i], prec[POP_ID[i]])
      True[i] ~ dcat(theta[1:nK,(Age[i]),POP_ID[i]])
    }
    
    for(i in 1:nPop)
    {
      for(k in 1:nAge){  
        theta[1:nK,k,i]~ ddirch(alpha[])          #prior for age model
      } 
      Linf[i] ~ dunif(0, 10000)
      k[i]    ~ dunif(0.001, 100) 
      t0[i]   ~ dunif(-1000, 2)
      prec[i] ~ dgamma(1.0E-4, 1.0E-4)
      sigma2[i]<-1/prec[i]    
    } 
  }
  
  # Run vbStarts for get reasonable initial values
  svTypical <<- as.vector(unlist(vbStarts( Lmm~Age )))
  
  print(svTypical)
  Inits <- function(svTypical.=svTypical){
    list(
      Linf=rep(unlist(svTypical[[1]]), times=nPop), 
      k=rep(unlist(svTypical[[2]]),times=nPop), 
      t0=rep(unlist(svTypical[[3]]),times=nPop), 
      prec=rep(1,times=nPop),
      theta=structure(.Data=rep(1/nK,length=nK*nK*nPop),.Dim=c(nK,nK,nPop))
    )
    list(
      Linf=rep(unlist(svTypical[[1]]), times=nPop), 
      k=rep(unlist(svTypical[[2]]),times=nPop), 
      t0=rep(unlist(svTypical[[3]]),times=nPop), 
      prec=rep(1,times=nPop),
      theta=structure(.Data=rep(1/nK,length=nK*nK*nPop),.Dim=c(nK,nK,nPop))
      #theta=array(apply(nAge,3, function(nK){structure(.Data=rep(1/nK,each=nK*nK),.Dim=c(nK,nK))}))
    )
  }
  
  # data <- c("N", "Age", "logLmm", "True", "nAge", "nK","alpha", "nPop", "POP_ID")
  
  data <- list(N=N, Age=Age, logLmm=logLmm, True=True, nAge=nAge, nK=nK, alpha=alpha, nPop=nPop, POP_ID=POP_ID)
  
  parameters = c("Linf","k","t0","theta","sigma2")  #True
  
  # using dclone
  cl <- makeCluster(n.chains, "SOCK")
  m2 <- jags.parfit(cl, data = data, params = parameters, model = LVBLOGPOPmodel,
                    inits = Inits, n.adapt =n.adapt, n.update = n.update,
                    n.iter = n.iter, thin = n.thin, n.chains = n.chains)
  stopCluster(cl)
  
  
  #jags.LVB = jags(data,Inits,parameters,LVBmodel,n.chains=n.chains,n.iter=n.iter,n.burnin=n.burnin,n.thin=n.thin,DIC=T) 
  #jags.LVB = do.call(jags.parallel, 
  #                   list(data=data,inits=Inits,parameters.to.save=parameters,model.file=LVBmodel,
  #                        n.chains=n.chains,n.iter=n.iter,n.burnin=n.burnin,n.thin=n.thin,DIC=T)
  #)
  detach(Data)
  # return(jags.LVB)
  return(m2)
} # End Jags code

# ======================= PREDICTION OF AGE CORRECTION ========================= ####

# ========================== Predict Normal LVB Model with stock index =======================

Run_Stock_LOG_LVB_Age_Pred <- function(Data,alpha, parameters){
  nPop = length(unique(Data$Pop))
  Data$logLmm <- log(Data$Lmm)
  N=length(Data$logLmm)
  
  nK <- nAge <-  6 # manual entry for max age classes
  attach(Data)
  
  
  # Model Specs
  n.chains=2
  # n.burnin=5000
  n.adapt=250
  n.update=250 
  n.iter=2000
  n.thin=4
  
  
  LVBLOGPOPpred <- function(){
    for (i in 1:N) {
      logLmm[i] ~ dnorm(eta[i], prec[POP_ID[i]])
      eta[i] <- log(Linf[POP_ID[i]] * (1-exp(-k[POP_ID[i]] * (Age[i] - t0[POP_ID[i]]))))
      Rep[i] ~ dnorm(eta[i], prec[POP_ID[i]])
      #True[i] ~ dcat(theta[1:nK,(Age[i]),POP_ID[i]])
     True[i] ~ dcat(alpha[1:nK,(Age[i]),POP_ID[i]])
    }
    
    for(i in 1:nPop)
    {
      #for(k in 1:nAge){  
      #  theta[1:nK,k,i]~ ddirch(alpha[1:nK,k,i])          #prior for age model
      #} 
      Linf[i] ~ dunif(0, 10000)
      k[i]    ~ dunif(0.001, 100) 
      t0[i]   ~ dunif(-1000, 2)
      prec[i] ~ dgamma(1.0E-4, 1.0E-4)
      sigma2[i]<-1/prec[i]    
    } 
  }
  
  # Run vbStarts for get reasonable initial values
  svTypical <<- as.vector(unlist(vbStarts( Lmm~Age )))
  
  print(svTypical)
  Inits <- function(svTypical.=svTypical){
    list(
      Linf=rep(unlist(svTypical[[1]]), times=nPop), 
      k=rep(unlist(svTypical[[2]]),times=nPop), 
      t0=rep(unlist(svTypical[[3]]),times=nPop), 
      prec=rep(1,times=nPop),
      True=rep(2, length=N)
    )
    list(
      Linf=rep(unlist(svTypical[[1]]), times=nPop), 
      k=rep(unlist(svTypical[[2]]),times=nPop), 
      t0=rep(unlist(svTypical[[3]]),times=nPop), 
      prec=rep(1,times=nPop),
      True=rep(3,length=N)
    )
  }
  
  # data <- c("N", "Age", "logLmm", "True", "nAge", "nK","alpha", "nPop", "POP_ID")
  
  data <- list(N=N, Age=Age, logLmm=logLmm, nK=nK, nPop=nPop,nAge=nAge, POP_ID=POP_ID, alpha=alpha)
  
  #parameters = c("Linf","k","t0","sigma2","True","theta")  #True
  #parameters = c("Linf","k","t0","sigma2","True","alpha")  #True
  # using dclone
  cl <- makeCluster(n.chains, "SOCK")
  m2 <- jags.parfit(cl, data = data, params = parameters, model = LVBLOGPOPpred,
                    inits = Inits, n.adapt =n.adapt, n.update = n.update,
                    n.iter = n.iter, thin = n.thin, n.chains = n.chains)
  stopCluster(cl)
  
  
  #jags.LVB = jags(data,Inits,parameters,LVBmodel,n.chains=n.chains,n.iter=n.iter,n.burnin=n.burnin,n.thin=n.thin,DIC=T) 
  #jags.LVB = do.call(jags.parallel, 
  #                   list(data=data,inits=Inits,parameters.to.save=parameters,model.file=LVBmodel,
  #                        n.chains=n.chains,n.iter=n.iter,n.burnin=n.burnin,n.thin=n.thin,DIC=T)
  #)
  detach(Data)
  # return(jags.LVB)
  return(m2)
} # End Jags code




# ============================ Age function =========================
AgePred <- function(Data, alpha){
  # Derive data
  nPop = length(unique(Data$REF_ID))
  N=length(Data[,1])
  nK <- nAge <-  6 #manual for now
  
  # Model Specs
  n.chains=2
  n.adapt=250
  n.update=250 
  n.iter=2000
  n.thin=4
  
  AGEmodel <- function(){ 
      for (i in 1:N){
        True[i] ~ dcat(theta[1:nK,(Age[i]),POP_ID[i]])
      }
      #Prior
      for(i in 1:nPop){
        for(k in 1:nAge){  
          theta[1:nK,k,i]~ ddirch(alpha[1:nK,k,i])          
        }   
      }  
  }# close AGEmodel
  
  inits <- function(){
    list(True = rep(2,length=N), theta=structure(.Data=rep(1/nK,length=nK*nK*nPop),.Dim=c(nK,nK,nPop)))
    list(True = rep(3,length=N), theta=structure(.Data=rep(1/nK,length=nK*nK*nPop),.Dim=c(nK,nK,nPop)))
  }
  
  data <- list(alpha=alpha, N=N, nK=nK, nPop=nPop, nAge=nAge, Age=Data$Age, POP_ID=Data$REF_ID) #************* REF_ID
  parameters <- c("True", "theta")
  cl <- makeCluster(n.chains, "SOCK")
  m2 <- jags.parfit(cl, data = data, params = parameters, model = AGEmodel,
                    inits = inits, n.adapt =n.adapt, n.update = n.update,
                    n.iter = n.iter, thin = n.thin, n.chains = n.chains)
  stopCluster(cl)
  return(m2)
}# close AgePred


# ============================ Age function =========================
AgePred_noprior <- function(Data, theta){
  # Derive data
  nPop = length(unique(Data$REF_ID))
  N=length(Data[,1])
  nK <- nAge <-  6 #manual for now
  
  # Model Specs
  n.chains=2
  n.adapt=250
  n.update=250 
  n.iter=2000
  n.thin=4
  
  AGEmodel <- function(){ 
    for (i in 1:N){
      True[i] ~ dcat(theta[1:nK,(Age[i]),POP_ID[i]])
    }
    
  }# close AGEmodel
  
  inits <- function(){
    list(True = rep(2,length=N))
    list(True = rep(3,length=N))
  }
  
  data <- list(theta=theta, N=N, nK=nK, nPop=nPop, nAge=nAge, Age=Data$Age, POP_ID=Data$REF_ID)
  parameters <- c("True")
  cl <- makeCluster(n.chains, "SOCK")
  m2 <- jags.parfit(cl, data = data, params = parameters, model = AGEmodel,
                    inits = inits, n.adapt =n.adapt, n.update = n.update,
                    n.iter = n.iter, thin = n.thin, n.chains = n.chains)
  stopCluster(cl)
  return(m2)
}# close AgePred


